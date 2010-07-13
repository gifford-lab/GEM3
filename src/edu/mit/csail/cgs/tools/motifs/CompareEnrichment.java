package edu.mit.csail.cgs.tools.motifs;

import java.util.*;
import java.io.*;
import java.sql.*;
import java.text.DecimalFormat;
import cern.jet.random.Binomial;
import cern.jet.random.engine.RandomEngine;

import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.io.parsing.FASTAStream;
import edu.mit.csail.cgs.utils.*;
//import edu.mit.csail.cgs.utils.probability.Binomial;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.motifs.*;
import edu.mit.csail.cgs.tools.motifs.*;
import edu.mit.csail.cgs.tools.utils.Args;

/** Compare the frequencies of a set of motifs between two FASTA files
 * usage:
 *   java edu.mit.csail.cgs.tools.motifs.CompareEnrichment --species "$SC;SGDv2" --first foo.fasta --second bar.fasta
 *
 * can also specify --accept to give a regex which the motif name must match
 * of --reject to specify a regex that the motif name must not match.  Remember the regex must match
 * the *entire* name, so use something like Hnf.*
 *
 * --cutoff .7 minimum percent (specify between 0 and 1) match to maximum motif score that counts as a match.
 * --filtersig .001 maximum pvalue for reporting an enrichment between the two files
 * --minfoldchange 1
 * --minfrac 0   minimum fraction of the sequences that must contain the motif (can be in either file)
 *
 * The comparison code will check all percent cutoffs between the value you specify as --cutoff and 1 (in increments of .05) 
 * to find the most significant threshold that also meets the other criteria.
 *
 * Output columns are
 * 1) log foldchange in frequency between the two files
 * 2) motif count in first set
 * 3) size of first set 
 * 4) frequency in first set
 * 5) motif count in second set
 * 6) size of second set
 * 7) frequency in second set
 * 8) pvalue of first count given second frequency
 * 9) motif name
 * 10) motif version
 * 11) percent cutoff used
 *  
 */

public class CompareEnrichment {

    public static final double step = .05;
    double cutoffpercent, minfrac, minfoldchange, filtersig, maxbackfrac;
    Collection<WeightMatrix> matrices;
    Map<String,char[]> foreground, background;
    Binomial binomial;

    public static Map<String,char[]> readFasta(BufferedReader reader) throws IOException {
        FASTAStream stream = new FASTAStream(reader);
        Map<String,char[]> output = new HashMap<String,char[]>();
        while (stream.hasNext()) {
            Pair<String,String> pair = stream.next();
            output.put(pair.car(), pair.cdr().toCharArray());
        }
        return output;
    }

    public static int countMeetsThreshold(Map<String,List<WMHit>> hits,
                                          double t) {
        int count = 0;
        for (String s : hits.keySet()) {
            List<WMHit> list = hits.get(s);
            for (int i = 0; i < list.size(); i++) {
                if (list.get(i).getScore() > t) {
                    count++;
                    break;
                }
            }
        }
        return count;
    }

    public   CEResult doScan(WeightMatrix matrix,
                             Map<String,char[]> fg, 
                             Map<String,char[]> bg) {
        Map<String,List<WMHit>> fghits = new HashMap<String,List<WMHit>>();
        Map<String,List<WMHit>> bghits = new HashMap<String,List<WMHit>>();
        double maxscore = matrix.getMaxScore();
        for (String s : fg.keySet()) {
            List<WMHit> hits = WeightMatrixScanner.scanSequence(matrix,
                                                                (float)(maxscore * cutoffpercent),
                                                                fg.get(s));
            fghits.put(s, hits);
        }
        for (String s : bg.keySet()) {
            List<WMHit> hits = WeightMatrixScanner.scanSequence(matrix,
                                                                (float)(maxscore * cutoffpercent),
                                                                bg.get(s));
            bghits.put(s, hits);
        }
        double percent = cutoffpercent;
        int fgsize = fg.size();
        int bgsize = bg.size();
        CEResult result = new CEResult();
        result.pval = 1.0;
        result.matrix = matrix;
        result.logfoldchange = 0;
        result.sizeone = fgsize;
        result.sizetwo = bgsize;
        while (percent <= 1.0) {
            percent += step;
            double t = maxscore * percent;
            int fgcount = countMeetsThreshold(fghits, t);
            int bgcount = countMeetsThreshold(bghits, t);
            if (bgcount == 0) {
                continue;
            }

            double thetaone = ((double)fgcount) / ((double)fgsize);
            double thetatwo = ((double)bgcount) / ((double)bgsize);
            //            double pval = Math.exp(Binomial.log_binomial_significance(fgcount, fgsize, thetatwo));
            if (fgsize <= 0 || thetatwo <= 0 || thetatwo >= 1) {
                continue;
            }

            //            System.err.println(String.format("Setting %d and %f for %d", fgsize, thetatwo, fgcount));
            binomial.setNandP(fgsize, thetatwo);
            double pval = 1 - binomial.cdf(fgcount);
            double fc = Math.log(thetaone / thetatwo);
            if (pval <= filtersig && 
                pval <= result.pval && 
                Math.abs(fc) >= Math.abs(result.logfoldchange) &&
                Math.abs(fc) >= Math.abs(Math.log(minfoldchange)) &&
                (thetatwo < maxbackfrac) && 
                (thetaone >= minfrac || thetatwo >= minfrac)) {
                result.pval = pval;
                result.percentString = Double.toString(percent);
                result.cutoffString = Double.toString(t);
                result.countone = fgcount;
                result.counttwo = bgcount;
                result.logfoldchange = fc;
                result.freqone = thetaone;
                result.freqtwo = thetatwo;
                //                System.err.println(String.format("Accepted %f, %d and %d give %f and %f. fc %f.  thresh %f ", pval, fgcount, bgcount, thetaone, thetatwo, fc, t));
            } else {
                //                System.err.println(String.format("Rejected %f, %d and %d give %f and %f. fc %f.  thresh %f ", pval, fgcount, bgcount, thetaone, thetatwo, fc, t));
            }
        }
        return result;
    }

    public void parseArgs(String args[]) throws Exception {
        binomial = new Binomial(100, .01, RandomEngine.makeDefault());
        String firstfname = null, secondfname = null;
        Collection<String> accept, reject;

        cutoffpercent = Args.parseDouble(args,"cutoff",.3);
        filtersig = Args.parseDouble(args,"filtersig",.001);
        minfoldchange = Args.parseDouble(args,"minfoldchange",1);
        minfrac = Args.parseDouble(args,"minfrac",0);
        maxbackfrac = Args.parseDouble(args,"maxbackfrac",1.0);
        firstfname = Args.parseString(args,"first",null);
        secondfname = Args.parseString(args,"second",null);
        MarkovBackgroundModel bgModel = null;
        String bgmodelname = Args.parseString(args,"bgmodel","whole genome zero order");
        BackgroundModelMetadata md = BackgroundModelLoader.getBackgroundModel(bgmodelname,
                                                                              1,
                                                                              "MARKOV",
                                                                              Args.parseGenome(args).cdr().getDBID());

        System.err.println("Need to reimplement minfrac");

        if (md != null) {
            bgModel = BackgroundModelLoader.getMarkovModel(md);
        } else {
            System.err.println("Couldn't get metadata for " + bgmodelname);
        }
                
        matrices = Args.parseWeightMatrices(args);
        if (bgModel == null) {
            for (WeightMatrix m : matrices) {
                m.toLogOdds();
            }            
        } else {
            for (WeightMatrix m : matrices) {
                m.toLogOdds(bgModel);
            }            
        }
        if (firstfname == null) {
            throw new RuntimeException("Must supply a --first");
        }
        if (secondfname == null) {
            throw new RuntimeException("Must supply a --second");
        }
        System.err.println("Going to scan for " + matrices.size() + " matrices");
        foreground = readFasta(new BufferedReader(new FileReader(firstfname)));
        background = readFasta(new BufferedReader(new FileReader(secondfname)));
    }
    
    public void doScan() {
        DecimalFormat nf = new DecimalFormat("0.000E000");
        for (WeightMatrix matrix : matrices) {            
            CEResult result = doScan(matrix,
                                     foreground, 
                                     background);
            if (result.pval <= filtersig && 
                Math.abs(result.logfoldchange) >= Math.abs(Math.log(minfoldchange)) &&
                (result.freqtwo <= maxbackfrac) && 
                (result.freqone >= minfrac || result.freqtwo >= minfrac)) {
                System.out.println(String.format("%2.1f",result.logfoldchange) + "\t" + 
                                   result.countone + "\t" + 
                                   foreground.size() + "\t" + 
                                   nf.format(result.freqone) + "\t" + 
                                   result.counttwo + "\t" + 
                                   background.size() + "\t" + 
                                   nf.format(result.freqtwo) + "\t" +
                                   nf.format(result.pval) + "\t" + 
                                   result.matrix.name + "\t" + 
                                   result.matrix.version + "\t" + result.percentString + "\t" + result.cutoffString);
            }
        }
    }

    public static void main(String args[]) throws Exception {
        CompareEnrichment ce = new CompareEnrichment();
        ce.parseArgs(args);
        ce.doScan();
    }

}

class CEResult {

    public double pval, logfoldchange, freqone, freqtwo;
    public WeightMatrix matrix;
    public int sizeone, sizetwo, countone, counttwo;
    public String percentString, cutoffString;
    
}
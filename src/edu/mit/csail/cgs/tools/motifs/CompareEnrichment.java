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
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.motifs.*;
import edu.mit.csail.cgs.tools.motifs.*;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.ewok.verbs.*;

/** Compare the frequencies of a set of motifs between two FASTA files
 * usage:
 *   java edu.mit.csail.cgs.tools.motifs.CompareEnrichment --species "$SC;SGDv2" [--first foo.fasta] [--second bar.fasta]
 *
 * can also specify --acceptwm to give a regex which the motif name must match
 * of --rejectwm to specify a regex that the motif name must not match.  Remember the regex must match
 * the *entire* name, so use something like Hnf.*
 *
 * If --first is not supplied, then regions are read from STDIN.  If --second is not supplied, then random 
 * regions are chosen according to --randombgcount and --randombgsize.  
 *
 * Input files ending in .fasta or .fa are parsed as fasta.  Otherwise, they're parsed as a list of regions.
 *
 * --cutoff .7 minimum percent (specify between 0 and 1) match to maximum motif score that counts as a match.
 * --filtersig .001 maximum pvalue for reporting an enrichment between the two files
 * --minfoldchange 1
 * --minfrac 0   minimum fraction of the sequences that must contain the motif (can be in either file)
 * --savedata motif_presence_file.txt
 * --savedatahits
 * --dumpfg output_fg.fasta
 * --dumpbg output_bg.fasta
 * --mask name;version;cutoff
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
    int randombgcount = 1000, randombgsize = 100, parsedregionexpand;
    Genome genome;
    double cutoffpercent, minfrac, minfoldchange, filtersig, maxbackfrac;
    Collection<WeightMatrix> matrices;
    Map<String,char[]> foreground, background;
    Map<WeightMatrix, Double> maskingMatrices;
    ArrayList<String> fgkeys;
    Binomial binomial;
    PrintWriter savedata = null, outfg = null, outbg = null;
    boolean savedatahits = false, matchedbg = false;

    public static void saveFasta(PrintWriter pw, Map<String, char[]> seqs) throws IOException {
        for (String s : seqs.keySet()) {
            pw.println(">" + s);
            int i = 0;
            char[] chars = seqs.get(s);
            while (i < chars.length) {
                int l = i + 60 < chars.length ? 60 : chars.length -i;
                pw.write(chars,i,l);
                i += 60;
                pw.println();
            }
        }
    }
    public static Map<String,char[]> readFasta(BufferedReader reader) throws IOException {
        FASTAStream stream = new FASTAStream(reader);
        Map<String,char[]> output = new HashMap<String,char[]>();
        while (stream.hasNext()) {
            Pair<String,String> pair = stream.next();
            output.put(pair.car(), pair.cdr().toCharArray());
        }
        return output;
    }
    public static Map<String,char[]> readRegions(Genome g, BufferedReader reader, int parsedregionexpand, Map<String,char[]> matchedRegions) throws IOException, NotFoundException {
        String line = null;
        Map<String,char[]> output = new HashMap<String,char[]>();
        SequenceGenerator seqgen = new SequenceGenerator();
        while ((line = reader.readLine()) != null) {
            Region region = Region.fromString(g, line);
            if (parsedregionexpand > 0) {
                region = region.expand(parsedregionexpand,parsedregionexpand);
            }
            output.put(line,
                       seqgen.execute(region).toCharArray());
            if (matchedRegions != null) {
                Region before = new Region(region.getGenome(),
                                           region.getChrom(),
                                           region.getStart() - region.getWidth(),
                                           region.getStart());
                Region after = new Region(region.getGenome(),
                                          region.getChrom(),
                                          region.getEnd(),
                                          region.getEnd() + region.getWidth());
                matchedRegions.put(before.toString(),
                                   seqgen.execute(before).toCharArray());
                matchedRegions.put(after.toString(),
                                   seqgen.execute(after).toCharArray());
            }                      
        }
        return output;
    }
    public static Map<String,char[]> randomRegions(Genome genome, int count, int size) {
        Map<String,char[]> output = new HashMap<String,char[]>();
        SequenceGenerator seqgen = new SequenceGenerator();
        Map<String,Integer> chromlengthmap = genome.getChromLengthMap();
        ArrayList<String> chromnames = new ArrayList<String>();
        ArrayList<Integer> chromlengths = new ArrayList<Integer>();
        chromnames.addAll(chromlengthmap.keySet());
        long totallength = 0;
        for (String s : chromnames) {
            totallength += chromlengthmap.get(s);
            chromlengths.add(chromlengthmap.get(s));
            
        }
        while (count-- > 0) {
            long target = (long)(Math.random() * totallength);
            for (int i = 0; i < chromlengths.size(); i++) {
                if (i == chromlengths.size() - 1 ||
                    target + size < chromlengths.get(i)) {
                    Region r = new Region(genome,
                                          chromnames.get(i),
                                          (int)target,
                                          (int)target+size);
                    String seq = seqgen.execute(r);
                    output.put(r.toString(),seq.toCharArray() );
                    break;
                } else {
                    target -= chromlengths.get(i);
                    if (target < 0) {
                        target = 10;
                    }
                }

            }
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
        for (String s : fgkeys) {
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
        double bestthresh = maxscore * percent;
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
            if (fgsize <= 0 || thetatwo <= 0 || thetatwo >= 1) {
                continue;
            }

            binomial.setNandP(fgsize, thetatwo);
            double pval = 1 - binomial.cdf(fgcount);
            double fc = Math.log(thetaone / thetatwo);
            if (pval <= filtersig && 
                pval <= result.pval && 
                Math.abs(fc) >= Math.abs(result.logfoldchange) &&
                Math.abs(fc) >= Math.abs(Math.log(minfoldchange)) &&
                (thetatwo < maxbackfrac) && 
                (thetaone >= minfrac || thetatwo >= minfrac)) {
                bestthresh = t;
                result.pval = pval;
                result.percentString = Double.toString(percent);
                result.cutoffString = Double.toString(t);
                result.countone = fgcount;
                result.counttwo = bgcount;
                result.logfoldchange = fc;
                result.freqone = thetaone;
                result.freqtwo = thetatwo;
            } 
        }
        if (savedata != null) {
            if (savedatahits) {
                ArrayList<WMHit> toprint = new ArrayList<WMHit>();
                for (String s : fgkeys) {
                    toprint.clear();
                    for (WMHit hit : fghits.get(s)) {
                        if (hit.getScore() >= bestthresh) {
                            toprint.add(hit);
                        }
                    }
                    savedata.print("\t" + toprint);
                }

            } else {
                for (String s : fgkeys) {
                    int count = 0;
                    for (WMHit hit : fghits.get(s)) {
                        if (hit.getScore() >= bestthresh) {
                            count++;
                        }
                    }
                    savedata.print("\t" + count);
                }
            }
        }  
        return result;
    }

    public void parseArgs(String args[]) throws Exception {
        binomial = new Binomial(100, .01, RandomEngine.makeDefault());
        String firstfname = null, secondfname = null;
        Collection<String> accept, reject;

        genome = Args.parseGenome(args).cdr();
        cutoffpercent = Args.parseDouble(args,"cutoff",.3);
        filtersig = Args.parseDouble(args,"filtersig",.001);
        minfoldchange = Args.parseDouble(args,"minfoldchange",1);
        minfrac = Args.parseDouble(args,"minfrac",0);
        maxbackfrac = Args.parseDouble(args,"maxbackfrac",1.0);
        firstfname = Args.parseString(args,"first",null);
        secondfname = Args.parseString(args,"second",null);
        randombgcount = Args.parseInteger(args,"randombgcount",1000);
        randombgsize = Args.parseInteger(args,"randombgsize",100);
        parsedregionexpand = Args.parseInteger(args,"parsedregionexpand",0);
        String savefile = Args.parseString(args,"savedata",null);
        savedatahits = Args.parseFlags(args).contains("savedatahits");
        matchedbg = Args.parseFlags(args).contains("matchedbg");
        String outfgfile = Args.parseString(args,"outfg",null);
        String outbgfile = Args.parseString(args,"outbg",null);
        if (savefile != null) {
            savedata = new PrintWriter(savefile);
        }
        if (outfgfile != null) {
            outfg = new PrintWriter(outfgfile);
        }
        if (outbgfile != null) {
            outbg = new PrintWriter(outbgfile);
        }
        maskingMatrices = new HashMap<WeightMatrix,Double>();

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

        WeightMatrixLoader wmloader = new WeightMatrixLoader();
        Collection<String> maskstrings = Args.parseStrings(args,"mask");
        for (String maskstring : maskstrings) {
            String[] pieces = maskstring.split(";");
            String name = pieces[0];
            String version = "";
            for (int i = 1; i < pieces.length - 1; i++) {
                version += (i > 1 ? ";" : "") + pieces[i];
            }

            Double threshold = Double.parseDouble(pieces[pieces.length - 1]);
            Collection<WeightMatrix> matrices = wmloader.query(name,version,null);
            for (WeightMatrix m : matrices) {
                if (threshold < 1) {
                    maskingMatrices.put(m, threshold * m.getMaxScore());
                } else {
                    maskingMatrices.put(m, threshold);
                }


            }


        }
        wmloader.close();

        System.err.println("Going to scan for " + matrices.size() + " matrices");
        if (matchedbg) {
            System.err.println("Using matched background");
            background = new HashMap<String,char[]>();
        }

        if (firstfname == null) {
            System.err.println("No --first specified.  Reading from stdin");
            foreground = readRegions(genome,
                                     new BufferedReader(new InputStreamReader(System.in)), 
                                     parsedregionexpand, 
                                     matchedbg ? background : null);
        } else {
            if (firstfname.matches(".*\\.fasta") ||
                firstfname.matches(".*\\.fa")) {
                foreground = readFasta(new BufferedReader(new FileReader(firstfname)));
            } else {
                foreground = readRegions(genome,
                                         new BufferedReader(new FileReader(firstfname)), 
                                         parsedregionexpand,
                                         matchedbg ? background : null);
            }
        }
        if (!matchedbg) {
            if (secondfname == null) {
                System.err.println("No background file given.  Generating " + randombgcount + " regions of size " + randombgsize);
                background = randomRegions(genome,
                                           randombgcount,
                                           randombgsize);
            } else {
                if (secondfname.matches(".*\\.fasta") ||
                    secondfname.matches(".*\\.fa")){
                    background = readFasta(new BufferedReader(new FileReader(secondfname)));
                } else {
                    background = readRegions(genome, new BufferedReader(new FileReader(secondfname)), parsedregionexpand,null);
                }
            }
        }
        if (outfg != null) {
            saveFasta(outfg, foreground);
            outfg.close();
            outfg = null;
        }
        if (outbg != null) {
            saveFasta(outbg, background);
            outbg.close();
            outbg = null;
        }
    }
    public void maskSequence() {
        for (WeightMatrix wm : maskingMatrices.keySet()) {
            double threshold = maskingMatrices.get(wm);
            for (String s : foreground.keySet()) {
                char[] seq = foreground.get(s);
                List<WMHit> hits = WeightMatrixScanner.scanSequence(wm,
                                                                    (float)threshold,
                                                                    seq);
                for (WMHit hit : hits) {
                    for (int i = hit.start; i < hit.end; i++) {
                        seq[i] = 'X';
                    }
                }
            }
            for (String s : background.keySet()) {
                char[] seq = background.get(s);
                List<WMHit> hits = WeightMatrixScanner.scanSequence(wm,
                                                                    (float)threshold,
                                                                    seq);
                for (WMHit hit : hits) {
                    for (int i = hit.start; i < hit.end; i++) {
                        seq[i] = 'X';
                    }
                }
            }

        }

    }
    public void doScan() {
        DecimalFormat nf = new DecimalFormat("0.000E000");
        fgkeys = new ArrayList<String>();
        fgkeys.addAll(foreground.keySet());
        Collections.sort(fgkeys);
        if (savedata != null) {
            savedata.print("Motif");
            for (String s : fgkeys) {
                savedata.print("\t" + s);
            }
            savedata.println();
        }

        for (WeightMatrix matrix : matrices) {            
            if (savedata != null) {
                savedata.print(matrix.toString());
            }

            CEResult result = doScan(matrix,
                                     foreground, 
                                     background);
            if (savedata != null) {
                savedata.println();
            }
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
        if (savedata != null) {
            savedata.close();
        }

    }

    public static void main(String args[]) throws Exception {
        CompareEnrichment ce = new CompareEnrichment();
        ce.parseArgs(args);
        ce.maskSequence();
        ce.doScan();
    }

}

class CEResult {

    public double pval, logfoldchange, freqone, freqtwo;
    public WeightMatrix matrix;
    public int sizeone, sizetwo, countone, counttwo;
    public String percentString, cutoffString;
    
}
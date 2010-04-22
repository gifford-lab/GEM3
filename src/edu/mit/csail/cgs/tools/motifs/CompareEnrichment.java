package edu.mit.csail.cgs.tools.motifs;

import java.util.*;
import java.io.*;
import java.sql.*;
import java.text.DecimalFormat;

import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.io.parsing.FASTAStream;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.probability.Binomial;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.motifs.*;
import edu.mit.csail.cgs.tools.motifs.*;
import edu.mit.csail.cgs.tools.utils.Args;

/** Compare the frequencies of a set of motifs between two FASTA files
 * usage:
 *   java edu.mit.csail.cgs.tools.motifs.CompareEnrichment --first foo.fasta --second bar.fasta
 *
 * can also specify --accept to give a regex which the motif name must match
 * of --reject to specify a regex that the motif name must not match.  Remember the regex must match
 * the *entire* name, so use something like Hnf.*
 *
 * --cutoff .7 minimum percent (specify between 0 and 1) match to maximum motif score that counts as a match.
 * --filtersig .001 maximum pvalue for reporting an enrichment between the two files
 * --perbase
 * --minfoldchange 1
 * --minfrac 0
 *
 * Output columns are
 * 1) foldchange in frequency
 * 2) motif count in first set
 * 3) size of first set 
 * 4) frequency in first set
 * 5) motif count in second set
 * 6) size of second set
 * 7) frequency in second set
 * 8) pvalue of first count given second frequency
 * 9) pvalue of second count given first frequency
 * 10) motif name
 * 11) motif version
 *  
 */

public class CompareEnrichment {
    static Hashtable<WeightMatrix,Integer> firstmotifcounts, secondmotifcounts;
    static int firstseqcount, secondseqcount;
    static Collection<WeightMatrix> matrices;

    public static Pair<Integer,Hashtable<WeightMatrix,Integer>> motifCountsFromFASTA(String fname,
                                                                                     Collection<WeightMatrix> matrices,
                                                                                     double cutoffpercent,
                                                                                     boolean perbasecounts) throws FileNotFoundException, IOException {
        Hashtable<WeightMatrix,Integer> filecounts = new Hashtable<WeightMatrix,Integer>();
        Set<char[]> seqs = new HashSet<char[]>();

        FASTAStream stream = new FASTAStream(new File(fname));
        int totalbases = 0;
        while (stream.hasNext()) {
            Pair<String,String> pair = stream.next();
            String name = pair.getFirst();
            String seq = pair.getLast();
            char[] aschars = seq.toCharArray();
            seqs.add(aschars);
            totalbases += aschars.length;
        }
        Pair<Integer,Hashtable<WeightMatrix,Integer>> out;
        if (perbasecounts) {
            out = new Pair<Integer,Hashtable<WeightMatrix,Integer>>(totalbases,
                                                                    filecounts);
        } else {
            out = new Pair<Integer,Hashtable<WeightMatrix,Integer>>(seqs.size(),
                                                                    filecounts);
        }


        int done = 0;
        for (WeightMatrix m : matrices) {
            int count = 0;
            for (char[] seq : seqs) {
                List<WMHit> hits = WeightMatrixScanner.scanSequence(m,
                                                                    (float)(m.getMaxScore() * cutoffpercent),
                                                                    seq);
                if (perbasecounts) {
                    count += hits.size();
                } else {
                    if (hits.size() > 0) {
                        count++;
                    }
                }
            }
            filecounts.put(m,count);
        }
        return out;
    }

    public static Collection<WeightMatrix> filterMatrices(Collection<String> accepts,
                                                          Collection<String> rejects,
                                                          Collection<WeightMatrix> matrices) {
        ArrayList<WeightMatrix> out = new ArrayList<WeightMatrix>();
        for (WeightMatrix wm : matrices) {
            boolean reject = false;
            for (String s : rejects) {
                if (wm.name.matches(s)) {
                    reject = true;
                }
            }
            if (reject) {
                continue;
            }
            reject = false;
            if (accepts.size() > 0) {
                boolean any = false;
                for (String s : accepts) {
                    if (wm.name.matches(s)) {
                        any = true;
                    }
                }
                reject = !any;
            }
            if (reject) {
                continue;
            }
            out.add(wm);
        }
        return out;
    }

    public static void main(String args[]) throws Exception {
        String firstfname = null, secondfname = null;
        Collection<String> accept, reject;
        matrices = new ArrayList<WeightMatrix>();
        matrices.addAll(WeightMatrix.getAllWeightMatrices());

        double cutoffpercent = Args.parseDouble(args,"cutoff",.7);
        double filtersig = Args.parseDouble(args,"filtersig",.001);
        double minfoldchange = Args.parseDouble(args,"minfoldchange",1);
        double minfrac = Args.parseDouble(args,"minfrac",0);
        boolean perbasecounts = Args.parseFlags(args).contains("perbase");
        accept = Args.parseStrings(args,"accept");
        reject = Args.parseStrings(args,"reject");
        firstfname = Args.parseString(args,"first",null);
        secondfname = Args.parseString(args,"second",null);
        MarkovBackgroundModel bgModel = null;
        String bgmodelname = Args.parseString(args,"bgmodel","whole genome zero order");
        BackgroundModelMetadata md = BackgroundModelLoader.getBackgroundModel(bgmodelname,
                                                                              1,
                                                                              "MARKOV",
                                                                              Args.parseGenome(args).cdr().getDBID());
        if (md != null) {
            bgModel = BackgroundModelLoader.getMarkovModel(md);
        } else {
            System.err.println("Couldn't get metadata for " + bgmodelname);
        }
                
        matrices = filterMatrices(accept,reject,matrices);
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

        Pair<Integer,Hashtable<WeightMatrix,Integer>> pair = motifCountsFromFASTA(firstfname, matrices, cutoffpercent,perbasecounts);
        firstseqcount = pair.getFirst();
        firstmotifcounts = pair.getLast();
        pair = motifCountsFromFASTA(secondfname, matrices, cutoffpercent,perbasecounts);
        secondseqcount = pair.getFirst();
        secondmotifcounts = pair.getLast();

        DecimalFormat nf = new DecimalFormat("0.000E000");
        for (WeightMatrix wm : matrices) {            
            int first = firstmotifcounts.get(wm);
            int second = secondmotifcounts.get(wm);
            double firstfreq = ((double)first) / firstseqcount;
            double secondfreq = ((double)second) / secondseqcount;
            System.err.println(String.format("WM %s,%s : %d/%d, %d/%d",
                                             wm.getName(), wm.getVersion(), first, firstseqcount, second, secondseqcount));
            if (firstfreq <= 0) {
                continue;
                //                throw new RuntimeException("firstfreq < 0 from " + first + "," + firstseqcount);
            }
            if (secondfreq <= 0) {
                continue;
                //                throw new RuntimeException("secondfreq < 0 from " + second + "," + secondseqcount);
            }
            if (firstfreq >= 1) {
                continue;
                //                throw new RuntimeException("firstfreq > 1 from " + first + "," + firstseqcount);
            }
            if (secondfreq >= 1) {
                continue;
                //                throw new RuntimeException("secondfreq > 1 from " + second + "," + secondseqcount);
            }
            double pfirst = Binomial.log_binomial_significance(first,firstseqcount,secondfreq);
            double sigfirst = Math.exp(pfirst);
            sigfirst = Math.min(sigfirst, 1-sigfirst);
            double psecond = Binomial.log_binomial_significance(second,secondseqcount,firstfreq);
            double sigsecond = Math.exp(psecond);
            sigsecond = Math.min(sigsecond, 1-sigsecond);
             System.err.println("sig : " + sigfirst + ", " + sigsecond + " :: " + filtersig);
             System.err.println("fold : " + (pfirst / psecond) + ", " + (psecond / pfirst) + " :: " + minfoldchange);
             System.err.println("frac : " + pfirst + ", " + psecond + " :: " + minfrac);
            if ((sigfirst < filtersig || sigsecond < filtersig) &&
                (minfoldchange < (firstfreq / secondfreq) || minfoldchange < (secondfreq / firstfreq)) &&
                (minfrac < firstfreq || minfrac < secondfreq)) {
                double foldchange = Math.max(firstfreq / secondfreq, secondfreq / firstfreq);
                System.out.println(nf.format(foldchange) + "\t" + first + "\t" + firstseqcount + "\t" + nf.format(firstfreq) + "\t" + 
                                   second + "\t" + secondseqcount + "\t" + nf.format(secondfreq) + "\t" +
                                   nf.format(sigfirst) + "\t" + nf.format(sigsecond) + "\t\t" + wm.name + "\t" + wm.version);
            }
        }

    }    
}

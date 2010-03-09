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
import edu.mit.csail.cgs.datasets.motifs.*;
import edu.mit.csail.cgs.tools.motifs.*;

/** usage:
 *   java edu.mit.csail.cgs.tools.motifs.CompareEnrichment --first foo.fasta --second bar.fasta
 *
 * can also specify --accept to give a regex which the motif name must match
 * of --reject to specify a regex that the motif name must not match
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
            char[] aschars = seq.toUpperCase().toCharArray();
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

    public static Collection<WeightMatrix> filterMatrices(List<String> accepts,
                                                          List<String> rejects,
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
        ArrayList<String> accept, reject;
        accept = new ArrayList<String>();
        reject = new ArrayList<String>();
        matrices = new ArrayList<WeightMatrix>();
        matrices.addAll(WeightMatrix.getAllWeightMatrices());

        double cutoffpercent = .7;
        double filtersig = .001;
        double minfoldchange = 1;
        double minfrac = 0;
        boolean perbasecounts = false;
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--first")) {
                firstfname = args[++i];
            }
            if (args[i].equals("--second")) {
                secondfname = args[++i];
            }
            if (args[i].equals("--accept")) {
                accept.add(args[++i]);
            }
            if (args[i].equals("--reject")) {
                reject.add(args[++i]);
            }
            if (args[i].equals("--cutoff")) {
                cutoffpercent = Double.parseDouble(args[++i]);
            }
            if (args[i].equals("--filtersig")) {
                filtersig = Double.parseDouble(args[++i]);
            }
            if (args[i].equals("--perbase")) {
                perbasecounts = true;
            }
            if (args[i].equals("--minfoldchange")) {
                minfoldchange = Double.parseDouble(args[++i]);
            }
            if (args[i].equals("--minfrac")) {
                minfrac = Double.parseDouble(args[++i]);
            }


        }
        matrices = filterMatrices(accept,reject,matrices);
        if (firstfname == null) {
            throw new RuntimeException("Must supply a --first");
        }
        if (secondfname == null) {
            throw new RuntimeException("Must supply a --second");
        }

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
//             System.err.println("sig : " + sigfirst + ", " + sigsecond + " :: " + filtersig);
//             System.err.println("fold : " + (pfirst / psecond) + ", " + (psecond / pfirst) + " :: " + minfoldchange);
//             System.err.println("frac : " + pfirst + ", " + psecond + " :: " + minfrac);
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

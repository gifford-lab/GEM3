package edu.mit.csail.cgs.tools.chipseq;

import java.util.*;
import java.io.*;
import cern.jet.random.ChiSquare;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.parsing.FASTAStream;
import edu.mit.csail.cgs.tools.motifs.CountKmers;


/**
 * Compare the frequency of kmers in two FASTA files.  
 *
 * java edu.mit.csail.cgs.tools.chipseq.CompareKmers --k 6 --one foo.fasta --two bar.fasta
 *
 */


public class CompareKmers {
    
    private FASTAStream fastaone, fastatwo;
    private int k;
    private CountKmers one, two;
    
    public static void main(String args[]) throws Exception {
        CompareKmers compare = new CompareKmers();
        compare.parseArgs(args);
        compare.read();
        compare.report();
    }

    public void parseArgs(String args[]) throws IOException {
        String fone = Args.parseString(args,"one",null);
        String ftwo = Args.parseString(args,"two",null);
        k = Args.parseInteger(args,"k",6);
        fastaone = new FASTAStream(new File(fone));
        fastatwo = new FASTAStream(new File(ftwo));
        one = new CountKmers();
        one.init(k,k);
        two = new CountKmers();
        two.init(k,k);
    }

    /* read the input files and file CountKmers one and two */
    public void read() {
        read(fastaone, one);
        read(fastatwo, two);
    }
    public void read(FASTAStream stream, CountKmers kmers) {
        while (stream.hasNext()) {
            Pair<String,String> pair = stream.next();
            kmers.addToCounts(pair.cdr());
        }
    }
    public void report() {
        Set<String> allKeys = one.getKeySet(k);
        allKeys.addAll(two.getKeySet(k));
        double stat = 0;
        int toosmall = 0;
        int zero = 0, inzero = 0;
        double countone = 0, counttwo = 0;
        for (String kmer : allKeys) {
            int obs = one.getCount(kmer, k);
            int exp = two.getCount(kmer,k);
            countone += obs;
            counttwo += exp;
        }
        double scale = counttwo / countone;
        System.out.println(String.format("Count one is %f and two is %f", countone, counttwo));
        for (String kmer : allKeys) {
            double obs = one.getCount(kmer, k) * scale;
            double exp = two.getCount(kmer, k);
            
            if (exp < 5) {
                toosmall++;
                if (exp == 0) {
                    zero++;
                    inzero += obs;
                    continue;
                }
            }
            stat += (obs - exp) * (obs - exp) / exp;
        }
        int n = allKeys.size() - zero;

        cern.jet.random.engine.RandomEngine engine = new cern.jet.random.engine.DRand();
        ChiSquare csquare = new ChiSquare(n - 1, engine);
        System.out.println(String.format("chi-squared value is %f with %d DOF.  CDF is %f",
                                         stat, n-1, csquare.cdf(stat)));
        System.out.println(String.format("There were %d below 5 and %d that were zero (%d observations)",
                                         toosmall, zero, inzero));

    }


}
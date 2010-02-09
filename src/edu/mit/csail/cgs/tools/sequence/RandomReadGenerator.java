package edu.mit.csail.cgs.tools.sequence;

import edu.mit.csail.cgs.utils.io.parsing.FASTAStream;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.tools.utils.Args;
import cern.jet.random.*;
import java.io.*;

public class RandomReadGenerator {

    public static final char[] letters = {'A','C','T','G'};

    public static void main(String args[]) throws IOException {
        String fastafile = Args.parseString(args,"fasta",null);
        /* uniform probability of error at each base */
        double errorP = Args.parseDouble(args,"uniformerrorp",0);
        /* in solexa model, error probability is
           base + multiplier * ((exp) ^ length)
           
           Based on the data in Dohm et al (2008), 
           base = .005, multiplier = .0001, and exp = 1.2
           their data only went to 32bp, so this may not be accurate for long reads
        */
        double solexaErrorBaseP = Args.parseDouble(args,"solexabase",0);
        double solexaErrorMultiplier = Args.parseDouble(args,"solexamultiplier",0);
        double solexaErrorExp = Args.parseDouble(args,"solexaexp",0);
        // read length
        int n = Args.parseInteger(args, "n", 25);
        /* can specify either absolute number of reads per sequence in the FASTA file
         * or an amount of coverage 
         */
        int nreads = Args.parseInteger(args,"nreads",100);
        double coverage = Args.parseDouble(args,"coverage",Double.NaN);
        
        BufferedReader reader;
        if (fastafile == null || fastafile.equals("-")) {
            reader = new BufferedReader(new InputStreamReader(System.in));
        } else {
            reader = new BufferedReader(new InputStreamReader(new FileInputStream(fastafile)));
        }

        cern.jet.random.engine.RandomEngine engine = new cern.jet.random.engine.DRand();
        cern.jet.random.Normal norm = new cern.jet.random.Normal(0,coverage,engine);
        cern.jet.random.Uniform uniform = new cern.jet.random.Uniform(0.0,1.0,engine);
        cern.jet.random.Uniform uniform4 = new cern.jet.random.Uniform(0.0,3,engine);

        FASTAStream stream = new FASTAStream(reader);
        int readcount = 0;
        while (stream.hasNext()) {
            String seq = stream.next().cdr();
            int reads;
            if (Double.isNaN(coverage)) {
                reads = nreads;
            } else {
                reads = (int) (seq.length() * norm.nextDouble());
            }

            for (int i = 0; i < reads; i++) {
                int start = (int)((seq.length() - (n+1)) * Math.random());
                String name = String.format("read_%d",readcount++);
                char[] read = seq.substring(start,start+n).toCharArray();
                if (errorP > 0) {
                    for (int j = 0; j < read.length; j++) {
                        if (uniform.nextDouble() < errorP) {
                            char old = read[j];
                            while (old == (read[j] = letters[uniform.nextIntFromTo(0,3)])) { }
                        }
                    }
                } else if (solexaErrorBaseP > 0) {
                    for (int j = 0; j < read.length; j++) {
                        double p = solexaErrorBaseP + solexaErrorMultiplier * Math.pow(solexaErrorExp, j);                    
                        //                        System.err.println(String.format("%d -> %f",j,p));
                        if (j > 0 && read[j-1] == 'G') {
                            p *= 2;
                        }
                        if (uniform.nextDouble() < p) {
                            //                            System.err.println("Mutating base " + j + " of " + i + " from " + read[j]);
                            char old = read[j];
                            while (old == (read[j] = letters[uniform.nextIntFromTo(0,3)])) { }
                        }
                    }
                }

                
                System.out.println(">" + name + "\n" + String.valueOf(read));
            }
        }
    }
}
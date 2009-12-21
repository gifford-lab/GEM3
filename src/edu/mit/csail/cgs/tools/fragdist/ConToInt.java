package edu.mit.csail.cgs.tools.fragdist;

import java.util.*;
import java.io.*;
import cern.jet.random.Gamma;
import cern.jet.random.engine.RandomEngine;
import edu.mit.csail.cgs.tools.utils.Args;

/* convert .con file (maps length to concentration of molecule) to
   array intensity (relative units):
   java edu.mit.csail.cgs.tools.fragdist --file foo.con > foo.int

   If multiple input files are specified then the output is the result of
   averaging across the input files.

   The conversion works in the following steps:
   1) read the input file of fragment concentrations.  Since the lengths
   in this file are floats, bin the lengths by floor(length) and average
   at each length.
   2) smooth the averaged concentrations
   3) subtract off the primer length from the fragment lengths
   4) recenter the data so the distribution starts at zero; call the 
   amount of the shift "min".  This is important
   because the gamma distribution starts at zero so you learn the wrong parameters
   if your distribution starts at a higher value
   5) learn alpha and beta for the gamma distribution
   6) alpha and beta give you p(l - min)
   7) Use p(l) to compute I(d) as described in the JBD paper.
*/

public class ConToInt {

    private boolean nogamma, nodeconv;
    private int maxdist;
    private double alpha, beta;
    private double smoothed[], deconv[], intensities[];
    private List<String> inputFiles;
    private int primerlength;

    public ConToInt () {
    }

    public void parseArgs(String args[]) {
        /* nogamma: don't model the fragment length distribution as a gamma
           distribution.  Use the actual observed distribution instead */
        nogamma = Args.parseFlags(args).contains("nogamma");
        /* don't deconvolve the data. Just read the input file, generate
           the smoothed version, and output that */
        nodeconv= Args.parseFlags(args).contains("nodeconv");
        /* instead of using an input file to learn the gamma distribution parameters,
           you can specify them directly with alpha and beta */
        alpha = Args.parseDouble(args, "alpha", Double.NaN);
        beta = Args.parseDouble(args, "beta", Double.NaN);
        /* the longest fragment you expect to see */
        maxdist = Args.parseInteger(args,"maxdist",1000);        
        /* the length of the primers.  This is subtracted off the observed fragment
           lengths during processing */
        primerlength = Args.parseInteger(args,"primerlen",50);
        inputFiles = Args.parseFile(args);
    }
    public void parseFiles() throws IOException {        
        if (inputFiles.size() == 0) {
            return;
        }
        ArrayList<Double> sum = new ArrayList<Double>();
        ArrayList<Integer> count = new ArrayList<Integer>();
        for (String fname : inputFiles) {
            BufferedReader reader = new BufferedReader(new FileReader(fname));
            String line;
            while ((line = reader.readLine()) != null) {
                String pieces[] = line.split("\\s+");
                int length = (int)Math.floor(Double.parseDouble(pieces[0]));
                Double value = Double.parseDouble(pieces[1]);
                if (length >= sum.size()) {
                    for (int i = sum.size(); i <= length; i++) {
                        sum.add(0.0);
                        count.add(0);
                    }
                }
                sum.set(length, sum.get(length) + value);
                count.set(length, count.get(length) + 1);
            }
        }        
        /* avg[] is the values averaged by floor(length) */
        double[] avg = new double[sum.size()];
        for (int i = 0; i < avg.length; i++) {
            avg[i] = sum.get(i) / count.get(i);
            if (Double.isNaN(avg[i])) {
                avg[i] = 0;
            }
            System.err.println(String.format("avg[%d] = %.20f", i, avg[i]));
        }
        smoothed = new double[sum.size()];
        smoothed[0] = avg[0];
        smoothed[1] = avg[1];
        smoothed[smoothed.length - 1] = avg[avg.length - 1];
        smoothed[smoothed.length - 2] = avg[avg.length - 2];
        for (int i = 2; i < smoothed.length - 2; i++) {
            smoothed[i] = .2 * avg[i-1] + 
                .2 * avg[i+1] + 
                .5 * avg[i] + 
                .1 * avg[i-2] + 
                .1 * avg[i+2];
        }

        /* subtract the primer length (50) */
        for (int i = 0; i < smoothed.length - primerlength; i++) {
            smoothed[i] = smoothed[i+primerlength];
        }
        for (int i = smoothed.length - primerlength; i < smoothed.length; i++) {
            smoothed[i] = 0;
        }
        
    }
    
    public void outputSmoothed() {
        double max = 0;
        for (int i = 0; i < smoothed.length; i++) {
            if (smoothed[i] > max) {
                max = smoothed[i];
            }
        }
        for (int i = 0; i < smoothed.length; i++) {
            System.out.print(String.format("%i\t%d\n",i,smoothed[i]/max));
        }
    }

    /* normalizes smoothed[] sum that the sum of the values is one.
       This converts it from a bunch of values on whatever scale the
       bioanalyzer uses to a probability distribution */
    public void normalizeSmoothed() {
        double sum = 0;
        for (int i = 0; i < smoothed.length; i++) {
            sum += smoothed[i];
        }
        for (int i = 0; i < smoothed.length; i++) {
            smoothed[i] /= sum;
        }
    }

    public void deconv() {
        if (nogamma) {
            /* deconvolve the actual data: convolution of two gammas is a gamma of half the size
             */
            deconv = new double[smoothed.length / 2 + 1];
            for (int i = 0; i < deconv.length; i++) {
                deconv[i] = 0;
            }
            for (int i = 0; i < smoothed.length; i++) {
                deconv[i / 2] += smoothed[i];
            }
            double max = 0;
            for (int i = 0; i < deconv.length; i++) {
                if (deconv[i] > max) {
                    max = deconv[i];
                }
            }
            for (int i = 0; i < deconv.length; i++) {
                deconv[i] /= max;
            }
        } else {
            double mean, var;
            int min = 0;
            if (!Double.isNaN(alpha) && !Double.isNaN(beta)) {
                mean = alpha * beta;
                var = alpha * beta * beta;
            } else {
                mean = 0;
                /* The gamma distribution isn't a good fit if it doesn't
                   start at ~0.  So find the first index at which there's much
                   probability mass and then effectively shift smoothed[] by
                   that much
                */
                for (int i = 0; i < smoothed.length; i++) {
                    if (smoothed[i] < .00000001) {
                        smoothed[i] = 0;
                        min = i;
                    } else {
                        break;
                    }
                }
                for (int i = min; i < smoothed.length; i++) {
                    System.err.println(String.format("smoothed[%d] = %.20f", i, smoothed[i]));
                    mean += smoothed[i] * (i - min);
                }
                var = 0;
                for (int i = min; i < smoothed.length; i++) {
                    var += smoothed[i] * (i - min - mean) * (i - min - mean);
                }
                beta = var / mean;
                alpha = mean / beta;
            }
            System.err.println(String.format("min %d, mean %f, var %f, alpha %f, beta %f",
                                             min, mean, var, alpha, beta));
            alpha = alpha / 2;  // do the deconvolution
            min /= 2;
            deconv = new double[maxdist];
            for (int i = 0; i < min; i++) {
                deconv[i] = 0;
            }
            RandomEngine engine = RandomEngine.makeDefault();
            Gamma gamma = new Gamma(alpha, beta, engine);
            for (int i = min; i < deconv.length; i++) {
                deconv[i] = gamma.pdf(i - min);
            }            
        }
    }
    
    /* requires deconv and maxdist.  creates intensities[]
     */
    public void generateIntensities() {
        intensities = new double[maxdist];
        double maxint = 0;
        for (int i = 0; i < maxdist; i++) {
            intensities[i] = 0;            
            for (int l = i; l < maxdist && l < deconv.length; l++) {
                double sum = 0;
                for (int x = i; x < l; x++) {
                    sum += deconv[x] * deconv[l - x];
                }
                intensities[i] += sum;
                if (intensities[i] > maxint) {
                    maxint = intensities[i];
                }
            }
        }
        for (int i = 0; i < maxdist; i++) {
            intensities[i] /= maxint;
        }
    }
    public void printIntensities() {
        for (int i = 0; i < intensities.length; i++) {
            System.out.println(String.format("%d\t%e",
                                             i,
                                             intensities[i]));
        }
    }

    public static void main(String args[]) throws Exception {
        ConToInt contoint = new ConToInt();
        contoint.parseArgs(args);
        contoint.parseFiles();
        if (contoint.nodeconv) {
            contoint.outputSmoothed();
            System.exit(0);
        } 
        if (contoint.inputFiles.size() > 0) {
            contoint.normalizeSmoothed();
        }
        contoint.deconv();
        contoint.generateIntensities();
        contoint.printIntensities();
    }

}
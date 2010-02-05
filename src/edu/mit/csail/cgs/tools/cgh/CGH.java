package edu.mit.csail.cgs.tools.cgh;

import java.util.*;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.ewok.verbs.cgh.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.chipchip.*;

/** Command line program to analyze CGH data.  Uses a three state HMM to determine when the log-ratio
 * in a region is below, equal to, or above zero.
 *
 *  command line params are
 *     --expt : experiment
 *     --species : eg 'Mus musculus;mm8'
 *  and descriptions of the three states, all optional.  All provide mean and stddev of the log-ratios
 *  for that state
 *      --low 'mean;stddev'
 *      --middle 'mean;stddev'
 *      --high 'mean;stddev'
 */
public class CGH {

    public static void main(String args[]) throws Exception {
        Genome genome = Args.parseGenome(args).cdr();
        ChipChipDataset dataset = new ChipChipDataset(genome);
        List<Region> regions = Args.parseRegionsOrDefault(args);

        double means[] = new double[CGHCallExpander.numStates], stddevs[] = new double[CGHCallExpander.numStates];
        String[] pieces = Args.parseString(args,"low",String.format("%f;%f",
                                                                    CGHCallExpander.defaultParams[0][0],
                                                                    Math.sqrt(CGHCallExpander.defaultParams[0][1]))).split(";");
        means[CGHCallExpander.LOW] = Double.parseDouble(pieces[0]);
        stddevs[CGHCallExpander.LOW] = Math.pow(Double.parseDouble(pieces[1]),2);

        pieces = Args.parseString(args,"middle",String.format("%f;%f",
                                                              CGHCallExpander.defaultParams[0][0],
                                                              Math.sqrt(CGHCallExpander.defaultParams[0][1]))).split(";");
        means[CGHCallExpander.ONE] = Double.parseDouble(pieces[0]);
        stddevs[CGHCallExpander.ONE] = Math.pow(Double.parseDouble(pieces[1]),2);

        pieces = Args.parseString(args,"high",String.format("%f;%f",
                                                            CGHCallExpander.defaultParams[0][0],
                                                            Math.sqrt(CGHCallExpander.defaultParams[0][1]))).split(";");
        means[CGHCallExpander.HIGH] = Double.parseDouble(pieces[0]);
        stddevs[CGHCallExpander.HIGH] = Math.pow(Double.parseDouble(pieces[1]),2);

        for (ExptNameVersion env : Args.parseENV(args)) {
            ChipChipData ccd = dataset.getData(env);
            CGHCallExpander cgh = new CGHCallExpander(ccd);
            cgh.setStateParameters(CGHCallExpander.LOW, 
                                   means[CGHCallExpander.LOW],
                                   stddevs[CGHCallExpander.LOW]);
            cgh.setStateParameters(CGHCallExpander.ONE, 
                                   means[CGHCallExpander.ONE],
                                   stddevs[CGHCallExpander.ONE]);
            cgh.setStateParameters(CGHCallExpander.HIGH, 
                                   means[CGHCallExpander.HIGH],
                                   stddevs[CGHCallExpander.HIGH]);
            cgh.resetHMM();
            for (Region r : regions) {
                List<ScoredRegion> calls = cgh.executeList(r);
                for (ScoredRegion call : calls) {
                    System.out.println(call.getScore() + "\t" + call.regionString());
                }
            }
        }
    }
}
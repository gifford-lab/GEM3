package edu.mit.csail.cgs.tools.cgh;

import java.util.*;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.ewok.verbs.cgh.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.chipchip.*;
import edu.mit.csail.cgs.utils.database.*;

/** Command line program to analyze a single channel of CGH data (that'd be GH data).
 *   First computes the median intensity in the data and assumes this is the intensity 
 *   from a single copy.  Then uses an HMM to determine whether the intensities are
 *   from zero copies (~0 * median), 1 copy (median), 2 copies (2 * median), etc
 *
 *  command line params are
 *   --expt : experiment as "name;version;replicate"
 *   --species : eg "Mus musculues;mm8"
 *   --cy5 or --cy3 (default is cy5) : which channel to analyze
 */

public class CGHIntensity {

    public static void main(String args[]) throws Exception {
        Genome genome = Args.parseGenome(args).cdr();
        ChipChipDataset dataset = new ChipChipDataset(genome);
        ChipChipMetadataLoader loader = new ChipChipMetadataLoader();
        List<Region> regions = Args.parseRegionsOrDefault(args);
        boolean cy5 = !Args.parseFlags(args).contains("cy3");

        for (ExptNameVersion env : Args.parseENV(args)) {
            ChipChipData ccd = dataset.getData(env);
            Collection<Experiment> expts = loader.loadExperiment(env);
            CGHIntensityCallExpander expander;
            if (expts.size() == 1) {
                expander = new CGHIntensityCallExpander(expts.iterator().next().getDBID(),
                                                        ccd,
                                                        cy5);
            } else {
                System.err.println("Can't resolve " + env + " to a single experiment.  Skipping!");
                continue;
            }
            for (Region r : regions) {
                List<ScoredRegion> calls = expander.executeList(r);
                for (ScoredRegion call : calls) {
                    System.out.println(call.getScore() + "\t" + call.regionString() + "\t" + (call.getWidth()) + "bp");
                }
            }
        }
    }
}
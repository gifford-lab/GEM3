package edu.mit.csail.cgs.tools.chipseq;

import java.util.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.datasets.alignments.*;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.*;



/**
 * Tries to determine one or more control experiments to use for
 * one or more chipseq experiments.
 *
 * Usage:
 * java edu.mit.csail.cgs.tools.chipseq.GuessControlExpt --species 'Mus musculus;mm9' --chipseq 'Sing ES CTCF E14;bowtie_unique'
 */
public class GuessControlExpt {

    private static ChipSeqLoader loader;
    private static Collection<String> controlFactorNames;
    private static Collection<Integer> controlFactorIDs;
    private static MetadataLoader core;
    private static Genome genome;
    public static void main(String args[]) throws Exception {
        genome = Args.parseGenome(args).cdr();
        List<ChipSeqLocator> locators = Args.parseChipSeq(args);
        loader = new ChipSeqLoader(false);
        core = new MetadataLoader();
        controlFactorNames = new ArrayList<String>();
        controlFactorNames.add("WCE");
        controlFactorNames.add("IgG");
        controlFactorNames.add("Input");
        controlFactorNames.add("GFP");        
        controlFactorIDs = new ArrayList<Integer>();

        boolean strict = Args.parseFlags(args).contains("strict");

        for (String f : controlFactorNames) {
            Factor factor = core.findFactor(f);
            if (factor != null) {
                controlFactorIDs.add(factor.getDBID());
            }
        }


        Set<String> output = new HashSet<String>();
        for (ChipSeqLocator l : locators) {
            Collection<ChipSeqAlignment> alignments = new ArrayList<ChipSeqAlignment>();
            if (l.getReplicates().size() == 0) {
                alignments.addAll(loader.loadAlignments(l.getExptName(),
                                                        null,
                                                        l.getAlignName(),
                                                        null,null,null, 
                                                        genome));
            } else {
                alignments.addAll(loader.loadAlignments(l, genome));
            }
            for (ChipSeqAlignment a : alignments) {
                Collection<ChipSeqAlignment> controls = controlsForAlignment(a);
                if (controls.size() == 0) {
                    System.err.println("No controls for " + a);
                    if (strict) {
                        System.exit(1);
                    }
                }
                for (ChipSeqAlignment c : controls) {
                    output.add(c.getExpt().getName() + ";" + c.getName());
                }
            }
        }
        for (String a : output) {
            System.out.println(a);
        }

    }
    public static Collection<ChipSeqAlignment> controlsForAlignment(ChipSeqAlignment a) throws Exception {
        ArrayList<ChipSeqAlignment> output = new ArrayList<ChipSeqAlignment>();
        ChipSeqExpt expt = a.getExpt();
        Cells cells = expt.getCells();
        Condition cond = expt.getCondition();
        for (Integer i : controlFactorIDs) {
            output.addAll(loader.loadAlignments(null,null,null,
                                                i, cells.getDBID(), cond.getDBID(), 
                                                genome));
        }
        return output;
    }
    

}
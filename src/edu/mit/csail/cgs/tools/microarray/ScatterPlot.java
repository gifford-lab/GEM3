package edu.mit.csail.cgs.tools.microarray;

import edu.mit.csail.cgs.viz.scatter.*;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.datasets.locators.*;
import edu.mit.csail.cgs.datasets.chipchip.*;


public class ScatterPlot {
    public static void main(String args[]) throws Exception {
        Pair<Organism, Genome> pair = Args.parseGenome(args);
        Genome g = pair.cdr();
        MicroarrayScatterSetupFrame frame = new MicroarrayScatterSetupFrame(g);
        if (g != null) {
            ChipChipDataset dataset = new ChipChipDataset(g);
            frame.pane.setLogScale(Args.parseFlags(args).contains("log"));
            frame.pane.setMvsA(Args.parseFlags(args).contains("mvsa"));
            String expt = Args.parseString(args, "exptone", null);
            if (expt != null) {
                String pieces[] = expt.split(";");
                if (pieces.length == 3) {
                    ExptLocator locator = new ChipChipLocator(dataset, 
                                                              pieces[0],
                                                              pieces[1],
                                                              pieces[2]);
                    System.err.println("Adding locator for expt one " + locator);
                    frame.pane.setExptOne(locator);
                }
            }
            expt = Args.parseString(args, "expttwo", null);
            if (expt != null) {
                String pieces[] = expt.split(";");
                if (pieces.length == 3) {
                    ExptLocator locator = new ChipChipLocator(dataset, 
                                                              pieces[0],
                                                              pieces[1],
                                                              pieces[2]);
                    frame.pane.setExptTwo(locator);
                }
            }
        }

    }
}
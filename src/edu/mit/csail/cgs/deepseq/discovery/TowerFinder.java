package edu.mit.csail.cgs.deepseq.discovery;

import java.sql.SQLException;
import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.deepseq.*;
import edu.mit.csail.cgs.deepseq.features.Feature;
import edu.mit.csail.cgs.deepseq.utilities.BindingModelGenerator;
import edu.mit.csail.cgs.ewok.verbs.ChromosomeGenerator;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;

/**
 * Find towers in data.  You should run this on control channels
 * from chipseq experiments.  
 *
 * java edu.mit.csail.cgs.deepseq.discovery.TowerFinder --species "$MM;mm9" --rdbexpt "YL ES WCE V6.5_129-C57BL6;1;bowtie_unique"
 *
 * This is based on ChipSeqPeakFinder and
 * takes all the same command line args plus:
 * --strandimbalance 2.4 : only return towers in which one strand contains this many
 *                         times more reads than the other.
 * --mintowerwidth 100 : only return towers that are at least this wide
 * --maxtowerwidth 1000 : only return towers that are less than this wide
 */

public class TowerFinder extends ChipSeqPeakFinder {

    private double strandImbalanceThreshold;
    private int minTowerWidth, maxTowerWidth;

	public TowerFinder(DeepSeqExpt signal) {
		super(signal, null);
	}
    public TowerFinder(String args[]) {
        super(args);
        parseArgs(args);
    }
    public void parseArgs(String args[]) {
        strandImbalanceThreshold = Math.log(Args.parseDouble(args,"strandimbalance",1.0));
        System.err.println("SIB " + strandImbalanceThreshold);
        minTowerWidth = Args.parseInteger(args,"mintowerwidth",1);
        maxTowerWidth = Args.parseInteger(args,"maxtowerwidth",100000);
    }
    public static void main(String[] args) throws SQLException, NotFoundException {
        TowerFinder tf = new TowerFinder(args);
        List<Feature> towers = tf.execute();

        towers = tf.filterStrandImbalance(towers);
        towers = tf.filterTowerWidth(towers);

        tf.signalFeatures.clear();
        tf.signalFeatures.addAll(towers);

        tf.printFeatures(true);

        tf.cleanup();
    }

    public List<Feature> filterStrandImbalance(List<Feature> features) {
        List<Feature> output = new ArrayList<Feature>();
        for (Feature feature : features) {
            double totalPlus = 0.0;
            double totalMinus = 0.0;
            for (Float f : signal.loadStrandedBaseCounts(feature.coords, '+').cdr()) {
                totalPlus += f;
            }
            for (Float f : signal.loadStrandedBaseCounts(feature.coords, '-').cdr()) {
                totalMinus += f;
            }
            double imb = Math.abs(Math.log(totalPlus / totalMinus));
            if (imb >= strandImbalanceThreshold) {
                output.add(feature);
            }
        }
        System.err.println("Kept " + output.size() + " after strand filtering");
        return output;
    }
    public List<Feature> filterTowerWidth(List<Feature> features) {
        List<Feature> output = new ArrayList<Feature>();
        for (Feature f : features) {
            int w = f.coords.getWidth();
            if (w >= minTowerWidth && w <= maxTowerWidth) {
                output.add(f);
            }
        }
        System.err.println("Kept " + output.size() + " after width filtering");
        return output;
    }


}
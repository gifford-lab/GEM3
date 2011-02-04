package edu.mit.csail.cgs.tools.chipseq;

import java.sql.SQLException;
import java.util.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.*;

/**
 * Just report on the overlap between the two experiments
 */

public class SimpleAnalysisOverlap extends CompareTwoAnalyses {

    public static void main(String args[]) throws Exception {
        SimpleAnalysisOverlap sao = new SimpleAnalysisOverlap();
        sao.parseArgs(args);
        sao.printReport();
    }
    public void printReport() throws SQLException {
        List<ChipSeqAnalysisResult> one = getResultsOne();
        List<ChipSeqAnalysisResult> two = getResultsTwo();
        Collections.sort(one);
        int onefound = 0, twofound = 0;
        for (ChipSeqAnalysisResult r : two) {
            if (containsMatch(one, r)) {
                twofound++;
            }
        }
        Collections.sort(two);
        for (ChipSeqAnalysisResult r : one) {
            if (containsMatch(two, r)) {
                onefound++;
            }
        }
        long genomeSize = getGenome().getGenomeSize();
        double bins = genomeSize / getMaxDistance();
        double pone = one.size() / bins;
        double ptwo  = two.size() / bins;
        //        System.err.println("Pone is " + pone + " and ptwo is " + ptwo);
        double expone = one.size() * ptwo;
        double exptwo = two.size() * pone;
        System.out.println(String.format("One: %d out of %d.  Expected %.2e.  Enrichment %.1f\nTwo: %d out of %d.  Expected %.2e.  Enrichment %.1f",
                                         onefound, one.size(), expone, onefound / expone,
                                         twofound, two.size(), exptwo, twofound / exptwo));
    }
    public List<ChipSeqAnalysisResult> getOutputEvents() { return null;}


}
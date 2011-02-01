package edu.mit.csail.cgs.tools.chipseq;

import java.sql.SQLException;
import java.util.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.*;

/**
 * Sorts analysis two based on some criteria (pvalue, fold change, etc).
 * Runs down list of events from two until less than some percent
 * is contained in analysis one.
 *
 * [--overlap .9]  go down list two until overlap with list one drops to this percent
 * [--sort {pvalue|foldchange}]  sort analysis two on this field
 * [--fcone 1.0] minimum fold change for events from analysis one
 * [--fctwo 1.0] minimum fold change for events from analysis two
 * [--pvalone .001] max pvalue for events from analysis one
 * [--pvaltwo .001] max pvalue for events from analysis two
 */

public class SortedAnalysisComparison extends CompareTwoAnalyses {

    private String sortType;
    private double minFoldChangeOne, minFoldChangeTwo;
    private double maxPvalueOne, maxPvalueTwo;
    private double overlapPercent = .9;

    private int firstCheckSize = 100;

    public SortedAnalysisComparison() {
        super();
    }
    public void parseArgs(String args[]) throws NotFoundException, SQLException {
        super.parseArgs(args);
        sortType = Args.parseString(args,"sort","pvalue");
        minFoldChangeOne = Args.parseDouble(args,"fcone",1.0);
        minFoldChangeTwo = Args.parseDouble(args,"fctwo",1.0);
        maxPvalueOne = Args.parseDouble(args,"pvalone",.001);
        maxPvalueTwo = Args.parseDouble(args,"pvaltwo",.001);
        overlapPercent = Args.parseDouble(args,"overlap",.9);
        if (!(sortType.equals("pvalue") || sortType.equals("foldchange"))) {
            throw new RuntimeException("Invalid sort type "+ sortType);
        }

    }
    public List<ChipSeqAnalysisResult> getResultsOne(Region region) throws SQLException {
        List<ChipSeqAnalysisResult> output = new ArrayList<ChipSeqAnalysisResult>();
        for (ChipSeqAnalysisResult r : super.getResultsOne(region)) {
            if (r.getPValue() <= maxPvalueOne && r.getFoldEnrichment() >= minFoldChangeOne) {
                output.add(r);
            }
        }
        return output;
    }
    public List<ChipSeqAnalysisResult> getResultsTwo(Region region) throws SQLException {
        List<ChipSeqAnalysisResult> output = new ArrayList<ChipSeqAnalysisResult>();
        for (ChipSeqAnalysisResult r : super.getResultsTwo(region)) {
            if (r.getPValue() <= maxPvalueTwo && r.getFoldEnrichment() >= minFoldChangeTwo) {
                output.add(r);
            }
        }
        return output;
    }
    public List<ChipSeqAnalysisResult> getOutputEvents() throws SQLException {
        List<ChipSeqAnalysisResult> listOne = getResultsOne();
        Collections.sort(listOne);

        List<ChipSeqAnalysisResult> listTwo = getResultsTwo();
        System.err.println("There are " + listOne.size() + " events in analysis one and " + listTwo.size() + " events in analysis two");
        if (sortType.equals("pvalue")) {
            Collections.sort(listTwo,new ChipSeqAnalysisResultPvalueComparator());
        } else if (sortType.equals("foldchange")) {
            Collections.sort(listTwo,new ChipSeqAnalysisResultEnrichmentComparator());
        } else {
            throw new RuntimeException("Invalid sort type " + sortType);
        }

        double overlap = 0;
        int i = 0;
        while (i < listTwo.size() &&
               (i < firstCheckSize || overlap/i > overlapPercent)) {
            if (containsMatch(listOne,listTwo.get(i))) {
                overlap++;
            }
            i++;
        }
        i--;
        System.err.println("Keeping through " + i);
        if (overlap/i > overlapPercent) {
            return listTwo.subList(0,i+1);
        } else {
            return new ArrayList<ChipSeqAnalysisResult>();
        }
    }
    
    public static void main(String args[]) throws Exception {
        SortedAnalysisComparison sac = new SortedAnalysisComparison();
        sac.parseArgs(args);
        sac.printOutputEvents();
    }
    

}
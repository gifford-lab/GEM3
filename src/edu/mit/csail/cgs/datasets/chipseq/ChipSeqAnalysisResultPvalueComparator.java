package edu.mit.csail.cgs.datasets.chipseq;

import java.util.Comparator;

public class ChipSeqAnalysisResultPvalueComparator implements Comparator<ChipSeqAnalysisResult> {

    public boolean equals(Object o) {
        return o instanceof ChipSeqAnalysisResultPvalueComparator;
    }
    public int compare(ChipSeqAnalysisResult a, ChipSeqAnalysisResult b) {
        return Double.compare(a.getPValue(),b.getPValue());
    }

}
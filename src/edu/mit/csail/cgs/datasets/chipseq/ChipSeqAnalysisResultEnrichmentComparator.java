package edu.mit.csail.cgs.datasets.chipseq;

import java.util.Comparator;

public class ChipSeqAnalysisResultEnrichmentComparator implements Comparator<ChipSeqAnalysisResult> {

    public boolean equals(Object o) {
        return o instanceof ChipSeqAnalysisResultEnrichmentComparator;
    }
    public int compare(ChipSeqAnalysisResult a, ChipSeqAnalysisResult b) {
        return Double.compare(b.getFoldEnrichment(), a.getFoldEnrichment());
    }

}
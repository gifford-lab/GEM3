package edu.mit.csail.cgs.datasets.chipseq;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;

public class ChipSeqAnalysisResult extends Region {
    public Integer position;
    public Double foregroundReadCount, backgroundReadCount, strength, shape, pvalue, foldEnrichment;
    public ChipSeqAnalysisResult(Genome g,
                                 String chrom,
                                 int start,
                                 int end, 
                                 Integer position,
                                 Double fgcount,
                                 Double bgcount,
                                 Double strength,
                                 Double shape,
                                 Double pvalue,
                                 Double foldEnrichment) {
        super(g,chrom,start,end);
        this.position = position;
        this.foregroundReadCount = fgcount;
        this.backgroundReadCount = bgcount;
        this.strength = strength;
        this.shape = shape;
        this.pvalue = pvalue;
        this.foldEnrichment = foldEnrichment;
    }
    
}
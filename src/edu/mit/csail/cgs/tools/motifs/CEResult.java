package edu.mit.csail.cgs.tools.motifs;

import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;

public class CEResult {

    public double pval, logfoldchange, freqone, freqtwo;
    public WeightMatrix matrix;
    public int sizeone, sizetwo, countone, counttwo;
    public String percentString, cutoffString;
    
    public String toString() {
        return String.format("%2.1f\t%d\t%d\t%.2f\t%d\t%d\t%.2f\t%.2e\t%s\t%s\t%s\t%s",
                             logfoldchange,
                             countone,sizeone,freqone,
                             counttwo,sizetwo,freqtwo,
                             pval,
                             matrix.name,matrix.version,
                             percentString, cutoffString);
    }
}
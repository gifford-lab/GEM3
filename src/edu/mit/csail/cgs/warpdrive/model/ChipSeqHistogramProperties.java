package edu.mit.csail.cgs.warpdrive.model;

public class ChipSeqHistogramProperties extends ModelProperties {

    public Integer BinWidth = 10;
    public Integer DeDuplicate = 2;
    public Boolean UseWeights = Boolean.TRUE;
    public Integer GaussianKernelWidth = 0;
    public Boolean ReadExtension = Boolean.TRUE;
    public Boolean ShowPairedReads = Boolean.FALSE;
    public Boolean ShowSingleReads = Boolean.TRUE;
    public Boolean ShowSelfLigationOverlap = Boolean.FALSE;
    public Integer SelfLigationCutoff = 10000;
    public Integer SmoothingWindowWidth = 0;
    public Boolean RightFlipped = Boolean.TRUE;
    public String TSS = "11:96164825";
    public String ReadDistribution = "";
    public String EventDistribution = "";
    public Boolean ShowInteractionProfile = Boolean.FALSE;
    
    private int totalReadCount = 0;

	public int getTotalReadCount() {
		return totalReadCount;
	}

	public void setTotalReadCount(int totalReadCount) {
		this.totalReadCount = totalReadCount;
	}
    
}
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
    public String Anchor = "3:34546470-34551382";
    public String ReadDistribution = "";
    public String EventDistribution = "";
    public String Kernel = "";
    public Integer TotalReads = 0;
    public Integer ChimericReads = 0;
    public Float PValueCutoff = 0.01f;
    public Boolean ShowInteractionProfile = Boolean.FALSE;
    public Boolean ShowInteractionHistogram = Boolean.FALSE;
    public Boolean ShowInteractionKernel = Boolean.FALSE;
    
    private int totalReadCount = 0;

	public int getTotalReadCount() {
		return totalReadCount;
	}

	public void setTotalReadCount(int totalReadCount) {
		this.totalReadCount = totalReadCount;
	}
    
}
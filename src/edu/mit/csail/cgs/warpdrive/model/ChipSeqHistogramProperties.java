package edu.mit.csail.cgs.warpdrive.model;

public class ChipSeqHistogramProperties extends ModelProperties {

    public Integer BinWidth = 10;
    public Boolean UseWeights = Boolean.TRUE;
    public Integer GaussianKernelWidth = 0;
    public Boolean ReadExtension = Boolean.TRUE;
    
    private int totalReadCount = 0;

	public int getTotalReadCount() {
		return totalReadCount;
	}

	public void setTotalReadCount(int totalReadCount) {
		this.totalReadCount = totalReadCount;
	}
    
}
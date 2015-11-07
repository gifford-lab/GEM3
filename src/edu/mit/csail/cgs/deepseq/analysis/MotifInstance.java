package edu.mit.csail.cgs.deepseq.analysis;

public class MotifInstance{
	public int motifID;
	public String motifName;
	public double score;
	/** Start position of the motif in the sequence */
	public int position;
	public char strand='*';
	public String matchSeq;
	public int seqID;
}

package edu.mit.csail.cgs.datasets.alignments.parsing;

public class BLAST8Hit {

	public BLAST8Hit() {
		
	}
	
	public String queryID;
	public String subjectID;
	public double percentIdentity;
	public int alignLength;
	public int mismatches;
	public int gapOpenings;
	public int queryStart;
	public int queryEnd;
	public int subjectStart;
	public int subjectEnd;
	public double eValue;
	public double bitScore;
	
}

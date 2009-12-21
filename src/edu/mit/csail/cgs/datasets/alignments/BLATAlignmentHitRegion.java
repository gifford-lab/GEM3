package edu.mit.csail.cgs.datasets.alignments;

import edu.mit.csail.cgs.datasets.general.NamedStrandedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;


public abstract class BLATAlignmentHitRegion extends NamedStrandedRegion {
	
	public double percentIdentity;
	public int alignmentLength;
	public int numMismatches;
	public int numGapOpenings;
	public int queryStart;
	public int queryEnd;
    
	
	/**
	 * Construct a new alignment hit region from an existing one
	 * @param r
	 */
    public BLATAlignmentHitRegion(BLATAlignmentHitRegion r) { 
    	super(r);
    	this.percentIdentity = r.getPercentIdentity();
    	this.alignmentLength = r.getAlignmentLength();
    	this.numMismatches = r.getNumMismatches();
    	this.numGapOpenings = r.getNumGapOpenings();
    	this.queryStart = r.getQueryStart();
    	this.queryEnd = r.getQueryEnd();
    }
    
    
    /**
     * Construct a new alignment hit region from an existing region and the
     * additional information that is described by this object
     * @param r : a Region that represents the Target's genome, chrom, start, and end
     * @param name : the Query name
     * @param strand : the strand of the alignment (in the Target)
     * @param percentIdentity : percent identity between sequences, as a fraction (0 <= percentIdentity <= 1).
     * @param alignmentLength : length of the alignment (in the Target) in bp
     * @param numMismatches 
     * @param numGapOpenings
     * @param queryStart
     * @param queryEnd
     */
    public BLATAlignmentHitRegion(Region r, String name, char strand, 
    		double percentIdentity, int alignmentLength, int numMismatches, 
    		int numGapOpenings, int queryStart, int queryEnd) {
    	super(r, name, strand);
    	this.percentIdentity = percentIdentity;
    	this.alignmentLength = alignmentLength;
    	this.numMismatches = numMismatches;
    	this.numGapOpenings = numGapOpenings;
    	this.queryStart = queryStart;
    	this.queryEnd = queryEnd;    	    
    }
    
    
    /**
     * Construct a new alignment hit region from scratch
     * @param genome
     * @param chrom
     * @param subjectStart
     * @param subjectEnd
     * @param name
     * @param strand
     * @param percentIdentity
     * @param alignmentLength
     * @param numMismatches
     * @param numGapOpenings
     * @param queryStart
     * @param queryEnd
     */
    public BLATAlignmentHitRegion(Genome genome, String chrom, 
    		int alignmentStart, int alignmentEnd, String name, char strand, 
    		double percentIdentity, int alignmentLength, int numMismatches, 
    		int numGapOpenings, int queryStart, int queryEnd) {
    	super(genome, chrom, alignmentStart, alignmentEnd, name, strand);
    	this.percentIdentity = percentIdentity;
    	this.alignmentLength = alignmentLength;
    	this.numMismatches = numMismatches;
    	this.numGapOpenings = numGapOpenings;
    	this.queryStart = queryStart;
    	this.queryEnd = queryEnd;
    }
    
    
    /**
     * 
     */
    public String toString() {
        return super.toString() + ",% identity: " + percentIdentity 
        + ", alignment length: " + alignmentLength + ", mismatches: " + numMismatches 
        + ", gap openings: " + numGapOpenings + ", query start: " + queryStart 
        + ", query end: " + queryEnd;
    }


	/**
	 * @return the percentIdentity (between 0 and 1, inclusive)
	 */
	public double getPercentIdentity() {
		return percentIdentity;
	}


	/**
	 * @return the alignmentLength
	 */
	public int getAlignmentLength() {
		return alignmentLength;
	}


	/**
	 * @return the numMismatches
	 */
	public int getNumMismatches() {
		return numMismatches;
	}


	/**
	 * @return the numGapOpenings
	 */
	public int getNumGapOpenings() {
		return numGapOpenings;
	}


	/**
	 * @return the queryStart
	 */
	public int getQueryStart() {
		return queryStart;
	}


	/**
	 * @return the queryEnd
	 */
	public int getQueryEnd() {
		return queryEnd;
	}


	/**
	 * @return the subjectStart
	 */
	public int getSubjectStart() {
		return super.getStart();
	}


	/**
	 * @return the subjectEnd
	 */
	public int getSubjectEnd() {
		return super.getEnd();
	}
}

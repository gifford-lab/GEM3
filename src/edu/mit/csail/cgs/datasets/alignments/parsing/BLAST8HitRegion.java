package edu.mit.csail.cgs.datasets.alignments.parsing;

import edu.mit.csail.cgs.datasets.alignments.BLATAlignmentHitRegion;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;

public class BLAST8HitRegion extends BLATAlignmentHitRegion {

	/**
	 * e-value The expected probability of the number of hits occurring by chance
	 * See: 
	 * http://www.ncbi.nlm.nih.gov/Education/BLASTinfo/glossary2.html
	 * http://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html
	 */
	private double eValue;
	
	/**
	 * Score obtained by normalizing raw alignment score with respect to the
	 * scoring system
	 * http://www.ncbi.nlm.nih.gov/Education/BLASTinfo/glossary2.html
	 */
	private double bitScore;
	

	/**
	 * Construct a BLAST8HitRegion from an existing BLAST8HitRegion
	 * @param r
	 */
    public BLAST8HitRegion(BLAST8HitRegion r) { 
    	super(r);
    	this.eValue = r.getEValue();
    	this.bitScore = r.getBitScore();
    }
    
    
    /**
     * Construct a BLAST8HitRegion from an existing BLATAlignmentHitRegion and
     * an e-value and bit score
     * @param r
     * @param eValue
     */
    public BLAST8HitRegion(BLATAlignmentHitRegion r, double eValue, double bitScore) { 
    	super(r);
    	this.eValue = eValue;
    	this.bitScore = bitScore;
    }
    
    
    /**
     * Construct a BLAST8HitRegion from scratch
     * @param genome
     * @param chrom
     * @param alignmentStart
     * @param alignmentEnd
     * @param name
     * @param strand
     * @param percentIdentity
     * @param alignmentLength
     * @param numMismatches
     * @param numGaps
     * @param queryStart
     * @param queryEnd
     * @param eValue
     * @param bitScore
     */
    public BLAST8HitRegion(Genome genome, String chrom, int alignmentStart, 
    		int alignmentEnd, String name, char strand, 
    		double percentIdentity, int alignmentLength, int numMismatches, 
    		int numGaps, int queryStart, int queryEnd, double eValue, double bitScore) {
    	super(genome, chrom, alignmentStart, alignmentEnd, name, strand, 
    			percentIdentity, alignmentLength, numMismatches, 
    			numGaps, queryStart, queryEnd);
    	this.eValue = eValue;
    	this.bitScore = bitScore;
    }
    
    
	/**
	 * @return the eValue
	 */
	public double getEValue() {
		return eValue;
	}
	
	
	/**
	 * @return the bitScore
	 */
	public double getBitScore() {
		return bitScore;
	}
	
	
	/**
	 * 
	 * @return
	 */
	public String getBLAST8String() {
		char delim = '\t';
		StringBuffer sb = new StringBuffer();
		sb.append(getName());
		sb.append(delim);
		sb.append(getChrom());
		sb.append(delim);
		sb.append(getPercentIdentity());
		sb.append(delim);
		sb.append(getAlignmentLength());
		sb.append(delim);
		sb.append(getNumMismatches());
		sb.append(delim);
		sb.append(getNumGapOpenings());
		sb.append(delim);
		sb.append(getQueryStart());
		sb.append(delim);
		sb.append(getQueryEnd());
		sb.append(delim);
		if (getStrand() == '+') {
			sb.append(getStart());
			sb.append(delim);
			sb.append(getEnd());
			sb.append(delim);
		}
		else {
			sb.append(getEnd());
			sb.append(delim);
			sb.append(getStart());
			sb.append(delim);			
		}
		sb.append(getEValue());
		sb.append(delim);
		sb.append(getBitScore());
		
		return sb.toString();
	}

    public String toString() {
        return getBLAST8String();
    }
}

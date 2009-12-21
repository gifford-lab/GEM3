package edu.mit.csail.cgs.datasets.alignments;

import java.util.List;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.strings.StringUtils;

public class PSLHitRegion extends BLATAlignmentHitRegion {

	int matches; //number of bases that match that aren't repeats
	int repMatches; //number of bases that match but are part of repeats
	int nCount; //number of 'N' bases
	int qNumInsert; //Number of inserts in query
	int qBaseInsert; //Number of bases inserted in query
	int tNumInsert; //Number of inserts in target
	int tBaseInsert; //Number of bases inserted in target
	int qSize; //Query Sequence size
    int tSize; //Target sequence size
    int blockCount; //Number of blocks in alignment. A block contains no gaps.
    List<Integer> blockSizes; //Size of each block
    List<Integer> qStarts; //Start of each block in query
    List<Integer> tStarts; //Start of each block in target
	
    /**
	 * Construct a PSLHitRegion from an existing PSLHitRegion
	 * @param r
	 */
    public PSLHitRegion(PSLHitRegion r) { 
    	super(r);
    	this.matches = r.getMatches();
    	this.repMatches = r.getRepMatches();
    	this.nCount = r.getNCount();
    	this.qNumInsert = r.getQNumInsert();
    	this.qBaseInsert = r.getQBaseInsert();
    	this.tNumInsert = r.getTNumInsert();
    	this.tBaseInsert = r.getTBaseInsert();
    	this.qSize = r.getQSize();
    	this.tSize = r.getTSize();
    	this.blockCount = r.getBlockCount();
    	this.blockSizes = r.getBlockSizes();
    	this.qStarts = r.getQStarts();
    	this.tStarts = r.getTStarts();
    }

    
    /**
     * Construct a PSLHitRegion from an existing BLATAlignmentHitRegion and
     * additional parameters
     * @param r
     * @param matches
     * @param repMatches
     * @param nCount
     * @param qNumInsert
     * @param qBaseInsert
     * @param tNumInsert
     * @param tBaseInsert
     * @param qSize
     * @param tSize
     * @param blockCount
     * @param blockSizes
     * @param qStarts
     * @param tStarts
     */
    public PSLHitRegion(BLATAlignmentHitRegion r, int matches, int repMatches,
    		int nCount, int qNumInsert, int qBaseInsert, int tNumInsert, int tBaseInsert, 
    		int qSize, int tSize, int blockCount, List<Integer> blockSizes,
    		List<Integer> qStarts, List<Integer> tStarts) { 
    	super(r);
    	this.matches = matches;
    	this.repMatches = repMatches;
    	this.nCount = nCount;
    	this.qNumInsert = qNumInsert;
    	this.qBaseInsert = qBaseInsert;
    	this.tNumInsert = tNumInsert;
    	this.tBaseInsert = tBaseInsert;
    	this.qSize = qSize;
    	this.tSize = tSize;
    	this.blockCount = blockCount;
    	this.blockSizes = blockSizes;
    	this.qStarts = qStarts;
    	this.tStarts = tStarts;
    }


    /**
     * Construct a PSLHitRegion from scratch
     * @param genome
     * @param chrom
     * @param alignmentStart
     * @param alignmentEnd
     * @param name
     * @param strand
     * @param percentIdentity
     * @param alignmentLength
     * @param numMismatches
     * @param numGapOpenings
     * @param queryStart
     * @param queryEnd
     * @param bitScore
     * @param matches
     * @param repMatches
     * @param nCount
     * @param qNumInsert
     * @param qBaseInsert
     * @param tNumInsert
     * @param tBaseInsert
     * @param qSize
     * @param tSize
     * @param blockCount
     * @param blockSizes
     * @param qStarts
     * @param tStarts
     */
    public PSLHitRegion(Genome genome, String chrom, 
    		int alignmentStart, int alignmentEnd, String name, char strand, 
    		double percentIdentity, int alignmentLength, int numMismatches, 
    		int numGapOpenings, int queryStart, int queryEnd,    		
    		int matches, int repMatches, int nCount, int qNumInsert, 
    		int qBaseInsert, int tNumInsert, int tBaseInsert, 
    		int qSize, int tSize, int blockCount, List<Integer> blockSizes,
    		List<Integer> qStarts, List<Integer> tStarts) { 
    	super(genome, chrom, alignmentStart, alignmentEnd, name, strand, 
    			percentIdentity, alignmentLength, numMismatches, numGapOpenings, 
    			queryStart, queryEnd);
    	this.matches = matches;
    	this.repMatches = repMatches;
    	this.nCount = nCount;
    	this.qNumInsert = qNumInsert;
    	this.qBaseInsert = qBaseInsert;
    	this.tNumInsert = tNumInsert;
    	this.tBaseInsert = tBaseInsert;
    	this.qSize = qSize;
    	this.tSize = tSize;
    	this.blockCount = blockCount;
    	this.blockSizes = blockSizes;
    	this.qStarts = qStarts;
    	this.tStarts = tStarts;
    }
    


	/**
	 * @return the matches
	 */
	public int getMatches() {
		return matches;
	}


	/**
	 * @return the repMatches
	 */
	public int getRepMatches() {
		return repMatches;
	}


	/**
	 * @return the nCount
	 */
	public int getNCount() {
		return nCount;
	}


	/**
	 * @return the qNumInsert
	 */
	public int getQNumInsert() {
		return qNumInsert;
	}


	/**
	 * @return the qBaseInsert
	 */
	public int getQBaseInsert() {
		return qBaseInsert;
	}


	/**
	 * @return the tNumInsert
	 */
	public int getTNumInsert() {
		return tNumInsert;
	}


	/**
	 * @return the tBaseInsert
	 */
	public int getTBaseInsert() {
		return tBaseInsert;
	}


	/**
	 * @return the qSize
	 */
	public int getQSize() {
		return qSize;
	}


	/**
	 * @return the tSize
	 */
	public int getTSize() {
		return tSize;
	}


	/**
	 * @return the blockCount
	 */
	public int getBlockCount() {
		return blockCount;
	}


	/**
	 * @return the blockSizes
	 */
	public List<Integer> getBlockSizes() {
		return blockSizes;
	}


	/**
	 * @return the qStarts
	 */
	public List<Integer> getQStarts() {
		return qStarts;
	}


	/**
	 * @return the tStarts
	 */
	public List<Integer> getTStarts() {
		return tStarts;
	}

    /*match   mis-    rep.    N's     Q gap   Q gap   T gap   T gap   strand  Q               Q       Q       Q       T               T       T       T       block   blockSizes      qStarts  tStarts
              match   match           count   bases   count   bases           name            size    start   end     name            size    start   end     count
         1     2       3     4        5        6        7      8       9      10              11       12     13      14               15     16      17        18       19            20       21     
    */

    public String toString() {
        return String.format("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s",
                             matches,
                             getWidth() - matches,
                             repMatches,
                             nCount,
                             qNumInsert,
                             qBaseInsert,
                             tNumInsert,
                             tBaseInsert,
                             getStrand(),
                             getName(),
                             qSize,
                             queryStart,
                             queryEnd,
                             getChrom(),
                             getGenome().getChromLength(getChrom()),
                             getStart(),
                             getEnd(),
                             blockSizes.size(),
                             StringUtils.join(blockSizes,","),
                             StringUtils.join(qStarts,","),
                             StringUtils.join(tStarts,","));
                             
    }
}

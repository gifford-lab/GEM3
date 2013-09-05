package edu.mit.csail.cgs.datasets.chipseq;

import edu.mit.csail.cgs.datasets.species.Genome;

public class ChiaPetHit extends ChipSeqHit {

	public ChiaPetHit(Genome g, String chrom, int start, int end, char strand,
			ChipSeqAlignment align, double weight) {
		super(g, chrom, start, end, strand, align, weight);
	}
	
	public ChiaPetHit extendHit(int ext) { 
	    if(getStrand() == '+') { 
	      return new ChiaPetHit(getGenome(), getChrom(), getStart() - ext, getEnd(), getStrand(), getAlignment(), getWeight());
	    } else { 
	      return new ChiaPetHit(getGenome(), getChrom(), getStart(), getEnd() + ext, getStrand(), getAlignment(), getWeight());
	    }
	  }

}

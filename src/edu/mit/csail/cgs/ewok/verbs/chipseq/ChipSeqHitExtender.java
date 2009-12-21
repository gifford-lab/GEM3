package edu.mit.csail.cgs.ewok.verbs.chipseq;

import edu.mit.csail.cgs.datasets.chipseq.ChipSeqHit;
import edu.mit.csail.cgs.ewok.verbs.Mapper;

public class ChipSeqHitExtender implements Mapper<ChipSeqHit,ChipSeqHit> {
	
	private int extension;
	
	public ChipSeqHitExtender(int ext) { 
		extension = ext;
	}

	public ChipSeqHit execute(ChipSeqHit a) {
        return a.extendHit(extension);
	}
}

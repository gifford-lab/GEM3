/*
 * Author: tdanford
 * Date: Dec 4, 2008
 */
package edu.mit.csail.cgs.utils.models.bns;

import edu.mit.csail.cgs.utils.graphs.DirectedGraph;

public interface BNModelScore {
	public Double graphScore(BN network);
}

class MDLGraphScore implements BNModelScore {
	
	private Double scale;
	
	public MDLGraphScore() { 
		this(1.0);
	}
	
	public MDLGraphScore(Double s) { 
		scale = s;
	}

	public Double graphScore(BN network) {
		Double mdl = Math.log((double)network.getData().size()) / 2.0;
		mdl *= (double)network.countParameters();
		return scale * mdl;
	} 
	
}

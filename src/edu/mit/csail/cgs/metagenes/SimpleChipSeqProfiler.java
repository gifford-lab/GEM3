/*
 * Author: tdanford
 * Date: Aug 19, 2008
 */
package edu.mit.csail.cgs.metagenes;

import java.util.*;

import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.ewok.verbs.chipseq.*;

public class SimpleChipSeqProfiler implements PointProfiler<Point,PointProfile> {
	
	private BinningParameters params;
	private List<ChipSeqExpander> expanders;
	private int extension; 
	
	public SimpleChipSeqProfiler(BinningParameters ps, ChipSeqExpander exp) { 
		this(ps, exp, 175);
	}
	public SimpleChipSeqProfiler(BinningParameters ps, ChipSeqExpander exp, int ext) {
		params = ps;
		expanders = new ArrayList<ChipSeqExpander>(); 
		expanders.add(exp);
		extension=ext;
	}
	public SimpleChipSeqProfiler(BinningParameters ps, List<ChipSeqExpander> exps, int ext) {
		params = ps;
		expanders = exps;
		extension=ext;
	}

	public BinningParameters getBinningParameters() {
		return params;
	}

	public PointProfile execute(Point a) {
		int window = params.getWindowSize();
		int left = window/2;
		int right = window-left-1;
		
		boolean strand = (a instanceof StrandedPoint) ? 
				((StrandedPoint)a).getStrand() == '+' : true;
		
		
		int start = Math.max(0, a.getLocation()-left);
		int end = Math.min(a.getLocation()+right, a.getGenome().getChromLength(a.getChrom())-1);
		
		Region query = new Region(a.getGenome(), a.getChrom(), start, end);
		
		double[] array = new double[params.getNumBins()];
		for(int i = 0; i < array.length; i++) { array[i] = 0; }
		
		for(ChipSeqExpander expander : expanders){
			Iterator<ChipSeqHit> hits = expander.execute(query);
			while(hits.hasNext()) {
				ChipSeqHit hit = hits.next().extendHit(extension);
				int startOffset = Math.max(0, hit.getStart()-query.getStart());
				int endOffset = Math.max(0, Math.min(query.getEnd(), hit.getEnd()-query.getStart()));
				
				if(!strand) { 
					int tmpEnd = window-startOffset;
					int tmpStart = window-endOffset;
					startOffset = tmpStart;
					endOffset = tmpEnd;
				}
				
				int startbin = params.findBin(startOffset);
				int endbin = params.findBin(endOffset);
				
				addToArray(startbin, endbin, array, 1.0);				
			}
		}		
		return new PointProfile(a, params, array, (a instanceof StrandedPoint));
	}

	private void addToArray(int i, int j, double[] array, double value) { 
		for(int k = i; k <= j; k++) { 
			array[k] += value;
		}
	}
	
	public void cleanup(){
		for(ChipSeqExpander e : expanders)
			e.close();
	}
}

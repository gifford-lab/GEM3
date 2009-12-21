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

public class ChipSeq5PrimeProfiler implements PointProfiler<Point,PointProfile> {
	
	private BinningParameters params;
	private List<ChipSeqExpander> expanders;
	private char strand;
	
	public ChipSeq5PrimeProfiler(BinningParameters ps, ChipSeqExpander exp, char strand) {
		params = ps;
		expanders = new ArrayList<ChipSeqExpander>(); 
		expanders.add(exp);
		this.strand = strand;
	}
	public ChipSeq5PrimeProfiler(BinningParameters ps, List<ChipSeqExpander> exps, char strand) {
		params = ps;
		expanders = exps;
		this.strand = strand;
	}

	public BinningParameters getBinningParameters() {
		return params;
	}

	public PointProfile execute(Point a) {
		int window = params.getWindowSize();
		int left = window/2;
		int right = window-left-1;
		
		int start = Math.max(0, a.getLocation()-left);
		int end = Math.min(a.getLocation()+right, a.getGenome().getChromLength(a.getChrom())-1);
		
		Region query = new Region(a.getGenome(), a.getChrom(), start, end);
		
		double[] array = new double[params.getNumBins()];
		for(int i = 0; i < array.length; i++) { array[i] = 0; }
		
		for(ChipSeqExpander expander : expanders){
			Iterator<ChipSeqHit> hits = expander.execute(query);
			List <ChipSeqHit> hitList =  filterDuplicateHits(hits);

			for(ChipSeqHit hit : hitList) {
				if (hit.getStrand()==this.strand){  //only count one strand
					if ((start<=hit.getFivePrime() && this.strand=='+')
							||(end>hit.getFivePrime() && this.strand=='-')){
						int hit5Prime = hit.getFivePrime()-start;
						array[params.findBin(hit5Prime)]++;
					}
				}				
			}
		}		
		return new PointProfile(a, params, array, (a instanceof StrandedPoint));
	}
	
	/*
	 * filter out duplicate reads (potential tower, needles, but could be real reads)
	 * assuming the reads are sorted
	 */
	public List<ChipSeqHit> filterDuplicateHits(Iterator<ChipSeqHit> hits){
		ChipSeqHit currentHit = hits.next();
		int count=1;
		List<ChipSeqHit> filteredReads = new ArrayList<ChipSeqHit>();
		while(hits.hasNext()) {
			ChipSeqHit hit = hits.next();
			// if read from a new position
			if (!(currentHit.getStart()==hit.getStart())){
				currentHit = hit;
				count=0;
				filteredReads.add(hit);
			}
			else {// if  duplicate
				count++;
				if (count<=3)
					filteredReads.add(hit);
			}
		}
		return filteredReads;
	}
	
	//No cleanup
	public void cleanup(){
		for(ChipSeqExpander e : expanders)
			e.close();
	}
}

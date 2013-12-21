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
	private int readShift=0; 
	private double perBaseMax=100;
	private boolean useFivePrime = false;
	private char readStrand ='/';
	
	public SimpleChipSeqProfiler(BinningParameters ps, ChipSeqExpander exp) { 
		this(ps, exp, 175, 100, '/');
	}
	public SimpleChipSeqProfiler(BinningParameters ps, ChipSeqExpander exp, int ext, double pbMax, char strand) {
		params = ps;
		expanders = new ArrayList<ChipSeqExpander>(); 
		expanders.add(exp);
		extension=ext;
		if(extension==-1){
			useFivePrime=true;
			extension = 0;
		}
		perBaseMax = pbMax;
		readStrand = strand;
	}
	public SimpleChipSeqProfiler(BinningParameters ps, List<ChipSeqExpander> exps, int ext, int shift, double pbMax, char strand) {
		params = ps;
		expanders = exps;
		extension=ext;
		if(extension==-1){
			useFivePrime=true;
			extension = 0;
		}
		readShift = shift;
		perBaseMax=pbMax;
		readStrand = strand;
	}

	public BinningParameters getBinningParameters() {
		return params;
	}
	public void setUseFivePrime(boolean ufp){useFivePrime = ufp;}

	public PointProfile execute(Point a) {
		int window = params.getWindowSize();
		int left = window/2;
		int right = window-left-1;
		
//		boolean isPlusStrand = (a instanceof StrandedPoint) ? 
//				((StrandedPoint)a).getStrand() == '+' : true;		// set to true if non-stranded
		
		int start = Math.max(0, a.getLocation()-left);
		int end = Math.min(a.getLocation()+right, a.getGenome().getChromLength(a.getChrom())-1);
		
		Region query = new Region(a.getGenome(), a.getChrom(), start, end);
		Region extQuery = new Region(a.getGenome(), a.getChrom(), start-extension>0 ? start-extension : 1, end+extension < a.getGenome().getChromLength(a.getChrom()) ? end+extension : a.getGenome().getChromLength(a.getChrom()) );
		
		double[] array = new double[params.getNumBins()];
		for(int i = 0; i < array.length; i++) { array[i] = 0; }
		
		for(ChipSeqExpander expander : expanders){
			Iterator<ChipSeqHit> hits = expander.execute(extQuery);
			HashMap<Region, Double> readFilter = new HashMap<Region, Double>();
			
			while(hits.hasNext()) {
				ChipSeqHit hit=null;		// hit is end exclusive
				if(useFivePrime)
					hit = hits.next().fivePrime().shiftExtendHit(0, readShift);
				else
					hit = hits.next().shiftExtendHit(extension, readShift);
				if(hit.overlaps(query) && (readStrand=='/' || hit.getStrand()==readStrand)){
					if(!readFilter.containsKey(hit))
						readFilter.put(hit, hit.getWeight());
					else
						readFilter.put(hit, readFilter.get(hit)+hit.getWeight());
					
					if(readFilter.get(hit)<=perBaseMax){			// skip higher count, not truncate
						int startOffset = Math.max(0, hit.getStart()-start);
						int endOffset = Math.max(0, Math.min(end, hit.getEnd()-1-start)); // -1: hit is end exclusive
					
						if(hit.getStrand()=='-') { 
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

package edu.mit.csail.cgs.deepseq;

import java.util.ArrayList;
import java.util.List;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.deepseq.features.ComponentFeature;

/**
 * It represents a stranded base position in the genome.
 * It has strand information.
 * Coordinate is the 5' end of the read. 
 * We do not store chromomosome info here because this class is always used
 * 		in the context of a chromosome or a region, no need to differentiate.
 * It records the number of reads mapped to this base position (or weighted count). 
 * 
 * @author yguo
 *
 */
public class StrandedBase implements Comparable<StrandedBase>{
	private char strand; 
	private int coordinate;
	private float count;
	
	public StrandedBase(char strand, int coord, float count){
		this.setStrand(strand);
		this.setCoordinate(coord);
		this.setCount(count);
	}

	public void setStrand(char strand) {
		this.strand = strand;
	}

	public char getStrand() {
		return strand;
	}

	public void setCoordinate(int coordinate) {
		this.coordinate = coordinate;
	}

	public int getCoordinate() {
		return coordinate;
	}

	public void setCount(float count) {
		this.count = count;
	}

	public float getCount() {
		return count;
	}
	// sort according to coordinate, considering strand
	public int compareTo(StrandedBase b) {
		double diff = coordinate-b.coordinate;
		return diff==0?
				(strand==b.strand?0:strand=='+'?1:-1):	// same coord, compare strand
					(diff<0?-1:1);						// diff coord
	}
	
	public String toString(){
		return coordinate+" "+strand+" "+count;
	}
	
	static public float countBaseHits(List<StrandedBase> bases) {
		float count = 0;
		for (StrandedBase b: bases){
			count += b.getCount();
		}
		return count;
	}
	static public List<StrandedBase> convertReadHits (List<ReadHit> reads){
		List<StrandedBase> result = new ArrayList<StrandedBase>();
		for (ReadHit h: reads){
			int fivePrime = h.getStrand()=='+'?h.getStart():h.getEnd();
			boolean found=false;
			for (StrandedBase b: result){
				if ((h.getStrand()=='+' && b.getStrand()=='+' && b.getCoordinate()==fivePrime)||
					(h.getStrand()=='-' && b.getStrand()=='-' && b.getCoordinate()==fivePrime)){
					found = true;
					b.setCount(b.getCount()+(float)h.getWeight());
					break;
				}
			}
			if (!found){
				StrandedBase base = new StrandedBase(h.getStrand(),fivePrime,(float)h.getWeight());
				result.add(base);
			}
		}
		return result;
	}
}

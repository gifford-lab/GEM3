package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.species.Genome;

public class RandomPointGenerator implements Iterator<Point> {

	private Genome genome;
	private double[] chromWeights;
	private String[] chroms;
	private int[] lengths;
	private Random rand;
	
	public RandomPointGenerator(Genome g) { 
		genome = g;
		rand = new Random();
		List<String> chromList = genome.getChromList();
		chroms = new String[chromList.size()];
		chromWeights = new double[chroms.length];
		lengths = new int[chroms.length];
		long totalLength = (long)0;
		
		int i =0;
		for(String chrom : chromList) { 
			chroms[i] = chrom;
			lengths[i] = genome.getChromLength(chroms[i]);
			totalLength += (long)lengths[i];
			i++;
		}
		
		for(i = 0; i < chromWeights.length; i++) { 
			chromWeights[i] = (double)lengths[i] / (double)totalLength;
		}
	}
	
	private int sampleChrom() { 
		double p = rand.nextDouble();
		for(int i = 0; i < chromWeights.length; i++) { 
			p -= chromWeights[i];
			if(p <= 0.0) { return i; }
		}
		return chromWeights.length-1;
	}
	
	private int sampleOffset(int chrom) { 
		return rand.nextInt(lengths[chrom]);
	}
	
	public boolean hasNext() {
		return true;
	}
	
	public Point next() {
		int chromIdx = sampleChrom();
		int offset = sampleOffset(chromIdx);
		return new Point(genome, chroms[chromIdx], offset);
	}
	
	public void remove() {
		throw new UnsupportedOperationException();
	}
	
	
}

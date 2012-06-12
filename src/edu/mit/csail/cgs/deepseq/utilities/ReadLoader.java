package edu.mit.csail.cgs.deepseq.utilities;

import java.util.List;
import java.util.TreeMap;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.ExtReadHit;
import edu.mit.csail.cgs.deepseq.ReadHit;
import edu.mit.csail.cgs.projects.readdb.PairedHit;

/**
 * Loads reads. Where those reads are sourced from is implementation-specific. 
 * This class contains methods that are required regardless of the source of the reads. For example, loading reads according to a region. 
 * 
 * @author shaun
 *
 */
public abstract class ReadLoader {
	protected Genome gen;
	protected double totalHits;
	protected double totalWeight;
	protected double totalForWeight=-1,totalRevWeight=-1;
	protected int readLength;
	protected boolean pairedEnd = false; //What to do with this flag is left to the subclass
	
	public ReadLoader(Genome g, int rLen){
		gen=g;
		readLength=rLen;
		totalHits=0;
		totalWeight=0;
	}
	//Accessors
	public Genome getGenome(){return(gen);}
	public int getReadLen(){return readLength;}
	public double getHitCount(){return(totalHits);}
	public double getTotalWeight(){return(totalWeight);}
	public double getStrandedWeight(char strand){
		//System.out.print(strand+" strand hits: ");
		if(strand=='+'){
			if(totalForWeight<=-1)
				totalForWeight=countStrandedWeight(strand);
			totalRevWeight=totalWeight-totalForWeight;
		//	System.out.print((int)totalForWeight+"\n");
			return(totalForWeight);
		}else{
			if(totalRevWeight<=-1)
				totalRevWeight=countStrandedWeight(strand);
			totalForWeight=totalWeight-totalRevWeight;
		//	System.out.print((int)totalRevWeight+"\n");
			return(totalRevWeight);
		}
	}
	
	//Abstract methods
	protected abstract double countHits();
	protected abstract double countStrandedWeight(char strand);
	public abstract List<ReadHit> loadHits(Region r);
	public abstract List<ExtReadHit> loadExtHits(Region r, int startShift, int fivePrimeExt, int threePrimeExt);
	public abstract List<ReadHit> loadPairs(Region r);
	public abstract List<PairedHit> loadPairsAsPairs(Region r);
	public abstract int countHits (Region r);
	public abstract double sumWeights (Region r);
	public abstract void cleanup();
	public abstract void setGenome(Genome g);
	public void setPairedEnd(boolean pe){pairedEnd=pe;} 
}

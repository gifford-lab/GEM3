package edu.mit.csail.cgs.deepseq;

import java.util.ArrayList;

import edu.mit.csail.cgs.datasets.general.Region;
//A container to hold diverse types of BackgroundModel and handle the various ways we may use them
public class BackgroundCollection {
	
	protected ArrayList<BackgroundModel> models = new ArrayList<BackgroundModel>();
	protected double highestHigh=0, lowestLow=0;
	
	public BackgroundCollection(){}
	
	//Accessor
	public int getGenomicModelThreshold(String tType){
		for(BackgroundModel m : models)
			if(m.isGenomeWide())
				return m.getThreshold();
		return -1;
	}
	public int getMaxThreshold(char str){
		int max=0;
		for(BackgroundModel m : models)
			if((m.strand==str ||str=='.' ) && m.getThreshold()>max)
				max = m.getThreshold();
		return max;
	}
	
	public void addBackgroundModel(BackgroundModel m){models.add(m);}
	
	//Update
	public void updateModels(Region currReg, int currOffset, double [] thisExptHitCounts, double [] otherExptHitCounts){
		for(BackgroundModel m : models)
			m.updateModel(currReg, currOffset, thisExptHitCounts, otherExptHitCounts);		
	}
	
	//Print
	public void printThresholds(){
		for(BackgroundModel m : models){
			System.out.println(m.modelType+"\t"+m.getThreshold()+"\t"+m.getStrand());
		}
	}
	
	//Passes?
	public boolean passesAllThresholds(int count){return passesAllThresholds(count, '.');}
	public boolean passesAllThresholds(int count, char strand){
		boolean pass=true;
		for(BackgroundModel m : models){
			if(m.getStrand()==strand)
				if(!m.passesThreshold(count))
					pass=false;
		}
		return pass;
	}
	public boolean passesGenomicThreshold(int count){return passesGenomicThreshold(count, '.');}
	public boolean passesGenomicThreshold(int count, char strand){
		boolean pass=true;
		for(BackgroundModel m : models){
			if(m.isGenomeWide() && m.getStrand()==strand)
				if(!m.passesThreshold(count))
					pass=false;
		}
		return pass;
	}
	
	//Under?
	public boolean underAllThresholds(int count){return underAllThresholds(count, '.');}
	public boolean underAllThresholds(int count, char strand){
		boolean pass=true;
		for(BackgroundModel m : models){
			if(m.getStrand()==strand)
				if(!m.underThreshold(count))
					pass=false;
		}
		return pass;
	}
	public boolean underGenomicThreshold(int count){return underGenomicThreshold(count, '.');}
	public boolean underGenomicThreshold(int count, char strand){
		boolean pass=true;
		for(BackgroundModel m : models){
			if(m.isGenomeWide() && m.getStrand()==strand)
				if(!m.underThreshold(count))
					pass=false;
		}
		return pass;
	}
}

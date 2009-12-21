package edu.mit.csail.cgs.ewok.nouns;

import java.util.Vector;

import edu.mit.csail.cgs.datasets.species.Gene;

public class GeneDomainTimeSeries {
	
	public static final int window = 30000;

	private Gene gene;
	private Vector<GeneDomainData> data;
	
	public GeneDomainTimeSeries(Gene g, int tps) { 
		gene = g;
		data = new Vector<GeneDomainData>();
		for(int i = 0; i < tps; i++) { 
			data.add(new GeneDomainData(gene, window));
		}
	}
	
	public void addDomain(SimpleDomain sd, int tp) { 
		data.get(tp).addDomain(sd);
	}
	
	public boolean isCovered(int tp) { 
		return data.get(tp).isTSSCovered();
	}
	
	public String getBitString() {
		String str = "";
		for(int i = 0; i < data.size(); i++) { 
			str += isCovered(i) ? "1" : "0";
		}
		return str;
	}
	
	public String toString() { return gene.getID() + " " + getBitString(); }
}

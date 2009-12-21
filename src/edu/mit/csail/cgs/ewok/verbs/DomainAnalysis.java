package edu.mit.csail.cgs.ewok.verbs;

import java.util.Iterator;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.ewok.nouns.GeneDomainData;

public class DomainAnalysis<DOMAIN extends Region> 
	implements Mapper<Gene,GeneDomainData> {

	private int window;
	private Expander<Region,DOMAIN> domainCaller;
	
	public DomainAnalysis(Expander<Region,DOMAIN> dc, int w) { 
		domainCaller = dc;
		window = w;
	}

	public GeneDomainData execute(Gene g) {
		GeneDomainData gdd = new GeneDomainData(g, window);
		Iterator<DOMAIN> itr = domainCaller.execute(gdd.getWindow());
		while(itr.hasNext()) { 
			gdd.addDomain(itr.next());
		}
		return gdd;
	}
	
	
}

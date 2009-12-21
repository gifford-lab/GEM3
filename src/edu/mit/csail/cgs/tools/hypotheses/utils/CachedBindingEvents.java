package edu.mit.csail.cgs.tools.hypotheses.utils;

import java.util.*;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.verbs.Expander;

public class CachedBindingEvents {
	
	private Genome genome;
	private Map<String,CachedChromosomeBindingEvents> chroms;

	public CachedBindingEvents(Genome g) { 
		genome = g;
		chroms = new HashMap<String,CachedChromosomeBindingEvents>();
		for(String chrom : genome.getChromList()) { 
			chroms.put(chrom, new CachedChromosomeBindingEvents(genome, chrom));
		}
	}
	
	public Genome getGenome() { return genome; }
	
	public void setEvents(Expander<Region,BindingEvent> exp) { 
		for(String chrom : chroms.keySet()) { 
			chroms.get(chrom).setEvents(exp);
		}
	}
	
	public int getNumEvents() { 
		int c = 0;
		for(String chrom : chroms.keySet()) { 
			c += chroms.get(chrom).getNumEvents(); 
		}
		return c;
	}
	
	public Collection<BindingEvent> findEvents(Region r) { 
		if(!chroms.containsKey(r.getChrom())) { 
			throw new IllegalArgumentException(r.getChrom());
		}
		
		return chroms.get(r.getChrom()).findEvents(r);
	}
}

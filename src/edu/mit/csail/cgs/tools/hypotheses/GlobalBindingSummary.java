package edu.mit.csail.cgs.tools.hypotheses;

import java.util.*;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.general.Factor;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.tools.hypotheses.utils.CachedBindingEvents;

public class GlobalBindingSummary {
	
	private Genome genome;
	private Map<Factor,CachedBindingEvents> events;

	public GlobalBindingSummary(Genome g) { 
		genome = g;
		events = new HashMap<Factor,CachedBindingEvents>();
	}
	
	public void addFactor(Factor f, Expander<Region,BindingEvent> exp) { 
		events.put(f, new CachedBindingEvents(genome));
		events.get(f).setEvents(exp);
	}
	
	public Set<Factor> getFactors() { 
		return new TreeSet<Factor>(events.keySet()); 
	}
	
	public Genome getGenome() { return genome; }
	
	public int getNumEvents(Factor f) { 
		return events.get(f).getNumEvents();
	}
	
	public Collection<BindingEvent> findEvents(Factor f, Region r) { 
		return events.get(f).findEvents(r);
	}
	
	public LocalBindingSummary createSummary(Region r) { 
		HashMap<Factor,Integer> counts = new HashMap<Factor,Integer>();
		for(Factor f : events.keySet()) { 
			Collection<BindingEvent> evts = events.get(f).findEvents(r);
			counts.put(f, evts.size());
		}
		return new LocalBindingSummary(r, counts);
	}
	
	public void printStatistics() { 
		for(Factor f : new TreeSet<Factor>(events.keySet())) {  
			System.out.println(f.getName());
			System.out.println("\t# Events: " + events.get(f).getNumEvents());
		}
	}
}

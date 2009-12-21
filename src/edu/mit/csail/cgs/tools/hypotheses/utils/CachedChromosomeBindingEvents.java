package edu.mit.csail.cgs.tools.hypotheses.utils;

import java.util.*;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.verbs.Expander;

public class CachedChromosomeBindingEvents {

	private Genome genome;
	private NamedRegion chrom;
	private BindingEvent[] eventArray;
	
	public CachedChromosomeBindingEvents(Genome g, String chromName) { 
		genome = g;
		Genome.ChromosomeInfo info = genome.getChrom(chromName);
		chrom = new NamedRegion(genome, info.getName(), 1, info.getLength(), info.getName());
		eventArray = null;
	}
	
	public void setEvents(Expander<Region,BindingEvent> exp) { 
		Vector<BindingEvent> events = new Vector<BindingEvent>();
		Iterator<BindingEvent> evts = exp.execute(chrom);
		while(evts.hasNext()) { 
			events.add(evts.next());
		}

		eventArray = events.toArray(new BindingEvent[events.size()]);
		Arrays.sort(eventArray);
	}
	
	public int getNumEvents() { 
		return eventArray != null ? eventArray.length : 0; 
	}
	
	public BindingEvent getEvent(int i) { 
		return eventArray != null ? eventArray[i] : null;
	}
	
	public Collection<BindingEvent> findEvents(Region r) { 
		if(eventArray == null) { throw new IllegalArgumentException("Uninitialized!"); }
		if(!chrom.overlaps(r)) { throw new IllegalArgumentException(r.toString()); }
		LinkedList<BindingEvent> found = new LinkedList<BindingEvent>();
		
		for(int i = 0; i < eventArray.length; i++) { 
			if(eventArray[i].getEnd() >= r.getStart()) { 
				for(int j = i; j < eventArray.length && eventArray[j].getStart() <= r.getEnd(); j++) { 
					if(eventArray[j].overlaps(r)) { 
						found.addLast(eventArray[j]);
					}
				}
				
				break;
			}
		}
		
		return found;
	}
}

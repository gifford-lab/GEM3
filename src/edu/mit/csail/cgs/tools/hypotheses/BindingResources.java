package edu.mit.csail.cgs.tools.hypotheses;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.general.Factor;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.utils.Closeable;

public class BindingResources {
	
	private Genome genome;
	private Map<Factor,Expander<Region,BindingEvent>> detectors;
	
	public BindingResources(Genome g) {
		genome = g;
		detectors = new HashMap<Factor,Expander<Region,BindingEvent>>();
	}
	
	public Genome getGenome() { return genome; }
	public Set<Factor> getFactors() { return new TreeSet<Factor>(detectors.keySet()); }

	public void addExpander(Factor f, Expander<Region,BindingEvent> exp) { 
		if(detectors.containsKey(f)) { 
			if(detectors.get(f) instanceof Closeable) { 
				Closeable c = (Closeable)detectors.get(f);
				c.close();
			}
		}
		detectors.put(f, exp);
	}
	
	public GlobalBindingSummary getSummary() { 
		GlobalBindingSummary gbs = new GlobalBindingSummary(genome);
		for(Factor f : detectors.keySet()) { 
			gbs.addFactor(f, detectors.get(f));
		}
		return gbs;
	}
}

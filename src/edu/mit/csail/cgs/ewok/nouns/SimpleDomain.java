package edu.mit.csail.cgs.ewok.nouns;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;

public class SimpleDomain extends Region {

	public SimpleDomain(Genome g, String c, int start, int end) {
		super(g, c, start, end);
	}	
	
	public SimpleDomain(Region r) { 
		super(r);
	}
	
	public BindingEvent getBindingEvent() { 
		return new BindingEvent(getGenome(), getChrom(), getStart(), getEnd(), 1.0, 1.0, "SimpleDomain");
	}
}

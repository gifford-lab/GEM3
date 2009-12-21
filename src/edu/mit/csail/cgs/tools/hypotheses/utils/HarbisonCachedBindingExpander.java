package edu.mit.csail.cgs.tools.hypotheses.utils;

import java.util.*;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.general.Condition;
import edu.mit.csail.cgs.datasets.general.Factor;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.nouns.HarbisonRegCodeProbe;
import edu.mit.csail.cgs.ewok.verbs.*;

public class HarbisonCachedBindingExpander implements Expander<Region,BindingEvent> {

	private Vector<HarbisonRegCodeProbe> probes;
	private Factor factor;
	private Condition cond;
	
	public HarbisonCachedBindingExpander(Collection<HarbisonRegCodeProbe> ps, Factor f, Condition c) { 
		probes = new Vector<HarbisonRegCodeProbe>(ps);
		factor = f;
		cond = c;
	}

	public Iterator<BindingEvent> execute(Region a) {
		LinkedList<BindingEvent> evts = new LinkedList<BindingEvent>();

		for(HarbisonRegCodeProbe probe : probes) { 
			if(probe.overlaps(a) && 
					probe.getFactors().contains(factor.getName())) {
				
				if(probe.getConditions(factor.getName()).contains(cond.getName())) { 
					int str = probe.getBindingStrength(factor.getName(), cond.getName());
					BindingEvent evt = 
						new BindingEvent(probe.getGenome(), probe.getChrom(),
								probe.getStart(), probe.getEnd(), 1.0, 
								(str == 0 ? 200.0 : 1000.0), 
								factor.getName() + ":" + cond.getName());
					evts.addLast(evt);
				}
			}
		}
		
		return evts.iterator();
	}

}

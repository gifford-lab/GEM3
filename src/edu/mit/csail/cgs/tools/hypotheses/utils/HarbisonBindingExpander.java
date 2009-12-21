package edu.mit.csail.cgs.tools.hypotheses.utils;

import java.util.*;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.general.Condition;
import edu.mit.csail.cgs.datasets.general.Factor;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.nouns.HarbisonRegCodeProbe;
import edu.mit.csail.cgs.ewok.verbs.*;

public class HarbisonBindingExpander implements Expander<Region,BindingEvent> {

	private HarbisonRegCodeProbeGenerator gen;
	private Factor factor;
	private Condition cond;
	
	public HarbisonBindingExpander(HarbisonRegCodeProbeGenerator g, Factor f, Condition c) { 
		gen = g;
		factor = f;
		cond = c;
	}

	public Iterator<BindingEvent> execute(Region a) {
		Iterator<HarbisonRegCodeProbe> probes = gen.execute(a);
		LinkedList<BindingEvent> evts = new LinkedList<BindingEvent>();
		
		while(probes.hasNext()) { 
			HarbisonRegCodeProbe probe = probes.next();
			
			if(probe.getFactors().contains(factor.getName())) {
				if(probe.getConditions(factor.getName()).contains(cond.getName())) { 
					int str = probe.getBindingStrength(factor.getName(), cond.getName());
					BindingEvent evt = 
						new BindingEvent(probe.getGenome(), probe.getChrom(),
								probe.getStart(), probe.getEnd(), 1.0, 
								(str == 0 ? 0.005 : 0.001), 
								factor.getName() + ":" + cond.getName());
					evts.addLast(evt);
				}
			}
		}
		
		return evts.iterator();
	}

}

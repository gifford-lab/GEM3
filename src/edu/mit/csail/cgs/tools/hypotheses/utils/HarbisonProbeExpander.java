package edu.mit.csail.cgs.tools.hypotheses.utils;

import java.sql.SQLException;
import java.util.*;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.general.Condition;
import edu.mit.csail.cgs.datasets.general.Factor;
import edu.mit.csail.cgs.datasets.general.MetadataLoader;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.ewok.nouns.HarbisonRegCodeProbe;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.utils.NotFoundException;

public class HarbisonProbeExpander implements Expander<Region,Region> {
	
	private HarbisonRegCodeProbeGenerator gen;
	
	public HarbisonProbeExpander(HarbisonRegCodeProbeGenerator g) { 
		gen = g;
	}

	public Iterator<Region> execute(Region a) {
		Iterator<HarbisonRegCodeProbe> probes = gen.execute(a);
		LinkedList<Region> evts = new LinkedList<Region>();
		
		while(probes.hasNext()) { 
			HarbisonRegCodeProbe probe = probes.next();
			evts.addLast(new Region(probe));
		}
		
		return evts.iterator();
	}

}

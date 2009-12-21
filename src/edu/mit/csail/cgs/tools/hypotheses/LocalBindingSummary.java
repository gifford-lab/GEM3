package edu.mit.csail.cgs.tools.hypotheses;

import java.util.*;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.general.Factor;
import edu.mit.csail.cgs.datasets.general.Region;

public class LocalBindingSummary {

	private Region region;
	private TreeMap<Factor,Integer> binding;
	
	public LocalBindingSummary(Region r, Map<Factor,Integer> cs) { 
		region = r;
		binding = new TreeMap<Factor,Integer>(cs);
	}
	
	public Region getRegion() { return region; }
	
	public Set<Factor> getBound() { 
		Set<Factor> b = new HashSet<Factor>();
		for(Factor k : binding.keySet()) { 
			if(binding.get(k) > 0) { 
				b.add(k);
			}
		}
		return b;
	}
	
	public Map<Factor,Integer> getCounts() { 
		return binding;
	}
	
	public String toString() { 
		StringBuilder sb = new StringBuilder();
		sb.append(region.getLocationString() + "\t");
		for(Factor f : binding.keySet()) { 
			sb.append(binding.get(f) > 0 ? "1" : "0");
		}
		return sb.toString();
	}
	
}

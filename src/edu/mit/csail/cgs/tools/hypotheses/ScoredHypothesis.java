package edu.mit.csail.cgs.tools.hypotheses;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import edu.mit.csail.cgs.datasets.general.Region;

public class ScoredHypothesis implements Comparable<ScoredHypothesis> {

	private BindingHypothesis hypothesis;
	private int inconsistent;
	private Set<Region> supportingRegions;
	
	public ScoredHypothesis(BindingHypothesis bh, int sc, Collection<Region> rs) { 
		hypothesis = bh;
		inconsistent = sc;
		supportingRegions = new HashSet<Region>(rs);
	}
	
	public ScoredHypothesis(BindingHypothesis bh) { 
		hypothesis = bh;
		inconsistent = 0;
		supportingRegions = new HashSet<Region>();
	}
	
	public void addRegion(Region r, boolean supports) { 
		if(supports) { 
			supportingRegions.add(r);
		} else { 
			inconsistent += 1;
		}
	}
	
	public void reset() { 
		inconsistent = 0;
		supportingRegions.clear();
	}
	
	public BindingHypothesis getHypothesis() { return hypothesis; }
	public int getInconsistentScore() { return inconsistent; }
	public Set<Region> getSupportingRegions() { return supportingRegions; }
	
	public String toString() { return "(" + inconsistent + ") " + hypothesis; }
	
	public int hashCode() { 
		int code = 17;
		code += hypothesis.hashCode(); code *= 37;
		code += inconsistent; code *= 37;
		return code;
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof ScoredHypothesis)) { return false; }
		ScoredHypothesis sh = (ScoredHypothesis)o;
		if(!hypothesis.equals(sh.hypothesis)) { return false; }
		if(inconsistent != sh.inconsistent) { return false; }
		return true;
	}

	public int compareTo(ScoredHypothesis sh) {
		if(inconsistent < sh.inconsistent) { return -1; }
		if(inconsistent > sh.inconsistent) { return 1; }
		String s1 = hypothesis.toString();
		String s2 = sh.hypothesis.toString();
		return s1.compareTo(s2);
	}
	
}

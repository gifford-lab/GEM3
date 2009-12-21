package edu.mit.csail.cgs.datasets.expression;

import edu.mit.csail.cgs.datasets.general.Region;

public class LocatedProbe extends Probe {
	
	private Region region;

	public LocatedProbe(Probe p, Region r) { 
		super(p);
		region = r;
	}
	
	public Region getRegion() { return region; }
	
	public int hashCode() { 
		int code = super.hashCode(); 
		code += region.hashCode(); code *= 37;
		return code;
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof LocatedProbe)) { return false; }
		LocatedProbe lp = (LocatedProbe)o;
		if(!super.equals(lp)) { return false; }
		return region.equals(lp.region);
	}
	
	public String toString() { 
		return super.toString() + " \"" + region.getLocationString() + "\"";
	}
}

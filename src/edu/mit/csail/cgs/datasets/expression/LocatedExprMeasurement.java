package edu.mit.csail.cgs.datasets.expression;

import edu.mit.csail.cgs.datasets.general.Region;

public class LocatedExprMeasurement extends ExprMeasurement {
	
	private Region region;
	
	public LocatedExprMeasurement(ExprMeasurement e, Region r) { 
		super(e);
		region = r;
	}
    
    public LocatedExprMeasurement(Experiment expt, Probe p, double val, Region r) { 
        super(expt, p, val);
        region = r;
    }

	public Region getRegion() { return region; }
	
	public String toString() { 
		return super.toString() + " \"" + region.getLocationString() + "\""; 
	}
	
	public int hashCode() { 
		int code = 17;
		code += region.hashCode(); 
		code *= 37;
		return code;
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof LocatedExprMeasurement)) { return false; }
		LocatedExprMeasurement lem = (LocatedExprMeasurement)o;
		if(!super.equals(lem)) { return false; }
		return region.equals(lem.region);
	}
}

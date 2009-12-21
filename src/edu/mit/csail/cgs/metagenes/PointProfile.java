/*
 * Author: tdanford
 * Date: Aug 12, 2008
 */
package edu.mit.csail.cgs.metagenes;

import java.util.*;

import edu.mit.csail.cgs.datasets.general.Point;

public class PointProfile extends ImmutableProfile {

	private BinningParameters params;
	private Point point;
	private String name;
	private double[] values;
	private Integer hashvalue;
	private double max, min, total;
	private boolean stranded=false;
	
	public PointProfile(Point p, BinningParameters ps, double[] vs) {
		this(p, ps, vs, false);
	}
	public PointProfile(Point p, BinningParameters ps, double[] vs, boolean str) {
		stranded=str;
		params = ps;
		point = p;
		name = point.toString();
		values = vs.clone();
		min = max = total = 0.0;
		
		if(values.length != params.getNumBins()) { 
			throw new IllegalArgumentException(String.format("# Bins (%d) must match array length (%d)", 
					params.getNumBins(), values.length));
		}
		
		for(int i = 0; i < values.length; i++) { 
			max = Math.max(max, values[i]);
			min = Math.min(min, values[i]);
			total += values[i];
		}
	}
	
	public PointProfile(Point p, BinningParameters ps, String n, double[] vs) {
		params = ps;
		point = p;
		name = n;
		values = vs.clone();
		min = max = total = 0.0;

		if(values.length != params.getNumBins()) { 
			throw new IllegalArgumentException(String.format("# Bins (%d) must match array length (%d)", 
					params.getNumBins(), values.length));
		}		

		for(int i = 0; i < values.length; i++) { 
			max = Math.max(max, values[i]);
			min = Math.min(min, values[i]);
			total += values[i];
		}
	}
	
	public BinningParameters getBinningParameters() { return params; }
	public double max() { return max; }
	public double min() { return min; }
	public double total() { return total; }
	public String getName() { return name; }
	public Point getPoint() { return point; }
	public double value(int i) { return values[i]; }
	public int length() { return values.length; }
	public void setStranded(boolean s){stranded = s;}
	public boolean isStranded(){return stranded;}

	public String toString() { return name; }
	
	public int hashCode() {
		if(hashvalue == null) { 
			int code = 17;
			code += point.hashCode(); code *= 37;
			code += name.hashCode(); code *= 37;
			hashvalue = code;
		}
		
		return hashvalue;
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof PointProfile)) { return false; }
		PointProfile p = (PointProfile)o;
		if(!name.equals(p.name)) { return false; }
		if(!point.equals(p.point)) { return false; }
		return true;
	}

	public int getNumProfiles() {
		return 1;
	}
}

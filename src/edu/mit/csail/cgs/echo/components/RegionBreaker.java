package edu.mit.csail.cgs.echo.components;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.ewok.verbs.Expander;

public class RegionBreaker implements Expander<Region,Region>, SelfDescribingVerb {
	
	private int window;
	
	public RegionBreaker(int w) { 
		window = w;
	}
	
	public RegionBreaker() { 
		window = 10000;
	}

	public Iterator<Region> execute(Region a) {
		LinkedList<Region> regs = new LinkedList<Region>();
		
		for(int start = a.getStart(); start <= a.getEnd(); start += window) { 
			int end = Math.min(a.getEnd(), start + window - 1);
			Region r = new Region(a.getGenome(), a.getChrom(), start, end);
			regs.addLast(r);
		}
		
		return regs.iterator();
	}
	
	public static final String[] inputNames = { "Regions" };
	public static final EchoType[] inputClasses = { new ClassType(Region.class) };
	public static final String[] paramNames = { "Window" };
	public static final EchoType[] paramClasses = { new ClassType(Integer.class) };

	public EchoType[] getInputClasses() {
		return inputClasses;
	}

	public String[] getInputNames() {
		return inputNames;
	}

	public EchoType getOutputClass() {
		return new ClassType(Region.class);
	}

	public EchoType[] getParameterClasses() {
		return paramClasses;
	}

	public String[] getParameterNames() {
		return paramNames;
	}

	public void init(Map<String, Object> params) {
		window = (Integer)params.get(paramNames[0]);
	}

}

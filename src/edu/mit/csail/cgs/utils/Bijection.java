/*
 * Author: tdanford
 * Date: Mar 21, 2009
 */
package edu.mit.csail.cgs.utils;

import java.util.*;

public class Bijection<X,Y> {
	
	private Map<X,Y> xy;
	private Map<Y,X> yx;

	public Bijection() { 
		xy = new HashMap<X,Y>();
		yx = new HashMap<Y,X>();
	}
	
	public Bijection(X[] x, Y[] y) { 
		this();
		if(x.length != y.length) { throw new IllegalArgumentException(); }
		for(int i = 0; i < x.length; i++) { 
			addPair(x[i], y[i]);
		}
	}
	
	public Iterator<X> x() { return xy.keySet().iterator(); }
	public Iterator<Y> y() { return yx.keySet().iterator(); }
	
	public void addPair(X xval, Y yval) { 
		if(xy.containsKey(xval)) { throw new IllegalArgumentException(); }
		if(yx.containsKey(yval)) { throw new IllegalArgumentException(); }
		
		xy.put(xval, yval);
		yx.put(yval, xval);		
	}
	
	public int size() { return xy.size(); }
	public Y findY(X v) { return xy.get(v); }
	public X findX(Y v) { return yx.get(v); }
	
	public boolean containsX(X v) { return xy.containsKey(v); }
	public boolean containsY(Y v) { return yx.containsKey(v); }
}

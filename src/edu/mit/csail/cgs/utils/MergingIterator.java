/*
 * Author: tdanford
 * Date: Apr 2, 2009
 */
package edu.mit.csail.cgs.utils;

import java.util.*;

public class MergingIterator<X extends Comparable<X>> implements Iterator<X> {
	
	private ArrayList<Iterator<X>> itrs;
	private ArrayList<X> pending;
	
	public MergingIterator(Iterator<X>... is) { 
		itrs = new ArrayList<Iterator<X>>();
		for(int i = 0; i < is.length; i++) { 
			itrs.add(is[i]);
		}
		
		pending = new ArrayList<X>();
		for(int i = 0; i < itrs.size(); i++) { 
			pending.add(itrs.get(i).hasNext() ? itrs.get(i).next() : null);
		}
	}

	public boolean hasNext() {
		for(int i = 0; i < itrs.size(); i++) { 
			if(pending.get(i) != null) { 
				return true;
			}
		}
		return false;
	}

	public X next() {
		int min = -1;
		for(int i = 0; i < pending.size(); i++) { 
			X val = pending.get(i);
			if(min == -1 || (val != null && val.compareTo(pending.get(min)) < 0)) { 
				min = i;
			}
		}
		X minVal = min != -1 ? pending.get(min) : null;
		
		if(min != -1) { 
			pending.set(min, itrs.get(min).hasNext() ? itrs.get(min).next() : null);
		}
		return minVal;
	}

	public void remove() {
		throw new UnsupportedOperationException();
	}
}

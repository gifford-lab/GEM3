/*
 * Author: tdanford
 * Date: Apr 9, 2009
 */
package edu.mit.csail.cgs.utils.iterators;

import java.util.*;

public class BacktrackingIterator<X> implements Iterator<X> {
	
	private LinkedList<X> pending;
	private Iterator<X> itr;
	
	public BacktrackingIterator(Iterator<X> i) { 
		itr = i;
		pending = new LinkedList<X>();
	}
	
	public void addNext(X v) { 
		pending.addFirst(v);
	}

	public boolean hasNext() {
		return !pending.isEmpty() || itr.hasNext();
	}

	public X next() {
		if(pending.isEmpty()) { 
			return itr.next(); 
		} else { 
			return pending.removeFirst();
		}
	}

	public void remove() {
		throw new UnsupportedOperationException();
	}
}

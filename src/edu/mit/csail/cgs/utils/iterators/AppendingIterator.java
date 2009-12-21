/*
 * Author: tdanford
 * Date: Apr 22, 2009
 */
package edu.mit.csail.cgs.utils.iterators;

import java.util.*;

/**
 * This is so not thread-safe, it's not even funny.  
 * The possibilities of using this class to shoot yourself in the foot are 
 * *NUMEROUS*.  
 * 
 * On the other hand, it's totally useful to have lying around. 
 * 
 * @author tdanford
 *
 * @param <X>
 */
public class AppendingIterator<X> implements Iterator<X> {
	
	private LinkedList<X> values;
	
	public AppendingIterator() { 
		values = new LinkedList<X>();
	}
	
	public void addValue(X v) { 
		values.addLast(v);
	}

	public boolean hasNext() {
		return !values.isEmpty();
	}

	public X next() {
		return values.removeFirst();
	}

	public void remove() {
		throw new UnsupportedOperationException();
	}

}

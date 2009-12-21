/*
 * Author: tdanford
 * Date: Apr 2, 2009
 */
package edu.mit.csail.cgs.utils.iterators;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;

public class SortingIterator<X extends Comparable<X>> implements Iterator<X> {
	
	private Iterator<X> itr;
	private Object[] array;
	private int idx;
	
	public SortingIterator(Iterator<X> itr) { 
		this.itr = itr;
		array = null;
		idx = -1;
	}
	
	private void loadAndSort() { 
		ArrayList<X> values = new ArrayList<X>();
		while(itr.hasNext()) { 
			values.add(itr.next());
		}
		array = values.toArray();
		Arrays.sort(array);
		idx = 0;
	}

	public boolean hasNext() {
		if(array == null) { loadAndSort(); }
		return idx < array.length;
	}

	public X next() {
		if(array == null) { loadAndSort(); }
		return idx < array.length ? (X)array[idx++] : null;
	}

	public void remove() {
		throw new UnsupportedOperationException();
	}
}

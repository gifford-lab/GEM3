/**
 * 
 */
package edu.mit.csail.cgs.viz.components;

import java.util.*;

/**
 * @author tdanford
 */
public class SelectionEvent<X> extends EventObject {
	
	private LinkedList<X> values;

	public SelectionEvent(Object src, Collection<X> vals) { 
		super(src);
		values = new LinkedList<X>(vals);
	}
	
	public SelectionEvent(Object src, X val) {
		super(src);
		values = new LinkedList<X>();
		values.addLast(val);
	}
	
	public Collection<X> getSelectedValues() { return values; }
    public X getFirstValue() { return values.getFirst(); }
}

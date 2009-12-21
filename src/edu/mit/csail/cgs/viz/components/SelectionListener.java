/*
 * Author: tdanford
 * Date: May 27, 2008
 */
package edu.mit.csail.cgs.viz.components;

public interface SelectionListener<X> {
	public void selectionMade(SelectionEvent<X> evt);
}

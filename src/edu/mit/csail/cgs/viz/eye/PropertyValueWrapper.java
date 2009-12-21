/*
 * Author: tdanford
 * Date: Sep 16, 2008
 */
package edu.mit.csail.cgs.viz.eye;

public class PropertyValueWrapper<X> extends ModelPaintableProperty<X> {

	private X value;
	
	public PropertyValueWrapper(String k, X v) { 
		super(k);
		value = v;
	}
	
	public X getValue() { return value; }
	
	public void setValue(X v) { 
		value = v;
		dispatchChangedEvent();
	}
}

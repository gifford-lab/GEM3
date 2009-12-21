/*
 * Author: tdanford
 * Date: Sep 16, 2008
 */
package edu.mit.csail.cgs.viz.eye;

public class SubPropertyWrapper<X> extends ModelPaintableProperty<X> {

	private String outerKey;
	private ModelPaintableProperty<X> inner;
	
	public SubPropertyWrapper(String k, ModelPaintableProperty<X> i) { 
		super(String.format("%s.%s", k, i.getKey()));
		outerKey = k;
		inner = i;
	}
	
	public X getValue() { return inner.getValue(); }
	public ModelPaintableProperty<X> getInnerProperty() { return inner; }
}

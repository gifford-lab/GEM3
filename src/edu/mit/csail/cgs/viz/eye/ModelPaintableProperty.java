/*
 * Author: tdanford
 * Date: Sep 16, 2008
 */
package edu.mit.csail.cgs.viz.eye;

import java.util.*;

public abstract class ModelPaintableProperty<X> {

	private String key;
	private LinkedList<ModelPaintablePropertyListener> listeners;
	
	public ModelPaintableProperty(String k) {
		key = k;
		listeners = new LinkedList<ModelPaintablePropertyListener>();
	}
	
	public String getKey() { return key; }
	public abstract X getValue();
	
	public void addListener(ModelPaintablePropertyListener l) { 
		listeners.add(l);
	}
	
	public void removeListener(ModelPaintablePropertyListener l) { 
		listeners.remove(l);
	}
	
	protected void dispatchChangedEvent() { 
		for(ModelPaintablePropertyListener list : listeners) { 
			list.propertyChanged(this);
		}
	}
	
	public int hashCode() { return key.hashCode(); }
	
	public boolean equals(Object o) { 
		if(!(o instanceof ModelPaintableProperty)) { return false; }
		ModelPaintableProperty p = (ModelPaintableProperty)o;
		return p.key.equals(key);
	}
	
	public String toString() { return key; }
}

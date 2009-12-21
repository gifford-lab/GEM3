/**
 * @author tdanford
 */
package edu.mit.csail.cgs.viz.utils;

import java.util.EventObject;

public class ProgressEvent extends EventObject {
	
	private String key;
	private int value;
	
	public ProgressEvent(Object src, String key, int val) { 
		super(src);
		this.key = key;
		value = val;
	}

	public String getKey() { return key; }
	public int getValue() { return value; }
}
        
/*
 * Author: tdanford
 * Date: Jan 19, 2009
 */
package edu.mit.csail.cgs.utils.models;

import java.util.Iterator;

public interface Timer {

	public void addTiming(Timing t);
	public void addTimings(Iterator<Timing> ts);
}

/*
 * Author: tdanford
 * Date: Jan 19, 2009
 */
package edu.mit.csail.cgs.utils.models;


public class Timing extends Model {
	
	public Integer size; 
	public Double seconds;
	
	public Timing() {}
	
	public Timing(Integer s, Double t) { 
		size = s;
		seconds = t;
	}
	
	public String toString() { 
		return String.format("%d -> %.3fs", size, seconds);
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof Timing)) { return false; }
		Timing t = (Timing)o;
		return t.size.equals(size) && t.seconds.equals(seconds);
	}

	public int hashCode() { 
		long secondbits = Double.doubleToLongBits(seconds);
		int code = 17;
		code += (int)(secondbits >> 32); code *= 37;
		code += size; code *= 37;
		return code;
	}
}



package edu.mit.csail.cgs.utils;

import java.io.*;
import java.util.Collection;
import java.util.LinkedList;
import java.util.Random;

public class Range implements Comparable, Saveable {
	
	protected long start, stop;
    protected int intStart, intStop;
	
	public Range(long start, long stop) { 
		this.start = start; this.stop = stop;
        intStart = (int)start;
        intStop = (int)stop;
	}
	
	public Range(DataInputStream dis) throws IOException { 
		start = dis.readLong();
		stop = dis.readLong();
        intStart = (int)start;
        intStop = (int)stop;
	}
	
	public int getStart() { return intStart; }
	public int getStop() { return intStop; }
    
    public long getLongStart() { return start; }
    public long getLongStop() { return stop; }
		
	public int compareTo(Object o) {
		Range b = (Range)o;
		if(start < b.start) { return -1; }
		if(start > b.start) { return 1; }
		if(stop < b.stop) { return -1; }
		if(stop > b.stop) { return 1; }
		return 0;
	}
	
	public boolean overlaps(Range b) { 
		if(start <= b.start && stop >= b.start) { return true; }
		if(b.start <= start && b.stop >= start) { return true; }
		return false;
	}
	
	public boolean contains(Range b) { 
		return start < b.start && stop > b.stop;
	}
	
	public long getWidth() { return stop - start + (long)1; }
	
	public boolean equals(Object o) { 
		if(!(o instanceof Range)) { return false; }
		Range b = (Range)o;
		if(start != b.start || stop != b.stop) { return false; }
		return true;
	}
	
	public int hashCode() { 
		int code = 17; 
		code += (int)(start >> 32); code *= 37;
		code += (int)(stop >> 32); code *= 37;
		return code;
	}
	
	public String toString() { 
		return start + "-" + stop;
	}

	public void save(DataOutputStream dos) throws IOException {
		dos.writeLong(start);
		dos.writeLong(stop);
	}
}

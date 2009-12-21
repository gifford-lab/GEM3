package edu.mit.csail.cgs.echo;

import java.util.*;

public class IteratorInputQueue<X> implements InputQueue<X> {
	
	private Iterator<X> itr;
	
	public IteratorInputQueue(Iterator<X> i) { itr = i; }

	public X getFirstValue() {
		return itr.next();
	}

	public boolean isEmpty() {
		return !itr.hasNext();
	}

	public boolean isFinished() {
		return !itr.hasNext();
	}
    
}

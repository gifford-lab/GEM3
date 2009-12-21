package edu.mit.csail.cgs.echo;

import java.util.*;

public class EchoQueue<X> implements InputQueue<X> {

	private LinkedList<X> internalList;
	private boolean finished;
	private boolean expanding;
	
	public EchoQueue() { 
		internalList = new LinkedList<X>();
		finished = false;
		expanding = true;
	}
	
	public EchoQueue(boolean exp) { 
		internalList = new LinkedList<X>();
		finished = false;
		expanding = exp;
	}
	
	public boolean isExpanding() { return expanding; }
	
	/* (non-Javadoc)
	 * @see edu.mit.csail.cgs.echo.InputQueue#isFinished()
	 */
	public synchronized boolean isFinished() { 
		return finished && internalList.isEmpty();
	}
	
	public synchronized void finish() { 
		finished = true;
	}
	
	public synchronized void addValue(X v) { 
		if(finished) { throw new IllegalStateException("Queue is finished."); }
		internalList.addLast(v);
	}
	
	/* (non-Javadoc)
	 * @see edu.mit.csail.cgs.echo.InputQueue#isEmpty()
	 */
	public synchronized boolean isEmpty() { 
		return internalList.isEmpty();
	}
	
	/* (non-Javadoc)
	 * @see edu.mit.csail.cgs.echo.InputQueue#getFirstValue()
	 */
	public synchronized X getFirstValue() { 
		X val = internalList.removeFirst();
		return val;
	}
	
	public synchronized void reset() { 
		internalList = new LinkedList<X>();
		finished = false;
	}
}

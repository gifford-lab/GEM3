package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;

/**
 * A Sink is a verb that consumes objects but produces no output
 */
public interface Sink<X> { 
	// this should basically be an alias to (a) calling init(), (b) calling
	// consume(val) on all the val's in the Iterator, and then (c) calling finish().
    public void consume(Iterator<X> itr);

	public void init();
	public void consume(X val);
	public void finish();
}

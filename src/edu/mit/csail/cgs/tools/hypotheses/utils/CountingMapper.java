package edu.mit.csail.cgs.tools.hypotheses.utils;

import java.util.Iterator;

import edu.mit.csail.cgs.ewok.verbs.Mapper;
import edu.mit.csail.cgs.utils.Closeable;

public class CountingMapper<X> implements Mapper<Iterator<X>, Integer> {
	
	public CountingMapper() { 
	}

	public Integer execute(Iterator<X> a) {
		int c = 0;
		while(a.hasNext()) { c += 1; }
		if(a instanceof Closeable) { 
			((Closeable)a).close();
		}
		return c;
	}

}

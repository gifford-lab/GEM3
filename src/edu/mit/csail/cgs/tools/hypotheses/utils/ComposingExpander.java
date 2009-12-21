package edu.mit.csail.cgs.tools.hypotheses.utils;

import java.util.Iterator;

import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.ewok.verbs.Mapper;
import edu.mit.csail.cgs.ewok.verbs.MapperIterator;

public class ComposingExpander<X,Y,Z> implements Expander<X,Z> {
	
	private Expander<X,Y> m1;
	private Mapper<Y,Z> m2;
	
	public ComposingExpander(Expander<X,Y> first, Mapper<Y,Z> second) { 
		m1 = first;
		m2 = second;
	}

	public Iterator<Z> execute(X a) {
		return new MapperIterator<Y,Z>(m2, m1.execute(a));
	}
}

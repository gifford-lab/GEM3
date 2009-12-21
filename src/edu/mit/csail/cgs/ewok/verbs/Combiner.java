package edu.mit.csail.cgs.ewok.verbs;

import java.util.Iterator;

public interface Combiner<A,B,C> {

    public C execute(A a, B b);

}

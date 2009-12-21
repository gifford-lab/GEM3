/*
 * Created on Feb 20, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.ewok.verbs;

import java.util.Iterator;

public class ExpanderCastWrapper<A,B,C extends B> implements Expander<A,B> {
    
    private Expander<A,C> internalExpander;
    
    public ExpanderCastWrapper(Expander<A,C> i) { internalExpander = i; }

    public Iterator<B> execute(A a) {
        return new CastingIterator<B,C>(internalExpander.execute(a));
    }
}

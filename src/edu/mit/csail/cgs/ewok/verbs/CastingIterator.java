/*
 * Created on Feb 20, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;

public class CastingIterator<A,B extends A> implements Iterator<A> {
    
    private Iterator<B> internal;
    
    public CastingIterator(Iterator<B> itr) { internal = itr; }

    public boolean hasNext() {
        return internal.hasNext();
    }

    public A next() {
        return (A)internal.next();
    }

    public void remove() {
        internal.remove();
    } 
}

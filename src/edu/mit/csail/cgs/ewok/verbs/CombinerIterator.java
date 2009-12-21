/*
 * Created on Mar 6, 2006
 */
package edu.mit.csail.cgs.ewok.verbs;

import java.util.Iterator;

/**
 * @author tdanford
 */
public class CombinerIterator<A,B,C> implements Iterator<C> {
    
    private Iterator<A> aitr;
    private Iterator<B> bitr;
    private Combiner<A,B,C> combiner;

    public CombinerIterator(Combiner<A,B,C> comb, Iterator<A> a, Iterator<B> b) {
        combiner = comb;
        aitr = a;
        bitr = b;
    }

    /* (non-Javadoc)
     * @see java.util.Iterator#hasNext()
     */
    public boolean hasNext() {
        return aitr.hasNext() || bitr.hasNext();
    }

    /* (non-Javadoc)
     * @see java.util.Iterator#next()
     */
    public C next() {
        A first = null;
        B second = null;
        
        if(aitr.hasNext()) { first = aitr.next(); }
        if(bitr.hasNext()) { second = bitr.next(); }
        
        return combiner.execute(first, second);
    }

    /* (non-Javadoc)
     * @see java.util.Iterator#remove()
     */
    public void remove() {
        throw new UnsupportedOperationException();
    }

}

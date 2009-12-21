/*
 * Created on Mar 8, 2006
 */
package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;

/**
 * @author tdanford
 */
public class MultiIterator<X> implements Iterator<X> {
    
    private LinkedList<Iterator<? extends X>> itrs;

    public MultiIterator(Collection<Iterator<? extends X>> c) { 
        itrs = new LinkedList<Iterator<? extends X>>(c);
        clearEmpty();
    }
    
    public MultiIterator(Iterator<? extends X> first, Iterator<? extends X> second) { 
        itrs = new LinkedList<Iterator<? extends X>>();
        itrs.addLast(first);
        itrs.addLast(second);
        clearEmpty();
    }

    /* (non-Javadoc)
     * @see java.util.Iterator#hasNext()
     */
    public boolean hasNext() {
        return !itrs.isEmpty() && itrs.getFirst().hasNext();
    }
    
    private void clearEmpty() { 
        while(!itrs.isEmpty() && !itrs.getFirst().hasNext()) { 
            itrs.removeFirst();
        }        
    }

    /* (non-Javadoc)
     * @see java.util.Iterator#next()
     */
    public X next() {
        Iterator<? extends X> nextItr = itrs.getFirst();
        X val = nextItr.next();
        clearEmpty();
        return val;
    }

    /* (non-Javadoc)
     * @see java.util.Iterator#remove()
     */
    public void remove() {
        throw new UnsupportedOperationException();
    }

}

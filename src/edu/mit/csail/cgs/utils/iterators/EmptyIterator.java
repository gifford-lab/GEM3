/*
 * Created on Mar 23, 2006
 */
package edu.mit.csail.cgs.utils.iterators;

import java.util.Iterator;

/**
 * @author tdanford
 */
public class EmptyIterator<X> implements Iterator<X> {

    public EmptyIterator() {
    }

    /* (non-Javadoc)
     * @see java.util.Iterator#hasNext()
     */
    public boolean hasNext() {
        return false;
    }

    /* (non-Javadoc)
     * @see java.util.Iterator#next()
     */
    public X next() {
        return null;
    }

    /* (non-Javadoc)
     * @see java.util.Iterator#remove()
     */
    public void remove() {
        throw new UnsupportedOperationException();
    }

}

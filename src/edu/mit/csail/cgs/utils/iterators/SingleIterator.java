package edu.mit.csail.cgs.utils.iterators;
import java.util.*;

public class SingleIterator<A> implements Iterator<A> {

    private A a;

    public SingleIterator(A a) {
        this.a = a;
    }
    public boolean hasNext() {
        return (a != null);
    }
    public A next() {
        A temp = a;
        a = null;
        return temp;
    }
    public void remove() {
        throw new UnsupportedOperationException("Can't remove from a SingleIterator");
    }

}

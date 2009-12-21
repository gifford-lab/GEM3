package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;

public class FilterIterator<A,B> implements Iterator<B> {
    private Filter<A,B> filter;
    private Iterator<A> input;
    private B next;

    public FilterIterator(Filter<A,B> filter, Iterator<A> input) {
        this.filter = filter;
        this.input = input;
        next = null;
        getNext();
    }

    private void getNext() {
        while (next == null && input != null) {
            if (input.hasNext()) {
                next = filter.execute(input.next());
            } else {
                next = null;
                input = null;
            }
        }
    }


    public boolean hasNext() {
        return next != null;
    }

    public B next() {
        B toreturn = next;
        next = null;
        getNext();
        return toreturn;
    }

    public void remove() {
        throw new UnsupportedOperationException("Can't remove from a FilterIterator");
    }


}

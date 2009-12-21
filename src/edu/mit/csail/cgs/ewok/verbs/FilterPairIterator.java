package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;
import edu.mit.csail.cgs.utils.Pair;

public class FilterPairIterator<A,B> implements Iterator<Pair<A,B>> {
    
    private Filter<A,B> filter;
    private Iterator<A> input;
    private Pair<A,B> next;

    public FilterPairIterator(Filter<A,B> filter, Iterator<A> input) {
        this.filter = filter;
        this.input = input;
        next = null;
        getNext();
    }

    private void getNext() {
        while (input != null && next.getLast() == null) {
            if (input.hasNext()) {
                A na = input.next();
                next = new Pair<A,B>(na, filter.execute(input.next()));
            } else {
                next = null;
                input = null;
            }
        }
    }


    public boolean hasNext() {
        return next != null;
    }

    public Pair<A,B> next() {
        Pair<A,B> toreturn = next;
        next = null;
        getNext();
        return toreturn;
    }

    public void remove() {
        throw new UnsupportedOperationException("Can't remove from a FilterIterator");
    }


}

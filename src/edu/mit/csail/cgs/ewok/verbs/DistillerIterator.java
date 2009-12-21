package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;

public class DistillerIterator<A,B> implements Iterator<B> {
    
    private Distiller<A,B> filter;
    private Iterator<A> input;
    private B next;
    private boolean finished;

    public DistillerIterator(Distiller<A,B> filter, Iterator<A> input) {
        this.filter = filter;
        this.input = input;
        next = null;
        finished = false;
        getNext();
    }

    private void getNext() {
        while(!finished && next == null && input != null) {
            if(input.hasNext()) {
                next = filter.execute(input.next());
            } else { 
                next = filter.getCurrent();
                finished = true;
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
        throw new UnsupportedOperationException("Can't remove from a DistillerIterator");
    }


}

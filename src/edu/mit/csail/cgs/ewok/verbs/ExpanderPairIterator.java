package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;
import edu.mit.csail.cgs.utils.Pair;

public class ExpanderPairIterator<A,B> implements Iterator<Pair<A,Iterator<B>>> {
    private Expander<A,B> expander;
    private Iterator<A> input;

    public ExpanderPairIterator(Expander<A,B> mapper, Iterator<A> input) {
        this.expander = mapper;
        this.input = input;
    }

    public boolean hasNext() {
        return input.hasNext();
    }

    public Pair<A,Iterator<B>> next() {
        A anext = input.next();
        Iterator<B> bnext = expander.execute(anext);
        return new Pair<A,Iterator<B>>(anext, bnext);
    }

    public void remove() {
        throw new UnsupportedOperationException("Can't remove from a MapperIterator");
    }


}

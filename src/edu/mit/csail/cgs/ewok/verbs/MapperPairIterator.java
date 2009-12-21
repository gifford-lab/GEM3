package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;
import edu.mit.csail.cgs.utils.Pair;

public class MapperPairIterator<A,B> implements Iterator<Pair<A,B>> {
    private Mapper<A,B> mapper;
    private Iterator<A> input;

    public MapperPairIterator(Mapper<A,B> mapper, Iterator<A> input) {
        this.mapper = mapper;
        this.input = input;
    }

    public boolean hasNext() {
        return input.hasNext();
    }

    public Pair<A,B> next() {
        A anext = input.next();
        B bnext = mapper.execute(anext);
        return new Pair<A,B>(anext, bnext);
    }

    public void remove() {
        throw new UnsupportedOperationException("Can't remove from a MapperIterator");
    }


}

package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;

/** 
 * A MapperIterator takes a Mapper<A,B> and an Iterator<A> and returns
 *  an Iterator<B> that is the result of applying the Mapper to each
 *  element of the input iterator.  
*/

public class MapperIterator<A,B> implements Iterator<B> {
    private Mapper<A,B> mapper;
    private Iterator<A> input;

    public MapperIterator(Mapper<A,B> mapper, Iterator<A> input) {
        this.mapper = mapper;
        this.input = input;
    }

    public boolean hasNext() {
        return input.hasNext();
    }

    public B next() {
        return mapper.execute(input.next());
    }

    public void remove() {
        throw new UnsupportedOperationException("Can't remove from a MapperIterator");
    }


}

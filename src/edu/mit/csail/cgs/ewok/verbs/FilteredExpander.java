package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;

public class FilteredExpander<A,B,C> implements Expander<A,C> {

    private Expander<A,B> expander;
    private Filter<B,C> filter;

    public FilteredExpander(Expander<A,B> e, Filter<B,C> f) {
        expander = e;
        filter = f;
    }

    public Iterator<C> execute(A input) {
        return new FilterIterator<B,C>(filter, expander.execute(input));
    }

}
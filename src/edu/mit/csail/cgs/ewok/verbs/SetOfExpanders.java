package edu.mit.csail.cgs.ewok.verbs;

import java.util.Iterator;
import java.util.Set;
import java.util.HashSet;

/* Holds a set of expanders and calls each of them when execute() is called.  The result
   contains the results of all the sub-expanders in an arbitrary order.  Duplicates
   are removed from the results 
*/

public class SetOfExpanders<A,B> implements Expander<A,B> {

    private Set<Expander<A,B>> expanders;

    public SetOfExpanders(Set<Expander<A,B>> initial) {
        expanders = new HashSet<Expander<A,B>>();
        expanders.addAll(initial);
    }
    public SetOfExpanders() {
        expanders = new HashSet<Expander<A,B>>();
    }
    public SetOfExpanders(Expander<A,B> initial) {
        expanders = new HashSet<Expander<A,B>>();
        expanders.add(initial);
    }

    public Iterator<B> execute(A a) {
        HashSet<B> results = new HashSet<B>();
        for (Expander<A,B> e : expanders) {
            Iterator<B> iter = e.execute(a);
            while (iter.hasNext()) {
                results.add(iter.next());
            }
        }
        return results.iterator();
    }


}

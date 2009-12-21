package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;

/* Maps an Iterator<X> to an Iterator<X> by only keeping the top N elements
   as defined by the natural ordering of the output of the mapper from X to comparable.
*/

public class TopNMapper<X> implements Mapper<Iterator<X>, Iterator<X>>, Comparator<X> {

    private Mapper<X,Comparable> mapper;
    private int n;

    public TopNMapper(Mapper<X,Comparable> m, int n) {
        mapper = m;
        this.n = n;
    }

    public Iterator<X> execute(Iterator<X> input) {
        List<X> list = new ArrayList<X>();
        while (input.hasNext()) {
            list.add(input.next());
        }
        Collections.sort(list,this);
        return list.subList(list.size() - n, list.size()).iterator();
    }

    public int compare(X a, X b) {
        return mapper.execute(a).compareTo(mapper.execute(b));
    }
}
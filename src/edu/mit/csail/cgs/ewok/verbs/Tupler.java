/*
 * Created on Aug 10, 2006
 */
package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;

/**
 * @author tdanford
 */
public class Tupler<X> implements Mapper<Iterator<X>,Iterator<Vector<X>>> {
    
    private int length;
    private LinkedList<X> current;

    public Tupler(int l) {
        length = l;
        current = new LinkedList<X>();
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.ewok.verbs.Distiller#execute(null)
     */
    public Iterator<Vector<X>> execute(Iterator<X> a) {
        current.clear();
        LinkedList<Vector<X>> tuples = new LinkedList<Vector<X>>();
        
        while(a.hasNext()) { 
            X val = a.next();
            current.addLast(val);
            if(current.size() > length) { 
                current.removeFirst();
                tuples.add(new Vector<X>(current));
            }
        }
        
        return tuples.iterator();
    }
}

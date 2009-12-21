/*
 * Created on Apr 1, 2006
 */
package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;

import edu.mit.csail.cgs.datasets.orthology.*;
import edu.mit.csail.cgs.datasets.species.Genome;

/**
 * @author tdanford
 */
public class OrthologyExpander implements Expander<String,String> {
    
    private Genome firstGenome, secondGenome;
    private OrthologyLoader loader;
    private OrthologyMapping mapping;

    public OrthologyExpander(OrthologyLoader l, OrthologyMapping m, Genome first, Genome second) {
        loader = l;
        firstGenome = first;
        secondGenome = second;
        mapping = m;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.ewok.verbs.Expander#execute(null)
     */
    public Iterator<String> execute(String input) {
        Collection<OrthologyPair> pairs = loader.getFirstNameGenomePairs(input, firstGenome);
        LinkedList<String> lst = new LinkedList<String>();
        
        for(OrthologyPair op : pairs) { 
            if(op.getGenome2().equals(secondGenome)) { 
                lst.addLast(op.getName2());
            }
        }
        
        return lst.iterator();
    }

}

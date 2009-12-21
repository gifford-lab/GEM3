/*
 * Created on Nov 28, 2006
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.ewok.verbs;

import java.util.Iterator;

import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;

public class RefGeneByNameGenerator implements Expander<String,Gene> { 

    private Genome genome;
    private RefGeneGenerator gen;
    
    public RefGeneByNameGenerator(Genome g) { 
        genome = g;
        gen = new RefGeneGenerator(g);
    }

    public Iterator<Gene> execute(String a) {
        return gen.byName(a);
    }
}

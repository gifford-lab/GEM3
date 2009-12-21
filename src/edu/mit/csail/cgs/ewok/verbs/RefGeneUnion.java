package edu.mit.csail.cgs.ewok.verbs;

import java.util.Iterator;
import java.util.Set;
import java.util.HashSet;
import java.util.Collections;
import java.util.ArrayList;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Gene;

/**
 * RefGeneUnion queries a set of RefGeneGenerators and presents the 
 * combined results.  Use this when you want a single generator but
 * want to pull genes from several sets of annotations.
 */
public class RefGeneUnion<X extends Region> implements Expander<X,Gene> {

    private Set<RefGeneGenerator<X>> generators;

    public RefGeneUnion() {
        generators = new HashSet<RefGeneGenerator<X>>();
    }

    public void addGenerator(RefGeneGenerator<X> gen) {
        generators.add(gen);
    }

    public Iterator<Gene> execute (X region) {
        ArrayList<Gene> output = new ArrayList<Gene>();
        for (RefGeneGenerator<X> gen : generators) {
            Iterator<Gene> iter = gen.execute(region);
            while (iter.hasNext()) {
                output.add(iter.next());
            }
        }
        Collections.sort(output);
        return output.iterator();
    }

}
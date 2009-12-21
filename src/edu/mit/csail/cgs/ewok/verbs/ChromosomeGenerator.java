package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;

/* gives an iterator over all of the chromosomes in a genome */

public class ChromosomeGenerator<X extends Genome> implements Expander<X,Region> {
    
    public Iterator<Region> execute (X genome) {
        List<String> names = genome.getChromList();
        List<Region> chroms = new ArrayList<Region>();
        for (int i = 0; i < names.size(); i++) {
            chroms.add(new Region(genome,
                                  names.get(i),
                                  0,
                                  genome.getChromLength(names.get(i))));
        }
        return chroms.iterator();
    }
}

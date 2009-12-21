package edu.mit.csail.cgs.ewok.verbs;

import org.biojava.bio.seq.*;
import java.util.*;
import edu.mit.csail.cgs.ewok.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;

public class FH2RegionIter implements Iterator<Region> {
    private Iterator<Feature> iter;
    private Genome g;
    public FH2RegionIter(FeatureHolder fh, Genome g) {
        iter = fh.features();
        this.g = g;
    }
    public boolean hasNext() {
        return iter.hasNext();
    }
    public Region next() {
        Feature f = iter.next();
        return new Region(g,
                          f.getSequence().getName(),
                          f.getLocation().getMin(),
                          f.getLocation().getMax());
    }
    public void remove() {
        throw new UnsupportedOperationException("Can't remove from a RH2RegionIter");
    }

}

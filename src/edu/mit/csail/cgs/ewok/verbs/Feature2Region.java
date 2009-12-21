package edu.mit.csail.cgs.ewok.verbs;

import org.biojava.bio.seq.*;
import java.util.*;
import edu.mit.csail.cgs.ewok.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.species.Genome;

public class Feature2Region implements Mapper<Feature,NamedRegion> {
    private Genome g;
    
    public Feature2Region(Genome g) {
        this.g = g;
    }
    
    public NamedRegion execute (Feature f) {
        return new NamedRegion(g,
                          f.getSequence().getName(),
                          f.getLocation().getMin(),
                          f.getLocation().getMax(), f.getType());
    }
}

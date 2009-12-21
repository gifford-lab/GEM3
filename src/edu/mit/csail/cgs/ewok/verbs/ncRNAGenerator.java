package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;

import org.biojava.bio.symbol.Location;
import org.biojava.bio.*;
import org.biojava.bio.symbol.*;
import org.biojava.bio.seq.*;



public class ncRNAGenerator implements Expander<Region,Gene> {
    private Genome genome;
    private FeatureFilter ncrnaFilter;
    public ncRNAGenerator(Genome g) {
    }
    public ncRNAGenerator (Genome g, String source) {
    }

    public Iterator<Gene> execute(Region region) {
        throw new UnsupportedOperationException("ncRNAGenerator needs to be updated for ucsc");
    }
}


package edu.mit.csail.cgs.ewok.verbs;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.ScoredRegion;
import edu.mit.csail.cgs.datasets.species.Genome;

public class PhastConsGenerator<X extends Region> extends ScoredRegionGenerator<X> {
    public PhastConsGenerator(Genome g, String tablename) {
        super(g,tablename);
    }
    public String getColumnsSQL() { return "chromStart, chromEnd, sumData / (chromEnd - chromStart + 1)";}
}
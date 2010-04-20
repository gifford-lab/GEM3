package edu.mit.csail.cgs.ewok.verbs;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.ScoredRegion;
import edu.mit.csail.cgs.datasets.species.Genome;

public class PhastConsGenerator<X extends Region> extends ScoredRegionGenerator<X> {
    public PhastConsGenerator(Genome g, String tablename) {
        super(g,tablename);
    }
    public String getSQL() {
        return "select chromStart, chromEnd, sumData / (chromEnd - chromStart + 1) from " + getTablename() + " where chrom = ? and " +
            "chromStart <= ? and chromEnd >= ? UNION " +
            "select chromStart, chromEnd, sumData / (chromEnd - chromStart + 1) from " + getTablename() + " where chrom = ? and " +
            "chromStart <= ?  and chromStart >= ?";
    }
}
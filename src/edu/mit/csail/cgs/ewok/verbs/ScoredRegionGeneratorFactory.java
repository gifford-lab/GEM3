package edu.mit.csail.cgs.ewok.verbs;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.ScoredRegion;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.RegionExpanderFactory;

public class ScoredRegionGeneratorFactory implements RegionExpanderFactory<ScoredRegion> {
    private String type;
    
    public ScoredRegionGeneratorFactory() {
        type = "ScoredRegion";
    }

    public void setType(String t) {type = t;}
    public String getType() {return type;}
    public String getProduct() {return "ScoredRegion";}
    public Expander<Region, ScoredRegion> getExpander(Genome g) {
        return getExpander(g,type);
    }

    public Expander<Region, ScoredRegion> getExpander(Genome g, String type) {
        return new ScoredRegionGenerator(g,type);
    }
}

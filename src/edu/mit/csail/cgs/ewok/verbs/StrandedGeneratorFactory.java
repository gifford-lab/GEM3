package edu.mit.csail.cgs.ewok.verbs;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.RegionExpanderFactory;

public class StrandedGeneratorFactory implements RegionExpanderFactory<StrandedRegion> {
    String type;

    public StrandedGeneratorFactory() {
        type = "StrandedRegion";
    }
    public StrandedGeneratorFactory(String t) {
        type = t;
    }

    public void setType(String t) {type = t;}
    public String getType() {return type;}
    public String getProduct() {return "StrandedRegion";}
    public Expander<Region, StrandedRegion> getExpander(Genome g) {
        return getExpander(g,type);
    }

    public Expander<Region, StrandedRegion> getExpander(Genome g, String type) {
        if (type == null) {
            throw new NullPointerException("StrandedGenerator must have a type");
        } else {
            return new StrandedRegionGenerator(g, type);
        }
    }
}

package edu.mit.csail.cgs.ewok.verbs;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.RegionExpanderFactory;

/**
 * @author tdanford
 */
public class HarbisonRegCodeGeneratorFactory implements RegionExpanderFactory<Gene> {
    private String type;

    public HarbisonRegCodeGeneratorFactory() {
    }

    public void setType(String t) {type = t;}
    public String getType() {return type;}
    public String getProduct() {return "HarbisonRegCodeRegion";}
    public Expander<Region, Gene> getExpander(Genome g) {
        return getExpander(g,type);
    }

    public Expander<Region, Gene> getExpander(Genome g, String type) {
        if (type == null) {
            return new HarbisonRegCodeGenerator(g);
        } else {
            return new HarbisonRegCodeGenerator(g, type);
        }
    }
}

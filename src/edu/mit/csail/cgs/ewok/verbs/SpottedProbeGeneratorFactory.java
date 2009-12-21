package edu.mit.csail.cgs.ewok.verbs;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.SpottedProbe;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.RegionExpanderFactory;

public class SpottedProbeGeneratorFactory implements RegionExpanderFactory<SpottedProbe> {
    String table;

    public SpottedProbeGeneratorFactory() {
        table = "spottedArray";
    }
    public SpottedProbeGeneratorFactory(String t) {
        table = t;
    }

    public void setType(String t) {table = t;}
    public String getType() {return table;}
    public String getProduct() {return "SpottedProbe";}
    
    public Expander<Region,SpottedProbe> getExpander(Genome g) {
        return getExpander(g,table);
    }

    public Expander<Region,SpottedProbe> getExpander(Genome g, String type) {
        if (type == null) {
            throw new NullPointerException("SpottedProbeGenerator must have a type");
        } else {
            return new SpottedProbeGenerator(g, table);
        }
    }
}

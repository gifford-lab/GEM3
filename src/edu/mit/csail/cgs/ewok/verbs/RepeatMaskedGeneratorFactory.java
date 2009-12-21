package edu.mit.csail.cgs.ewok.verbs;

import edu.mit.csail.cgs.datasets.general.NamedTypedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.RegionExpanderFactory;

public class RepeatMaskedGeneratorFactory implements RegionExpanderFactory<NamedTypedRegion> {
    private String type;
    
    public RepeatMaskedGeneratorFactory() {
        type = "RepeatMasked";
    }

    public void setType(String t) {type = t;}
    public String getType() {return type;}
    public String getProduct() {return "NamedTypedRegion";}
    public Expander<Region, NamedTypedRegion> getExpander(Genome g) {
        return getExpander(g,type);
    }

    public Expander<Region, NamedTypedRegion> getExpander(Genome g, String type) {
        return new RepeatMaskedGenerator(g);
    }
}

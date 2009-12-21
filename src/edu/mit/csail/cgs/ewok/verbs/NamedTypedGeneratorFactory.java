/*
 * Created on Sep 28, 2006
 */
package edu.mit.csail.cgs.ewok.verbs;

import edu.mit.csail.cgs.datasets.general.NamedTypedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.RegionExpanderFactory;

public class NamedTypedGeneratorFactory implements RegionExpanderFactory<NamedTypedRegion> {
    String type;

    public NamedTypedGeneratorFactory() {
        type = "NamedTypedRegion";
    }
    public NamedTypedGeneratorFactory(String t) {
        type = t;
    }

    public void setType(String t) {type = t;}
    public String getType() {return type;}
    public String getProduct() {return "NamedTypedRegion";}
    public Expander<Region, NamedTypedRegion> getExpander(Genome g) {
        return getExpander(g,type);
    }

    public Expander<Region, NamedTypedRegion> getExpander(Genome g, String type) {
        if (type == null) {
            throw new NullPointerException("NamedTypedGenerator must have a type");
        } else {
            return new NamedTypedGenerator(g, type);
        }
    }
}

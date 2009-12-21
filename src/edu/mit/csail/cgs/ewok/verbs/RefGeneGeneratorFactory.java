/*
 * Created on Sep 28, 2006
 */
package edu.mit.csail.cgs.ewok.verbs;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.GeneFactory;
import edu.mit.csail.cgs.ewok.RegionExpanderFactory;

/**
 * @author tdanford
 */
public class RefGeneGeneratorFactory implements RegionExpanderFactory<Gene>, GeneFactory {
    private String type;

    public RefGeneGeneratorFactory() {
    }

    public void setType(String t) {type = t;}
    public String getType() {return type;}
    public String getProduct() {return "Gene";}
    public Expander<Region, Gene> getExpander(Genome g) {
        return getExpander(g,type);
    }

    public Expander<Region, Gene> getExpander(Genome g, String type) {
        if (type == null) {
            return new RefGeneGenerator(g);
        } else {
            RefGeneGenerator gg = new RefGeneGenerator(g, type);
            return gg;
        }
    }
}

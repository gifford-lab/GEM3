package edu.mit.csail.cgs.ewok;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.verbs.Expander;

public interface RegionExpanderFactory<PRODUCT> {
    public void setType(String type);    
    public String getType();
    public String getProduct();
    public Expander<Region,PRODUCT> getExpander(Genome g);
    public Expander<Region,PRODUCT> getExpander(Genome g, String type);
}

package edu.mit.csail.cgs.ewok;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.verbs.Mapper;

public interface RegionMapperFactory<PRODUCT> {
    public void setType(String type);    
    public String getType();
    public String getProduct();
    public Mapper<Region,PRODUCT> getMapper(Genome g);
    public Mapper<Region,PRODUCT> getMapper(Genome g, String type);
}

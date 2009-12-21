/*
 * Created on Sep 28, 2006
 */
package edu.mit.csail.cgs.ewok.verbs;

import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.RegionExpanderFactory;

public class NamedGeneratorFactory implements RegionExpanderFactory<NamedRegion> {
    
    private String tableName;
    
    public NamedGeneratorFactory() {
        tableName = "sgdOther";
    }
    
    public NamedGeneratorFactory(String t) { 
        tableName = t;
    }

    public String getProduct() {return "NamedRegion";}

    public Expander<Region, NamedRegion> getExpander(Genome g, String table) {
        if (table == null) {
            throw new NullPointerException("NamedGenerator must have a tablename");
        } else {
            return new NamedRegionGenerator<Region>(g, table);
        }
    }

    public Expander<Region, NamedRegion> getExpander(Genome g) {
        return getExpander(g, tableName);
    }

    public String getType() {
        return tableName;
    }

    public void setType(String type) {
        tableName = type;
    }
}

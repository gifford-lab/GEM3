package edu.mit.csail.cgs.ewok.verbs.chippet;

import java.sql.SQLException;

import edu.mit.csail.cgs.datasets.chippet.ChipPetDatum;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.ScoredRegion;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.RegionExpanderFactory;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.utils.Interval;

public class ChipPetExpanderFactory implements RegionExpanderFactory<ChipPetDatum> {
    
    private String name;
    
    public ChipPetExpanderFactory() {
        name = "";
    }

    public void setType(String t) {name = t;}
    
    public String getType() {return name;}
    
    public String getProduct() {return "ChipPet";}
    
    public Expander<Region,ChipPetDatum> getExpander(Genome g) {
        return getExpander(g, name);
    }

    public Expander<Region,ChipPetDatum> getExpander(Genome g, String name) {
        try {
            return new ChipPetExpander(name);
        } catch (SQLException e) {
            e.printStackTrace();
            throw new IllegalArgumentException();
        }
    }

}

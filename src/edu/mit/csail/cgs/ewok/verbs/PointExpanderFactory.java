package edu.mit.csail.cgs.ewok.verbs;

import java.sql.SQLException;
import java.util.Iterator;

import edu.mit.csail.cgs.datasets.chipseq.ChipSeqHit;
import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.ScoredRegion;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.RegionExpanderFactory;
import edu.mit.csail.cgs.utils.Closeable;
import edu.mit.csail.cgs.utils.Interval;

public class PointExpanderFactory implements RegionExpanderFactory<Point> {
    
    private String name;
    
    public PointExpanderFactory() {
        name = "";
    }

    public void setType(String t) {name = t;}
    
    public String getType() {return name;}
    
    public String getProduct() {return "Point";}
    
    public Expander<Region,Point> getExpander(Genome g) {
        return getExpander(g, name);
    }

    public Expander<Region,Point> getExpander(Genome g, String name) {
        try {
            if(name.startsWith("albert")) { 
                return new H2AZMidpointExpander(exptName, alignName);
            } else { 
                throw new IllegalArgumentException();
            }
        } catch (SQLException e) {
            e.printStackTrace();
            throw new IllegalArgumentException();
        }
    }
    public static String exptName = "Albert H2A.Z ChIP-Seq";
    public static String alignName = "0.90 Albert Alignment";    
}

package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;
import java.sql.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.database.*;

/**
 *  Maps a region to a set of scored regions that describe the
 *  aggregate MultiZ alignment score from the UCSC annotations.
*/

public class MultiZGenerator<X extends Region> extends ScoredRegionGenerator {
    
    private static Map<String,String> genomeToTables;
    
    static { 
        genomeToTables = new HashMap<String,String>();
        ResourceBundle res = ResourceBundle.getBundle("edu.mit.csail.cgs.ewok.multiz");
        Enumeration<String> keys = res.getKeys();
        while(keys.hasMoreElements()) { 
            String key = keys.nextElement();
            String value = res.getString(key);
            genomeToTables.put(key, value);
            System.err.println("G2T " + key + " , " + value);
        }
    }
    
    public MultiZGenerator(Genome g) {
        super(g,null);
        String tablename;
        String version = g.getVersion();
        
        if(genomeToTables.containsKey(g.getVersion())) { 
            tablename = genomeToTables.get(g.getVersion());
        } else {
            throw new DatabaseException("Don't know the name of the multiz table for " + version);
        }
        setTablename(tablename);
    }
}


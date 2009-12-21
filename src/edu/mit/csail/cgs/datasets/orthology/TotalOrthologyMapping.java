/*
 * Created on Oct 20, 2005
 */
package edu.mit.csail.cgs.datasets.orthology;

import java.util.*;
import java.io.*;
import java.sql.*;

import edu.mit.csail.cgs.utils.*;

/**
 * @author tdanford
 */
public class TotalOrthologyMapping {
    
    private OrthologyMapping mapping;
    private Set<OrthologyPair> pairs;

    public TotalOrthologyMapping(OrthologyMapping m) {
        mapping = m;
        pairs = new HashSet<OrthologyPair>();
    }

    public OrthologyMapping getMapping() { return mapping; }
    public Collection<OrthologyPair> getPairs() { return pairs; }
    public int size() { return pairs.size(); }
    
    public void addPair(OrthologyPair op) { 
        if(!op.getMapping().equals(mapping)) { throw new IllegalArgumentException("Mapping objects don't match."); }
        if(!pairs.contains(op)) {
            pairs.add(op);
        }
    }
    
    public void insertIntoDB(PreparedStatement mapInsert, PreparedStatement pairInsert, 
            int nextMapID, int nextPairID) throws SQLException { 

        mapping.insertIntoDB(mapInsert, nextMapID);
        for(OrthologyPair p : pairs) { 
            p.insertIntoDB(pairInsert, nextPairID++);
        }
    }    
}

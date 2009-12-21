/*
 * Created on Mar 13, 2008
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;
import java.sql.*;

import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.iterators.EmptyIterator;

public class SpottedProbeGenerator implements Expander<Region,SpottedProbe> {
    
    private Genome genome;
    private String table;
    
    public SpottedProbeGenerator(Genome g, String t) { 
        genome = g;
        table = t;
    }
    
    public Iterator<SpottedProbe> execute(Region a) {
        try {
            LinkedList<SpottedProbe> ps = new LinkedList<SpottedProbe>();
            java.sql.Connection cxn = genome.getUcscConnection();
            Statement s = cxn.createStatement();
            
            String condition = "where chrom='%s' and ((chromStart <= %d and chromEnd >= %d) or " +
                    "(chromStart >= %d and chromStart <= %d))";
            condition = String.format(condition, a.getChrom(), a.getStart(), a.getStart(), a.getStart(), a.getEnd());
            String query = 
                "select chrom, chromStart, chromEnd, name, binding from spottedArray " + condition;

            ResultSet rs = s.executeQuery(query);
            
            while(rs.next()) { 
                SpottedProbe sp = new SpottedProbe(genome, rs);
                ps.addLast(sp);
            }
            
            rs.close();
            s.close();
            DatabaseFactory.freeConnection(cxn);
            
            return ps.iterator();
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return new EmptyIterator<SpottedProbe>();
    }

}

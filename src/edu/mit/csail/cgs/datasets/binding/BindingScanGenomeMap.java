/*
 * Created on Nov 3, 2006
 */
package edu.mit.csail.cgs.datasets.binding;

import java.util.*;
import java.io.*;
import java.sql.*;

import edu.mit.csail.cgs.utils.database.*;

/**
 * @author tdanford
 */
public class BindingScanGenomeMap {
    
    private Map<Integer,Set<Integer>> scan2Genomes, genome2Scans;

    public BindingScanGenomeMap(java.sql.Connection c) throws SQLException {
        scan2Genomes = new HashMap<Integer,Set<Integer>>();
        genome2Scans = new HashMap<Integer,Set<Integer>>();
        
        Statement s = c.createStatement();
        ResultSet rs = s.executeQuery("select scan, genome from bindingScanToGenome");
        
        while(rs.next()) { 
            int scan = rs.getInt(1);
            int genome = rs.getInt(2);
            addMapping(scan, genome);
        }
        
        rs.close();
        s.close();
    }
    
    public boolean containsScan(int scan) { return scan2Genomes.containsKey(scan); }
    public boolean containsGenome(int genome) { return genome2Scans.containsKey(genome); }
    
    public Set<Integer> getAllScans() { return scan2Genomes.keySet(); }
    public Set<Integer> getAllGenomes() { return genome2Scans.keySet(); }

    public Set<Integer> getGenomes(int scan) { return scan2Genomes.get(scan); }
    public Set<Integer> getScans(int genome) { return genome2Scans.get(genome); }
    
    // intended only for use from the BindingScanLoader.deleteScan() method.
    void removeScan(int dbid) { 
        if(scan2Genomes.containsKey(dbid)) {
            Set<Integer> empty = new HashSet<Integer>();
            for(int g : scan2Genomes.get(dbid)) { 
                genome2Scans.get(g).remove(dbid);
                if(genome2Scans.get(g).isEmpty()) { empty.add(g); }
            }
            for(int g : empty) { genome2Scans.remove(g); }
            scan2Genomes.remove(dbid); 
        }
    }
    
    private void addMapping(int scan, int genome) { 
        if(!scan2Genomes.containsKey(scan)) { scan2Genomes.put(scan, new HashSet<Integer>()); }
        if(!genome2Scans.containsKey(genome)) { genome2Scans.put(genome, new HashSet<Integer>()); }
        genome2Scans.get(genome).add(scan);
        scan2Genomes.get(scan).add(genome);        
    }
    
    public void insertNewMapping(java.sql.Connection c, int scan, int genome) throws SQLException {
        if(scan2Genomes.containsKey(scan) && scan2Genomes.get(scan).contains(genome)) { return; }
        Statement s = c.createStatement();
        String insert = "insert into bindingscanToGenome (scan, genome) values (" + scan + "," + genome + ")";
        s.executeUpdate(insert);
        s.close();
        
        addMapping(scan, genome);
    }
}

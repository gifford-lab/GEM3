/*
 * Created on Sep 12, 2006
 */
package edu.mit.csail.cgs.datasets.chippet;

import java.util.*;
import java.sql.*;

import edu.mit.csail.cgs.utils.database.*;

/**
 * @author tdanford
 */
public class ChipPetExptToGenome {
    
    private Map<Integer,Set<Integer>> expt2Genomes, genome2Expts;

    public ChipPetExptToGenome(java.sql.Connection c) throws SQLException {
        expt2Genomes = new HashMap<Integer,Set<Integer>>();
        genome2Expts = new HashMap<Integer,Set<Integer>>();
        
        Statement s = c.createStatement();
        ResultSet rs = s.executeQuery("select expt, genome from chippetToGenome");
        
        while(rs.next()) { 
            int expt = rs.getInt(1), genome = rs.getInt(2);
            addPairing(expt, genome);
            System.out.println(expt + " --> " + genome);
        }
        
        rs.close();
        s.close();
    }
    
    public void addPairing(int expt, int genome) { 
        if(!expt2Genomes.containsKey(expt)) { expt2Genomes.put(expt, new HashSet<Integer>()); }
        if(!genome2Expts.containsKey(genome)) { genome2Expts.put(genome, new HashSet<Integer>()); }
        expt2Genomes.get(expt).add(genome);
        genome2Expts.get(genome).add(expt);
    }
    
    public Set<Integer> getGenomes(int expt) { return expt2Genomes.get(expt); }
    public Set<Integer> getExpts(int genome) { return genome2Expts.get(genome); }

}

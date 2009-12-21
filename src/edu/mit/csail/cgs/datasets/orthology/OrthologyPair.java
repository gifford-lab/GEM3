/*
 * Created on Oct 19, 2005
 */
package edu.mit.csail.cgs.datasets.orthology;

import java.util.*;
import java.io.*;
import java.sql.*;

import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;

/**
 * @author tdanford
 */
public class OrthologyPair {

    /**
     * orth_pair
     * 
     * Name                                      Null?    Type
     * ----------------------------------------- -------- ----------------------------
     * ID                                        NOT NULL NUMBER(10)
     * MAPPING                                   NOT NULL NUMBER(10)
     * NAME1                                     NOT NULL VARCHAR2(100)
     * GENOME1                                   NOT NULL NUMBER(10)
     * NAME2                                     NOT NULL VARCHAR2(100)
     * GENOME2                                   NOT NULL NUMBER(10)
     * 
     */
    
    private int dbid;
    private OrthologyMapping mapping;
    private String name1, name2;
    private Genome genome1, genome2;
    
    public OrthologyPair(OrthologyMapping m, String n1, Genome g1, String n2, Genome g2) { 
        dbid = -1;
        mapping = m;
        name1 = n1; genome1 = g1;
        name2 = n2; genome2 = g2;
    }
    
    public OrthologyPair(ResultSet rs, OrthologyLoader loader) throws NotFoundException {
        try { 
            dbid = rs.getInt(1);
            int mappingID = rs.getInt(2);
            mapping = loader.loadMapping(mappingID);
            
            name1 = rs.getString(3);
            int gid1 = rs.getInt(4);
            
            name2 = rs.getString(5);
            int gid2 = rs.getInt(6);
            
            genome1 = Organism.findGenome(gid1);
            genome2 = Organism.findGenome(gid2);
        } catch(SQLException se) { 
            throw new DatabaseException("Error: " + se.getMessage(), se);
        }
    }
    
    public void insertIntoDB(PreparedStatement ps, int id) throws SQLException {
        if(dbid != -1) { throw new IllegalArgumentException("Already in DB!"); }
        ps.setInt(1, id);
        ps.setInt(2, mapping.getDBID());
        ps.setString(3, name1);
        ps.setInt(4, genome1.getDBID());
        ps.setString(5, name2);
        ps.setInt(6, genome2.getDBID());
        ps.executeUpdate();
        dbid = id;
    }
    
    public void deleteFromDB(PreparedStatement ps) throws SQLException { 
        if(dbid == -1) { throw new IllegalArgumentException("Can't delete non-DB OrthologyPair from DB!"); }
        ps.setInt(1, dbid);
        ps.executeUpdate();
        dbid = -1;
    }
    
    public int getDBID() { return dbid; }
    public OrthologyMapping getMapping() { return mapping; }
    public String getName1() { return name1; }
    public String getName2() { return name2; }
    public Genome getGenome1() { return genome1; }
    public Genome getGenome2() { return genome2; }
    
    public String toString() { 
        return name1 + " (" + genome1.getVersion() + ") --> " + name2 + " (" + genome2.getVersion() + ")";
    }
    
    public int hashCode() { 
        int code = 17;
        code += dbid; code *= 37;
        code += mapping.hashCode(); code *= 37;
        code += name1.hashCode(); code *= 37; 
        code += name2.hashCode(); code *= 37;
        code += genome1.getVersion().hashCode(); code *= 37;
        code += genome2.getVersion().hashCode(); code *= 37;
        return code;
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof OrthologyPair)) { return false; }
        OrthologyPair p = (OrthologyPair)o;
        if(dbid != p.dbid) { return false; }
        if(!name1.equals(p.name1)) { return false; }
        if(!name2.equals(p.name2)) { return false; }
        if(!genome1.getVersion().equals(p.genome1.getVersion())) { return false; }
        if(!genome2.getVersion().equals(p.genome2.getVersion())) { return false; }
        return true;
    }
    
    public static PreparedStatement prepareDeleteStatement(Connection cxn) throws SQLException { 
        return cxn.prepareStatement("delete from orth_pair where id=?");
    }
    
    public static PreparedStatement prepareInsertStatement(Connection cxn) throws SQLException { 
        return cxn.prepareStatement("insert into orth_pair (id, mapping, name1, genome1, name2, genome2) values " +
                "(?, ?, ?, ?, ?, ?)");
    }

    public static PreparedStatement prepareLoadByID(Connection cxn) throws SQLException { 
        return cxn.prepareStatement("select id, mapping, name1, genome1, name2, genome2 from orth_pair " +
                "where id=?");
    }
    
    public static PreparedStatement prepareLoadByMapping(Connection cxn) throws SQLException { 
        return cxn.prepareStatement("select id, mapping, name1, genome1, name2, genome2 from orth_pair " +
                "where mapping=?");
    }
    
    public static PreparedStatement prepareLoadByFirstNameGenome(Connection cxn) throws SQLException { 
        return cxn.prepareStatement("select id, mapping, name1, genome1, name2, genome2 from orth_pair " +
                "where name1=? and genome1=?");
    }

    public static PreparedStatement prepareLoadBySecondNameGenome(Connection cxn) throws SQLException { 
        return cxn.prepareStatement("select id, mapping, name1, genome1, name2, genome2 from orth_pair " +
                "where name2=? and genome2=?");
    }

    public static PreparedStatement prepareLoadByFirstNameGenomeWithMapping(Connection cxn) throws SQLException { 
        return cxn.prepareStatement("select id, mapping, name1, genome1, name2, genome2 from orth_pair " +
                "where name1=? and genome1=? and mapping=?");
    }

    public static PreparedStatement prepareLoadBySecondNameGenomeWithMapping(Connection cxn) throws SQLException { 
        return cxn.prepareStatement("select id, mapping, name1, genome1, name2, genome2 from orth_pair " +
                "where name2=? and genome2=? and mapping=?");
    }
}

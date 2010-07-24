/*
 * Created on Sep 8, 2006
 */
package edu.mit.csail.cgs.datasets.chipseq;

import java.util.*;
import java.io.*;
import java.sql.*;

import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.*;
import java.sql.Connection;

/**
 * @author tdanford
 * 
 * Represents a row in the chipseqexpts table, from the *chipseq schema.
 * 
 * create sequence chipseqexpt_id;
 * 
 * create table chipseqexpts (
 * 
 * 	id number(10) unique not null,
 * 	name varchar2(100),
 * 	replicate varchar2(100),
 *  species number(10) constraint cse_species not null,
 * 	readlength number(5) constraint ces_readlen not null,
 *  cells number(10) constraint cse_cellsone not null,
 *  condition number(10) constraint cse_conditionone not null,
 *  factor number(10) constraint cse_factorone not null,
 *  
 * 	constraint chipseqexpts_pk primary key (name, replicate)
 * );
 */
public class ChipSeqExpt {
    
    private int dbid, readlength;
    private String name, replicate;
    private Factor factor;
    private Condition condition;
    private Cells cells;
    private Organism species;
    
    public ChipSeqExpt(ResultSet rs, ChipSeqLoader loader) throws SQLException { 
        dbid = rs.getInt(1);
        name = rs.getString(2);
        replicate = rs.getString(3);
        int speciesID = rs.getInt(4);
        try {
            species = Organism.getOrganism(speciesID);
        } catch (NotFoundException e) {
            species = null;
            e.printStackTrace();
        }
        readlength = rs.getInt(5);
        int cellsID = rs.getInt(6); 
        int condID = rs.getInt(7);
        int factorID = rs.getInt(8);
        
        MetadataLoader mloader = loader == null ? new MetadataLoader() : loader.getMetadataLoader();

        factor = mloader.loadFactor(factorID);
        cells = mloader.loadCells(cellsID);
        condition = mloader.loadCondition(condID);
    }
    
    public int getDBID() { return dbid; }
    public String getName() { return name; }
    public String getReplicate() { return replicate; }
    public Organism getOrganism() { return species; }
    public Factor getFactor() { return factor; }
    public Cells getCells() { return cells; }
    public Condition getCondition() { return condition; }
    public int getReadLength() {return readlength;}
    
    public String toString() { 
    	return String.format("%s (%s)", name, replicate);
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof ChipSeqExpt)) { return false; }
        ChipSeqExpt c = (ChipSeqExpt)o;
        if(dbid != c.dbid) { return false; }
        return true;
    }
    
    public int hashCode() { 
        int code = 17;
        code += dbid; code *= 37;
        return code;
    }
    
    public static PreparedStatement createLoadAll(java.sql.Connection c) throws SQLException { 
        return c.prepareStatement("select id, name, replicate, species, readlength, cells, condition, factor from chipseqexpts");
    }
    
    public static PreparedStatement createLoadByDBID(java.sql.Connection c) throws SQLException { 
        return c.prepareStatement("select id, name, replicate, species, readlength, cells, condition, factor from chipseqexpts where id=?");
    }
    
    public static PreparedStatement createLoadByName(java.sql.Connection c) throws SQLException { 
        return c.prepareStatement(
        		"select id, name, replicate, species, readlength, cells, condition, factor " +
        		"from chipseqexpts where name=?");
    }
    
    public static PreparedStatement createLoadByNameReplicate(java.sql.Connection c) throws SQLException { 
        return c.prepareStatement(
        		"select id, name, replicate, species, readlength, cells, condition, factor " +
        		"from chipseqexpts where name=? and replicate=?");
    }
    
    public static PreparedStatement createInsert(java.sql.Connection c) throws SQLException { 
    	String query = String.format(
                "insert into chipseqexpts (id, name, replicate, species, readlength, cells, condition, factor) " +
                "values (%s, ?, ?, ?, ?, ?, ?, ?)",
                edu.mit.csail.cgs.utils.database.Sequence.getInsertSQL(c, "chipseqexpt_id"));
    	return c.prepareStatement(query);
    }
}

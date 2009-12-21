/*
 * Created on Sep 8, 2006
 */
package edu.mit.csail.cgs.datasets.chippet;

import java.util.*;
import java.io.*;
import java.sql.*;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.*;
import java.sql.Connection;

/**
 * @author tdanford
 */
public class ChipPetExpt {
    
    public static Vector<ChipPetExpt> loadAllExpts(java.sql.Connection c) throws SQLException { 
        Vector<ChipPetExpt> expts = new Vector<ChipPetExpt>();
        PreparedStatement ps = createLoadAll(c);
        ResultSet rs = ps.executeQuery();
        while(rs.next()) { 
            int id = rs.getInt(1);
            ChipPetExpt expt = new ChipPetExpt(id, rs.getString(2));
            expts.add(expt);
        }
        rs.close();
        ps.close();
        return expts;
    }
    
    public static ChipPetExpt loadExpt(String name) throws SQLException { 
        try {
            java.sql.Connection cxn = DatabaseFactory.getConnection("chippet");
            ChipPetExpt expt = null;
            PreparedStatement ps = ChipPetExpt.createLoadByName(cxn);
            expt = new ChipPetExpt(name, ps);
            if(expt.getDBID() == -1) { expt = null; }
            ps.close();
            DatabaseFactory.freeConnection(cxn);
            return expt;
        } catch (UnknownRoleException e) {
            e.printStackTrace();
            throw new IllegalArgumentException();
        }
    }
    
    private int dbid;
    private String name;
    
    public ChipPetExpt(ResultSet rs) throws SQLException { 
        dbid = rs.getInt(1);
        name = rs.getString(2);
    }
    
    public ChipPetExpt(int id, String n) { 
        dbid = id;
        name = n;
    }
    
    public ChipPetExpt(String n) { 
        dbid = -1;
        name = n;
    }
    
    public ChipPetExpt(String name, PreparedStatement ps) throws SQLException { 
        ps.setString(1, name);
        ResultSet rs = ps.executeQuery();
        dbid = -1;
        if(rs.next()) { 
            dbid = rs.getInt(1);
            name = rs.getString(2);
        }
        rs.close();        
    }

    public ChipPetExpt(int id, PreparedStatement ps) throws SQLException { 
        ps.setInt(1, id);
        ResultSet rs = ps.executeQuery();
        if(rs.next()) { 
            dbid = id;
            name = rs.getString(2);
        }
        rs.close();
    }

    public int getDBID() { return dbid; }
    public String getName() { return name; }
    
    public String toString() { 
        return "(" + dbid + ") " + name;
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof ChipPetExpt)) { return false; }
        ChipPetExpt c = (ChipPetExpt)o;
        if(dbid != c.dbid) { return false; }
        if(!name.equals(c.name)) { return false; }
        return true;
    }
    
    public int hashCode() { 
        int code = 17;
        code += dbid; code *= 37;
        code += name.hashCode(); code *= 37;
        return code;
    }
    
    public static PreparedStatement createLoadAll(java.sql.Connection c) throws SQLException { 
        return c.prepareStatement("select id, name from chippetexpt");
    }
    
    public static PreparedStatement createLoadByDBID(java.sql.Connection c) throws SQLException { 
        return c.prepareStatement("select id, name from chippetexpt where id=?");
    }
    
    public static PreparedStatement createLoadByName(java.sql.Connection c) throws SQLException { 
        return c.prepareStatement("select id, name from chippetexpt where name=?");
    }
    
    public static PreparedStatement createInsert(java.sql.Connection c) throws SQLException { 
        return c.prepareStatement("insert into chippetexpt (id, name) values (?, ?)");
    }
}

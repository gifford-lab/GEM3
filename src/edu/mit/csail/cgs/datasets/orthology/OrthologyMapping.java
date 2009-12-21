/*
 * Created on Oct 19, 2005
 */
package edu.mit.csail.cgs.datasets.orthology;

import java.util.*;
import java.io.*;
import java.sql.*;

import edu.mit.csail.cgs.utils.*;

/**
 * @author tdanford
 */
public class OrthologyMapping {
    
    /**
     * orth_mapping 
     * 
     * Name                                      Null?    Type
     * ----------------------------------------- -------- ----------------------------
     * ID                                        NOT NULL NUMBER(10)
     * NAME                                      NOT NULL VARCHAR2(1000)
     * VERSION                                   NOT NULL VARCHAR2(1000)
     */
    
    private int id;
    private String name, version;
    
    public OrthologyMapping(String n, String v) { 
        id = -1;
        name = n; version = v;
    }

    public OrthologyMapping(ResultSet rs) throws SQLException {
        id = rs.getInt(1);
        name = rs.getString(2);
        version = rs.getString(3);
    }
    
    public void insertIntoDB(PreparedStatement ps, int dbid) throws SQLException { 
        if(id != -1) { throw new IllegalArgumentException("Already in DB!"); }
        ps.setInt(1, dbid);
        ps.setString(2, name);
        ps.setString(3, version);
        ps.executeUpdate();
        id = dbid;
    }
    
    public void deleteFromDB(PreparedStatement ps) throws SQLException {
        if(id == -1) { throw new IllegalArgumentException("Can't delete a non-DB OrthologyMapping from DB."); }
        ps.setInt(1, id);
        ps.executeUpdate();
        id = -1;
    }
    
    public String getName() { return name; }
    public String getVersion() { return version; }
    public int getDBID() { return id; }
    
    public int hashCode() { 
        int code = 17;
        code += id; code *= 37;
        code += name.hashCode(); code *= 37;
        code += version.hashCode(); code *= 37;
        return code;
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof OrthologyMapping)) { return false; }
        OrthologyMapping m = (OrthologyMapping)o;
        if(id != m.id) { return false; }
        if(!name.equals(m.name) || !version.equals(m.version)) { return false; }
        return true;
    }
    
    public String toString() { return "Orthology Mapping: " + name + " (" + version + ")"; }
    
    public static PreparedStatement preparePairDeleteStatement(Connection cxn) throws SQLException { 
        return cxn.prepareStatement("delete from orth_pair where mapping=?");
    }
    
    public static PreparedStatement prepareDeleteStatement(Connection cxn) throws SQLException { 
        return cxn.prepareStatement("delete from orth_mapping where id=?");
    }
    
    public static PreparedStatement prepareLoadAllStatement(Connection cxn) throws SQLException { 
        return cxn.prepareStatement("select id, name, version from orth_mapping");
    }
    
    public static PreparedStatement prepareInsertStatement(Connection cxn) throws SQLException { 
        return cxn.prepareStatement("insert into orth_mapping (id, name, version) values " +
                "(?, ?, ?)");
    }

    public static PreparedStatement prepareLoadByID(Connection cxn) throws SQLException { 
        PreparedStatement ps = cxn.prepareStatement("select id, name, version from orth_mapping where id=?");
        return ps;
    }

    public static PreparedStatement prepareLoadByNameVersion(Connection cxn) throws SQLException { 
        PreparedStatement ps = cxn.prepareStatement("select id, name, version from orth_mapping where name=? and " +
                "version=?");
        return ps;
    }

}

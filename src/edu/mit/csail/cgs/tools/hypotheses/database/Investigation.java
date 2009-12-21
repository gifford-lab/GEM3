/*
 * Created on Jul 16, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.tools.hypotheses.database;

import java.sql.*;

public class Investigation {

    private int dbid;
    private String name;
    
    public Investigation(ResultSet rs) throws SQLException { 
        dbid = rs.getInt(1);
        name = rs.getString(2);
    }
    
    public int getDBID() { return dbid; }
    public String getName() { return name; }
    
    public int hashCode() { 
        int code = 17;
        code += name.hashCode(); code *= 37;
        return code;
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof Investigation)) { return false; }
        Investigation in = (Investigation)o;
        if(!(name.equals(in.name))) { return false; }
        return true;
    }
    
    public String toString() { return "[Investigation: \"" + name + "\"]"; }
}

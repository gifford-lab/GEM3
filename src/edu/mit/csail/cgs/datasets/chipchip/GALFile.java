/**
 * 
 */
package edu.mit.csail.cgs.datasets.chipchip;

import java.util.*;
import java.io.*;
import java.sql.*;

import edu.mit.csail.cgs.datasets.locators.ExptLocator;
import edu.mit.csail.cgs.utils.database.*;

/**
 * @author tdanford
 * 
 * Each object of the GALFile class describes a row (entry) in the "galfiles"
 * table in the DB.
 */
public class GALFile {
    
    public static void main(String[] args) { 
        try {
            Collection<GALFile> values = loadDBGALFiles();
            for(GALFile v : values) { System.out.println(v.toString()); }
            
        } catch (SQLException e) {
            e.printStackTrace();
        } catch (UnknownRoleException e) {
            e.printStackTrace();
        }        
    }
	
	public static Collection<GALFile> loadDBGALFiles()
		throws SQLException, UnknownRoleException { 
		
        LinkedList<GALFile> cells = new LinkedList<GALFile>();
		java.sql.Connection c = 
			DatabaseFactory.getConnection(ExptLocator.dbRole);
		Statement s = c.createStatement();
		ResultSet rs = s.executeQuery("select id, name from galfiles");
		
		while(rs.next()) { 
			GALFile cell = new GALFile(rs);
			cells.addLast(cell);
		}
		
		rs.close();
		s.close();
		DatabaseFactory.freeConnection(c);
		return cells;
	}

	private int dbid;
	private String name;
    
    public GALFile(ResultSet rs) throws SQLException { 
        dbid = rs.getInt(1);
        name = rs.getString(2);
    }
	
	public String getName() { return name; }
	public int getDBID() { return dbid; }
    
    public String toString() { return name + " (#" + dbid + ")"; }
	
	public int hashCode() { 
		int code = 17;
		code += dbid; code *= 37;
		return code;
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof GALFile)) { return false; }
		GALFile c = (GALFile)o;
		if(dbid != c.dbid) { return false; }
		return true;
	}
}

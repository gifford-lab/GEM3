/**
 * 
 */
package edu.mit.csail.cgs.datasets.function;

import java.io.*;
import java.util.*;
import java.sql.*;

import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.Sequence;

/**
 * @author Timothy Danford
 */
public class FunctionVersion {

	/**
	 * func_version:
	 * id number(10)
	 * name varchar(500)
	 */
	
	private int id;
	private String name;
	
	public FunctionVersion(ResultSet rs) throws SQLException { 
		id = rs.getInt(1);
		name = rs.getString(2);
	}
    
    FunctionVersion(String n) { 
        id = -1;
        name = n;
    }
	
	public int getID() { return id; }
	public String getName() { return name; }
    
	public boolean equals(Object o) { 
		if(!(o instanceof FunctionVersion)) { return false; }
		FunctionVersion fv = (FunctionVersion)o;
		return id == fv.id && name.equals(fv.name);
	}
	
	public int hashCode() { 
		int code = 17;
		code += id; code *= 37;
        code += name.hashCode(); code *= 37;
		return code;
	}
	
	public String toString() { 
		return "Function Version: (#" + id + ") " + name;
	}
    
    public static PreparedStatement prepareInsertStatement(java.sql.Connection cxn) 
        throws SQLException { 
        return cxn.prepareStatement("insert into func_version values (" + 
        		Sequence.getInsertSQL(cxn, "id") + ", ?)"); 
    }
	
	public static PreparedStatement prepareReadAll(Connection cxn) throws SQLException { 
		return cxn.prepareStatement("select id, name from func_version");
	}
	
	public static PreparedStatement prepareReadByID(Connection cxn) throws SQLException { 
		return cxn.prepareStatement("select id, name from func_version where id=?");
	}

	public static PreparedStatement prepareReadByName(Connection cxn) throws SQLException { 
		return cxn.prepareStatement("select id, name from func_version where name=?");
	}
}

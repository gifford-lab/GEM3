/**
 * 
 */
package edu.mit.csail.cgs.datasets.function;

import java.io.*;
import java.util.*;
import java.sql.*;

import edu.mit.csail.cgs.utils.*;

/**
 * @author Timothy Danford
 *
 */
public class Assignment {
	
	/**
	 * func_assignment:
	 * object varchar(100)
	 * category number(10)
	 */
	
	private String object;
	private Category category;
	
	public Assignment(ResultSet rs, DatabaseFunctionLoader loader) throws SQLException { 
		object = rs.getString(1);
		int cid = rs.getInt(2);
		category = loader.getCategory(cid);
	}
	public Assignment(ResultSet rs, GOFunctionLoader loader) throws SQLException { 
		object = rs.getString(1);
		int cid = rs.getInt(2);
		category = loader.getCategory(cid);
	}
    
    Assignment(String o, Category c) { 
        object = o;
        category = c;
    }
	
	public Category getCategory() { return category; }
	public String getObject() { return object; }
	
	public String toString() { 
		return "\"" + object + "\" --> " + category.toString();
	}
	
	public int hashCode() { 
		int code = 17;
		code += object.hashCode(); code *= 37;
		code += category.hashCode(); code *= 37;
		return code;
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof Assignment)) { return false; }
		Assignment a = (Assignment)o;
		if(!category.equals(a.category)) { return false; }
		if(!object.equals(a.object)) { return false; }
		return true;
	}
    
    public static PreparedStatement prepareInsert(Connection cxn) throws SQLException { 
        return cxn.prepareStatement("insert into func_assignment (object, category) " +
        		"values (?, ?)");
    }

	public static PreparedStatement createReadByCategory(Connection cxn) throws SQLException { 
		return cxn.prepareStatement("select object, category from func_assignment " +
				"where category=?");
	}
	
    public static PreparedStatement createReadIDByObjectVersion(Connection cxn) throws SQLException { 
		return cxn.prepareStatement("select fa.object, fa.category from func_assignment fa, " +
		"func_category fc where fa.object=? and fa.category=fc.id and fc.version=?");
	}
    
    public static PreparedStatement createReadByVersion(Connection cxn) throws SQLException { 
		return cxn.prepareStatement("select fa.object, fa.category from func_assignment fa, " +
			"func_category fc where fa.category=fc.id and fc.version=?");    	
    }
}

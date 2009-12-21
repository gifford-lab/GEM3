/**
 * 
 */
package edu.mit.csail.cgs.datasets.function;

import java.io.*;
import java.util.*;
import java.sql.*;

import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.Sequence;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

/**
 * @author Timothy Danford
 *
 */
public class Category implements Comparable<Category> {
	
	/**
	 * func_category:
	 * id number(10)
	 * version number(10)
	 * name varchar(100)
     * description varchar(1000)
	 */
	
	/**
	 * func_subcategory:
	 * child_id number(10)
	 * parent_id number(10)
     * version number(10)
	 */

	private int id;
	private FunctionVersion version;
	private String name, description;
	private SortedSet<Category> parents;
    
	public Category(ResultSet rs, DatabaseFunctionLoader loader) throws SQLException { 
		id = rs.getInt(1);
		int versionID = rs.getInt(2);
		version = loader.getVersion(versionID);
		name = rs.getString(3);
        description = rs.getString(4);
		parents = new TreeSet<Category>();		
	}
	public Category(ResultSet rs, GOFunctionLoader loader) throws SQLException { 
		id = rs.getInt(1);
        version = null;
		name = rs.getString(2);
        description = rs.getString(2);
		parents = new TreeSet<Category>();		
	}
    
    Category(FunctionVersion v, String n, String d) {  
        id = -1;
        version = v;
        name = n;
        description = d;
        parents = new TreeSet<Category>();
    }
    Category(FunctionVersion v, String n, String d, int id) {
        this.id = id;
        version = v;
        name = n;
        description = d;
        parents = new TreeSet<Category>();
    }
	
    public Collection<Category> getAllDescendants(Collection<Category> total) { 
        HashSet<Category> cs = new HashSet<Category>();
        for(Category c : total) { 
            if(c.hasAncestor(this)) { 
                cs.add(c);
            }
        }
        return cs;
    }
    
	void addParents(ResultSet parentRS, DatabaseFunctionLoader loader) throws SQLException { 
		while(parentRS.next()) { 
			int child = parentRS.getInt(1);
			int parent = parentRS.getInt(2);
			if(child != id) { throw new IllegalArgumentException(child + " != " + id); }
			parents.add(loader.getCategory(parent));
		}		
	}

	void addParent(Category c) { 
		parents.add(c);
	}
	
    public boolean containsAParent(Set<Category> c) { 
    	for(Category p : parents) { 
    		if(p.equals(c)) { return true; }
    	}
    	return false;
    }
    
    public boolean isChildCategory(Category parent) { 
        for(Category p : parents) { 
            if(p.equals(parent)) { return true; }
            if(p.isChildCategory(parent)) { return true; }
        }
        return false;
    }
    
	public int getID() { return id; }
	public String getName() { return name; }
    public String getDescription() { return description; }
	public FunctionVersion getVersion() { return version; }
	public Collection<Category> getParents() { return parents; }
    
    public Collection<Category> getAncestors() { 
        Set<Category> ancestors = new HashSet<Category>();
        ancestors.addAll(parents);
        for(Category p : parents) { ancestors.addAll(p.getAncestors()); }
        return ancestors;
    }
    
    public boolean hasAncestor(Category a) { 
        if(equals(a)) { return true; }
        if(parents.contains(a)) { return true; }
        for(Category p : parents) { if(p.hasAncestor(a)) { return true; } }
        return false;
    }
	
	public boolean equals(Object o) { 
		if(!(o instanceof Category)) { return false; }
		Category c = (Category)o;
		return id == c.id && c.name.equals(name);
	}
	
	public int hashCode() { 
		int code = 17;
		code += id; code *= 37;
        code += name.hashCode(); code *= 37;
		return code;
	}
	
	public String toString() { 
        StringBuilder sb = new StringBuilder();
        sb.append(name);
        sb.append(" [");
        for(Category p : parents) { 
            sb.append(" " + p.getName());
        }
        sb.append(" ] --> \"" + description + "\"");
        return sb.toString();
	}
    
	public static PreparedStatement prepareInsert(java.sql.Connection cxn) throws SQLException { 
	    return cxn.prepareStatement("insert into func_category (id, version, name, description)" +
	    		" values (" + Sequence.getInsertSQL(cxn, "id")+ ", ?, ?, ?)");
	}
    
    public static PreparedStatement prepareInsertParents(Connection cxn) 
        throws SQLException { 
        return cxn.prepareStatement("insert into func_subcategory (child_id, parent_id) values (?, ?)");
    }
	
	public static PreparedStatement prepareReadByID(Connection cxn) throws SQLException { 
		return cxn.prepareStatement("select id, version, name, description from func_category " +
				"where id=?");
	}
	
	public static PreparedStatement prepareReadIDByVersion(Connection cxn) throws SQLException { 
		return cxn.prepareStatement("select id, version, name, description from func_category " +
				"where version=?");
	}
	
	public static PreparedStatement prepareReadByNameVersion(Connection cxn) throws SQLException { 
		return cxn.prepareStatement("select id, version, name, description from func_category " +
				"where name=? and version=?");
	}
	
	public static PreparedStatement prepareReadParentsByID(Connection cxn) throws SQLException { 
		return cxn.prepareStatement("select child_id, parent_id from func_subcategory" +
				" where child_id = ?");
	}

    public static PreparedStatement prepareReadChildrenByID(Connection cxn) throws SQLException { 
        return cxn.prepareStatement("select child_id, parent_id from func_subcategory" +
                " where parent_id = ?");
    }

    /* (non-Javadoc)
     * @see java.lang.Comparable#compareTo(java.lang.Object)
     */
    public int compareTo(Category c) {
        return name.compareTo(c.name);
    }
}

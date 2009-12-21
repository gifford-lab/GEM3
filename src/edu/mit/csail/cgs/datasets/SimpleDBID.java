/*
 * Author: tdanford
 * Date: Aug 25, 2008
 */
package edu.mit.csail.cgs.datasets;

public class SimpleDBID implements DBID { 
	private int dbid;
	
	public SimpleDBID(int id) { dbid = id; }
	
	public String toString() { return String.format("DBID:%d", dbid); }
	
	public int getID() { return dbid; }
	
	public boolean equals(Object o) { 
		if(!(o instanceof SimpleDBID)) { return false; }
		return dbid == ((SimpleDBID)o).dbid;
	}
	
	public int hashCode() { 
		int code = 17; code += dbid; code *= 37;
		return code;
	}
}


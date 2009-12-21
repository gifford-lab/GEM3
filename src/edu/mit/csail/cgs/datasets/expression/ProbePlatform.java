package edu.mit.csail.cgs.datasets.expression;

import java.sql.ResultSet;
import java.sql.PreparedStatement;
import java.sql.SQLException;

import java.sql.Connection;
import edu.mit.csail.cgs.utils.database.Sequence;

public class ProbePlatform {

    public static final int UNKNOWN_TYPE = -1;
    public static final int PROBE_TYPE = 0;
    public static final int GENE_TYPE = 1;

	private int dbid;
	private String name;
	private int type;

	public ProbePlatform(ResultSet rs) throws SQLException { 
		dbid = rs.getInt(1);
		name = rs.getString(2);
		type = rs.getInt(3);
	}
	
	public int getDBID() { return dbid; }
	public String getName() { return name; }
	public int getType() { return type; }
	
	public boolean equals(Object o) { 
		if(!(o instanceof ProbePlatform)) { return false; }
		ProbePlatform p = (ProbePlatform)o;
		if(dbid != p.dbid) { return false; }
		return true;
	}
	
	public int hashCode() { 
		int code = 17;
		code += dbid; code *= 37;
		return code;
	}
	
	public String toString() { return "Probe Platform \"" + name + "\""; }
	
	public static PreparedStatement prepareLoadByID(Connection cxn) throws SQLException { 
		String query = "select id, name, type from probe_platform where id=?";
		return cxn.prepareStatement(query);
	}
	
	public static PreparedStatement prepareInsert(Connection cxn) throws SQLException { 
        String nextID = Sequence.getInsertSQL(cxn, "probe_platform_id");
		String insert = "insert into probe_platform (id, name, type) values (" + nextID + ", ?, ?)";
		return cxn.prepareStatement(insert);
	}
}

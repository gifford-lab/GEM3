package edu.mit.csail.cgs.datasets.alignments;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

import edu.mit.csail.cgs.utils.database.Sequence;

public class InsertableAlignmentVersion { 
	
	public Integer id;
	public String name;
	
	public InsertableAlignmentVersion(String n) { 
		id = null;
		name = n;
	}
	
	public void setPreparedStatement(PreparedStatement ps) throws SQLException { 
		ps.setString(1, name);
	}
	
	public static PreparedStatement prepareStatement(Connection cxn) throws SQLException { 
		String insert = String.format("insert into alignment_version (id, name) values (%s, ?)", 
				Sequence.getInsertSQL(cxn, "alignment_version_id"));
		return cxn.prepareStatement(insert);
	}

	public int getInsertedID(Connection cxn, Statement s) throws SQLException { 
		ResultSet rs = s.executeQuery(Sequence.getLastSQLStatement(cxn, "alignment_version_id"));
		rs.next();
		id = rs.getInt(1);
		rs.close();
		return id;
	}

	public boolean check() {
		System.out.println("Version Name: " + name);
		return name != null && name.length() > 0 && id == null;
	}
}
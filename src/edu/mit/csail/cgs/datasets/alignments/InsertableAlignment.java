package edu.mit.csail.cgs.datasets.alignments;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

import edu.mit.csail.cgs.utils.database.Sequence;

public class InsertableAlignment { 
	
	public Integer id;
	public String params;
	public InsertableAlignmentVersion version;
	public Double score;
	
	public InsertableAlignment(String p, InsertableAlignmentVersion av, double sc) { 
		id = null;
		params = p;
		version = av;
		score = sc;
	}
	
	public void setPreparedStatement(PreparedStatement ps) throws SQLException { 
		ps.setString(1, params);
		ps.setInt(2, version.id);
		ps.setDouble(3, score);
	}

	public static PreparedStatement prepareStatement(Connection cxn) throws SQLException { 
		String insert = String.format("insert into alignment (id, params, version, score) " +
				"values (%s, ?, ?, ?)", 
				Sequence.getInsertSQL(cxn, "alignment_id"));
		return cxn.prepareStatement(insert);
	}
	
	public int getInsertedID(Connection cxn, Statement s) throws SQLException { 
		ResultSet rs = s.executeQuery(Sequence.getLastSQLStatement(cxn, "alignment_id"));
		rs.next();
		id = rs.getInt(1);
		rs.close();
		return id;
	}

	public boolean check() {
		System.out.println(String.format("Alignment Params: \"%s\"", params));
		System.out.println(String.format("Alignment Score: \"%.3f\"", score));
		return id == null && params != null && score != null;
	}
}


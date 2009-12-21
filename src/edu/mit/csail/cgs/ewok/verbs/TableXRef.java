/**
 */
package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;
import java.sql.*;
import java.io.*;

import edu.mit.csail.cgs.utils.iterators.EmptyIterator;

/**
 * @author tdanford
 */
public class TableXRef implements Expander<String,String>, edu.mit.csail.cgs.utils.Closeable {
	
	private Connection cxn;
	private PreparedStatement ps;
	private String tableName, matchField, returnField;
	
	public TableXRef(Connection c, String t, String match, String ret) 
		throws SQLException {
		
		cxn = c;
		tableName = t;
		matchField = match;
		returnField = ret;
		
		ps = cxn.prepareStatement("select " + returnField + " from " + tableName + " " +
				"where " + matchField + "=?");
	}

	public Iterator<String> execute(String input) {
		try {
			ps.setString(1, input);
			ResultSet rs = ps.executeQuery();
			LinkedList<String> results = new LinkedList<String>();
			while(rs.next()) { 
				results.add(rs.getString(1));
			}
			rs.close();
			return results.iterator();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		return new EmptyIterator<String>();
	}

	public void close() {
		try {
			ps.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		ps = null;
		cxn = null;
	}

	public boolean isClosed() {
		return cxn==null;
	}
}

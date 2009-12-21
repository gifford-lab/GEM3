package edu.mit.csail.cgs.datasets.expression;

import java.sql.*;
import java.util.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.utils.Closeable;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

public class ExpressionMetadataLoader implements Closeable {

	private MetadataLoader loader;
	private boolean shouldCloseLoader;
	
	private java.sql.Connection cxn;
	
	public ExpressionMetadataLoader() throws SQLException { 
		cxn = DatabaseFactory.getConnection("expression");
		loader = new MetadataLoader();
		shouldCloseLoader = true;
	}
	
	public ExpressionMetadataLoader(MetadataLoader l) throws SQLException { 
		cxn = DatabaseFactory.getConnection("expression");
		loader = l;
		shouldCloseLoader = false;
	}

	public void close() {
		if(shouldCloseLoader && !loader.isClosed()) { 
			loader.close();
		}
		DatabaseFactory.freeConnection(cxn);
		cxn = null;
		loader = null;
	}

	public boolean isClosed() {
		return loader == null || loader.isClosed();
	}
	
	public Collection<Cells> getAllCells() throws SQLException { 
		LinkedList<Integer> cells = new LinkedList<Integer>();
		Statement s = cxn.createStatement();
		ResultSet rs = s.executeQuery("select distinct(e.cells) from experiment e");
		
		while(rs.next()) { 
			cells.add(rs.getInt(1));
		}
		
		rs.close();
		s.close();
		
		return loader.loadAllCells(cells);
	}

	public Collection<Condition> getAllConditions() throws SQLException { 
		LinkedList<Integer> conds = new LinkedList<Integer>();
		Statement s = cxn.createStatement();
		ResultSet rs = s.executeQuery("select distinct(e.condition) from experiment e");
		
		while(rs.next()) { 
			conds.add(rs.getInt(1));
		}
		
		rs.close();
		s.close();
		
		return loader.loadAllConditions(conds);
	}
}

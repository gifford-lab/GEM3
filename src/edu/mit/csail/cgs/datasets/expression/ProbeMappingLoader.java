package edu.mit.csail.cgs.datasets.expression;

import java.sql.*;
import java.util.*;

import edu.mit.csail.cgs.utils.Closeable;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.database.Sequence;

public class ProbeMappingLoader implements Closeable {
	
	private static String role = "expression";

	private java.sql.Connection cxn;
	private ExpressionLoader exprLoader;
	
	private PreparedStatement loadMappingByID, loadMappingPairsByID;
	private PreparedStatement insertMapping, insertMappingPair;
	
	public ProbeMappingLoader(ExpressionLoader el) throws SQLException {
		exprLoader = el;
		
		try { 
			cxn = DatabaseFactory.getConnection(role);
		} catch(UnknownRoleException uke) { 
			throw new IllegalArgumentException("Unknown Role: " + role, uke);
		}
		
		loadMappingByID = ProbeMapping.prepareLoadByID(cxn);
		loadMappingPairsByID = ProbeMapping.prepareLoadPairsByID(cxn);
		
		insertMapping = ProbeMapping.prepareInsert(cxn);
		insertMappingPair = ProbeMapping.prepareInsertPair(cxn);
	}

	public void close() {
		exprLoader = null;
		
		try {
			loadMappingByID.close(); loadMappingByID = null;
			loadMappingPairsByID.close(); loadMappingPairsByID = null;
			
			insertMapping.close(); insertMapping = null;
			insertMappingPair.close(); insertMappingPair = null;
		} catch (SQLException e) {
			e.printStackTrace();
		}
		
		DatabaseFactory.freeConnection(cxn);
		
		cxn = null;
	}

	public boolean isClosed() {
		return cxn == null;
	}
	
	public Collection<ProbeMapping> loadAllProbeMappings() throws SQLException {
		LinkedList<ProbeMapping> mappings = new LinkedList<ProbeMapping>();
		Statement s = cxn.createStatement();
		ResultSet rs = s.executeQuery("select id, name from probe_mapping order by name");
		while(rs.next()) { 
			ProbeMapping pm = new ProbeMapping(rs, exprLoader);
			mappings.addLast(pm);
		}
		rs.close();
		s.close();
		return mappings;
	}
	
	public ProbeMapping loadProbeMapping(int dbid) throws SQLException { 
        ProbeMapping pm;
        synchronized (loadMappingByID) {
            loadMappingByID.setInt(1, dbid);
            loadMappingPairsByID.setInt(1, dbid);
            
            ResultSet rs = loadMappingByID.executeQuery();
            ResultSet prs = loadMappingPairsByID.executeQuery();
            
            pm = new ProbeMapping(rs, prs, exprLoader);
            
            rs.close();
            prs.close();
        }
		
		return pm;
	}
	
	public void fillProbeMapping(ProbeMapping pm) throws SQLException { 
        synchronized (loadMappingPairsByID) {
            loadMappingPairsByID.setInt(1, pm.getDBID());
            ResultSet prs = loadMappingPairsByID.executeQuery();
            pm.fillMapping(prs, exprLoader);
            prs.close();
        }
	}
	
	private int getLastProbeMappingID() throws SQLException { 
		int id = -1;
		Statement s = cxn.createStatement();
		ResultSet rs = s.executeQuery(Sequence.getLastSQLStatement(cxn, "probe_mapoing_id"));
		if(rs.next()) { 
			id = rs.getInt(1);
		} else { 
			throw new IllegalArgumentException("No probe_mapping_id sequence defined.");
		}
		rs.close();
		s.close();
		return id;
	}
	
	public int insertProbeMapping(String name, ProbePlatform from, ProbePlatform to) throws SQLException { 
        int id;
        synchronized (insertMapping) {
            insertMapping.setString(1, name);
            insertMapping.setInt(2, from.getDBID());
            insertMapping.setInt(3, to.getDBID());
            
            insertMapping.executeUpdate();
            id = getLastProbeMappingID();
        }
		return id;
	}
	
	public void insertProbeMappingPair(ProbeMapping mapping, Probe from, Probe to) throws SQLException { 
        synchronized (insertMappingPair) {
            insertMappingPair.setInt(1, mapping.getDBID());
            insertMappingPair.setInt(2, from.getDBID());
            insertMappingPair.setInt(3, to.getDBID());
            
            insertMappingPair.executeUpdate();
            
            mapping.addMapping(from, to);
        }
	}
}

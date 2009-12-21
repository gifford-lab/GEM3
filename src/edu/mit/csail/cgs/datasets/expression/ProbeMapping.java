/*
 * Created on Mar 14, 2007
 */
package edu.mit.csail.cgs.datasets.expression;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import java.sql.Connection;
import edu.mit.csail.cgs.utils.database.Sequence;

public class ProbeMapping {
	
	private int dbid;
	private String name;
	private Map<Probe,Set<Probe>> forwardMapping;
	
	private ProbePlatform fromPlatform, toPlatform;
	
	public ProbeMapping(ResultSet rs, ResultSet prs, ExpressionLoader loader) 
		throws SQLException {
		
		dbid = rs.getInt(1);
		name = rs.getString(2);
		
		int fromID = rs.getInt(3);
		int toID = rs.getInt(4);
		
		fromPlatform = loader.loadPlatform(fromID);
		toPlatform = loader.loadPlatform(toID);
		
		forwardMapping = new HashMap<Probe,Set<Probe>>();
		fillMapping(prs, loader);
	}
	
	public ProbeMapping(ResultSet rs, ExpressionLoader loader) throws SQLException { 
		dbid = rs.getInt(1);
		name = rs.getString(2);

		int fromID = rs.getInt(3);
		int toID = rs.getInt(4);
		
		fromPlatform = loader.loadPlatform(fromID);
		toPlatform = loader.loadPlatform(toID);
		
		forwardMapping = new HashMap<Probe,Set<Probe>>();		
	}
	
	void fillMapping(ResultSet prs, ExpressionLoader loader) throws SQLException { 
		while(prs.next()) { 
			int fromID = prs.getInt(1);
			int toID = prs.getInt(2);
			Probe from = loader.loadProbe(fromID);
			Probe to = loader.loadProbe(toID);
		}		
	}
	
	void addMapping(Probe f, Probe t) { 
		if(!forwardMapping.containsKey(f)) { forwardMapping.put(f, new HashSet<Probe>()); }
		forwardMapping.get(f).add(t);
	}
	
	public int getDBID() { return dbid; }
	public String getName() { return name; }
	public ProbePlatform getFromPlatform() { return fromPlatform; }
	public ProbePlatform getToPlatform() { return toPlatform; }
	
	public int size() { return forwardMapping.size(); }
	
	public boolean hasMapFrom(Probe f) { return forwardMapping.containsKey(f); }
	public Set<Probe> getMappingFrom(Probe f) { return forwardMapping.get(f); }
	
	public String toString() { 
		return "Probe Mapping \"" + name + "\"";
	}
	
	public int hashCode() { 
		int code = 17;
		code += dbid; code *= 37;
		return code;
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof ProbeMapping)) { return false; }
		ProbeMapping pm = (ProbeMapping)o;
		if(dbid != pm.dbid) { return false; }
		return true;
	}
	
	public static PreparedStatement prepareLoadByID(Connection cxn) throws SQLException { 
		String query = "select id, name, from_platform, to_platform from probe_mapping where id=?";
		return cxn.prepareStatement(query);
	}
	
	public static PreparedStatement prepareLoadPairsByID(Connection cxn) throws SQLException { 
		String query = "select from_probe, to_probe from probe_mapping_pair where mapping=?";
		return cxn.prepareStatement(query);
	}

	public static PreparedStatement prepareInsert(Connection cxn) throws SQLException {
        String nextID = Sequence.getInsertSQL(cxn, "probe_mapping_id");
		String insert = "insert into probe_mapping (id,name,from_platform,to_platform) " +
				"values (" + nextID + ", ?, ?, ?)";
		return cxn.prepareStatement(insert);
	}
	
	public static PreparedStatement prepareInsertPair(Connection cxn) throws SQLException { 
		String insert = "insert into probe_mapping_pair (mapping, from_probe, to_probe) " +
				"values (?, ?, ?)";
		return cxn.prepareStatement(insert);
	}
}

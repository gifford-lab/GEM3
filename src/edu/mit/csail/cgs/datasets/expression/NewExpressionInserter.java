/*
 * Author: tdanford
 * Date: Sep 11, 2008
 */
package edu.mit.csail.cgs.datasets.expression;

import java.util.*;
import java.sql.*;

import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;

import edu.mit.csail.cgs.utils.database.*;

public class NewExpressionInserter {

	
}

class InsertableExpressionDataset { 
	
}

/*
create sequence probe_id;
create table probe (
	id number(10) constraint probe_id unique not null,
	platform number(10) constraint fk_probe_platform references probe_platform(id) not null,
	name varchar2(50) constraint probe_name unique not null
);
 */
class InsertableProbe { 
	public Integer id;
	public InsertableProbePlatform platform;
	public String name;
	
	public InsertableProbe(InsertableProbePlatform p, String n) { 
		id = null;
		platform = p;
		name = n;
	}
	
	public void setPreparedStatement(PreparedStatement ps) throws SQLException { 
		ps.setInt(1, platform.id);
		ps.setString(2, name);
	}
	
	public void loadLastID(PreparedStatement ps) throws SQLException { 
		ResultSet rs = ps.executeQuery();
		id = rs.getInt(1);
		rs.close();		
	}
	
	public static PreparedStatement prepareStatement(Connection cxn) throws SQLException { 
		String insert = "INSERT INTO probe (id, platform, name) VALUES (%s, ?, ?)";
		insert = String.format(insert, Sequence.getInsertSQL(cxn, "probe_id"));
		return cxn.prepareStatement(insert);
	}
}

/*
create sequence probe_platform_id;
create table probe_platform (
	id number(10) constraint pk_platform_id unique not null,
	name varchar2(100) constraint cst_platform_name not null,
	type number(10) constraint cst_platform_type not null
);
 */
class InsertableProbePlatform { 
	public Integer id;
	public String name;
	public Integer type;
	
	public InsertableProbePlatform(String n, Integer t) { name = n; type = t; id = null; }
	
	public void setPreparedStatement(PreparedStatement ps) throws SQLException { 
		ps.setString(1, name);
		ps.setInt(2, type);
	}
	
	public void loadLastID(PreparedStatement ps) throws SQLException { 
		ResultSet rs = ps.executeQuery();
		id = rs.getInt(1);
		rs.close();
	}
	
	public static PreparedStatement prepareStatement(Connection cxn) throws SQLException { 
		String insert = "INSERT INTO probe_platform (id, name, type) VALUES (%s, ?, ?)";
		insert = String.format(insert, Sequence.getInsertSQL(cxn, "probe_platform_id"));
		return cxn.prepareStatement(insert);
	}
}

class InsertableProbePlatformDataset { 
	
	private InsertableProbePlatform platform;
	private Map<String,InsertableProbe> probes; 
	
	public InsertableProbePlatformDataset(InsertableProbePlatform pp, Collection<InsertableProbe> ps) { 
		platform = pp;
		probes = new HashMap<String,InsertableProbe>();
		for(InsertableProbe p : ps) { 
			probes.put(p.name, p);
		}
	}
	
	public void insert(Connection cxn) throws SQLException { 
		PreparedStatement platformInsert = InsertableProbePlatform.prepareStatement(cxn);
		PreparedStatement probeInsert = InsertableProbe.prepareStatement(cxn);
	}
}

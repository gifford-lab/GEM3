package edu.mit.csail.cgs.utils.database;

import java.util.*;
import java.sql.SQLException;

public abstract class RoleConfiguration {

	private UserTokens user;
	private String schema;
	
	public RoleConfiguration(UserTokens ut, String s) { 
		user = ut;
		schema = s;
	}
	
	public RoleConfiguration(Properties p) { 
		user = new UserTokens(p);
		schema = p.getProperty("schema");
	}
	
	public void setUser(UserTokens ut) { user = ut; }
	public void setSchema(String s) { schema = s; }
	
}

class MySQLRoleConfiguration extends RoleConfiguration {
	
	public MySQLRoleConfiguration(Properties props) { 
		super(props);
	}

}

class OracleRoleConfiguration extends RoleConfiguration {
	
	public OracleRoleConfiguration(Properties props) { 
		super(props);
	}

}

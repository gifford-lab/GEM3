package edu.mit.csail.cgs.utils.database;

import java.util.*;

public abstract class RoleConfiguration {

	public RoleConfiguration(UserTokens ut, String s) {
	}
	
	public RoleConfiguration(Properties p) { 
		new UserTokens(p);
		p.getProperty("schema");
	}
	
	public void setUser(UserTokens ut) { }
	public void setSchema(String s) { }
	
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

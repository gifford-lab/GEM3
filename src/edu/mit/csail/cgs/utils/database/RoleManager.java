package edu.mit.csail.cgs.utils.database;

import java.util.*;
import java.sql.*;

public class RoleManager {

	private Map<String,RoleConfiguration> roles;
	
	public RoleManager() { 
		roles = new HashMap<String,RoleConfiguration>();
	}
	
	public void setRoleConfiguration(String role, RoleConfiguration config) { 
		roles.put(role, config);
	}
	
}

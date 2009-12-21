package edu.mit.csail.cgs.utils.database;

import java.util.Properties;

public class UserTokens {

	private String username;
	private String password;
	
	public UserTokens(String u, String p) { 
		username = u;
		password = p;
	}
	
	public UserTokens(Properties props) { 
		username = props.getProperty("user");
		password = props.getProperty("password");
	}
	
	public String getUserName() { return username; }
	public String getPassword() { return password; }
}

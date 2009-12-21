package edu.mit.csail.cgs.utils.database;

import java.sql.Connection;
import java.sql.SQLException;
import java.util.Map;
import java.util.HashMap;
import java.util.Properties;
import java.net.URL;
import java.io.*;
import java.util.regex.*;

/**
 * <code>DatabaseFactory</code> is the class that all code should use to obtain
 * database connections (don't use CxnPool directly).  DatabaseFactory
 * provides database connections for *roles* and manages the connection pool for you.
 * A role is a database resource that you want to access, eg, <b>ucsc_SGDv1</b> for 
 * annotations for SGDv1 or <b>chipchip</b> for ChIP-Chip data.
 *
 * @author <a href="mailto:arolfe@mit.edu">Alex Rolfe</a>
 * @version 1.0
 */
public abstract class DatabaseFactory {
    
    static Pattern oraclePattern = Pattern.compile(".*oracle.*",Pattern.CASE_INSENSITIVE);
    static Pattern mysqlPattern = Pattern.compile(".*mysql.*",Pattern.CASE_INSENSITIVE);
    public static int ORACLE = 1, MYSQL = 2, UNKNOWN;

    /*
     * Static Methods
     */
    private static Map<String,CxnPool>pools = new HashMap<String,CxnPool>();
    private static Map<String,String>defaultUsers = new HashMap<String,String>();
    private static Map<String,String>defaultSchemas = new HashMap<String,String>();
    private static Map<Connection, CxnPool> cxnSource = new HashMap<Connection, CxnPool>();

    /**
     * returns the default read-only connection for this role 
     */
    public static Connection getConnection(String role) throws SQLException, UnknownRoleException {
        Properties props = null;
        try {
            role = getRealRole(role);
            String user, schema;
            if (defaultUsers.containsKey(role) && 
                defaultSchemas.containsKey(role)) {
                user = defaultUsers.get(role);
                schema = defaultSchemas.get(role);
            } else {
                props = getPropertiesForRole(role);
                user = props.getProperty("user");
                schema = props.getProperty("schema");
                defaultUsers.put(role,user);
                defaultSchemas.put(role,schema);
            }
            String key = role + user + schema;
            if (!pools.containsKey(key)) {
                if (props == null) {
                    props = getPropertiesForRole(role);
                }
                addPool(key,props);
            }
            CxnPool pool = pools.get(key);            
            Connection cxn = pool.getConnection();
            cxnSource.put(cxn,pool);
            return cxn;
        } catch (IOException ex) {
            throw new RuntimeException("Couldn't read properties for " + role,ex);
        }
    }

    /**
     * Returns a RW connection for the specified role using the specified username and password 
     */
    public static Connection getConnection(String role, String uname, String passwd) throws SQLException, UnknownRoleException {
        Properties props;
        try {
            role = getRealRole(role);
            props = getPropertiesForRole(role);
            props.setProperty("user",uname);
            props.setProperty("passwd",passwd);
            String key = role + props.getProperty("user") + props.getProperty("schema");
            if (!pools.containsKey(key)) {
                addPool(key,props);
            }
            CxnPool pool = pools.get(key);
            Connection cxn = pool.getConnection();
            cxnSource.put(cxn,pool);
            return cxn;
        } catch (IOException ex) {
            throw new RuntimeException("Couldn't read properties for " + role,ex);
        }
    }

    /** for a given role name, eg chipchpi, look up
     * CHIPCHIPROLE in the environment.  If it's set,
     * use its value as the return value which is the
     * actual role to use.
     *
     * For example, you might set CHIPCHIPROLE=finkchipchip
     * to connect to a chipchip schema that isn't the one
     * in your .chipchip_passwd
     */
    public static String getRealRole(String role) {
        String envkey = role.toUpperCase() + "ROLE";
        String val = System.getenv(envkey);
        if (val == null) {
            return role;
        } else {
            return val;
        }
    }

    public static void addPool(String key, Properties props) throws SQLException{
        String cs = props.getProperty("jdbcconnectstring");
        Matcher o = oraclePattern.matcher(cs);
        Matcher m = mysqlPattern.matcher(cs);
        if (o.matches()) {
            //            System.err.println("Creating oracle Connection pool for " + key);
            pools.put(key,new OracleCxnPool(props));
        } else if (m.matches()) {
            //            System.err.println("Creating mysql Connection pool for " + key);
            pools.put(key,new MySQLCxnPool(props));
        } else {
            throw new SQLException("Unknown database type in " + cs);
        }
    }
    /**
     * Returns the Connection to the appropriate pool
     * or frees it if it didn't come from any CxnPool.
     */
    public static void freeConnection(Connection c) {
        CxnPool p = cxnSource.get(c);
        if (p != null) {
            p.freeConnection(c);
        } else {
            try {
                c.close();
            } catch (SQLException e) {
                e.printStackTrace();
            }
        }
    }

    private static Properties getPropertiesForRole(String role) throws UnknownRoleException, IOException {
        role = role.replaceAll("\\W+","_");
        if (System.getenv(role.toUpperCase() + "ROLE") != null) {
        	String oldrole = role;
            role = System.getenv(role.toUpperCase() + "ROLE");
            //System.err.println(String.format("Found ENV for role %s = \"%s\"", oldrole, role));
        }
        String homedir = System.getenv("HOME");
        //System.err.println(String.format("ENV variable HOME=%s", homedir));
        String basename = role + "_passwd";
        String fname = homedir + "/." + basename;
        File propfile = new File(fname);
        //        System.err.println("Trying to connect as role " + role);
        if (propfile.exists() && propfile.canRead()) {
        	//System.err.println(String.format("Reading HOME/. properties file %s: \"%s\"", propfile.getName(), propfile.getAbsolutePath()));
            return readPasswd(propfile);
        }
        fname = homedir + "/" + basename;
        propfile = new File(fname);
        if (propfile.exists() && propfile.canRead()) {
        	//System.err.println(String.format("Reading HOME properties file %s: \"%s\"", propfile.getName(), propfile.getAbsolutePath()));
            return readPasswd(propfile);
        } else {
        	//System.err.println("Looking in classpath for properties file...");
            try {
                ClassLoader cl = ClassLoader.getSystemClassLoader();
                URL url = cl.getResource(basename);
                if (url != null) {
                	//System.err.println(String.format("Found properties URL: %s", url.toString()));
                    return readPasswdStream(url.openStream());
                }
            } catch (Exception ex) {
                throw new UnknownRoleException("Couldn't find properties file for role " + role,ex);
            }
            throw new UnknownRoleException("Couldn't find properties file for role " + role);
        }         
    }

    private static Properties readPasswd(File propfile) throws IOException {
        return readPasswdStream(new FileInputStream(propfile));
    }

    private static Properties readPasswdStream(InputStream propstream) throws IOException {
        String line;
        Properties props = new Properties();
        BufferedReader reader = new BufferedReader(new InputStreamReader(propstream));        
        
        /**
         * TODO
         * What the hell?  Why is this parsing the thing line-by-line, and not using the PropertyResourceBundle implementation
         * provided by Java itself?  Did *I* do this?  (If so, sorry.)  -Tim
         */
        
        while ((line = reader.readLine()) != null) {
        	try { 
        		int p = line.indexOf('=');
        		if(p < 0) { 
        			continue;
        		}
        		String key = line.substring(0,p);
        		String value = line.substring(p+1);
        		props.setProperty(key,value);
        	} catch(RuntimeException e) { 
        		System.err.println(line);
        		throw e;
        	}
        }
        reader.close();
        return props;
    }
    /**
     * returns the type (ORACLE, MYSQL, UNKNOWN) of this
     * Connection
     */
    public static int getType(Connection c) {
        CxnPool p = cxnSource.get(c);
        if (p == null) {
            return UNKNOWN;
        } else {
            return p.getType();
        }
    }
    /**
     * Returns true iff this is a connection to an Oracle DB
     */
    public static boolean isOracle(Connection c){
        return getType(c) == ORACLE;
    }
    public static boolean isMySQL(Connection c){
        return getType(c) == MYSQL;
    }
    public static boolean isPostGres(Connection c){
        return false;
    }
    public static boolean isSQLLite(Connection c){
        return false;
    }

}


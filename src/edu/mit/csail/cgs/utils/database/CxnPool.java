/*
 * Created on Feb 28, 2006
 */
package edu.mit.csail.cgs.utils.database;

import java.util.*;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.sql.Connection;

/**
 * @author tdanford
 */
public abstract class CxnPool {

    /* availPool is the set of available connections.  fullPool
       is all of the connections we've created.  This freeConnection()
       check that the connection can be returned ot the availPool
       before it does so */
    private LinkedList<Connection> availPool, fullPool;
    private int poolSize;
    protected Properties props;
    
    public CxnPool (Properties p) {
        availPool = new LinkedList<Connection>(); 
        fullPool = new LinkedList<Connection>();
        poolSize = 1;
        props = p;
    }
    
    private Connection createConnection() { 
        try {
            Connection cxn = DriverManager.getConnection(connectString(), username(), password());
            initializeConnection(cxn);
            availPool.add(cxn);
            fullPool.add(cxn);
            return cxn;
        } catch (SQLException ex) {
            System.err.println("Couldn't create a database connection " + ex.toString());
            ex.printStackTrace();
            throw new DatabaseException("Couldn't create database connection",ex);
        }
    }

    protected void initializeConnection(java.sql.Connection cxn) throws SQLException {
    }

    public synchronized Connection getConnection () {
        if (availPool.size() == 0) {
            for (int i = 0; i < poolSize; i++) {
                createConnection();
            }   
        }
        Connection c = null;
        
        /*
        c = availPool.removeFirst();
        availPool.remove(c);
        */
        c = fullPool.getFirst();
        
        return c;
    }

    public synchronized void freeConnection(Connection n) {
        if(fullPool.contains(n)) {
            availPool.add(n);            
        }
    }

    public synchronized void remove(Connection n) {
        (new RuntimeException("removing connection " + n)).printStackTrace();
        fullPool.remove(n);
        availPool.remove(n);
    }

    private static void reportClosed (Connection c) {
        
    }

    public abstract int getType();

    public String connectString() {return props.getProperty("jdbcconnectstring");}
    public String username() {return props.getProperty("user");}
    public String password() {return props.getProperty("passwd");}    
}

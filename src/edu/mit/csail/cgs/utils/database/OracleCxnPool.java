package edu.mit.csail.cgs.utils.database;

import java.util.*;
import java.sql.*;


public class OracleCxnPool extends CxnPool {
    
    public OracleCxnPool(Properties props) {
        super(props);
        try {
            Class.forName("oracle.jdbc.OracleDriver").newInstance();
        } catch (Exception ex) {
            System.err.println("Couldn't load Oracle JDBC driver");
        }
        try {
            DriverManager.registerDriver(new oracle.jdbc.driver.OracleDriver());
        } catch (SQLException ex) {
            throw new DatabaseException("Can't register driver",ex);
        }

    }    

    protected void initializeConnection(java.sql.Connection cxn) throws SQLException {
        super.initializeConnection(cxn);        
        java.sql.Statement s = cxn.createStatement();
        s.execute("alter session set current_schema=" + props.getProperty("schema"));
        s.close();
    }
    public int getType(){ return DatabaseFactory.ORACLE;}
}

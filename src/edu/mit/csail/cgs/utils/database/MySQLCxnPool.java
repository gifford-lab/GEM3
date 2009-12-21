package edu.mit.csail.cgs.utils.database;

import java.util.*;
import java.sql.*;

//         connectString = "jdbc:mysql://opteron.csail.mit.edu/" + dbname + "?user=" + username +
//-             "&password=" + password;



public class MySQLCxnPool extends CxnPool {

    public MySQLCxnPool(Properties props) {
        super(props);
        try {
            Class.forName("com.mysql.jdbc.Driver").newInstance();
        } catch (Exception ex) {
            System.err.println("Couldn't load MySQL JDBC driver");
        }
    }
    public int getType() {return DatabaseFactory.MYSQL;}
}
        

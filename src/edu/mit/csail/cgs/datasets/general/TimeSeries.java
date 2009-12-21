/*
 * Created on Mar 14, 2007
 * 
 * @author Timothy Danford
 */
package edu.mit.csail.cgs.datasets.general;

import java.sql.*;
import edu.mit.csail.cgs.utils.database.Sequence;

/*
create sequence timeseries_id; 

create table timeseries (
    id number(10) constraint timeseries_id unique not null,
    name varchar2(100) constraint timeseries_name not null
);
 */

public class TimeSeries {

    private int dbid;
    private String name;
    
    public TimeSeries(ResultSet rs) throws SQLException { 
        dbid = rs.getInt(1);
        name = rs.getString(2);
    }
    
    public String getName() { return name; }
    public int getDBID() { return dbid; }

    public int hashCode() { 
        int code = 17;
        code += dbid; code *= 37;
        code += name.hashCode(); code *= 37;
        return code;
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof TimeSeries)) { return false; }
        TimeSeries ts = (TimeSeries)o;
        if(dbid != ts.dbid) { return false; }
        if(!name.equals(ts.name)) { return false; }
        return true;
    }
    
    public String toString() { return "Time-Series \"" + name + "\""; }
    
    public static PreparedStatement prepareLoadByID(java.sql.Connection cxn) throws SQLException { 
        String query = "select id, name from timeseries where id=?";
        return cxn.prepareStatement(query);
    }    
    
    public static PreparedStatement prepareInsert(java.sql.Connection cxn) throws SQLException {
        String nextID = Sequence.getInsertSQL(cxn, "timeseries_id");
        String query = "insert into timeseries (id, name) values (" + nextID + ", ?)";
        return cxn.prepareStatement(query);
    }
}

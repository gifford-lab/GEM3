/*
 * Created on Mar 14, 2007
 * 
 * @author Timothy Danford
 */
package edu.mit.csail.cgs.datasets.general;

import java.sql.*;
import edu.mit.csail.cgs.utils.database.Sequence;

/*
create sequence timepoint_id;

create table timepoint (
    id number(10),
    time_series number(10),
    series_order number(10),
    name varchar2(100)
);
 */

public class TimePoint {

    private int dbid;
    private TimeSeries series;
    private int order;
    private String name;
    
    public TimePoint(ResultSet rs, TimeSeriesLoader loader) throws SQLException { 
        dbid = rs.getInt(1);
        series = loader.loadTimeSeries(rs.getInt(2));
        order = rs.getInt(3);
        name = rs.getString(4);
    }
    
    public int getDBID() { return dbid; }
    public TimeSeries getSeries() { return series; }
    public int getOrder() { return order; }
    public String getName() { return name; }
    
    public int hashCode() { 
        int code = 17;
        code += dbid; code *= 37;        
        return code;
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof TimePoint)) { return false; }
        TimePoint tp = (TimePoint)o;
        if(dbid != tp.dbid) { return false; }
        return true;
    }
    
    public String toString() { 
        return "TimePoint #" + order + " in \"" + series.getName() + "\" (" + name + ")"; 
    }
    
    public static PreparedStatement prepareLoadByID(java.sql.Connection cxn) throws SQLException { 
        String query = "select id, time_series, series_order, name from timepoint where id=?"; 
        return cxn.prepareStatement(query);
    }
    
    public static PreparedStatement prepareLoadBySeriesID(java.sql.Connection cxn) throws SQLException { 
        String query = "select id, time_series, series_order, name from timepoint where series=? order by series_order";
        return cxn.prepareStatement(query);
    }
    
    public static PreparedStatement prepareInsert(java.sql.Connection cxn) throws SQLException { 
        String nextID = Sequence.getInsertSQL(cxn, "timepoint_id");
        String query = "insert into timepoint (id, time_series, series_order, name) values (" + nextID + ", ?, ?, ?)";
        return cxn.prepareStatement(query);
    }
}

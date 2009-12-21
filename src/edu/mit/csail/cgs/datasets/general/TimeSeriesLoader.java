/*
 * Created on Mar 15, 2007
 * 
 * @author Timothy Danford
 */
package edu.mit.csail.cgs.datasets.general;

import java.sql.*;
import java.util.*;

import edu.mit.csail.cgs.utils.Closeable;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.database.Sequence;

public class TimeSeriesLoader implements Closeable {
    
    public static final String role = "core";
    
    private java.sql.Connection cxn;
    private PreparedStatement loadSeriesByID, insertSeries;
    private PreparedStatement loadPointByID, loadPointsBySeries, insertPoint;
    
    private Map<Integer,TimeSeries> cachedSeries;
    private Map<Integer,TimePoint> cachedPoints;
    
    public TimeSeriesLoader() throws SQLException {
        try {
            cxn = DatabaseFactory.getConnection(role);
        } catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        }
        
        loadSeriesByID = TimeSeries.prepareLoadByID(cxn);
        loadPointByID = TimePoint.prepareLoadByID(cxn);
        loadPointsBySeries = TimePoint.prepareLoadBySeriesID(cxn);
        
        insertSeries = TimeSeries.prepareInsert(cxn);
        insertPoint = TimePoint.prepareInsert(cxn);
        
        cachedSeries = new HashMap<Integer,TimeSeries>();
        cachedPoints = new HashMap<Integer,TimePoint>();
    }

	public void beginTransaction() throws SQLException { 
		cxn.setAutoCommit(false);
	}

	public void commitTransaction() throws SQLException { 
		cxn.commit();
		cxn.setAutoCommit(true);
	}
    
    public void close() {
        cachedSeries.clear();
        cachedPoints.clear();
        
        try {
            loadSeriesByID.close(); loadSeriesByID = null;
            loadPointByID.close(); loadPointByID = null;
            loadPointsBySeries.close(); loadPointsBySeries = null;
            
            insertSeries.close(); insertSeries = null;
            insertPoint.close(); insertPoint = null;
            
        } catch (SQLException e) {
            e.printStackTrace();
        }
        
        DatabaseFactory.freeConnection(cxn);
        cxn = null;
    }
    
    public boolean isClosed() { return cxn==null; }
    
    public int insertTimePoint(TimeSeries series, String name, int order) throws SQLException {
        synchronized(insertPoint) {
            insertPoint.setInt(1, series.getDBID());
            insertPoint.setInt(2, order);
            insertPoint.setString(3, name);
            
            insertPoint.executeUpdate();
        }
        
        return getLastTimePointID();
    }
    
    public int insertTimeSeries(String name) throws SQLException {
        synchronized(insertSeries) {
            insertSeries.setString(1, name);        
            insertSeries.executeUpdate();
        }
        
        return getLastTimeSeriesID();
    }
    
    private int getLastID(String tableName) throws SQLException { 
        int id = -1;
        Statement s = cxn.createStatement();
        ResultSet rs = s.executeQuery(Sequence.getLastSQLStatement(cxn, tableName));
        if(rs.next()) { 
            id = rs.getInt(1);
        } else { 
            throw new IllegalStateException("Couldn't get last value of " + tableName + "?");
        }
        s.close();
        
        return id;              
    }

	private int getLastTimeSeriesID() throws SQLException { return getLastID("timeseries_id"); }
	private int getLastTimePointID() throws SQLException { return getLastID("timepoint_id"); }

	public Collection<TimeSeries> loadAllTimeSeries() throws SQLException { 
		LinkedList<TimeSeries> series = new LinkedList<TimeSeries>();
        String query = "select id, name from timeseries order by name";
		PreparedStatement ps = cxn.prepareStatement(query);

		ResultSet rs = ps.executeQuery();
		while(rs.next()) { 
			int id = rs.getInt(1);
			if(cachedSeries.containsKey(id)) { 
				series.add(cachedSeries.get(id));
			} else { 
				TimeSeries ts = new TimeSeries(rs);
				cachedSeries.put(id, ts);
				series.add(ts);
			}
		}
		rs.close();

		ps.close();
		return series;
	}

    public TimeSeries loadTimeSeries(int dbid) throws SQLException {
        if(cachedSeries.containsKey(dbid)) { return cachedSeries.get(dbid); }
        
        TimeSeries ts = null;
        synchronized (loadSeriesByID) {
            loadSeriesByID.setInt(1, dbid);
            ResultSet rs = loadSeriesByID.executeQuery();
            if(rs.next()) { ts = new TimeSeries(rs); }
            rs.close();
        }
        
        cachedSeries.put(ts.getDBID(), ts);
        return ts;
    }
    
    public TimePoint loadTimePoint(int dbid) throws SQLException { 
        if(cachedPoints.containsKey(dbid)) { return cachedPoints.get(dbid); }
        
        TimePoint tp = null;
        synchronized (loadPointByID) {
            loadPointByID.setInt(1, dbid);
            ResultSet rs = loadPointByID.executeQuery();
            if(rs.next()) { tp = new TimePoint(rs, this); }
            rs.close();
        }
        cachedPoints.put(tp.getDBID(), tp);
        
        return tp;
    }
    
    public Collection<TimePoint> loadTimePoints(TimeSeries ts) throws SQLException { 
        LinkedList<TimePoint> points = new LinkedList<TimePoint>();
        
        synchronized (loadPointsBySeries) {
            loadPointsBySeries.setInt(1, ts.getDBID());
            ResultSet rs = loadPointsBySeries.executeQuery();
            while(rs.next()) { 
                int dbid = rs.getInt(1);
                TimePoint tp = null;
                if(cachedPoints.containsKey(dbid)) { 
                    tp = cachedPoints.get(dbid);
                } else { 
                    tp = new TimePoint(rs, this);
                    cachedPoints.put(tp.getDBID(), tp);
                }
                points.addLast(tp);
            }
            rs.close();
        }
        
        return points;
    }
}

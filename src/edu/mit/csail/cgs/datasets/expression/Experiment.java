/*
 * Created on Mar 13, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.datasets.expression;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import edu.mit.csail.cgs.datasets.general.Cells;
import edu.mit.csail.cgs.datasets.general.Condition;
import edu.mit.csail.cgs.datasets.general.MetadataLoader;
import edu.mit.csail.cgs.datasets.general.TimePoint;
import edu.mit.csail.cgs.datasets.general.TimeSeriesLoader;

import java.sql.Connection;
import edu.mit.csail.cgs.utils.database.Sequence;

/*
 *  These are copied from the expression.oracle file, for reference. 
 * 

create sequence expr_experiment_id;

create table expr_experiment (
    id number(10) constraint expr_experiment_id unique not null,
    name varchar2(100) constraint expr_experiment_name not null,
    value_type number(5) constraint expr_experiment_value_type not null,
    log_scale number(5) constraint expr_experiment_log_scale not null,
    cells number(10) constraint expr_experiment_cells references cells(id) not null,
    condition number(10) constraint expr_experiment_condition references condition(id) not null,
    timepoint number(10) constraint expr_experiment_timepoint references timepoint(id)
);

create table expr_experiment_params (
    experiment number(10) constraint expr_experiment_params_experiment references expr_experiment(id) not null,
    key varchar2(50) constraint expr_experiment_params_key not null,
    value varchar2(200)
);

 */

public class Experiment {
    
    public final static int UNKNOWN_TYPE = -1;
    public final static int CONTINUOUS_TYPE = 0;
    public final static int DISCRETE_TYPE = 1;

    // The fact that the database's ID field is also a field of this class means that
	// objects of this class should *always* correspond to pre-existing rows in the 
	// expr_experiment table.  The ExpressionLoader class will handle inserting new rows
	// into that table, through a different interface.
    private int dbid;
    
    private String name;
    private int value;
    private boolean logScale;
    private Map<String,String> params;
        
    private Cells cells;
    private Condition condition;
    private ProbePlatform platform;
    private TimePoint timePoint;
    
    public Experiment(ResultSet rs, ResultSet paramsRs, 
    		ExpressionLoader exprLoader,
    		TimeSeriesLoader timeLoader,
    		MetadataLoader chipLoader) 
    	throws SQLException {
    	
    	/*
    	 * Notice that the order in which we read out the values from the two result-set objects
    	 * in this constructor matches the order in which the PreparedStatements (returned by 
    	 * the static methods at the bottom of this class) retrieve the values.  
    	 * 
    	 * The ExpressionLoader class will hold those PreparedStatements.  When asked for 
    	 * an experiment with a given ID (or by some other criteria), it will create the 
    	 * corresponding ResultSet objects and pass them along to this constructor.  
    	 */
    	
        dbid = rs.getInt(1);
        name = rs.getString(2);
        value = rs.getInt(3);
        logScale = rs.getInt(4) == 1;
        
        int cellsID = rs.getInt(5);
        int condID = rs.getInt(6);        
        int timePointID = rs.getInt(7);
        int platformID = rs.getInt(8);
        
        // These next two lines, using the ChipChipLoader object, implicitly do another 
        // table-lookup against the cells and conditions tables 
        cells = chipLoader.loadCells(cellsID);
        condition = chipLoader.loadCondition(condID);
        
        params = new HashMap<String,String>();
        
        while(paramsRs.next()) { 
            String key = paramsRs.getString(1), value = paramsRs.getString(2);
            params.put(key, value);
        }
        
        platform = exprLoader.loadPlatform(platformID);

        if(timePointID != 0) { 
        	timePoint = timeLoader.loadTimePoint(timePointID); 
        } else { 
        	timePoint = null;
        }
    }
    
    public int getDBID() { return dbid; }
    public String getName() { return name; }
    public int getValueType() { return value; }
    public boolean isLogScale() { return logScale; }
    public Cells getCells() { return cells; }
    public Condition getCondition() { return condition; }
    public ProbePlatform getPlatform() { return platform; }
    public TimePoint getTimePoint() { return timePoint; }
    
    public boolean hasParam(String k) { return params.containsKey(k); }
    public Set<String> getParams() { return params.keySet(); }
    public String getParamValue(String k) { return params.get(k); }

    // this method has package visibility, because it will be called appropriately 
    // from ExpressionLoader.
    void setParam(String k, String v) { 
    	params.put(k, v);
    }
    
    public int hashCode() { 
        int hash = 17;
        hash += dbid; hash *= 37;
        return hash;
    }
    
    public boolean equals(Object o) { 
    	if(!(o instanceof Experiment)) { return false; }
    	Experiment e = (Experiment)o;
    	if(dbid != e.dbid) { return false; }
    	return true;
    }
    
    public String toString() { return "Expression Experiment \"" + name + "\""; }
    
    public static PreparedStatement prepareLoadByID(Connection cxn) throws SQLException { 
    	String query = "select id, name, value_type, log_scale, cells, condition, timepoint, platform " +
                "from experiment where id=?";
    	return cxn.prepareStatement(query);
    }
    
    public static PreparedStatement prepareLoadParamsByID(Connection cxn) throws SQLException { 
    	String query = "select key, value from experiment_params where experiment=?";
    	return cxn.prepareStatement(query);
    }
    
    public static PreparedStatement prepareInsert(Connection cxn) throws SQLException {
        String nextID = Sequence.getInsertSQL(cxn, "experiment_id");
        String insert = "insert into experiment " +
        		"(id, name, value_type, log_scale, cells, condition, timepoint, platform) " +
                " values (" + nextID + ", ?, ?, ?, ?, ?, ?, ?)";
        return cxn.prepareStatement(insert);
    }
    
    public static PreparedStatement prepareInsertParam(Connection cxn) throws SQLException { 
        String insert = "insert into experiment_params (experiment, key, value) " +
        		"values (?, ?, ?)";
        return cxn.prepareStatement(insert);
    }
}

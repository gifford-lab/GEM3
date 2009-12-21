/*
 * Created on Mar 14, 2007
 * @author Timothy Danford
 */
package edu.mit.csail.cgs.datasets.expression;

import java.util.*;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import edu.mit.csail.cgs.utils.database.Sequence;
import java.sql.Connection;

/*
create sequence expr_processing_id;                                                    
create table expr_processing (                                                         
    id number(10) constraint cst_expr_processsing_id unique not null,
    type varchar2(100) constraint cst_expr_processing_type not null
);                                                                                     
                                                                                       
create table expr_processing_params (                                                  
    processing number(10) constraint fk_expr_proc_params_processing references expr_processing(id) not null,
    key varchar2(50) constraint cst_expr_proc_params_key not null,                     
    value varchar2(200),                                                               
    constraint expr_proc_params_pk primary key (processing,key)                        
);                                                                                     

create table expr_processing_inputpair (                                               
    processing number(10) constraint fk_expr_inputpair_processing references expr_processing(id) not null,
    input number(10) constraint fk_expr_inputpair_input references expr_experiment(id) not null,     
    output number(10) constraint fk_expr_inputpair_output references exper_experiment(id) not null,  
    constraint expr_proc_inputpair_pk primary key (processing,input,output)
);                                                                                     
 */

public class Processing {

	private int dbid;
	private String type;
	
	private Map<String,String> params;
	private Map<Experiment,Set<Experiment>> output2inputs;
	
	public Processing(ResultSet rs, ResultSet paramrs, ResultSet inputrs, ExpressionLoader loader) 
		throws SQLException { 
		
		dbid = rs.getInt(1);
		type = rs.getString(2);
		
		params = new HashMap<String,String>();
		while(paramrs.next()) { 
			String key = paramrs.getString(1), value = paramrs.getString(2);
			params.put(key, value);
		}
		
		output2inputs = new HashMap<Experiment,Set<Experiment>>();
		while(inputrs.next()) { 
			int inputID = inputrs.getInt(1);
			int outputID = inputrs.getInt(2);
			Experiment input = loader.loadExperiment(inputID);
			Experiment output = loader.loadExperiment(outputID);
			
			addInputPair(input, output);
		}
	}
	
	public int getDBID() { return dbid; }
	public String getType() { return type; }
	
	public Set<Experiment> getInputs(Experiment output) { 
		return output2inputs.get(output);
	}
	
	public boolean hasOutput(Experiment output) { return output2inputs.containsKey(output); }
	public boolean hasInputOutputPair(Experiment input, Experiment output) { 
		return hasOutput(output) && output2inputs.get(output).contains(input);
	}
	
	public Set<String> getParams() { return params.keySet(); }
	public boolean hasParam(String key) { return params.containsKey(key); }
	public String getParamValue(String key) { return params.get(key); } 
	
	// package accessibility, because this will be called from ExpressionLoader 
	// as appropriate.
	void setParam(String k, String v) { params.put(k, v); }
	
	void addInputPair(Experiment input, Experiment output) { 
		if(!(output2inputs.containsKey(output))) { 
			output2inputs.put(output, new HashSet<Experiment>());
		}
		output2inputs.get(output).add(input);
	}
	
	public int hashCode() { 
		int code = 17;
		code += dbid; code *= 37;
		return code;
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof Processing)) { return false; }
		Processing p = (Processing)o;
		if(dbid != p.dbid) { return false; }
		return true;
	}
	
	public static PreparedStatement prepareLoadByID(Connection cxn) throws SQLException { 
		String query = "select id, type from processing where id=?";
		return cxn.prepareStatement(query);
	}
	
	public static PreparedStatement prepareLoadParamsByID(Connection cxn) throws SQLException { 
		String query = "select key, value from processing_params where processing=?";
		return cxn.prepareStatement(query);
	}
	
	public static PreparedStatement prepareLoadInputsByID(Connection cxn) throws SQLException { 
		String query = "select input, output from processing_inputpair where processing=?";
		return cxn.prepareStatement(query);
	}
	
	public static PreparedStatement prepareInsert(Connection cxn) throws SQLException {
        String nextID = Sequence.getInsertSQL(cxn, "processing_id");
		String insert = "insert into processing (id, type) values (" + nextID + ", ?)";
		return cxn.prepareStatement(insert);
	}
	
	public static PreparedStatement prepareInsertParam(Connection cxn) throws SQLException { 
		String insert = "insert into processing_params (processing, key, value) values (?, ?, ?)";
		return cxn.prepareStatement(insert);
	}
	
	public static PreparedStatement prepareInsertInput(Connection cxn) throws SQLException { 
		String insert = "insert into processing_inputpair (processing, input, output) values (?, ?, ?)";
		return cxn.prepareStatement(insert);
	}
}

/*
 * Created on Mar 18, 2007
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
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Set;

import java.sql.Connection;
import edu.mit.csail.cgs.utils.database.Sequence;

/*
create sequence expr_experiment_set_id;
create table expr_experiment_set (
    id number(10) constraint expr_experiment_set_pk unique not null,
    name varchar2(100) constraint cst_expr_experiment_set_name not null,
    type varchar2(100) constraint cst_expr_experiment_set_type not null
);
create index ix_expr_experiment_set_id on expr_experiment_set(id);

create table expr_experiment_set_member (
    experiment_set number(10) constraint fk_expr_experiment_set_member_exptset references expr_experiment_set(id),
    experiment number(10) constraint fk_expr_experiment_set_member_expt references expr_experiment(id),
    constraint pk_expr_experiment_set_member primary key(experiment_set, experiment)
);
create index ix_expr_experiment_set_member on expr_experiment_set_member(experiment_set, experiment);
*/

public class ExperimentSet {

    private int dbid;
    private String name; 
    private int type;
    private Set<Experiment> experiments;
    
    public ExperimentSet(ResultSet rs, ResultSet mrs, ExpressionLoader loader) throws SQLException { 
        
        dbid = rs.getInt(1);
        name = rs.getString(2);
        type = rs.getInt(3);
        
        experiments = new HashSet<Experiment>();
        
        while(mrs.next()) { 
            int exptID = mrs.getInt(1);
            experiments.add(loader.loadExperiment(exptID));
        }
    }
    
    public int getDBID() { return dbid; }
    public String getName() { return name; }
    public int getType() { return type; }
    
    public int size() { return experiments.size(); }
    public boolean containsExperiment(Experiment e) { return experiments.contains(e); }
    public Collection<Experiment> getExperiments() { return new LinkedList<Experiment>(experiments); }
    
    void addMember(Experiment e) { 
        experiments.add(e);
    }
    
    public int hashCode() { 
        int code = 17;
        code += dbid; code *= 37;
        return code;
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof ExperimentSet)) { return false; }
        ExperimentSet es = (ExperimentSet)o;
        if(dbid != es.dbid) { return false; }
        return true;
    }
    
    public String toString() { 
        return "Experiment Set \"" + name + "\" (" + type + ")";
    }
    
    public static PreparedStatement prepareLoadByID(Connection cxn) throws SQLException { 
        String query = "select id, name, type from expr_experiment_set where id=?";
        return cxn.prepareStatement(query);
    }
    
    public static PreparedStatement prepareLoadMembersByID(Connection cxn) throws SQLException { 
        String query = "select experiment_set from expr_experiment_set_member where experiment=?";
        return cxn.prepareStatement(query);
    }
    
    public static PreparedStatement prepareInsert(Connection cxn) throws SQLException { 
        String nextID = Sequence.getInsertSQL(cxn, "experiment_set_id");
        String insert = "insert into experiment_set (id, name, type) values (" + nextID + ", ?, ?)";
        return cxn.prepareStatement(insert);
    }
    
    public static PreparedStatement prepareInsertMember(Connection cxn) throws SQLException { 
        String insert = "insert int experiment_set_member (experiment_set, experiment) values (?, ?)";
        return cxn.prepareStatement(insert);
    }
}

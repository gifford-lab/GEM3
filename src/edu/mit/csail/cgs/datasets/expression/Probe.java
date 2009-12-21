/*
 * Created on Mar 14, 2007
 */
package edu.mit.csail.cgs.datasets.expression;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;

import java.sql.Connection;
import edu.mit.csail.cgs.utils.database.Sequence;

/*
create sequence probe_id;
create table probe (
    id number(10) constraint expr_probe_id unique not null,
    name varchar2(50) constraint expr_probe_name not null
);
 */

public class Probe {

    private int dbid;
    private String name;
    private ProbePlatform platform;
    
    protected Probe(Probe p) { 
    	dbid = p.dbid;
    	name = p.name;
    	platform = p.platform;
    }
    
    public Probe(ResultSet rs, ExpressionLoader loader) throws SQLException { 
        dbid = rs.getInt(1);
        name = rs.getString(2);
        int platID = rs.getInt(3);
        platform = loader.loadPlatform(platID);
    }
    
    public int getDBID() { return dbid; }
    public String getName() { return name; }
    public ProbePlatform getPlatform() { return platform; }
    
    public int hashCode() { 
        int code = 17;
        code += dbid; code *= 37;
        return code;
    }
    
    public String toString() { return "Probe \"" + name + "\""; }
    
    public boolean equals(Object o) { 
        if(!(o instanceof Probe)) { return false; }
        Probe p = (Probe)o;
        if(dbid != p.dbid) { return false; }
        return true;
    }
    
    public static PreparedStatement prepareLoadByID(Connection cxn) throws SQLException { 
        String query = "select id, name, platform from probe where id=?";
        return cxn.prepareStatement(query);
    }
    
    public static PreparedStatement prepareLoadByExperimentID(Connection cxn) throws SQLException { 
        String query = "select ep.id, ep.name, ep.platform from probe ep, measurement em " +
                "where ep.id=em.probe and em.experiment=?";
        return cxn.prepareStatement(query);
    }
    
    public static PreparedStatement prepareInsert(Connection cxn) throws SQLException {
        String nextID = Sequence.getInsertSQL(cxn, "probe_id");
        String insert = "insert into probe (id, name, platform) values (" + nextID + ", ?, ?)";
        return cxn.prepareStatement(insert);
    }
}

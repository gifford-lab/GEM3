/*
 * Created on Mar 15, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.datasets.expression;

import java.sql.*;

/*
create table expr_measurement (
    probe number(10) constraint expr_measurement_probe references probe(id) not null,
    experiment number(10) constraint expr_measurement_experiment references expr_experiment(id) not null,
    value binary_float constraint expr_measurement_value not null
);
 */

public class ExprMeasurement {

    private Probe probe;
    private Experiment expt;
    private double value;
    
    // package-scope to the constructor.
    ExprMeasurement(Experiment e, Probe p, double v) { 
        probe = p;
        expt = e;
        value = v;
    }
    
    protected ExprMeasurement(ExprMeasurement em) { 
    	probe = em.probe;
    	expt = em.expt;
    	value = em.value;
    }
    
    public Probe getProbe() { return probe; }
    public Experiment getExperiment() { return expt; }
    public double getValue() { return value; }
    
    public boolean equals(Object o) { 
        if(!(o instanceof ExprMeasurement)) { return false; }
        ExprMeasurement em = (ExprMeasurement)o;
        if(!probe.equals(em.probe)) { return false; }
        if(!expt.equals(em.expt)) { return false; }
        return true;
    }
    
    public int hashCode() { 
        int code = 17;
        code += probe.hashCode(); code *= 37;
        code += expt.hashCode(); code *= 37;
        return code;
    }
    
    public String toString() { return value + " @ " + probe.getName() + " (" + expt.getName() + ")"; }
    
    public static PreparedStatement prepareLoadByBothIDs(Connection cxn) throws SQLException { 
        String query = "select probe, experiment, value from measurement where experiment=? and probe=?";
        return cxn.prepareStatement(query);
    }
    
    // this one's a little tricky, because we need to initialize *both* an ExprMeasurement and a Probe
    // object from its results.  Notice that the corresponding ResultSet could be used to initialize a Probe 
    // object itself, since it has the probe ID and name in the #1 and #2 locations.
    public static PreparedStatement prepareLoadByExptID(Connection cxn) throws SQLException { 
        String query = "select em.probe, ep.name, ep.platform, em.value from measurement em," +
                " probe ep where em.probe=ep.id and em.experiment=?";
        return cxn.prepareStatement(query);
    }

    public static PreparedStatement prepareLoadByProbeID(Connection cxn) throws SQLException { 
        String query = "select probe, experiment, value from measurement em where probe=?";
        return cxn.prepareStatement(query);
    }
    
    public static PreparedStatement prepareInsert(Connection cxn) throws SQLException { 
        String insert = "insert into measurement (experiment, probe, value) values (?, ?, ?)";
        return cxn.prepareStatement(insert);
    }
}

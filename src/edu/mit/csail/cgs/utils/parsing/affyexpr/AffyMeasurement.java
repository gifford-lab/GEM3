/*
 * Created on Sep 6, 2006
 */
package edu.mit.csail.cgs.utils.parsing.affyexpr;

/**
 * @author tdanford
 */
public class AffyMeasurement {
    
    private AffyProbe probe;
    private double value;
    private String call;

    public AffyMeasurement(AffyProbeSet probeSet, String line) {
        String[] array = line.split("\t");
        String pn = array[0];
        value = Double.parseDouble(array[1]);
        call = array[2].trim();
        probe = probeSet.getAffyProbe(pn);
        //if(probe == null) { throw new IllegalArgumentException(line); }
    }
    
    public boolean isPresent() { return call.equals("P"); }
    public boolean isMarginal() { return call.equals("M"); }
    public boolean isAbsent() { return call.equals("A"); }
    public AffyProbe getProbe() { return probe; }
    public String getCall() { return call; }
    public double getValue() { return value; }
    
    public boolean equals(Object o) { 
        if(!(o instanceof AffyMeasurement)) { return false; }
        AffyMeasurement am = (AffyMeasurement)o;
        return probe.equals(am.probe);
    }
    
    public int hashCode() { return probe.hashCode(); }
    
    public String toString() { 
        String appended = isPresent() ? "(P)" : "   ";
        return appended + " " + probe.toString() + " --> " + value;
    }
}

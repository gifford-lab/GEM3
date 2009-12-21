/*
 * Created on Sep 19, 2005
 */
package edu.mit.csail.cgs.utils.parsing.expression;

/**
 * @author tdanford
 */
public class ExprProbeValue {
    
    private String exptName;
    private String probeName;
    private double value;

    public ExprProbeValue(String e, String n, double v) {
        exptName = e;
        probeName = n; 
        value = v;
    }

    public String getExptName() { return exptName; }
    public String getProbeName() { return probeName; }
    public double getValue() { return value; }
    
    public String toString() { return exptName + "(" + probeName + ":" + value + ")"; }
    
    public int hashCode() { 
        int code = 17;
        code += exptName.hashCode(); code *= 37;
        code += probeName.hashCode(); code *= 37;
        long bits = Double.doubleToLongBits(value);
        code += (int)(bits >> 32);
        code *= 37;
        return code;
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof ExprProbeValue)) { return false; } 
        ExprProbeValue epv = (ExprProbeValue)o;
        
        return epv.exptName.equals(exptName) && 
            epv.probeName.equals(probeName) && 
            epv.value == value;
    }
}

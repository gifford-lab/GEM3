/*
 * Created on Aug 21, 2005
 */
package edu.mit.csail.cgs.utils;

/**
 * @author tdanford
 */
public class NamedValue<Y> { 
    private Y data;
    private String name;
    
    public NamedValue(String n, Y d) { 
        name = n; 
        data = d;
    }
    
    public NamedValue(Y d) { 
        name = d.toString();
        data = d;
    }
    
    public String toString() { 
        return name;
    }
    
    public String getName() { return name; }
    public Y getData() { return data; }
    
    public int hashCode() { 
        int code = 17;
        code += name.hashCode(); code *= 37;
        code += data.hashCode(); code *= 37;
        return code;
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof NamedValue)) { return false; }
        NamedValue nv = (NamedValue)o;
        return name.equals(nv.getName()) && data.equals(nv.getData()); 
    }
}
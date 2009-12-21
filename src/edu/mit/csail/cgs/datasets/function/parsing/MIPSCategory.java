/*
 * Created on Dec 19, 2005
 */
package edu.mit.csail.cgs.datasets.function.parsing;

/**
 * @author tdanford
 */
public class MIPSCategory { 
    private String name, desc;
    
    public MIPSCategory(String inputLine) { 
        int splitIndex = inputLine.indexOf(" ");
        name = inputLine.substring(0, splitIndex);
        desc = inputLine.substring(splitIndex+1, inputLine.length());
    }
    
    public boolean isParent(MIPSCategory cat) { 
        return cat.name.matches(name + "\\.[^\\.]+");
    }
    
    public boolean isChild(MIPSCategory cat) { return cat.isParent(this); }
    
    public String getName() { return name; }
    public String getDesc() { return desc; }
    
    public boolean equals(Object o) { 
        if(!(o instanceof MIPSCategory)) { return false; }
        MIPSCategory cat = (MIPSCategory)o;
        return name.equals(cat.name) && desc.equals(cat.desc);
    }
    
    public int hashCode() { 
        int code = 17;
        code += name.hashCode(); code *= 37;
        code += desc.hashCode(); code *= 37;
        return code;
    }
}

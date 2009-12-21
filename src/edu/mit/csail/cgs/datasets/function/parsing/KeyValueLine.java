/*
 * Created on Nov 19, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.datasets.function.parsing;

public class KeyValueLine {
    
    public String key, value;
    
    public KeyValueLine(String line) { 
        int si = line.indexOf(":");
        if(si == -1) { throw new IllegalArgumentException(line); }
        key = line.substring(0, si).trim();
        value = line.substring(si+1, line.length()).trim();
    }
    
    public String getKey() { return key; }
    public String getValue() { return value; }
}
/*
 * Created on Apr 13, 2007
 */
package edu.mit.csail.cgs.echo;

import java.util.*;

/**
 * @author Timothy Danford
 *
 * Pulls out "key=value" pairs, from a given string, with different pairs
 * split by ";" tokens.  Currently used in the ConstantParser code.
 */
public class OptionsParser {
    
    private Map<String,String> valueMap;

    public OptionsParser(String t) { 
        valueMap = new HashMap<String,String>();
        String[] array = t.split(";");
        for(int i = 0; i < array.length; i++) { 
            if(array[i].indexOf("=") != -1) { 
                int idx = array[i].indexOf("=");
                String key = array[i].substring(0, idx);
                String value = array[i].substring(idx+1, array[i].length());
                valueMap.put(key, value);
            } else {
                String key = array[i];
                valueMap.put(key, "true");
            }
        }
    }
    
    public Map<String,String> getValueMap() { return valueMap; }
}

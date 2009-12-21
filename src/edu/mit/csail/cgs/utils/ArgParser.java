/*
 * Created on Mar 24, 2005
 */
package edu.mit.csail.cgs.utils;

import java.util.*;

/**
 * @author tdanford
 */
public class ArgParser {
    
    private String[] originalArgs;
    private Map<String,String> keyedValues;
    private LinkedList<String> unkeyedValues;
    private LinkedList<String> orderedKeys;
    
    public boolean hasKey(String k) { return keyedValues.containsKey(k); }
    
    public List<String> getUnkeyedValues() { return new LinkedList<String>(unkeyedValues); }
    public Set<String> getKeys() { return keyedValues.keySet(); }
    public String getKeyValue(String k) { return keyedValues.get(k); }
    public String getOriginalArg(int i) { return originalArgs[i]; }
    public List<String> getOrderedKeys() { return new LinkedList<String>(orderedKeys); }

    public int getNumKeyedValues() { return keyedValues.size(); }
    public int getNumUnkeyedValues() { return unkeyedValues.size(); }
    public int getNumOriginalArgs() { return originalArgs.length; }
    
    public ArgParser(String[] args) {
        
        originalArgs = (String[])args.clone();
        SortedSet<Integer> unkeyedInds = new TreeSet<Integer>();
        Map<Integer, Integer> keyMap = new HashMap<Integer, Integer>();
        orderedKeys = new LinkedList<String>();

		for(int i = 0; i < args.length; i++) { unkeyedInds.add(i); }

        for(int i = 0; i < args.length; i++) { 
            if(args[i].startsWith("--")) {
                orderedKeys.addLast(args[i].substring(2, args[i].length()));
                unkeyedInds.remove(i);
                if(i+1 < args.length) { 
                    if(!args[i+1].startsWith("--")) {
                        keyMap.put(i, i+1);
						unkeyedInds.remove(i+1);
                    } else { 
                        keyMap.put(i, -1);
                    }
                } else { 
                    keyMap.put(i, -1);
                }
            }
        }

        unkeyedValues = new LinkedList<String>();
        keyedValues = new HashMap<String,String>();
        
        for(int i : keyMap.keySet()) { 
            String k = args[i].substring(2, args[i].length());
            String v = null;
            if(keyMap.get(i) != -1) { v = args[keyMap.get(i)]; }
            keyedValues.put(k, v);
        }
        
        for(int i : unkeyedInds) { 
            unkeyedValues.addLast(args[i]);
        }
    }
   
    public String toString() { 
        StringBuilder sb = new StringBuilder();
        SortedSet<String> keys = new TreeSet<String>(keyedValues.keySet());
        for(String k : keys) { 
            String v = keyedValues.get(k);
            sb.append(" --" + k);
            if(v != null) { sb.append(" " + v); }
        }
        
        for(String v : unkeyedValues) { 
            sb.append(" " + v);
        }
        return sb.toString();
    }
}

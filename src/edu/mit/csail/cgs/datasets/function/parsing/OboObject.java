/*
 * Created on Nov 19, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.datasets.function.parsing;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.Vector;

public class OboObject { 

    private String termTag;
    private boolean initialized;
    private Vector<KeyValueLine> lines;
    
    public OboObject(BufferedReader br) throws IOException {
        String line = null;
        initialized = false;
        boolean continueReading = true;
        lines = new Vector<KeyValueLine>();
        
        line = br.readLine();
        if(line != null) {
            
            line = line.trim();
            if(line.startsWith("[") && line.endsWith("]")) { 
                termTag = line.substring(1, line.length()-1);
            } else {
                termTag = null;
                lines.add(new KeyValueLine(line));
            }
            initialized = true;
            
            while(continueReading && ((line = br.readLine()) != null)) { 
                line = line.trim();
                if(line.length() == 0) { 
                    continueReading = false; 
                } else { 
                    lines.add(new KeyValueLine(line));
                }
            }
        }
    }
    
    public Map<String,LinkedList<String>> getKeyValueMap() { 
        Map<String,LinkedList<String>> map = new HashMap<String,LinkedList<String>>();
        for(KeyValueLine line : lines) { 
            String key = line.key;
            if(!map.containsKey(key)) { 
                map.put(key, new LinkedList<String>());
            }
            
            LinkedList<String> lst = map.get(key);
            lst.addLast(line.value);
        }
        return map;
    }
    
    public void debugPrint(PrintStream ps) { 
        if(termTag != null) { ps.println("[" + termTag + "]"); }
        for(KeyValueLine l : lines) { 
            ps.println(l.getKey() + ": " + l.getValue());
        }
    }
    
    public String getTermTag() { return termTag; }
    public int size() { return lines.size(); }
    public KeyValueLine getLine(int i) { return lines.get(i); }
    public boolean isInitialized() { return initialized; }
    
    public int countValues(String key) { 
        int c = 0;
        for(KeyValueLine l : lines) { 
            if(l.getKey().equals(key)) { 
                c += 1; 
            }
        }
        return c;
    }
    
    public Collection<String> getValues(String key) { 
        LinkedList<String> values = new LinkedList<String>();
        for(KeyValueLine l : lines) { 
            if(l.getKey().equals(key)) { 
                values.addLast(l.getValue());
            }
        }
        return values;
    }
}
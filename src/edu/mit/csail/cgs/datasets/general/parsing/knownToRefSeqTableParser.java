/*
 * Created on May 24, 2006
 */
package edu.mit.csail.cgs.datasets.general.parsing;

import java.io.*;
import java.sql.SQLException;
import java.util.*;

import edu.mit.csail.cgs.datasets.function.DatabaseFunctionLoader;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

/**
 * @author tdanford
 */
public class knownToRefSeqTableParser {
    
    private static class Entry { 
        private String knownID, refSeqID;
        
        public Entry(String line) { 
            String[] array = line.split("\\s+");
            knownID = array[0];
            refSeqID = array[1];
        }
        
        public String getKnownID() { return knownID; }
        public String getRefSeqID() { return refSeqID; }
    }

    private LinkedList<Entry> entries;

    public knownToRefSeqTableParser(File f) throws IOException { 
        BufferedReader br = new BufferedReader(new FileReader(f));
        entries = new LinkedList<Entry>();
        String line;
        
        while((line = br.readLine()) != null) { 
            line = line.trim();
            if(line.length() > 0) { 
                Entry e = new Entry(line);
                entries.add(e);
            }
        }
        
        br.close();
    }
    
    public Map<String,Set<String>> createKnownToRefSeqMap() { 
        Map<String,Set<String>> map = new HashMap<String,Set<String>>();
        for(Entry e : entries) { 
            String k = e.getKnownID(), r = e.getRefSeqID();
            if(!map.containsKey(k)) { map.put(k, new HashSet<String>()); }
            map.get(k).add(r);
        }
        return map;
    }
    
    public Map<String,Set<String>> createSwissProtToRefSeqMap(kgXrefTableParser xrefParser) { 
        Map<String,Set<String>> map = new HashMap<String,Set<String>>();
        Map<String, ? extends Collection<String>> kg2spMap = xrefParser.build_kgID2spID_map();

        for(Entry e : entries) { 
            String kg = e.getKnownID();
            String refseq = e.getRefSeqID();
            
            if(kg2spMap.containsKey(kg)) { 
                for(String spID : kg2spMap.get(kg)) { 
                    if(!map.containsKey(spID)) { map.put(spID, new HashSet<String>()); }
                    map.get(spID).add(refseq);
                }
            }
        }
        
        return map;
    }
}

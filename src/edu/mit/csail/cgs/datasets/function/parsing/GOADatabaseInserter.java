/*
 * Created on Feb 12, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.datasets.function.parsing;

import java.util.*;
import java.io.*;
import java.sql.*;

import edu.mit.csail.cgs.datasets.function.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.*;

public class GOADatabaseInserter {
    
    private SetTools<String> tools;
    private File oboFile, goaFile, xrefFile;
    
    public GOADatabaseInserter(File obo, File goa, File xref) {
        oboFile = obo;
        goaFile = goa;
        xrefFile = xref;
        tools = new SetTools<String>();
    }
    
    public Map<String,Set<String>> buildAssignmentMap() throws IOException { 
        Map<String,Set<String>> assigns = new HashMap<String,Set<String>>();
        Map<String,Set<String>> revAssigns = new HashMap<String,Set<String>>();
        
        GOAIterator itr = new GOAIterator(goaFile);
        while(itr.hasNext()) { 
            GOALine line = itr.next();
            String go = line.getGOID();
            String id = line.getObjectID();
            String symbol = line.getObjectSymbol();
            
            if(!assigns.containsKey(go)) { assigns.put(go, new HashSet<String>()); }
            assigns.get(go).add(id);
            assigns.get(go).add(symbol);
            
            if(!revAssigns.containsKey(id)) { revAssigns.put(id, new HashSet<String>()); }
            if(!revAssigns.containsKey(symbol)) { revAssigns.put(symbol, new HashSet<String>()); }
            revAssigns.get(id).add(go);
            revAssigns.get(symbol).add(go);
        }
        
        if(xrefFile != null) { 
            GOAXRefIterator xitr = new GOAXRefIterator(xrefFile);
            while(xitr.hasNext()) { 
                GOAXRefLine xref = xitr.next();
                Collection<String> names = xref.getAllNames();
                for(String name : names) { 
                    if(revAssigns.containsKey(name)) { 
                        for(String go : revAssigns.get(name)) { 
                            assigns.get(go).add(name);
                        }
                    }
                }
            }
        }
        
        return assigns;
    }

}

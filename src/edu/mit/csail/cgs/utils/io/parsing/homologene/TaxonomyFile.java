/*
 * Created on Aug 14, 2006
 */
package edu.mit.csail.cgs.utils.io.parsing.homologene;

import java.util.*;
import java.io.*;

/**
 * @author tdanford
 */
public class TaxonomyFile {
    
    private Map<String,String> id2name, name2id; 

    public TaxonomyFile(File f) throws IOException {
        parse(f);
    }

    public TaxonomyFile() throws IOException { 
        parse(new File("taxid_taxname"));
    }
    
    private void parse(File f) throws IOException { 
        id2name = new HashMap<String,String>();
        name2id = new HashMap<String,String>();
        
        String line;
        BufferedReader br = new BufferedReader(new FileReader(f));
        while((line = br.readLine()) != null) { 
            line = line.trim();
            if(line.length() > 0) { 
                String[] array = line.split("\t");
                String id = array[0], name = array[1];
                if(id2name.containsKey(id) || name2id.containsKey(name)) {
                    throw new IllegalArgumentException("Duplicate: " + line);
                }
                
                id2name.put(id, name);
                name2id.put(name, id);
            }
        }
        br.close();
    }
    
    public String getName(String id) { return id2name.get(id); }
    public String getID(String name) { return name2id.get(name); }
    public Set<String> getNames() { return name2id.keySet(); }
    public Set<String> getIDs() { return id2name.keySet(); }
    public int size() { return id2name.size(); }
    
}

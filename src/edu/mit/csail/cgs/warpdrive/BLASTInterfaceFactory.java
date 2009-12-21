/*
 * Created on Sep 28, 2006
 */
package edu.mit.csail.cgs.warpdrive;

import java.util.*;
import java.io.*;
import java.net.MalformedURLException;
import java.net.URL;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

/**
 * @author tdanford
 */
public class BLASTInterfaceFactory {
    
    public static BLASTInterfaceFactory defaultFactory;
    
    static { 
        defaultFactory= new BLASTInterfaceFactory();
    }
    
    private ResourceBundle gfRes;
    private Map<String,BLASTInterface> interfaces;

    public BLASTInterfaceFactory() {
        interfaces = new HashMap<String,BLASTInterface>();
        try {
            gfRes = ResourceBundle.getBundle("blast_servers");
        
            Enumeration<String> keys = gfRes.getKeys();
            while(keys.hasMoreElements()) { 
                String key = keys.nextElement();
                String genome = key; 
                String url = gfRes.getString(key);
                addInterface(genome, url);
            }
        } catch (Exception ex) {
            System.err.println("Couldn't find blast_servers resource");
        }
    }

    
    public int size() { return interfaces.size(); }
    
    public BLASTInterface getInterface(Genome g) { 
        if(!interfaces.containsKey(g.getVersion())) { return null; }
        return interfaces.get(g.getVersion());
    }
    
    public BLASTInterface getInterface(String genomeName) { 
        if(!interfaces.containsKey(genomeName)) { return null; }
        return interfaces.get(genomeName);
    }
    
    public Set<String> getGenomes() { return new TreeSet<String>(interfaces.keySet()); }
    
    private void addInterface(String genome, String url) {
        try {
            interfaces.put(genome, new BLASTInterface(genome, url));
        } catch (MalformedURLException e) {
            e.printStackTrace();
        }
    }
}

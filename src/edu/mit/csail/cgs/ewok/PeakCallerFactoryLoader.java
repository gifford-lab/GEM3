/*
 * Created on Sep 28, 2006
 */
package edu.mit.csail.cgs.ewok;

import java.util.*;
import java.io.*;

import edu.mit.csail.cgs.datasets.species.Genome;

/**
 * @author tdanford
 */
public class PeakCallerFactoryLoader {
    
    private ResourceBundle gfRes;
    private Map<String,PeakCallerFactory> factories;

    public PeakCallerFactoryLoader() {
        gfRes = ResourceBundle.getBundle("edu.mit.csail.cgs.ewok.peak_callers");
        factories = new HashMap<String,PeakCallerFactory>();
        Enumeration<String> keys = gfRes.getKeys();
        while(keys.hasMoreElements()) { 
            String key = keys.nextElement();
            addFactory(key, gfRes.getString(key));
        }
    }
    
    public PeakCallerFactory getFactory(String type) { 
        if(!factories.containsKey(type)) { return null; }
        return factories.get(type);
    }
    
    public Set<String> getTypes() { 
        return factories.keySet();
    }
    
    private void addFactory(String type, String pathName) {
        ClassLoader loader = ClassLoader.getSystemClassLoader();
        try {
            Class gfClass = loader.loadClass(pathName);
            PeakCallerFactory gf = (PeakCallerFactory)gfClass.newInstance();
            factories.put(type, gf);
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        } catch (InstantiationException e) {
            e.printStackTrace();
        } catch (IllegalAccessException e) {
            e.printStackTrace();
        }        
    }
}

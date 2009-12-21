package edu.mit.csail.cgs.ewok;

import java.util.*;
import java.io.*;

import edu.mit.csail.cgs.datasets.species.Genome;


public class RegionExpanderFactoryLoader<PRODUCT> {
    private ResourceBundle res;
    private Map<String,Map<String,RegionExpanderFactory<PRODUCT>>> factories;

    public RegionExpanderFactoryLoader(String ftype) {
        res = ResourceBundle.getBundle("edu.mit.csail.cgs.ewok."+ftype+"_factories");
        factories = new HashMap<String,Map<String,RegionExpanderFactory<PRODUCT>>>();
        Enumeration<String> keys = res.getKeys();
        while(keys.hasMoreElements()) { 
            String key = keys.nextElement();
            String[] keyArray = key.split(",");
            String genome = keyArray[0], type = keyArray[1];
            addFactory(genome, type, res.getString(key));
        }
    }
    
    public RegionExpanderFactory<PRODUCT> getFactory(Genome g, String type) { 
        if(!factories.containsKey(g.getVersion())) { return null; }
        if(!factories.get(g.getVersion()).containsKey(type)) { return null; }
        return factories.get(g.getVersion()).get(type);
    }
    
    public Set<String> getGenomes() { return factories.keySet(); }
    
    public Set<String> getTypes(Genome g) { 
        if(!factories.containsKey(g.getVersion())) { return new HashSet<String>(); }
        return factories.get(g.getVersion()).keySet();
    }
    
    public Set<String> getTypes(String gversion) { 
        if(!factories.containsKey(gversion)) { return new HashSet<String>(); }
        return factories.get(gversion).keySet();
    }
    
    private void addFactory(String genome, String type, String pathName) {
        ClassLoader loader = ClassLoader.getSystemClassLoader();
        try {
            Class fClass = loader.loadClass(pathName);
            RegionExpanderFactory<PRODUCT> f = (RegionExpanderFactory<PRODUCT>)fClass.newInstance();
            f.setType(type);
            if(!factories.containsKey(genome)) { factories.put(genome, new HashMap<String,RegionExpanderFactory<PRODUCT>>()); }
            Map<String,RegionExpanderFactory<PRODUCT>> genomeMap = factories.get(genome);
            genomeMap.put(type, f);
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        } catch (InstantiationException e) {
            e.printStackTrace();
        } catch (IllegalAccessException e) {
            e.printStackTrace();
        }        
    }   
}

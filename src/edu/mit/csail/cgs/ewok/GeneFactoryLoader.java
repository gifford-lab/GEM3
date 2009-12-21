/*
 * Created on Sep 28, 2006
 */
package edu.mit.csail.cgs.ewok;

import java.util.*;
import java.io.*;

import edu.mit.csail.cgs.datasets.species.Genome;

/**
 * @author tdanford
 *
 * <code>GeneFactoryLoader</code> is basically a Singleton.  It reads
 * a properties file in the classpath and creates GeneFactory objects.
 * Creating a GeneFactory with <code>getFactory</code> requires a
 * Genome and a type; the type generally indicates the set of gene
 * annotations to be accessed.
 *
 * The indirection through GeneFactoryLoader and GeneFactory allows
 * the properties file to specify more parameters than just the type
 */
public class GeneFactoryLoader {
    
    public static void main(String[] args) { 
        GeneFactoryLoader loader = new GeneFactoryLoader();
        for(String genome : loader.getGenomes()) { 
            System.out.println(genome);
            for(String type : loader.getTypes(genome)) { 
                System.out.println("\t" + type);
            }
        }
    }
    
    private ResourceBundle gfRes;
    private Map<String,Map<String,GeneFactory>> factories;

    public GeneFactoryLoader() {
        gfRes = ResourceBundle.getBundle("edu.mit.csail.cgs.ewok.gene_factories");
        factories = new HashMap<String,Map<String,GeneFactory>>();
        Enumeration<String> keys = gfRes.getKeys();
        while(keys.hasMoreElements()) { 
            String key = keys.nextElement();
            String[] keyArray = key.split(",");
            String genome = keyArray[0], type = keyArray[1];
            addFactory(genome, type, gfRes.getString(key));
        }
    }
    
    public GeneFactory getFactory(Genome g, String type) { 
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
            Class gfClass = loader.loadClass(pathName);
            GeneFactory gf = (GeneFactory)gfClass.newInstance();
            if(!factories.containsKey(genome)) { factories.put(genome, new HashMap<String,GeneFactory>()); }
            Map<String,GeneFactory> genomeMap = factories.get(genome);
            genomeMap.put(type, gf);
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        } catch (InstantiationException e) {
            e.printStackTrace();
        } catch (IllegalAccessException e) {
            e.printStackTrace();
        }        
    }
}

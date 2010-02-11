/*
 * Created on Sep 6, 2006
 */
package edu.mit.csail.cgs.utils.io.parsing.affyexpr;

import java.util.*;
import java.io.*;

/**
 * @author tdanford
 */
public class AffyProbeSet {
    
    private Map<String,AffyProbe> probes;
    private Map<String,Set<AffyProbe>> symbol, unigene, locus, name;

    public AffyProbeSet(File f) throws IOException { 
        probes = new HashMap<String,AffyProbe>();
        symbol = new HashMap<String,Set<AffyProbe>>();
        unigene = new HashMap<String,Set<AffyProbe>>();
        locus = new HashMap<String,Set<AffyProbe>>();
        name = new HashMap<String,Set<AffyProbe>>();
        
        BufferedReader br = new BufferedReader(new FileReader(f));
        String line = null;
        while((line = br.readLine()) != null) { 
            line = line.trim();
            if(line.length() > 0) { 
                AffyProbe ap = new AffyProbe(line);
                addProbe(ap);
            }
        }
        br.close();
    }
    
    private void addProbe(AffyProbe ap) { 
        if(probes.containsKey(ap.getProbeName())) { throw new IllegalArgumentException(ap.getProbeName()); }
        probes.put(ap.getProbeName(), ap);
        
        if(!symbol.containsKey(ap.getGeneSymbol())) { 
            symbol.put(ap.getGeneSymbol(), new HashSet<AffyProbe>()); 
        }
        symbol.get(ap.getGeneSymbol()).add(ap);
        
        if(!unigene.containsKey(ap.getUnigene())) { 
            unigene.put(ap.getUnigene(), new HashSet<AffyProbe>()); 
        }
        unigene.get(ap.getUnigene()).add(ap);
        
        if(!locus.containsKey(ap.getLocusID())) { 
            locus.put(ap.getLocusID(), new HashSet<AffyProbe>()); 
        }
        locus.get(ap.getLocusID()).add(ap);
        
        if(!name.containsKey(ap.getGeneName())) { 
            name.put(ap.getGeneName(), new HashSet<AffyProbe>()); 
        }
        name.get(ap.getGeneName()).add(ap);
    }
    
    public AffyProbe getAffyProbe(String pn) { return probes.get(pn); }
    public int size() { return probes.size(); }
    public Set<AffyProbe> getAllAffyProbes() { return new HashSet<AffyProbe>(probes.values()); }
    
    public Set<AffyProbe> lookupGeneSymbol(String v) { return symbol.get(v); }
    public Set<AffyProbe> lookupUnigene(String v) { return unigene.get(v); }
    public Set<AffyProbe> lookupLocusID(String v) { return locus.get(v); }
    public Set<AffyProbe> lookupGeneName(String v) { return name.get(v); }
}

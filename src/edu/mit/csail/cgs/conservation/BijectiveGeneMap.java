/*
 * Created on Aug 4, 2006
 */
package edu.mit.csail.cgs.conservation;

import java.util.*;
import java.io.*;
import edu.mit.csail.cgs.utils.Saveable;

/**
 * @author tdanford
 * 
 * A BijectiveGeneMap is a SpeciesGeneMap which returns one-to-one GeneMaps 
 * for two particular species.
 */
public class BijectiveGeneMap implements GeneMap, Saveable {
	
	/**
	 * SpeciesMapWrapper "wraps" a SpeciesGeneMap object, so that every map
	 * returned by the original species-map is itself "bijective."
	 * 
	 * @author tdanford
	 */
	public static class SpeciesMapWrapper implements SpeciesGeneMap {
		
		private SpeciesGeneMap innerMap;
		
		public SpeciesMapWrapper(SpeciesGeneMap m) { 
			innerMap = m;
		}

		public GeneMap getMap(String startSpecies, String targetSpecies) {
			GeneMap inner = innerMap.getMap(startSpecies, targetSpecies);
			return new BijectiveGeneMap(inner);
		}

		public boolean hasMap(String startSpecies, String targetSpecies) {
			return innerMap.hasMap(startSpecies, targetSpecies);
		} 
		
	}
    
    private Map<String,String> forward, reverse; 

    public BijectiveGeneMap() {
        forward = new HashMap<String,String>();
        reverse = new HashMap<String,String>();
    }
    
    /**
     * Creates a bijective map out of the "uniquely-mapped" elements of another GeneMap.
     * @param gm
     */
    public BijectiveGeneMap(GeneMap gm) { 
        forward = new HashMap<String,String>();
        reverse = new HashMap<String,String>();
        
        for(String id : gm.getDomainIDs()) { 
            Set<String> targets = gm.mapID(id);
            if(targets.size() == 1) { 
                Iterator<String> itr = targets.iterator();
                String target = itr.next();
                addMapping(id, target);
            }
        }
    }
    
    public BijectiveGeneMap(File f) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(f));
        String line = null;
        while((line = br.readLine()) != null) { 
            line = line.trim();
            if(line.length() > 0 && line.charAt(0) != '#') { 
                String[] array = line.split("\\s+");
                addMapping(array[0], array[1]);
            }
        }
        br.close();
    }
    
    public BijectiveGeneMap(DataInputStream dis) throws IOException {
        int s = dis.readInt();
        forward = new HashMap<String,String>();
        reverse = new HashMap<String,String>();        
        for(int i =0; i < s; i++) { 
            String id1 = dis.readUTF();
            String id2 = dis.readUTF();
            addMapping(id1, id2);
        }
    }
    
    public void outputPairs(PrintStream ps) { 
    	for(String k : forward.keySet()) { 
    		String v = forward.get(k);
    		ps.println(k + "\t" + v);
    	}
    }
    
    public void save(DataOutputStream dos) throws IOException {
        dos.writeInt(forward.size());
        for(String id : forward.keySet()) { 
            dos.writeUTF(id); 
            dos.writeUTF(forward.get(id));
        }
    }
    
    public void addMapping(String s1, String s2) { 
        if(forward.containsKey(s1)) { throw new IllegalArgumentException(s1); }
        if(reverse.containsKey(s2)) { throw new IllegalArgumentException(s2); }
        
        forward.put(s1, s2); 
        reverse.put(s1, s2);
    }

    public Set<String> getDomainIDs() {
        return forward.keySet();
    }

    public Set<String> getRangeIDs() {
        return reverse.keySet();
    }

    public Set<String> mapID(String id) {
        HashSet<String> ids = new HashSet<String>();
        if(forward.containsKey(id)) { ids.add(forward.get(id)); }
        return ids;
    }

    public Set<String> mapIDs(Collection<String> ids) {
        HashSet<String> targets = new HashSet<String>();
        for(String id : ids) { 
            if(forward.containsKey(id)) { 
                targets.add(forward.get(id));
            }
        }
        return targets;        
    }

}

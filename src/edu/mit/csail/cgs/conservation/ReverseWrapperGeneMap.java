/**
 * 
 */
package edu.mit.csail.cgs.conservation;

import java.util.*;
import java.io.*;

/**
 * @author tdanford
 * 
 * WrapperGeneMap composes a String->Stringset map with another, "internal" 
 * GeneMap to create a single composite GeneMap -- right now, I use this to 
 * compose the RefSeq -> Entrez map with the Entrez->Entrez Homologene-based
 * GeneMap.
 */
public class ReverseWrapperGeneMap implements GeneMap {
	
	private Map<String,Set<String>> layerMap;
	private GeneMap innerMap;
	private Set<String> range;
	
	public void outputPairs(PrintStream ps) { 
		for(String k : innerMap.getDomainIDs()) { 
			for(String v1 : innerMap.mapID(k)) { 
				if(layerMap.containsKey(v1)) { 
					for(String v2 : layerMap.get(v1)) { 
						ps.println(k + "\t" + v2);
					}
				}
			}
		}
	}
	
	public ReverseWrapperGeneMap(Map<String,Set<String>> layer, GeneMap gm) { 
		layerMap = layer;
		innerMap = gm;
		range = new HashSet<String>();
		for(String id : layerMap.keySet()) { 
			range.addAll(layerMap.get(id));
		}
	}

	public Set<String> getDomainIDs() {
		return innerMap.getDomainIDs();
	}

	public Set<String> getRangeIDs() {
		return range;
	}

	public Set<String> mapID(String id) {
		HashSet<String> targets = new HashSet<String>();
		for(String middle : innerMap.mapID(id)) { 
			if(layerMap.containsKey(middle)) { 
				targets.addAll(layerMap.get(middle));
			}
		}
		return targets;
	}

	public Set<String> mapIDs(Collection<String> ids) {
		HashSet<String> total = new HashSet<String>();
		for(String id : ids) { total.addAll(mapID(id)); }
		return total;
	}
}

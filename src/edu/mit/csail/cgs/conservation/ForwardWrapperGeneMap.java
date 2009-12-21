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
public class ForwardWrapperGeneMap implements GeneMap {
	
	private Map<String,Set<String>> layerMap;
	private GeneMap innerMap;
	
	public ForwardWrapperGeneMap(Map<String,Set<String>> layer, GeneMap gm) { 
		layerMap = layer;
		innerMap = gm;
	}
	
	public void outputPairs(PrintStream ps) { 
		for(String k : layerMap.keySet()) { 
			for(String v1 : layerMap.get(k)) { 
				for(String v2 : innerMap.mapID(v1)) { 
					ps.println(k + "\t" + v2);
				}
			}
		}
	}

	public Set<String> getDomainIDs() {
		return layerMap.keySet();
	}

	public Set<String> getRangeIDs() {
		return innerMap.getRangeIDs();
	}

	public Set<String> mapID(String id) {
		HashSet<String> targets = new HashSet<String>();
		if(layerMap.containsKey(id)) { 
			for(String middle : layerMap.get(id)) { 
				if(innerMap.getDomainIDs().contains(middle)) { 
					targets.addAll(innerMap.mapID(middle));
				}
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

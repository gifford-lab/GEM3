/**
 * 
 */
package edu.mit.csail.cgs.conservation;

import java.util.*;
import java.io.*;

import edu.mit.csail.cgs.utils.parsing.homologene.*;

/**
 * @author tdanford
 * 
 * HomologeneEntrezMap is a helper class, used to build an Entrez -> Entrez ID map
 * for a particular species, by using the Homologene identifiers for the 
 * Entrez genes, and mapping identifiers to each other if they share the same
 * Homologene identifiers.
 */
public class HomologeneEntrezMap implements GeneMap {
	
	private Set<String> range;
	private Map<String,Set<String>> mapping;
	
	public HomologeneEntrezMap(HomoloGeneAssignments assigns, 
			String tax1, String tax2) { 

		range = new HashSet<String>();
		mapping = new HashMap<String,Set<String>>();
		
		Map<String,Set<String>> tax1rev, tax2rev;
		tax1rev = new HashMap<String,Set<String>>();
		tax2rev = new HashMap<String,Set<String>>();
		
		Set<HomoloGeneAssignments.Entry> tax1Entries = 
			assigns.selectEntries(null, tax1, null);
		Set<HomoloGeneAssignments.Entry> tax2Entries = 
			assigns.selectEntries(null, tax2, null);
		
		for(HomoloGeneAssignments.Entry e : tax1Entries) { 
			if(!tax1rev.containsKey(e.homologeneID)) { 
				tax1rev.put(e.homologeneID, new HashSet<String>());
			}
			tax1rev.get(e.homologeneID).add(e.geneID);
		}
		for(HomoloGeneAssignments.Entry e : tax2Entries) { 
			if(!tax2rev.containsKey(e.homologeneID)) { 
				tax2rev.put(e.homologeneID, new HashSet<String>());
			}
			tax2rev.get(e.homologeneID).add(e.geneID);
		}
		
		for(String hgID : tax1rev.keySet()) { 
			for(String e1 : tax1rev.get(hgID)) { 
				if(tax2rev.containsKey(hgID)) { 
					for(String e2 : tax2rev.get(hgID)) { 
						if(!mapping.containsKey(e1)) { 
							mapping.put(e1, new HashSet<String>()); 
						}
						mapping.get(e1).add(e2);
					}
				}
			}
		}
	}

	public void outputPairs(PrintStream ps) { 
		for(String k : mapping.keySet()) { 
			for(String v : mapping.get(k)) { 
				ps.println(k + "\t" + v);
			}
		}
	}
	
	public Set<String> getDomainIDs() {
		return mapping.keySet();
	}

	public Set<String> getRangeIDs() {
		return range;
	}

	public Set<String> mapID(String id) {
		return mapping.get(id);
	}

	public Set<String> mapIDs(Collection<String> ids) {
		HashSet<String> total = new HashSet<String>();
		for(String id : ids) { total.addAll(mapID(id)); }
		return total;
	}

}

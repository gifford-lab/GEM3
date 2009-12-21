/*
 * Created on Aug 14, 2006
 */
package edu.mit.csail.cgs.conservation;

import java.io.*;
import java.util.*;

import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.utils.parsing.homologene.*;
import edu.mit.csail.cgs.utils.parsing.ncbi.*;
import edu.mit.csail.cgs.utils.parsing.ncbi.Gene2RefSeqParser.Entry;

/**
 * @author tdanford
 */
public class HomologeneEntrezSpeciesGeneMap implements SpeciesGeneMap {
    
    public static void main(String[] args) { 
        try { 
        	File tf = new File(args[0]);
        	File hf = new File(args[1]);
            
            Set<String> admitSpecies = new HashSet<String>();
            admitSpecies.add("9606");
            admitSpecies.add("10090");
            
        	HomologeneEntrezSpeciesGeneMap hsgm = new HomologeneEntrezSpeciesGeneMap(tf, hf, admitSpecies);
        	GeneMap gm = hsgm.getMap("9606", "10090");
        	
        	String line;
        	System.out.print(">"); System.out.flush();
        	BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
        	while((line = br.readLine()) != null) { 
        		line = line.trim();
        		if(line.length() > 0) { 
					if(line.startsWith("dom")) { 
						String[] array = line.split("\\s+");
						int min = array.length > 0 ? Integer.parseInt(array[1]) : 0;
						Set<String> domain = gm.getDomainIDs();
						for(String id : domain) { 
							Set<String> mapped = gm.mapID(id);
							if(mapped.size() > min) { 
								System.out.print(id + "(" + mapped.size() + ") ");
							}
						}
						System.out.println(" total: " + domain.size());
					} else { 
						Set<String> mapped = gm.mapID(line);
						System.out.print(line + " --> ");
						for(String id : mapped) { System.out.print(" " + id); }
						System.out.println();
					}
        		} else { 
        			return;
        		}
            	System.out.print(">"); System.out.flush();
        	}
        	
        } catch(IOException ie) { 
            ie.printStackTrace(System.err);
        }
    }
    
    private Set<String> speciesSet;
    private HomoloGeneAssignments assignments;

    public HomologeneEntrezSpeciesGeneMap(File taxFile, File homologeneDataFile, Set<String> admitSpecies) 
        throws IOException {
    	
    	System.out.println("Building HomoloGene-based Map...");
    	assignments = new HomoloGeneAssignments(homologeneDataFile);
    	System.out.println("\tBuilt HomologeneAssignments data-object.");
    	
    	TaxonomyFile tax_parser = new TaxonomyFile(taxFile);
        speciesSet = new HashSet<String>(tax_parser.getIDs());
        
        if(admitSpecies != null) { 
            Iterator<String> itr = speciesSet.iterator();
            while(itr.hasNext()) { 
                String sp = itr.next();
                if(!admitSpecies.contains(sp)) { 
                    itr.remove();
                }
            }
        }
    	
    	System.out.println("Done.");
    	
    }

    public boolean hasMap(String startSpecies, String targetSpecies) {
        return speciesSet.contains(startSpecies) && speciesSet.contains(targetSpecies);
    }

    public GeneMap getMap(String startSpecies, String targetSpecies) {
    	if(!hasMap(startSpecies,targetSpecies)) { return null; }
    	HomologeneEntrezMap inner = 
    		new HomologeneEntrezMap(assignments, startSpecies, targetSpecies);
        return inner;
    }
    
    public static Map<String,Set<String>> reverseMap(Map<String,Set<String>> original) { 
        HashMap<String,Set<String>> newMap = new HashMap<String,Set<String>>();
        
        for(String k : original.keySet()) { 
            for(String v : original.get(k)) { 
                if(!newMap.containsKey(v)) { newMap.put(v, new HashSet<String>()); }
                newMap.get(v).add(k);
            }
        }
        
        return newMap;
    }    
}

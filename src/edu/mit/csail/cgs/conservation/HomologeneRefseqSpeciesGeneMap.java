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
public class HomologeneRefseqSpeciesGeneMap implements SpeciesGeneMap {
    
    public static void main(String[] args) { 
        try { 
        	File tf = new File(args[0]);
        	File hf = new File(args[1]);
        	File nf = new File(args[2]);
            
            Set<String> admitSpecies = new HashSet<String>();
            admitSpecies.add("9606");
            admitSpecies.add("10090");
            
        	HomologeneRefseqSpeciesGeneMap hsgm = new HomologeneRefseqSpeciesGeneMap(tf, hf, nf, admitSpecies);
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
    private Map<String,Map<String,Set<String>>> refseqMaps;
    private HomoloGeneAssignments assignments;

    public HomologeneRefseqSpeciesGeneMap(File taxFile, 
    		File homologeneDataFile, File ncbi_gene2refseq, Set<String> admitSpecies) throws IOException {
    	
    	System.out.println("Building HomoloGene-based Map...");
    	assignments = new HomoloGeneAssignments(homologeneDataFile);
    	System.out.println("\tBuilt HomologeneAssignments data-object.");
    	
    	TaxonomyFile tax_parser = new TaxonomyFile(taxFile);
    	refseqMaps = new HashMap<String,Map<String,Set<String>>>();
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
    	
    	for(String species : speciesSet) { 
    		refseqMaps.put(species,buildRefSeq2EntrezMap(ncbi_gene2refseq,species));
    		System.out.println("\t++ " + species);
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
    	GeneMap forward = new ForwardWrapperGeneMap(refseqMaps.get(startSpecies),inner);
    	GeneMap reverse = new ReverseWrapperGeneMap(reverseMap(refseqMaps.get(targetSpecies)),forward);
    	return reverse;
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
    
    public static Map<String,Set<String>> buildEntrez2RefSeqMap(File gene2refseq, String taxid) throws IOException { 
        HashMap<String,Set<String>> map = new HashMap<String,Set<String>>();
        
        Iterator<Gene2RefSeqParser.Entry> entries = new Gene2RefSeqParser.DynamicLoader(gene2refseq);
        Iterator<Gene2RefSeqParser.Entry> taxEntries = 
            new FilterIterator<Gene2RefSeqParser.Entry,Gene2RefSeqParser.Entry>(new TaxFilter(taxid), entries);

        while(taxEntries.hasNext()) { 
            Gene2RefSeqParser.Entry entry = taxEntries.next();
            if(!entry.rnaAcc.equals("-")) {
                if(!map.containsKey(entry.geneID)) { map.put(entry.geneID, new HashSet<String>()); }
                map.get(entry.geneID).add(entry.rnaAcc);
            }
        }
        
        return map;
    }
    
    /**
     * This uses the NCBI-parsing code in utils.parsing.ncbi to build a map which 
     * has, as keys, ref-seq IDs and as values, the sets of Entrez gene identifiers
     * which those refseq IDs map to.  These sets usually (pretty much) only contain
     * one element.
     * 
     * @param gene2refseq  The file for the "gene2refseq" index from NCBI
     * @param taxid The particular taxID to search, or NULL for all taxonomies.
     * @return The map of refseq->entrez identifiers.
     * @throws IOException 
     */
    public static Map<String,Set<String>> buildRefSeq2EntrezMap(File gene2refseq, String taxid) throws IOException { 
        HashMap<String,Set<String>> map = new HashMap<String,Set<String>>();
        
        Iterator<Gene2RefSeqParser.Entry> entries = new Gene2RefSeqParser.DynamicLoader(gene2refseq);
        Iterator<Gene2RefSeqParser.Entry> taxEntries = 
            new FilterIterator<Gene2RefSeqParser.Entry,Gene2RefSeqParser.Entry>(new TaxFilter(taxid), entries);

        while(taxEntries.hasNext()) { 
            Gene2RefSeqParser.Entry entry = taxEntries.next();
            if(!entry.rnaAcc.equals("-")) {
                if(!map.containsKey(entry.rnaAcc)) { map.put(entry.rnaAcc, new HashSet<String>()); }
                map.get(entry.rnaAcc).add(entry.geneID);
            }
        }
        
        return map;
    }
    
    private static class TaxFilter implements Filter<Gene2RefSeqParser.Entry,Gene2RefSeqParser.Entry> {
        
        private String tax;
        
        public TaxFilter(String t) { tax = t; }

        public Gene2RefSeqParser.Entry execute(Gene2RefSeqParser.Entry a) {
            if(a.taxID.equals(tax)) { 
                return a;
            } else { 
                return null;
            }
        } 
        
    }
}

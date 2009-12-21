package edu.mit.csail.cgs.conservation;

import java.sql.SQLException;
import java.util.*;
import java.io.*;

import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.general.NamedStrandedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.*;

import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.parsing.ncbi.*;

public class LocationMapping implements Saveable {
    
    public static void main(String[] args) { 
        try {
            Genome genome = Organism.findGenome(args[0]);
            String geneType = args[1];
            int upstream = Integer.parseInt(args[2]), downstream = Integer.parseInt(args[3]);
            LocationMapping mapping = new LocationMapping(args[0], genome, geneType, upstream, downstream);
            
            if(args.length > 4) { 
                File g2rs = new File(args[4]);
                HashSet<String> speciesSet = new HashSet<String>();
                speciesSet.add(args[5]);
                Map<String,Set<String>> name_transfer = Gene2RefSeqParser.buildEntrez2RefSeqMap(g2rs, speciesSet);
                mapping = new LocationMapping(mapping, name_transfer);
				System.out.println("Built Entrez Gene mapping.");
            }
            
            File output = new File(args[0] + "_" + upstream + "_" + downstream + ".fasta");
            mapping.outputFASTA(output);
            
        } catch (NotFoundException e) {
            e.printStackTrace();
        } catch (SQLException e) {
            e.printStackTrace();
        } catch (UnknownRoleException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

	private String species;
	public Map<String,Set<LocationMappedID>> map;
	
	public LocationMapping(String s) { 
		species = s;
        map = new HashMap<String,Set<LocationMappedID>>();
	}
    
    public LocationMapping(LocationMapping base, Map<String,Set<String>> newKeyMap) { 
        species = base.species;
        map = new HashMap<String,Set<LocationMappedID>>();
        
        for(String k : newKeyMap.keySet()) {
            HashSet<LocationMappedID> ids = new HashSet<LocationMappedID>();
            for(String ok : newKeyMap.get(k)) {
                if(base.map.containsKey(ok)) { 
                    ids.addAll(base.map.get(ok));
                }
            }
            if(ids.size() > 0) { map.put(k, ids); }
        }
    }
	
    public LocationMapping(String s, Genome g, String geneSrc, int tssUp, int tssDown) throws SQLException, UnknownRoleException { 
		species = s;
		map = new HashMap<String,Set<LocationMappedID>>();

		Expander<NamedRegion,Gene> gen = new RefGeneGenerator<NamedRegion>(g, geneSrc);
		Iterator<NamedRegion> chroms = new ChromRegionIterator(g);
		Iterator<Gene> genes = new ExpanderIterator<NamedRegion,Gene>(gen, chroms);
		Iterator<NamedStrandedRegion> proms = 
			new MapperIterator<Gene,NamedStrandedRegion>(new GeneToPromoter(tssUp, tssDown), genes);
		
		while(proms.hasNext()) { 
			NamedStrandedRegion reg = proms.next();
			LocationMappedID lmi = new LocationMappedID(reg);
			addMapping(lmi);
		}
    }
	
	public LocationMapping(String s, File input, Genome g) throws IOException { 
		species = s;
		map = new HashMap<String,Set<LocationMappedID>>();
		BufferedReader br = new BufferedReader(new FileReader(input));
		String line = null;
		while((line = br.readLine()) != null) {
			line = line.trim();
			if(line.length() > 0 && line.charAt(0) != '#') { 
				String[] array = line.split("\\s+");
				String id = array[0], chrom = array[1];
				int start = Integer.parseInt(array[2]),
					end = Integer.parseInt(array[3]);
				Region r = new Region(g, chrom, start, end);
				LocationMappedID lmi = new LocationMappedID(id, r);
				addMapping(lmi);
			}
		}
		br.close();
	}
	
	public LocationMapping(Genome g, DataInputStream dis) throws IOException { 
		species = dis.readUTF();
		int size = dis.readInt();
		map = new HashMap<String,Set<LocationMappedID>>();
		for(int i = 0; i < size; i++) { 
			String id = dis.readUTF();
			map.put(id, new HashSet<LocationMappedID>());
			int ms = dis.readInt();
			for(int j = 0; j < ms; j++) { 
				LocationMappedID lmi = new LocationMappedID(g, dis);
				map.get(id).add(lmi);
			}
		}
	}
    
	public void save(DataOutputStream dos) throws IOException { 
		dos.writeUTF(species);
		dos.writeInt(map.size());
		for(String id : map.keySet()) { 
			dos.writeUTF(id);
			dos.writeInt(map.get(id).size());
			for(LocationMappedID lmi : map.get(id)) { 
				lmi.save(dos);
			}
		}
	}
	
    public void outputFASTA(File f) throws IOException { 
        PrintStream ps = new PrintStream(new FileOutputStream(f));
        for(String id : map.keySet()) { 
            for(LocationMappedID lmi : map.get(id)) { 
                lmi.outputFASTA(id, ps);
            }
        }
        ps.close();
    }
    
	public void addMapping(LocationMappedID lmi) { 
		if(!map.containsKey(lmi.getID())) { 
			map.put(lmi.getID(), new HashSet<LocationMappedID>()); 
		}
		map.get(lmi.getID()).add(lmi);
	}

	public String getSpecies() { return species; }
	public Set<String> getMappedIDs() { return map.keySet(); }
	
	public Set<LocationMappedID> getMappings(String id) { 
		return map.get(id);
	}
	
	public void debugPrint(PrintStream ps) { 
		for(String key : map.keySet()) { 
			ps.println(key);
			for(LocationMappedID id : map.get(key)) { 
				ps.println("\t" + id.getLocation().getLocationString());
			}
		}
	}

	/**
	 * Removes all the location-mappings which *aren't* in the given set of 
	 * "admissible" IDs.
	 * 
	 * @param admitIDs  The set of admissible IDs.
	 * @return The total number of mapped locations removed.
	 */
	public int removeLocations(Set<String> admitIDs) {
		int count = 0;
		Set<String> remove = new HashSet<String>();
		for(String k : map.keySet()) { 
			if(!admitIDs.contains(k)) { 
				count += map.get(k).size(); 
				remove.add(k);
			}
		}
		for(String k : remove) { 
			map.remove(k);
		}
		return count;
	}
	
	/**
	 * Given a Filter on Region objects, this removes all LocationMappedID's
	 * from the LocationMapping whose internal region does *not* satisfy the 
	 * filter.  Keys which, after this removal, have no associated Locations 
	 * left are also removed.  This is used, in practice, to remove Locations
	 * which contain no probes, when the LocationMapping is built from a 
	 * gene-set.
	 * 
	 * @param filter
	 */
	public int removeLocations(Filter<Region,Region> filter) {
		HashSet<String> emptyKeys = new HashSet<String>();
		int count = 0;
		for(String k : map.keySet()) { 
			Iterator<LocationMappedID> ids = map.get(k).iterator();
			while(ids.hasNext()) { 
				LocationMappedID id = ids.next();
				Region r = id.getLocation();
				if(filter.execute(r) == null) { 
					ids.remove();
					count += 1;
				}
			}
			
			if(map.get(k).size() == 0) { emptyKeys.add(k); }
		}
		
		for(String ek : emptyKeys) { 
			map.remove(ek);
		}
		return count;
	}
	
	// returns those IDs where the filter accepts *any* of the mapped regions
	// for that ID -- this is, in effect, a logical-OR over the mapped
	// regions.
	public Set<String> getOrAcceptedIDs(Filter<Region,Region> filter) { 
		HashSet<String> ids = new HashSet<String>();
		for(String id : map.keySet()) { 
			Set<LocationMappedID> lmis = map.get(id);
			Iterator<LocationMappedID> itr = lmis.iterator();
			boolean sat = false;
			while(itr.hasNext() && !sat) { 
				LocationMappedID lmi = itr.next();
				Region r = lmi.getLocation();
				if(filter.execute(r) != null) { 
					sat = true;
				}
			}
			
			if(sat) { 
				ids.add(id);
			}
		}
		return ids;
	}
}

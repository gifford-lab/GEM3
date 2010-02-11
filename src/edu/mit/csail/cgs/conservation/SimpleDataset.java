/**
 * 
 */
package edu.mit.csail.cgs.conservation;

import java.sql.SQLException;
import java.util.*;
import java.io.*;

import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.verbs.binding.*;
import edu.mit.csail.cgs.ewok.verbs.probers.*;
import edu.mit.csail.cgs.ewok.nouns.*;

import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.io.parsing.ncbi.*;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.chipchip.MSPProbe;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.locators.MSPLocator;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;

/**
 * @author tdanford
 *
 */
public class SimpleDataset implements ConservationDataset, Saveable {
	
	public static void main(String[] args) { 
		try {
			Organism mouse = Organism.getOrganism("Mus musculus");
			Organism human = Organism.getOrganism("Homo sapiens");
			Genome mm6 = mouse.getGenome("mm6");
			Genome hg17 = human.getGenome("hg17");
            
            String hgSpecies = "9606";
            String mmSpecies = "10090";
            
            File g2rs = new File("gene2refseq");
            Map<String,Set<String>> hg_entrez2refseq = 
                HomologeneRefseqSpeciesGeneMap.buildEntrez2RefSeqMap(g2rs, hgSpecies);
            Map<String,Set<String>> mm_entrez2refseq = 
                HomologeneRefseqSpeciesGeneMap.buildEntrez2RefSeqMap(g2rs, mmSpecies);

			System.out.println(hg_entrez2refseq.size() + " Human Entrez IDs.");
			System.out.println(mm_entrez2refseq.size() + " Mouse Entrez IDs.");
            
			int up = 8000;
			int down = 2000;
			
			LocationMapping humanMapping = 
				new LocationMapping(hgSpecies, hg17, "refGene", up, down);
            humanMapping = new LocationMapping(humanMapping, hg_entrez2refseq);
            
			LocationMapping mouseMapping = 
				new LocationMapping(mmSpecies, mm6, "refGene", up, down);
            mouseMapping = new LocationMapping(mouseMapping, mm_entrez2refseq);

			System.out.println(humanMapping.getMappedIDs().size() + " Mapped Human IDs.");
			System.out.println(mouseMapping.getMappedIDs().size() + " Mapped Mouse IDs.");
            
            Vector<String> expts = new Vector<String>();
            expts.add("Hs Oct4:hES9:culture vs WCE:hES9:culture;H10_MERGE_1.2_2.2_x");
            expts.add("Hs Nanog:hES9:culture vs WCE:hES9:culture;H10_MERGE_1.1_2.1");
            expts.add("Mm OCT4:MES:culture vs WCE:MES:culture;rep1");
            expts.add("Mm OCT4:MES:culture vs WCE:MES:culture;rep2");
            expts.add("Mm NANOG:MES:culture vs WCE:MES:culture;rep1");
            expts.add("Mm NANOG:MES:culture vs WCE:MES:culture;rep2");
			
			int humanRemoved = 0, mouseRemoved = 0;
            for(String expt : expts) { 
				boolean isHuman = expt.startsWith("Hs");
				Genome g = isHuman ? hg17 : mm6;
				
				ProbedRegionFilter filter = createProbedRegionFilter(g, expt);
				if(isHuman) { 
					humanRemoved += humanMapping.removeLocations(filter);
				} else { 
					mouseRemoved += mouseMapping.removeLocations(filter);
				}
				filter.close();
			}
			
			System.out.println("Removed " + humanRemoved + " human regions.");
			System.out.println("Removed " + mouseRemoved + " mouse regions.");
			
			SimpleDataset mouseData, humanData;
			mouseData = new SimpleDataset(mouseMapping);
			humanData = new SimpleDataset(humanMapping);
			
			System.out.println("# Human IDs: " + humanData.getIDs().size());
			System.out.println("# Mouse IDs: " + mouseData.getIDs().size());
			
            for(String expt : expts) { 
				boolean isHuman = expt.startsWith("Hs");
				Genome g = isHuman ? hg17 : mm6;
				PeakCaller caller = createCaller(g, expt);
				ExptDescriptor ed = new ExptDescriptor(expt);
				ExptBindingEvents events = 
					new ExptBindingEvents((isHuman ? hg17 : mm6), ed, caller);
				int boundIDs = 0;
				if(isHuman) { 
					boundIDs = humanData.addExptBindingEvents(events);
				} else { 
					boundIDs = mouseData.addExptBindingEvents(events);
				}
				
				System.out.println(boundIDs + " --> " + expt);
			}
            
            ExptDescriptor andED = new ExptDescriptor("Mm Oct4 Combined");
            ExptDescriptor a1 = new ExptDescriptor(expts.get(2));
            ExptDescriptor a2 = new ExptDescriptor(expts.get(3));
            int andCount = mouseData.addANDExpt(andED, a1, a2);
            System.out.println(andED.getName() + " --> " + andCount);

            andED = new ExptDescriptor("Mm Nanog Combined");
            a1 = new ExptDescriptor(expts.get(4));
            a2 = new ExptDescriptor(expts.get(5));
            andCount = mouseData.addANDExpt(andED, a1, a2);
            System.out.println(andED.getName() + " --> " + andCount);
            
            DataOutputStream dos;
            
            String hgFname = "human_data_"+ up + "_" + down + ".dat";
            dos = new DataOutputStream(new FileOutputStream(new File(hgFname)));
            humanData.save(dos);
            dos.close();

            String mmFname = "mouse_data_" + up + "_" + down + ".dat";
            dos = new DataOutputStream(new FileOutputStream(new File(mmFname)));
            mouseData.save(dos);
            dos.close();
			
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
    	
	public static ProbedRegionFilter createProbedRegionFilter(Genome g, String expt) { 
		String[ ]array = expt.split(";");
		MSPLocator loc = new MSPLocator(g, array[0], array[1]);
		MSPImmediateProbeGenerator probeGen = new MSPImmediateProbeGenerator(g, loc);
		RegionProber<MSPProbe> prober = new RegionProber.Wrapper<MSPProbe>(probeGen);
		return new ProbedRegionFilter<MSPProbe>(prober);
	}
	
	public static PeakCaller createCaller(Genome g, String expt) { 
		PeakCaller caller = null;
		String[ ]array = expt.split(";");
		MSPLocator loc = new MSPLocator(g, array[0], array[1]);
		MSPImmediateProbeGenerator probeGen = new MSPImmediateProbeGenerator(g, loc);
		RegionProber<MSPProbe> prober = new RegionProber.Wrapper<MSPProbe>(probeGen);
		YoungLabRosettaTuplePeakFinder finder = 
			new YoungLabRosettaTuplePeakFinder(g, array[0] + "," + array[1]);
		
		caller = new PeakCaller.FromTupler<MSPProbe>(prober, finder, 3);		
		return caller;
	}

	private LocationMapping mapping;
	private Map<ExptDescriptor,ExptBindingEvents> eventMap;
	private Map<ExptDescriptor,Filter<Region,Region>> bindingFilters;
	private Map<ExptDescriptor,Set<String>> boundIDs;
	
	public SimpleDataset(LocationMapping map) {
		mapping = map;
		eventMap = new HashMap<ExptDescriptor,ExptBindingEvents>();
		bindingFilters = new HashMap<ExptDescriptor,Filter<Region,Region>>();
		boundIDs = new HashMap<ExptDescriptor,Set<String>>();
	}
	
	public SimpleDataset(Genome g, DataInputStream dis) throws IOException { 
		mapping = new LocationMapping(g, dis);
		eventMap = new HashMap<ExptDescriptor,ExptBindingEvents>();
		bindingFilters = new HashMap<ExptDescriptor,Filter<Region,Region>>();
		boundIDs = new HashMap<ExptDescriptor,Set<String>>();
		
		int s = dis.readInt();
		for(int i = 0; i < s; i++) { 
			ExptDescriptor ed = new ExptDescriptor(dis);
			ExptBindingEvents events = new ExptBindingEvents(g, dis);
			Set<String> ids = new HashSet<String>();
			int ss = dis.readInt();
			for(int ii = 0; ii < ss; ii++) { 
				ids.add(dis.readUTF());
			}
			
			eventMap.put(ed, events);
			boundIDs.put(ed, ids);
			bindingFilters.put(ed, eventMap.get(ed).getBindingFilter());
		}
	}
    
    public void outputPythonBindingData(File f) throws IOException { 
        PrintStream ps = new PrintStream(new FileOutputStream(f));
        for(String key : mapping.getMappedIDs()) { 
            for(LocationMappedID id : mapping.getMappings(key)) { 
                String regionName = key + ";" + id.getID() + ";" + id.getLocation().getLocationString();
                ps.println(">" + regionName);
                for(ExptDescriptor ed : eventMap.keySet()) { 
                    ps.print(ed.getName());
                    Vector<BindingEvent> events = eventMap.get(ed).getEventSubset(id.getLocation());
                    for(BindingEvent evt : events) { 
                        int start = evt.getStart() - id.getLocation().getStart();
                        int end = evt.getEnd() - id.getLocation().getStart();
                        ps.print("\t" + start + "," + end);
                    }
                    ps.println();
                }
            }
        }
        ps.close();
    }
	
	public void constrainIDs(Set<String> ids) { 
		mapping.removeLocations(ids);
		for(ExptDescriptor ed : boundIDs.keySet()) { 
			Iterator<String> itr = boundIDs.get(ed).iterator();
			while(itr.hasNext()) { 
				String key = itr.next();
				if(!ids.contains(key)) { 
					itr.remove();
				}
			}
		}
	}
	
	public void save(DataOutputStream dos) throws IOException { 
		mapping.save(dos);
		dos.writeInt(eventMap.size());
		for(ExptDescriptor ed : eventMap.keySet()) { 
			ed.save(dos);
			eventMap.get(ed).save(dos);
			dos.writeInt(boundIDs.get(ed).size());
			for(String id : boundIDs.get(ed)) { 
				dos.writeUTF(id);
			}
		}
	}
	
	public int addANDExpt(ExptDescriptor newED, 
			ExptDescriptor ed1, ExptDescriptor ed2) { 
		
		if(!boundIDs.containsKey(ed1) || !boundIDs.containsKey(ed2)) { 
			throw new IllegalArgumentException();
		}
		
		Filter<Region,Region> f1 = bindingFilters.get(ed1);
		Filter<Region,Region> f2 = bindingFilters.get(ed2);
		Filter<Region,Region> filter = 
			new Filter.Compose<Region, Region, Region>(f1, f2);
		
		SetTools<String> strtools = new SetTools<String>();
		
		bindingFilters.put(newED, filter);
		eventMap.put(newED, new ExptBindingEvents(newED));
		boundIDs.put(newED, strtools.intersection(boundIDs.get(ed1), 
				boundIDs.get(ed2)));
		
		return boundIDs.get(newED).size();
	}
	
	public int addExptBindingEvents(ExptBindingEvents ebe) { 
		if(!eventMap.containsKey(ebe.getExpt())) { 
			eventMap.put(ebe.getExpt(), ebe);
			bindingFilters.put(ebe.getExpt(), ebe.getBindingFilter());
			boundIDs.put(ebe.getExpt(), 
					mapping.getOrAcceptedIDs(bindingFilters.get(ebe.getExpt())));
			return boundIDs.get(ebe.getExpt()).size();
			
		} else { 
			throw new IllegalArgumentException("Duplicate descriptor.");
		}
	}
	
	public Set<String> getBound(ExptDescriptor ed, BindingOptions bo) {
		return boundIDs.containsKey(ed) ? boundIDs.get(ed) : new HashSet<String>();
	}

	public Vector<ExptDescriptor> getExpts() {
		return new Vector<ExptDescriptor>(eventMap.keySet());
	}

	public Set<String> getIDs() {
		return mapping.getMappedIDs();
	}

	public boolean isBound(String id, ExptDescriptor ed, BindingOptions bo) {
		return boundIDs.containsKey(ed) && boundIDs.get(ed).contains(id);
	}
}

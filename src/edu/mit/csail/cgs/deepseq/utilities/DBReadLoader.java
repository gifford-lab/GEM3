package edu.mit.csail.cgs.deepseq.utilities;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.Vector;

import edu.mit.csail.cgs.datasets.chipseq.ChipSeqAlignment;
import edu.mit.csail.cgs.datasets.chipseq.ChipSeqExpt;
import edu.mit.csail.cgs.datasets.chipseq.ChipSeqHit;
import edu.mit.csail.cgs.datasets.chipseq.ChipSeqLoader;
import edu.mit.csail.cgs.datasets.chipseq.ChipSeqLocator;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.ExtReadHit;
import edu.mit.csail.cgs.deepseq.ReadHit;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

/**
 * Load reads from a ChipSeq experiment in the database
 * @author shaun
 *
 */
public class DBReadLoader extends ReadLoader{

	private List<String> exptNames =new ArrayList<String>();
	private List<ChipSeqLoader> loaders =new ArrayList<ChipSeqLoader>();
	private HashMap<ChipSeqLoader, List<ChipSeqAlignment>> loaderAligns = new HashMap<ChipSeqLoader, List<ChipSeqAlignment>>();
	
	
	public DBReadLoader(Genome g, List<ChipSeqLocator> locs, int rLen) throws NotFoundException, SQLException {
		super(g, rLen);
		try {
			
			//Initialize
			for(ChipSeqLocator locator : locs){
				String exptName = locator.getExptName(); exptNames.add(exptName);
				ChipSeqLoader loader = new ChipSeqLoader();
				 
				LinkedList<ChipSeqAlignment> alignments = new LinkedList<ChipSeqAlignment>();
		        if (locator.getAlignName() == null) {
		            if(locator.getReplicates().isEmpty()) { //No alignment name, no replicate names
		            	Collection<ChipSeqExpt> expts = loader.loadExperiments(locator.getExptName());
		        		for(ChipSeqExpt expt : expts) { 
		                	Collection<ChipSeqAlignment> aligns;
							aligns = loader.loadAllAlignments(expt);
							for (ChipSeqAlignment currentAlign : aligns) {
		            			if (currentAlign.getGenome().equals(g)) { 
		            				ChipSeqLocator currentLoc = new ChipSeqLocator(expt.getName(), 
		                                    expt.getReplicate(), currentAlign.getName());
		            				alignments.add(currentAlign);
		    						break;
		    					}
		            		}
		    			}
		    			List<ChipSeqLocator> collapsedLocs = new Vector<ChipSeqLocator>(this.collapseLocatorsByName(locs));
		    			if (collapsedLocs.size() != 1) {
		    				System.err.println(collapsedLocs.size() + " collapsed locators");
		    				System.exit(0);
		    			}
		    			locator = collapsedLocs.get(0);
		            } else { //No alignment name, given replicate names
		                for(String repName : locator.getReplicates()) { 
		                    ChipSeqExpt expt = loader.loadExperiment(locator.getExptName(), repName);
		                    ChipSeqAlignment alignment = 
		                        loader.loadAlignment(expt, locator.getAlignName(), g);
		                    if(alignment != null) { 
		                        locator = new ChipSeqLocator(locator.getExptName(),
		                                                     locator.getReplicates(),
		                                                     alignment.getName());
		                        alignments.add(alignment);
		                        break;
		                    }
		                }
		            }
		        } else {
		        	if(locator.getReplicates().isEmpty()) {//Given alignment name, no replicate names
		        		Collection<ChipSeqExpt> expts = loader.loadExperiments(locator.getExptName());
		        		for(ChipSeqExpt expt : expts) { 
		                	Collection<ChipSeqAlignment> aligns;
							aligns = loader.loadAllAlignments(expt);
							for (ChipSeqAlignment currentAlign : aligns) {
		            			if (currentAlign.getGenome().equals(g) && currentAlign.getName().equals(locator.getAlignName())) { 
		            				ChipSeqLocator currentLoc = new ChipSeqLocator(expt.getName(), 
		                                    expt.getReplicate(), currentAlign.getName());
		            				alignments.add(currentAlign);
		    						break;
		    					}
		            		}
		    			}
		    			List<ChipSeqLocator> collapsedLocs = new Vector<ChipSeqLocator>(this.collapseLocatorsByName(locs));
		    			if (collapsedLocs.size() != 1) {
		    				System.err.println(collapsedLocs.size() + " collapsed locators");
		    				System.exit(0);
		    			}
		    			locator = collapsedLocs.get(0);
		            }else{
		            	for (String replicate : locator.getReplicates()) {//Given alignment name, given replicate names
		        			alignments.add(loader.loadAlignment(loader.loadExperiment(locator.getExptName(),
                                                                                      replicate), 
                                                                locator.getAlignName(),
                                                                g));
		        		}
		            }
		        }
		        loaders.add(loader);
		        loaderAligns.put(loader, alignments);
			}countHits();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	protected double countHits(){
		totalHits=0;
		totalWeight=0;
		try {
			for(ChipSeqLoader loader : loaders){
				List<ChipSeqAlignment> alignments = loaderAligns.get(loader);
				for(ChipSeqAlignment alignment : alignments) { 
					double currHits = (double)loader.countAllHits(alignment);
					totalHits+=currHits;
					double currWeight = (double)loader.weighAllHits(alignment);
					totalWeight +=currWeight;
				}
			}
		}catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}			
		return totalHits;
	}
	
	//Load reads in a region
	public List<ReadHit> loadHits(Region r) {
		try {
			ArrayList<ReadHit> total = new ArrayList<ReadHit>();
			for(ChipSeqLoader loader : loaders){
				Collection<ChipSeqHit> hits = loader.loadByRegion(loaderAligns.get(loader), r);
				total.addAll(hits2reads(hits));		
			}
			return total;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return new ArrayList<ReadHit>();
		}
	}
	
	//Load extended reads in a region
	public List<ExtReadHit> loadExtHits(Region r, int startShift, int fivePrimeExt, int threePrimeExt) {
		try {
			ArrayList<ExtReadHit> total = new ArrayList<ExtReadHit>();
			for(ChipSeqLoader loader : loaders){
				Collection<ChipSeqHit> hits = loader.loadByRegion(loaderAligns.get(loader), r);
				total.addAll(hits2extreads(hits, startShift, fivePrimeExt, threePrimeExt));		
			}
			return total;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return new ArrayList<ExtReadHit>();
		}
	}
	/**
	 * Get sorted start positions of all reads (regardless of strand) in one chrom
	 */
	public int[] getStartCoords(String chrom) {
        try {
            ChipSeqLoader loader = loaders.get(0);
            List<Integer> list = loader.positionsByRegion(loaderAligns.get(loader).get(0),
                                                          new Region(gen, chrom, 0, gen.getChromLength(chrom)));
            int[] allReads = new int[list.size()];
            for (int i=0;i<allReads.length;i++){
                allReads[i] = list.get(i);
            }		
            return allReads ;
        } catch (IOException e) {
            e.printStackTrace();
            return null;
        }
	}
	
	// count number of reads in region
	public int countHits (Region r) {
		try {
			int count = 0;
			for(ChipSeqLoader loader : loaders){
				count += loader.countByRegion(loaderAligns.get(loader), r);
			}
			return count;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return 0;
		}
	}
	// sum weight of reads in region
	public double sumWeights (Region r) {
		try {
			double sum = 0;
			for(ChipSeqLoader loader : loaders){
				sum += loader.weightByRegion(loaderAligns.get(loader), r);
			}
			return sum;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return 0;
		}
	}
	
	//Collapse locators
	private Collection<ChipSeqLocator> collapseLocatorsByName(Collection<ChipSeqLocator> locs) { 
        LinkedHashMap<String,Map<String,Set<String>>> map = 
            new LinkedHashMap<String,Map<String,Set<String>>>();
        
        for(ChipSeqLocator loc : locs) { 
            String exptName = loc.getExptName();
            String alignName = loc.getAlignName();
            if(!map.containsKey(exptName)) { map.put(exptName, new LinkedHashMap<String,Set<String>>()); }
            if(!map.get(exptName).containsKey(alignName)) { map.get(exptName).put(alignName, new TreeSet<String>()); }
            map.get(exptName).get(alignName).addAll(loc.getReplicates());
        }
        
        LinkedList<ChipSeqLocator> collapsed = new LinkedList<ChipSeqLocator>();
        
        for(String exptName : map.keySet()) { 
            for(String alignName : map.get(exptName).keySet()) { 
                ChipSeqLocator newloc = new ChipSeqLocator(exptName, map.get(exptName).get(alignName), alignName);
                collapsed.add(newloc);
            }
        }        
        return collapsed;
    }
	
	//Convert ChipSeqHits to Reads
	private Collection<ReadHit> hits2reads(Collection<ChipSeqHit> hits){
		ArrayList<ReadHit> r = new ArrayList<ReadHit>();
		for(ChipSeqHit h : hits){
			r.add(new ReadHit(h.getGenome(), -1, h.getChrom(), h.getStart(), h.getEnd(), h.getStrand(), h.getWeight()));
		}return(r);
	}
	//Convert ChipSeqHits to ExtReads
	private Collection<ExtReadHit> hits2extreads(Collection<ChipSeqHit> hits, int startShift, int fivePrimeExt, int threePrimeExt){
		ArrayList<ExtReadHit> r = new ArrayList<ExtReadHit>();
		for(ChipSeqHit h : hits){
			r.add(new ExtReadHit(h.getGenome(), -1,h.getChrom(), h.getStart(), h.getEnd(), h.getStrand(), h.getWeight(), startShift, fivePrimeExt, threePrimeExt));
		}return(r);
	}

	//Count total read hits on one strand
	protected double countStrandedWeight(char strand) {
		double count=0;
		try {
			for(ChipSeqLoader loader : loaders){
				List<ChipSeqAlignment> alignments = loaderAligns.get(loader);
				for(ChipSeqAlignment alignment : alignments) {
					Pair<Long,Double> cw=loader.getAlignmentStrandedCountWeight(alignment,strand);	
					count+= cw.cdr();
				}				
			}
			return count;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return -1;
		}
	}
	//Close the loaders
	public void cleanup(){
		for(ChipSeqLoader l : loaders){
			l.close();
		}
	}
}

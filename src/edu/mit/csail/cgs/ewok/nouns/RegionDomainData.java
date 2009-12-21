package edu.mit.csail.cgs.ewok.nouns;

import java.util.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;

public class RegionDomainData {

	protected Vector<Region> domains;
	private Region region;
	protected Vector<Region> uncoveredRegions;
	
	public RegionDomainData(Region r) { 
		domains = new Vector<Region>();
		region = r;
		uncoveredRegions = new Vector<Region>();
		uncoveredRegions.add(region);
	}
	
	public Region getRegion() { return region; }
	public int getNumDomains() { return domains.size(); }
	
	public int getUncoveredBP() { 
		int l = 0;
		for(Region uc : uncoveredRegions) { 
			l += uc.getWidth();
		}
		return l;
	}
	
	public int getCoveredBP() { 
		return region.getWidth() - getUncoveredBP(); 
	} 
	
	public int getBP() { return region.getWidth(); }
	
	public void addDomain(Region dom) { 
		if(region.overlaps(dom)) { 
			domains.add(dom);
			
			Iterator<Region> itr = uncoveredRegions.iterator();
			LinkedList<Region> nurs = new LinkedList<Region>();
			while(itr.hasNext()) { 
				Region uc = itr.next();
				if(uc.overlaps(dom)) { 
					itr.remove();
					nurs.addAll(subtractRegion(uc, dom)); 
				}
			}
			
			uncoveredRegions.addAll(nurs);
		}
	}
	
	public Collection<Region> subtractRegion(Region r1, Region r2) { 
		LinkedList<Region> regs = new LinkedList<Region>();
		Genome g = r1.getGenome();
		String chrom = r1.getChrom();
		if(r1.overlaps(r2)) { 
			if (r2.contains(r1)) { 
				// do nothing
			} else { 
				if(r1.getStart() < r2.getStart()) { 
					regs.add(new Region(g, chrom, r1.getStart(), 
							r2.getStart()-1));
				}

				if(r1.getEnd() > r2.getEnd()) { 
					regs.add(new Region(g, chrom, r2.getEnd()+1,
							r1.getEnd()));
				}
			}
		}
		return regs;
	}
}

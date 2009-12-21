/*
 * Author: tdanford
 * Date: Mar 20, 2009
 */
package edu.mit.csail.cgs.ewok.verbs.sequence;

import java.util.*;

import edu.mit.csail.cgs.datasets.chippet.RunningOverlapSum;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.utils.iterators.EmptyIterator;

public class RegionMatching<RBase extends Region,RCover extends Region> {
	
	private Expander<Region,RBase> baseFinder;
	private Map<RBase,Set<RCover>> cover;
	
	public RegionMatching(Expander<Region,RBase> base) {
		baseFinder = base;
		cover = new HashMap<RBase,Set<RCover>>();
	}
	
	public void addToCoverage(RCover c) { 
		Iterator<RBase> covered = baseFinder.execute(c);
		while(covered.hasNext()) { 
			RBase b = covered.next();
			if(!cover.containsKey(b)) { 
				cover.put(b, new HashSet<RCover>());
			}
			cover.get(b).add(c);
		}
	}
	
	public int numCovered() { return cover.size(); }
	
	public Iterator<RBase> covered() { return cover.keySet().iterator(); }
	
	public Iterator<RCover> coverage(RBase b) { 
		return cover.containsKey(b) ? cover.get(b).iterator() : 
			new EmptyIterator<RCover>();
	}
	
	public Iterator<RCover> allCoverage() { 
		return new ExpanderIterator<RBase,RCover>(
				new Expander<RBase,RCover>() { 
					public Iterator<RCover> execute(RBase b) { 
						return coverage(b); 
					}
				}, covered());
	}
	
	public int totalContiguousCoverage(RBase b) { 
		int s = 0;
		Iterator<Region> cov = contiguousCoverage(b);
		while(cov.hasNext()) { 
			s += cov.next().getWidth();
		}
		return s;
	}
	
	public Iterator<Region> contiguousCoverage(RBase b) { 
		ArrayList<Region> contig = new ArrayList<Region>();
		RunningOverlapSum sum = new RunningOverlapSum(b.getGenome(), b.getChrom());
		Iterator<RCover> cov = coverage(b);
		while(cov.hasNext()) { sum.addRegion(cov.next()); }
		
		Collection<Region> rs = sum.collectRegions(1);
		for(Region r : rs) { 
			if(b.contains(r)) { 
				contig.add(r);
			} else { 
				int s = Math.max(b.getStart(), r.getStart());
				int e = Math.min(b.getEnd(), r.getEnd());
				contig.add(new Region(b.getGenome(), b.getChrom(), s, e));
			}
		}
		
		return contig.iterator();
	}
	
	public Iterator<RCover> pointCoverage(RBase b, int offset) { 
		Point p = new Point(b.getGenome(), b.getChrom(), b.getStart()+offset);
		return new FilterIterator<RCover,RCover>(
				new ContainsFilter(p), coverage(b));
	}
	
	private class ContainsFilter implements Filter<RCover,RCover> { 
		private Point p;
		public ContainsFilter(Point pp) { p = pp; }
		public RCover execute(RCover c) { 
			return c.contains(p) ? c : null;
		}
	}
}

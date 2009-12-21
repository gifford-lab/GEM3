/*
 * Author: tdanford
 * Date: Mar 30, 2009
 */
package edu.mit.csail.cgs.utils;

import java.util.*;

import edu.mit.csail.cgs.utils.iterators.EmptyIterator;
import edu.mit.csail.cgs.utils.iterators.SerialIterator;
import edu.mit.csail.cgs.utils.models.Model;

/**
 * 
 * Collects a series of "UnitCoverage" objects, keyed by arbitrary strings.
 * This is most useful when we're looking at the coverage of an entire genome, 
 * in which case each UnitCoverage corresponds to a Chromosome and the Coverage
 * object itself is the genome.  
 * 
 * See also: ewok.verbs.sequence.RegionCoverage, which extends the functionality
 * of this class by appending a Genome object and allowing the user to deal 
 * directly with Region objects instead of abstract intervals (as here.)  
 * 
 * @author tdanford
 */
public class Coverage {
	
	private Map<String,UnitCoverage> covermap;
	
	public Coverage() { 
		covermap = new TreeMap<String,UnitCoverage>();
	}
	
	public Coverage(CoverageModel m) { 
		this();
		for(int i = 0; i < m.keys.length; i++) { 
			covermap.put(m.keys[i], new UnitCoverage(m.coverages[i]));
		}
	}
	
	public Coverage(Coverage c) { 
		this(c.covermap);
	}
	
	public Coverage(Map<String,UnitCoverage> kmap) { 
		covermap = new TreeMap<String,UnitCoverage>();
		for(String k : kmap.keySet()) { 
			covermap.put(k, new UnitCoverage(kmap.get(k)));
		}
	}
	
	public Integer[] rightNearest(String key, int pos) { 
		return covermap.containsKey(key) ? covermap.get(key).rightNearest(pos) : null;
	}
	
	public Integer[] leftNearest(String key, int pos) { 
		return covermap.containsKey(key) ? covermap.get(key).leftNearest(pos) : null;
	}
	
	public boolean hasCoverage(String key, int start, int end) { 
		return covermap.containsKey(key) ? 
				covermap.get(key).hasOverlap(start, end) : false;
	}
	
	public boolean isContained(String key, int start, int end) { 
		return covermap.containsKey(key) ? 
				covermap.get(key).isContained(start, end) : false;
	}
	
	public Iterator<String> keys() { 
		return covermap.keySet().iterator();
	}
	
	public Iterator<Integer[]> units() { 
		ArrayList<Iterator<Integer[]>> units = new ArrayList<Iterator<Integer[]>>();
		for(String k : covermap.keySet()) { units.add(covermap.get(k).covered()); }
		return new SerialIterator<Integer[]>(units.iterator());
	}
	
	public Iterator<Integer[]> units(String k) { 
		return covermap.get(k).units();
	}
	
	public int size() { 
		int s = 0;
		for(String k : covermap.keySet()) { 
			s += covermap.get(k).size();
		}
		return s;
	}
	
	public Coverage union(Coverage c) { 
		Map<String,UnitCoverage> kmap = new TreeMap<String,UnitCoverage>();
		for(String k : covermap.keySet()) {
			if(c.covermap.containsKey(k)) { 
				kmap.put(k, covermap.get(k).union(c.covermap.get(k)));
			} else { 
				kmap.put(k, new UnitCoverage(covermap.get(k)));
			}
		}
		for(String k : c.covermap.keySet()) { 
			if(!covermap.containsKey(k)) { 
				kmap.put(k, new UnitCoverage(c.covermap.get(k)));
			}
		}
		return new Coverage(kmap);
	}

	public Coverage subtract(Coverage c) { 
		Map<String,UnitCoverage> kmap = new TreeMap<String,UnitCoverage>();
		for(String k : covermap.keySet()) {
			if(c.covermap.containsKey(k)) { 
				kmap.put(k, covermap.get(k).subtract(c.covermap.get(k)));
			} else { 
				kmap.put(k, new UnitCoverage(covermap.get(k)));
			}
		}
		return new Coverage(kmap);
	}
	
	public Coverage intersection(Coverage c) { 
		Map<String,UnitCoverage> kmap = new TreeMap<String,UnitCoverage>();
		for(String k : covermap.keySet()) {
			if(c.covermap.containsKey(k)) { 
				kmap.put(k, covermap.get(k).intersection(c.covermap.get(k)));
			}
		}		
		return new Coverage(kmap);
	}
	
	public void addInterval(String key, int start, int end) { 
		if(!covermap.containsKey(key)) { covermap.put(key, new UnitCoverage()); }
		covermap.get(key).addInterval(start, end);
	}
	
	public int coverage(String key, int start, int end) { 
		return covermap.containsKey(key) ? covermap.get(key).coverage(start, end) : 0;
	}
	
	public Collection<Integer[]> covered(String key, int start, int end) { 
		return covermap.containsKey(key) ? covermap.get(key).covered(start, end) : 
			new ArrayList<Integer[]>();
	}
	
	public boolean hasUnitCoverage(String k) { 
		return covermap.containsKey(k);
	}
	
	public UnitCoverage getUnitCoverage(String k) { 
		return covermap.containsKey(k) ? covermap.get(k) : null;
	}
	
	public int area() { 
		int a = 0;
		for(String key : covermap.keySet()) { 
			a += covermap.get(key).area();
		}
		return a;
	}
	
	public CoverageModel asModel() { 
		return new CoverageModel(this);
	}
	
	public static class CoverageModel extends Model {
		
		public String[] keys;
		public UnitCoverage.UnitCoverageModel[] coverages;
		
		public CoverageModel() {}
		
		public CoverageModel(Coverage c) {
			keys = c.covermap.keySet().toArray(new String[0]); 
			coverages = new UnitCoverage.UnitCoverageModel[keys.length];
			for(int i = 0; i < keys.length; i++) { 
				coverages[i] = c.covermap.get(keys[i]).asModel();
			}
		}
	}

	public Integer[] findNearestRight(String key, int location) {
		return covermap.containsKey(key) ? covermap.get(key).rightNearest(location) : null; 
	}
	
	public Integer[] findNearestLeft(String key, int location) {
		return covermap.containsKey(key) ? covermap.get(key).leftNearest(location) : null; 
	}

	public Collection<Pair<Integer[],Integer[]>> 
		findRightPairs(String key, int maxDist, Coverage coverage) {

		if(covermap.containsKey(key) && coverage.covermap.containsKey(key)) { 
			return covermap.get(key).findRightPairs(maxDist, coverage.covermap.get(key));
		} else { 
			return new ArrayList<Pair<Integer[],Integer[]>>();
		}
	}	

	public Collection<Pair<Integer[],Integer[]>> 
		findLeftPairs(String key, int maxDist, Coverage coverage) {

		if(covermap.containsKey(key) && coverage.covermap.containsKey(key)) { 
			return covermap.get(key).findLeftPairs(maxDist, coverage.covermap.get(key));
		} else { 
			return new ArrayList<Pair<Integer[],Integer[]>>();
		}
	}
}



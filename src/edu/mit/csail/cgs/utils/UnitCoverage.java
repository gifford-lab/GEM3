/*
 * Author: tdanford
 * Date: Mar 31, 2009
 */
package edu.mit.csail.cgs.utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.Iterator;
import java.util.TreeSet;

import edu.mit.csail.cgs.utils.iterators.BacktrackingIterator;
import edu.mit.csail.cgs.utils.models.Model;

/*
 * UnitCoverage provides an interface to a set of disjoint intervals, 
 * by maintaining an array of their bounds sorted by the start-coordinate. 
 *
 */
public class UnitCoverage { 
	
	public static void main(String[] args) { 
		UnitCoverage a = new UnitCoverage(), b = new UnitCoverage();
		
		a.addInterval(35, 45);
		a.addInterval(12, 20);
		a.addInterval(60, 80);
		a.addInterval(0, 10);

		b.addInterval(30, 50);
		b.addInterval(8, 16);
		b.addInterval(70, 75);
		
		System.out.println("A: " + a.toString());
		System.out.println("B: " + b.toString());
		System.out.println("A-B: " + a.subtract(b).toString());
		System.out.println("B-A: " + b.subtract(a).toString());
		System.out.println("A|B: " + a.union(b).toString());
		System.out.println("A&B: " + a.intersection(b).toString());
	}

	private ArrayList<Integer[]> units;
	
	public UnitCoverage() { 
		units = new ArrayList<Integer[]>();
	}
	
	public UnitCoverage(UnitCoverageModel m) { 
		this();
		
		for(int i = 0; i < m.unitpairs.length; i+=2) { 
			units.add(new Integer[] { m.unitpairs[i], m.unitpairs[i+1] });
		}
		signalErrorOnUnsort();
	}
	
	public void signalErrorOnUnsort() { 
		if(!checkSorted()) { 
			System.err.flush();
			System.out.flush();
			System.err.println(String.format("ERROR! not sorted!"));
			for(int i = 0; i < units.size(); i++) {
				System.err.print(String.format("%d \t", i));
				if(i < units.size()) { 
					Integer[] u = units.get(i);
					Integer[] v = i > 0 ? units.get(i-1) : null;
					boolean err = i > 0 && v[1] >= u[0];
					String errstr = err ? "*" : " ";
					System.err.print(String.format(" %s %d,%d", errstr, u[0], u[1]));
				}
				System.err.println();
			}
			System.err.println();
			throw new IllegalStateException();
		}		
	}
	
	UnitCoverage(Collection<Integer[]> us) { 
		units = new ArrayList<Integer[]>(us);
		signalErrorOnUnsort();
	}
	
	public UnitCoverage(UnitCoverage c) { 
		units = new ArrayList<Integer[]>();
		for(int i = 0; i < c.units.size(); i++) { 
			units.add(c.units.get(i).clone());
		}
		signalErrorOnUnsort();
	}
	
	public Collection<Interval> compareToOverlapSum(OverlapSum s, int thresh) { 
		Collection<Interval> intvs = s.collect(thresh);
		ArrayList<Interval> failing = new ArrayList<Interval>();
		TreeSet<Integer> matched = new TreeSet<Integer>();
		
		for(Interval intv : intvs) { 
			int istart = intv.start, iend = intv.end;
			Integer[] rquery = range(istart, iend);
			
			if(rquery == null || rquery[0] < rquery[1] || 
				units.get(rquery[0])[0] != istart || 
				units.get(rquery[1])[1] != iend) { 
			
				failing.add(intv);
			} else { 
				matched.add(rquery[0]);
			}
		}
		
		for(int i = 0; i < units.size(); i++) { 
			if(!matched.contains(i)) { 
				failing.add(new Interval(units.get(i)[0], units.get(i)[1]));
			}
		}
		
		return failing;
	}
	
	public int size() { return units.size(); }
	
	public Iterator<Integer[]> units() { return units.iterator(); }
	
	public UnitCoverage subtract(UnitCoverage c) { 
		ArrayList<Integer[]> newUnits = new ArrayList<Integer[]>();
		
		BacktrackingIterator<Integer[]> these = 
			new BacktrackingIterator<Integer[]>(units.iterator());
		BacktrackingIterator<Integer[]> those = 
			new BacktrackingIterator<Integer[]>(c.units.iterator());

		while(these.hasNext()) { 
			if (!those.hasNext()) { 
				newUnits.add(these.next());
			} else { 
				Integer[] u1 = these.next(), u2 = those.next();
				
				if(overlaps(u1, u2)) {
					Integer[] lft = new Integer[] { u1[0], u2[0]-1 };
					Integer[] rgt = new Integer[] { u2[1]+1, u1[1] };
					
					if(rgt[1] >= rgt[0]) { these.addNext(rgt); }
					if(lft[1] >= lft[0]) { these.addNext(lft); }

					if(u2[1] > u1[1]) { 
						those.addNext(u2);
					}
					
				} else { 
					if(u1[1] < u2[0]) {
						newUnits.add(u1);
						those.addNext(u2);
					} else { 
						these.addNext(u1);
					}
				}
			}
		}
		
		return new UnitCoverage(newUnits);		
	}
	
	public UnitCoverage union(UnitCoverage c) { 
		ArrayList<Integer[]> newUnits = new ArrayList<Integer[]>();
		
		BacktrackingIterator<Integer[]> these = 
			new BacktrackingIterator<Integer[]>(units.iterator());
		BacktrackingIterator<Integer[]> those = 
			new BacktrackingIterator<Integer[]>(c.units.iterator());
			
		while(these.hasNext() || those.hasNext()) { 
			if(!these.hasNext()) { 
				newUnits.add(those.next()); 
			} else if (!those.hasNext()) { 
				newUnits.add(these.next());
			} else { 
				Integer[] u1 = these.next(), u2 = those.next();
				
				if(overlaps(u1, u2)) { 
					if(u1[1] <= u2[1]) { 
						those.addNext(union(u1, u2));
					} else { 
						these.addNext(union(u1, u2));
					}
				} else { 
					if(u1[1] < u2[0]) { 
						newUnits.add(u1);
						those.addNext(u2);
					} else { 
						newUnits.add(u2);
						these.addNext(u1);
					}
				}
			}
		}
		
		return new UnitCoverage(newUnits);		
	}

	public UnitCoverage intersection(UnitCoverage c) { 
		ArrayList<Integer[]> newUnits = new ArrayList<Integer[]>();
		
		BacktrackingIterator<Integer[]> these = 
			new BacktrackingIterator<Integer[]>(units.iterator());
		BacktrackingIterator<Integer[]> those = 
			new BacktrackingIterator<Integer[]>(c.units.iterator());

		while(these.hasNext() && those.hasNext()) { 
			Integer[] u1 = these.next(), u2 = those.next();
			if(overlaps(u1, u2)) {
				Integer[] inter = intersection(u1, u2);
				newUnits.add(inter);
				
				if (u1[1] >= u2[1]) {
					these.addNext(u1);
				} else { 
					those.addNext(u2);
				}
			} else if (u1[1] < u2[0]) { 
				those.addNext(u2);
			} else { 
				these.addNext(u1);
			}
		}
		
		return new UnitCoverage(newUnits);		
	}
	
	public String toString() { 
		StringBuilder sb = new StringBuilder();
		sb.append("[ ");
		for(Integer[] u : units) { 
			sb.append(String.format("%d,%d ", u[0], u[1]));
		}
		sb.append("]");
		return sb.toString();
	}
	
	public int addInterval(int start, int end) {
		Integer[] r = range(start, end);
		//System.out.println(String.format("%d,%d -> %s", start, end, toString()));

		if(r == null) { 
			//System.out.println("\tNo overlap");
			insertInterval(start, end);
		} else { 
			//System.out.println(String.format("\nelmt: %d,%d", start, end));
			//System.out.println(String.format("\tover: r[0]=%d, r[1]=%d", r[0], r[1]));
			Integer[] nunit = new Integer[] { start, end };
			for(int i = r[1]; i > r[0]; i--) {
				if(overlaps(nunit, units.get(i))) { 
					Integer[] u = units.get(i);
					//System.out.println(String.format("\t\t%d - u: %d-%d", i, u[0], u[1]));
					nunit = union(u, nunit);
					units.remove(i);
					//System.out.println(String.format("\t\tnunit: %d-%d", nunit[0], nunit[1]));
				}
			}
			//System.out.println(String.format("\t\t%d - u: %d-%d", r[0], units.get(r[0])[0], units.get(r[0])[1]));
			nunit = union(units.get(r[0]), nunit);
			//System.out.println(String.format("\t\tnunit: %d-%d", nunit[0], nunit[1]));
			units.set(r[0], nunit);
		}

		signalErrorOnUnsort();
		return 1;
	}
	
	private void insertInterval(int start, int end) { 
		Integer left = leftBinarySearch(start);
		if(left == null) { 
			//System.out.println(String.format("\t%d,%d -> ?", start, end));
			units.add(0, new Integer[] { start, end });
		} else { 
			//System.out.println(String.format("\t%d,%d -> %d", start, end, left));
			units.add(left+1, new Integer[] { start, end });
		}
		signalErrorOnUnsort();
	}
	
	/*
	 * Finds either the interval that contains the given location, or the 
	 * interval that occurs immediately to that location's *left*.  (Or null, 
	 * if no such interval exists.)  
	 */
	private Integer leftBinarySearch(int loc) { 
		if(units.isEmpty()) { return null; }
		int lower = 0, upper = units.size()-1;
		if(units.get(upper)[1] < loc) { return upper; }
		if(units.get(lower)[0] > loc) { return null; }
		if(contains(units.get(lower), loc)) { return lower; }
		if(contains(units.get(upper), loc)) { return upper; }
		
		while(upper-lower > 1) { 
			//System.out.println(String.format("<(%d %d) [%d-%d / %d-%d]", lower, upper, units.get(lower)[0], units.get(lower)[1], units.get(upper)[0], units.get(upper)[1]));
			int middle = (upper+lower)/2;
			if(contains(units.get(middle), loc)) { 
				return middle;
			} else if(units.get(middle)[0] < loc) { 
				lower = middle;
			} else { 
				upper = middle;
			}
		}
		
		//System.out.println(String.format("left: %d -> %d", loc, lower));
		
		return lower;
	}

	private Integer rightBinarySearch(int loc) { 
		if(units.isEmpty()) { return null; }
		int lower = 0, upper = units.size()-1;
		if(units.get(lower)[0] > loc) { return lower; }
		if(units.get(upper)[1] < loc) { return null; }
		if(contains(units.get(lower), loc)) { return lower; }
		if(contains(units.get(upper), loc)) { return upper; }
		
		while(upper-lower > 1) {
			//System.out.println(String.format(">(%d %d) [%d-%d / %d-%d]", lower, upper, units.get(lower)[0], units.get(lower)[1], units.get(upper)[0], units.get(upper)[1]));
			int middle = (upper+lower)/2;
			if(contains(units.get(middle), loc)) { 
				return middle;
			} else if(units.get(middle)[1] > loc) { 
				upper = middle;
			} else { 
				lower = middle;
			}
		}

		//System.out.println(String.format("right: %d -> %d", loc, upper));

		return upper;
	}
	
	/*
	 *  Returns a pair of indices, into the unit array: [i, j]
	 *  i and j define, inclusively, a range of Intervals that are overlapped
	 *  by the given start and end query.  
	 *  
	 *  If no intervals are overlapped, a 'null' is returned.   
	 */
	private Integer[] range(int start, int end) {
		if(units.isEmpty() || start > units.get(units.size()-1)[1] || end < units.get(0)[0]) {
			System.out.println("Returning from range() early.");
			return null;
		}
		
		Integer lower = rightBinarySearch(start);
		Integer upper = lower;
		while(upper < units.size() && units.get(upper)[0] <= end) { 
			upper += 1;
		}
		System.out.println(String.format("%d,%d -> [%d,%d]", start, end, lower, upper));
		
		if(lower >= upper) { 
			return null;
		}
		
		return new Integer[] { lower, upper-1 };
	}
	
	public Integer[] rightNearest(Integer pos) { 
		Integer upper = rightBinarySearch(pos);
		return upper != null ? units.get(upper) : null;
	}
	
	public Integer[] leftNearest(Integer pos) { 
		Integer lower = leftBinarySearch(pos);
		return lower != null ? units.get(lower) : null;
	}
	
	/**
	 * Calculates the total area covered by all units.
	 */
	public int area() { 
		int a = 0;
		for(Integer[] unit : units) { 
			a += (unit[1] - unit[0] + 1);
		}
		return a;
	}

	/**
	 * Calculates the covered area within the given query.  
	 * @param start the left-most coordinate of the query
	 * @param end the right-most coordinate of the query.  
	 * @return
	 */
	public int coverage(int start, int end) { 
		int coverage = 0;
		Integer[] query = new Integer[] { start, end };
		Integer[] range = range(start, end);
		if(range != null) { 
			for(int i = range[0]; i<= range[1]; i++) { 
				Integer[] unit = units.get(i);
				if(!overlaps(query, unit)) {
					if(!checkSorted()) { 
						System.err.println(String.format("ERROR! units not sorted."));
					}
					System.err.println(String.format("ERROR: %d,%d not in %d,%d", unit[0], unit[1], query[0], query[1]));
				}
				int s = Math.max(unit[0], start);
				int e = Math.min(unit[1], end);
				if(e >= s) { 
					int w = e - s + 1;
					coverage += w;
				}
			}
		}
		return coverage;
	}
	
	public boolean checkSorted() { 
		for(int i = 1; i < units.size(); i++) { 
			Integer[] u = units.get(i-1), v = units.get(i);
			if(u[1] >= v[0]) { 
				return false;
			}
		}
		return true;
	}
	
	/**
	 * Returns a list of all regions *within the given query* 
	 * which are covered.
	 * 
	 * @param start
	 * @param end
	 * @return
	 */
	public Collection<Integer[]> covered(int start, int end) { 
		ArrayList<Integer[]> lst = new ArrayList<Integer[]>();
		Integer[] query = new Integer[] { start, end };
		Integer[] range = range(start, end);
		if(range != null) { 
			for(int i = range[0]; i<= range[1]; i++) { 
				lst.add(intersection(units.get(i), query));
			}
		}
		return lst;
	}
	
	public Iterator<Integer[]> covered() { 
		return units.iterator();
	}
	
	/**
	 * Returns true if there is any overlap with the given query.  
	 * 
	 * @param start
	 * @param end
	 * @return
	 */
	public boolean hasOverlap(int start, int end) { 
		return range(start, end) != null;
	}

	public boolean isContained(int start, int end) {
		Integer[] r = range(start, end);
		return r != null && r[0] == r[1] && 
			contains(units.get(r[0]), new Integer[] { start, end });
	}

	private static Integer[] intersection(Integer[] r1, Integer[] r2) { 
		return new Integer[] { Math.max(r1[0], r2[0]), Math.min(r1[1], r2[1]) };		
	}
	
	private static Integer[] union(Integer[] r1, Integer[] r2) { 
		return new Integer[] { Math.min(r1[0], r2[0]), Math.max(r1[1], r2[1]) };
	}
	
	private static boolean equal(Integer[] r1, Integer[] r2) { 
		return r1[0].equals(r2[0]) && r1[1].equals(r2[1]);
	}
	
	private static String string(Integer[] r) { 
		return String.format("%d,%d", r[0], r[1]);
	}
	
	private static Integer[] subtract(Integer[] r1, Integer[] r2) {
		if(equal(r1, r2)) { 
			return null; 
		} else if(!overlaps(r1, r2)) {
			return r1; 
		} else { 
			if(r1[0] < r2[0]) { 
				return new Integer[] { r1[0], Math.min(r1[1], r2[1])-1 };				
			} else { 
				return new Integer[] { Math.max(r1[0], r2[0])+1, r1[1] };								
			}
		}
	}
	
	private static Integer[] union(Integer[]... rs) { 
		if(rs.length == 0) { 
			throw new IllegalArgumentException();
		} else if (rs.length == 1) { 
			return rs[0];
		} else { 
			return union(rs[0], union(ArrayUtils.tail(rs)));
		}
	}
	
	private static boolean overlaps(Integer[] r1, Integer[] r2) { 
		return contains(r1, r2[0]) || 
			contains(r2, r1[0]);
	}
	
	private static boolean contains(Integer[] r1, Integer[] r2) { 
		return r1[0] <= r2[0] && r1[1] >= r2[1];
	}
	
	private static boolean contains(Integer[] r, Integer l) { 
		return r[0] <= l && r[1] >= l;
	}

	public static class UnitCoverageModel extends Model { 
		public Integer[] unitpairs; 
		
		public UnitCoverageModel() {}
		
		private UnitCoverageModel(Collection<Integer[]> us) { 
			unitpairs = new Integer[us.size()*2];
			int i = 0; 
			for(Integer[] p : us) { 
				if(p.length != 2) { throw new IllegalArgumentException(); }
				unitpairs[i] = p[0];
				unitpairs[i+1] = p[1];
				i += 2;
			}
		}
	}
	
	public UnitCoverageModel asModel() { 
		signalErrorOnUnsort();
		return new UnitCoverageModel(units);
	}

	public Collection<Pair<Integer[], Integer[]>> findLeftPairs(int maxDist, UnitCoverage coverage) {
		ArrayList<Pair<Integer[],Integer[]>> pairs = new ArrayList<Pair<Integer[],Integer[]>>();
		for(Integer[] right : units) { 
			Integer[] left = coverage.leftNearest(right[0]-1);
			if(left != null && (contains(left, right[0]) || 
					(!hasOverlap(left[1]+1, right[0]-1) && right[0] - left[1] <= maxDist))) { 
				pairs.add(new Pair<Integer[],Integer[]>(right, left));
			}
		}
		return pairs;
	}

	public Collection<Pair<Integer[], Integer[]>> findRightPairs(int maxDist, UnitCoverage coverage) {
		ArrayList<Pair<Integer[],Integer[]>> pairs = new ArrayList<Pair<Integer[],Integer[]>>();
		for(Integer[] left : units) { 
			Integer[] right = coverage.rightNearest(left[1]+1);
			if(right != null && (contains(right, left[1]) || 
					(!hasOverlap(left[1]+1, right[0]-1) && right[0] - left[1] <= maxDist))) { 
				pairs.add(new Pair<Integer[],Integer[]>(left, right));
			}
		}
		return pairs;
	}
}

class PairComparator implements Comparator<Integer[]> { 
	public int compare(Integer[] f, Integer[] s) { 
		if(f[0] < s[0]) { return -1; }
		if(f[0] > s[1]) { return 1; }
		return f[1].compareTo(s[1]);
	}
}

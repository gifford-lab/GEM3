/*
 * Author: tdanford
 * Date: Nov 16, 2008
 */
package edu.mit.csail.cgs.utils;

import java.util.*;

/**
 * Collapses a set of Interval objects into non-overlapping (contiguous) regions.
 * 
 * @author tdanford
 *
 * @param <X>
 */
public class NonOverlappingIntervalSet<X> {

	private ArrayList<Interval<Set<X>>> blocks;
	
	public NonOverlappingIntervalSet() { 
		blocks = new ArrayList<Interval<Set<X>>>();
	}
	
	public Collection<Interval<Set<X>>> blocks() { return blocks; }

	public Collection<Interval<Set<X>>> overlapping(Interval<X> intv) { 
		LinkedList<Interval<Set<X>>> oblocks = new LinkedList<Interval<Set<X>>>();
		for(int i = 0; i < blocks.size(); i++) {  
			if(blocks.get(i).overlaps(intv)) { 
				oblocks.add(blocks.get(i));
			}
		}
		return oblocks;
	}
	
	public void addIntervals(Iterator<Interval<X>> itr) { 
		while(itr.hasNext()) {
			 addInterval(itr.next());
		}
	}

	public void addInterval(Interval<X> intv) { 
		LinkedList<Integer> overlapping = overlappingIndices(intv);
		if(overlapping.isEmpty()) { 
			HashSet<X> set = new HashSet<X>();
			set.add(intv.data);
			blocks.add(new Interval<Set<X>>(intv.start, intv.end, set));
		}
		else if (overlapping.size() == 1) { 
			Integer overIdx = overlapping.getFirst();
			Interval<Set<X>> overlap = blocks.get(overIdx);
			overlap.data.add(intv.data);
			overlap.start = Math.min(intv.start, overlap.start);
			overlap.end= Math.max(intv.end, overlap.end);
		}
		else { 
			Set<X> total = new HashSet<X>();
			int start = intv.start;
			int end = intv.end;
			total.add(intv.data);
			
			for(Integer idx : overlapping) { 
				Interval<Set<X>> overlap = blocks.get(idx);
				start = Math.min(start, overlap.start);
				end = Math.max(end, overlap.end);
				total.addAll(overlap.data);
			}
			
			Iterator<Interval<Set<X>>> itr = blocks.iterator();
			while(itr.hasNext()) { 
				Interval<Set<X>> over = itr.next();
				if(over.overlaps(intv)) { 
					itr.remove();
				}
			}
			
			blocks.add(new Interval<Set<X>>(start, end, total));
		}
	}
	
	private LinkedList<Integer> overlappingIndices(Interval<X> intv) { 
		LinkedList<Integer> oblocks = new LinkedList<Integer>();
		for(int i = 0; i < blocks.size(); i++) {  
			if(blocks.get(i).overlaps(intv)) { 
				oblocks.add(i);
			}
		}
		return oblocks;
	}
}

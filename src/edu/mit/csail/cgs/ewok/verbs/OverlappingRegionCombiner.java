package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;
import edu.mit.csail.cgs.datasets.general.Region;

public class OverlappingRegionCombiner implements Iterator<Region> {

    private Iterator<Region> iter;
    
    public OverlappingRegionCombiner(Collection<? extends Region> rgs, int dummy) { 
        Vector<Region> regions = new Vector<Region>();
        for(Region r : rgs) {
            boolean foundOverlap = false;
            for(int i = 0; i < regions.size(); i++) { 
                if(r.overlaps(regions.get(i))) { 
                    int start = Math.min(r.getStart(), regions.get(i).getStart());
                    int end = Math.max(r.getEnd(), regions.get(i).getEnd());
                    regions.set(i, new Region(r.getGenome(), r.getChrom(), start, end));
                    foundOverlap = true;
                    break;
                }
            }
            
            if(!foundOverlap) { 
                regions.add(r);
            }
        }
        
        iter = regions.iterator();
    }

    public OverlappingRegionCombiner(Collection<Iterator<? extends Region>> collection) {
        ArrayList<Region> regions = new ArrayList<Region>();
        for (Iterator<? extends Region> iter : collection) {
            while (iter.hasNext()) {
                Region r = iter.next();
                boolean found = false;
                for (int i = 0; i < regions.size(); i++) {
                    if (r.overlaps(regions.get(i))) {
                        regions.set(i,r.combine(regions.get(i)));
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    regions.add(r);
                }
            }
        }
        iter = regions.iterator();
    }
    public boolean hasNext() {
        return iter.hasNext();
    }
    public Region next() {
        return iter.next();
    }

    public void remove() {
    	throw new UnsupportedOperationException();
	}

}

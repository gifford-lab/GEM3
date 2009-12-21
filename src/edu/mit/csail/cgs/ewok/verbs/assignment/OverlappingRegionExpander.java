/*
 * Created on Jan 11, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.ewok.verbs.assignment;

import java.util.*;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.*;

public class OverlappingRegionExpander<X extends Region> implements Expander<Region,X> {

    private HashMap<String,Set<X>> regions;
    
    public OverlappingRegionExpander(Collection<X> regs) { 
        regions = new HashMap<String,Set<X>>();
        for(X r : regs) { 
        	if(!regions.containsKey(r.getChrom())) { 
        		regions.put(r.getChrom(), new HashSet<X>());
        	}
        	regions.get(r.getChrom()).add(r);
        }
    }
    
    public void addRegion(X r) { 
        String chrom = r.getChrom();
        if(!regions.containsKey(chrom)) { regions.put(chrom, new HashSet<X>()); }
        regions.get(chrom).add(r);
    }

    public Iterator<X> execute(Region a) {
        LinkedList<X> overlaps = new LinkedList<X>();
        
        if(regions.containsKey(a.getChrom())) { 
        	for(X r : regions.get(a.getChrom())) { 
        		if(a.overlaps(r)) { 
        			overlaps.addLast(r);
        		}
        	}
        }
        
        return overlaps.iterator();
    }
}

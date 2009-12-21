/*
 * Created on Oct 25, 2006
 */
package edu.mit.csail.cgs.conservation;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.verbs.Filter;

/**
 * @author tdanford
 */
public class OverlappingRegionFilter implements Filter<Region,Region> {
    
    private Region target;

    public OverlappingRegionFilter(Region t) {
        target = t;
    }
    
    public void setTarget(Region t) { target = t; }

    public Region execute(Region a) {
        return (target.overlaps(a)) ? a : null;
    }

}

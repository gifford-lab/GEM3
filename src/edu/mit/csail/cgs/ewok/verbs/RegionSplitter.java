/*
 * Created on Mar 13, 2006
 */
package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.nouns.*;

/**
 * @author tdanford
 */
public class RegionSplitter implements Expander<Region,Region> {
    
    private int window, overlap;
    
    public RegionSplitter(int w) {
        if(w <= 0) { throw new IllegalArgumentException(); }
        window = w;
        overlap = Math.max(0, window-1);
    }

    public RegionSplitter(int w, int over) {
        if(w <= 0) { throw new IllegalArgumentException(); }
        if(over <= 0) { throw new IllegalArgumentException(); }
        
        window = w;
        overlap = over;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.ewok.verbs.Filter#execute(null)
     */
    public Iterator<Region> execute(Region region) {
        return new LazySplittingItr(region);
    }

    private class LazySplittingItr implements Iterator<Region> {
        
        private Region region;
        private int location;
        
        public LazySplittingItr(Region r) { 
            region = r;
            location = region.getStart();
        }

        /* (non-Javadoc)
         * @see java.util.Iterator#hasNext()
         */
        public boolean hasNext() {
            return location <= region.getEnd();
        }

        /* (non-Javadoc)
         * @see java.util.Iterator#next()
         */
        public Region next() {
            Region nextR = new Region(region.getGenome(),region.getChrom(), location, location+window-1);
            location += (window - overlap);
            return nextR;
        }

        /* (non-Javadoc)
         * @see java.util.Iterator#remove()
         */
        public void remove() {
            throw new UnsupportedOperationException();
        }  
        
    }

}

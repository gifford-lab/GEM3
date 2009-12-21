/*
 * Created on Mar 27, 2006
 */
package edu.mit.csail.cgs.viz;

import java.util.*;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.utils.*;

/**
 * @author tdanford
 * 
 * NonOverlappingLayout is a non-painting utility class, which is used to 
 * annotate a set of (possibly overlapping) Region objects with "tracks."  If 
 * upon painting each of the "tracks" is drawn on a separate line, then none of the 
 * Region objects will (visually) overlap.  This is used (primarily) in the 
 * GenePainter class.
 * 
 * The two main interface functions are setRegions() and getTrack().  
 * To use the NonOverlappingLayout, create an object of the class.  Then, every
 * time you want to layout a set of Regions, call setRegions() with that 
 * Region collection.  Then, for each Region in the collection, you can query to 
 * ask the track that it was laid out on.
 */
public class NonOverlappingLayout<X extends Region> {
    
    private static class RegionComparator implements Comparator<Region> { 
        public int compare(Region r1, Region r2) { 
            if(!r1.getChrom().equals(r2.getChrom())) { return r1.getChrom().compareTo(r2.getChrom()); }
            if(r1.getStart() < r2.getStart()) { return -1; }
            if(r1.getStart() > r2.getStart()) { return 1; }
            if(r1.getEnd() < r2.getEnd()) { return -1; }
            if(r2.getEnd() > r2.getEnd()) { return 1; }
            return 0;
        }
    }
    
    private static class LayoutTrack<Y extends Region> { 
        
        private Vector<Y> regions;
        private int index;
        
        public LayoutTrack(int ind) { 
            regions = new Vector<Y>();
            index = ind;
        }
        
        public int getIndex() { return index; }
        public Vector<Y> getRegions() { return regions; }
        
        public void addRegion(Y r) { regions.add(r); }
        
        public boolean acceptsRegion(Y r) { 
            for(Y rt : regions) { 
                if(rt.overlaps(r)) { return false; }
            }
            return true;
        }
    }
    
    private static Comparator<Region> comp;
    
    static { 
        comp = new RegionComparator();
    }
    
    private Region[] regions;
    private Vector<LayoutTrack<X>> tracks;
    private Map<Region,LayoutTrack> trackMap;

    public NonOverlappingLayout() {
        regions = null;
        tracks = new Vector<LayoutTrack<X>>();
        trackMap = new HashMap<Region,LayoutTrack>();
    }
    
    public void setRegions(Iterator<X> rs) { 
    	ArrayList<X> rrs = new ArrayList<X>();
    	while(rs.hasNext()) { 
    		rrs.add(rs.next());
    	}
    	setRegions(rrs);
    }
    
    public void setRegions(Collection<X> rs) {
        tracks.clear();
        trackMap.clear();

        regions = rs.toArray(new Region[rs.size()]);
        Arrays.sort(regions, comp);
        
        doLayout();
    }
    
    public Collection<X> getRegions() { 
    	ArrayList<X> rrs = new ArrayList<X>();
    	for(LayoutTrack<X> track : tracks) {
    		rrs.addAll(track.getRegions());
    	}
    	return rrs;
    }
    
    private void doLayout() {
        for(int i = 0; i < regions.length; i++) { 
            X r = (X)(regions[i]);
            int currentTrack = 0;
            
            while(currentTrack < tracks.size() && 
            		!tracks.get(currentTrack).acceptsRegion(r)) { 
                currentTrack += 1;
            }
            
            if(currentTrack >= tracks.size()) { 
                tracks.add(new LayoutTrack(tracks.size()));
            }
            
            tracks.get(currentTrack).addRegion(r);
            trackMap.put(r, tracks.get(currentTrack));
        }
    }
    
    public boolean hasTrack(Region r) { return trackMap.containsKey(r); }
    public int getNumTracks() { return tracks.size(); }
    public int getTrack(Region r) { return trackMap.get(r).getIndex(); }
}

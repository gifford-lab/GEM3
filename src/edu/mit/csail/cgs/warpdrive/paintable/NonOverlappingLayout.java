/*
 * Created on Mar 27, 2006
 */
package edu.mit.csail.cgs.warpdrive.paintable;

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
    
    private static class LayoutTrack { 
        
        private Vector<Region> regions;
        private int index;
        
        public LayoutTrack(int ind) { 
            regions = new Vector<Region>();
            index = ind;
        }
        
        public int getIndex() { return index; }
        public Vector<Region> getRegions() { return regions; }
        
        public void addRegion(Region r) { regions.add(r); }
        
        public boolean acceptsRegion(Region r) { 
            for(Region rt : regions) { 
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
    private Vector<LayoutTrack> tracks;
    private Map<Region,LayoutTrack> trackMap;
    /* never return more than this many tracks.  If the nonoverlapping layout would
       have needed more, then everything past this just gets dumped into the
       last track.  That's ugly, but it lets you avoid memory explosions.
    */
    private int maxTracks;

    public NonOverlappingLayout() {
        regions = null;
        tracks = new Vector<LayoutTrack>();
        trackMap = new HashMap<Region,LayoutTrack>();
        maxTracks = 1000;
    }
    public int getMaxTracks() {return maxTracks;}
    public void setMaxTracks(int m) {
        maxTracks = m;
        if (m < tracks.size()) {
            doLayout();
        }

    }
    public void setRegions(Collection<X> rs) {
        tracks.clear();
        trackMap.clear();

        regions = rs.toArray(new Region[rs.size()]);
        boolean sorted = true;
        for (int i = 0; i < regions.length - 1; i++) {
            if (regions[i].compareTo(regions[i+1]) > 0) {
                sorted = false;
                break;
            }
        }

        if (!sorted) {
            Arrays.sort(regions, comp);     
        }
        doLayout();
    }
    
    private void doLayout() {
        ArrayList<Integer> lastEnd = new ArrayList<Integer>();
        for(int i = 0; i < regions.length; i++) { 
            Region r = regions[i];
            int currentTrack = 0;
            
            while(currentTrack < tracks.size() && 
                  (lastEnd.get(currentTrack) > r.getStart()) &&
                  currentTrack < maxTracks) {
                currentTrack += 1;
            }
        
            if(currentTrack >= tracks.size()) { 
                tracks.add(new LayoutTrack(tracks.size()));
                lastEnd.add(0);
            }
            
            tracks.get(currentTrack).addRegion(r);
            trackMap.put(r, tracks.get(currentTrack));
            lastEnd.set(currentTrack, r.getEnd());
        }
    }
    
    public boolean hasTrack(Region r) { return trackMap.containsKey(r); }
    public int getNumTracks() { return tracks.size(); }
    public int getTrack(Region r) { return trackMap.get(r).getIndex(); }
}

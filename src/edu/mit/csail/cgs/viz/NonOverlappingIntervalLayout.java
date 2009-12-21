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
 */
public class NonOverlappingIntervalLayout<X> {
    
    private static class LayoutTrack<Y> { 
        
        private Vector<Interval<Y>> intervals;
        private int index;
        
        public LayoutTrack(int ind) { 
            intervals = new Vector<Interval<Y>>();
            index = ind;
        }
        
        public int getIndex() { return index; }
        public Vector<Interval<Y>> getIntervals() { return intervals; }
        
        public void addInterval(Interval<Y> i) { intervals.add(i); }
        
        public boolean acceptsInterval(Interval<Y> r) { 
            for(Interval<Y> rt : intervals) { 
                if(rt.overlaps(r)) { return false; }
            }
            return true;
        }
    }
    
    private Interval[] regions;
    private Vector<LayoutTrack<X>> tracks;
    private Map<Interval<X>,LayoutTrack<X>> trackMap;

    public NonOverlappingIntervalLayout() {
        regions = null;
        tracks = new Vector<LayoutTrack<X>>();
        trackMap = new HashMap<Interval<X>,LayoutTrack<X>>();
    }
    
    public NonOverlappingIntervalLayout(Iterator<Interval<X>> intvs) { 
    	this();
    	ArrayList<Interval<X>> intervals = new ArrayList<Interval<X>>();
    	while(intvs.hasNext()) { intervals.add(intvs.next()); }
    	regions = intervals.toArray(new Interval[0]); 
    	Arrays.sort(regions);
    	
    	doLayout();
    }
    
    private void clearTracks() { 
    	tracks.clear();
    	trackMap.clear();
    }
    
    public void setRegions(Collection<Interval<X>> rs) {
    	clearTracks();

        regions = rs.toArray(new Interval[rs.size()]);
        Arrays.sort(regions);
        
        doLayout();
    }
    
    public void addInterval(Interval<X> intv) {
    	clearTracks();

    	int oldLength = regions != null ? regions.length : 0;
    	Interval[] newr = new Interval[oldLength+1];
    	for(int i = 0; i < oldLength; i++) { newr[i] = regions[i]; }
    	newr[oldLength] = intv;
    	regions = newr;
    	Arrays.sort(regions);
    	
    	doLayout();
    }
    
    private void doLayout() {
        for(int i = 0; i < regions.length; i++) { 
            Interval<X> r = regions[i];
            int currentTrack = 0;
            
            while(currentTrack < tracks.size() && !tracks.get(currentTrack).acceptsInterval(r)) { 
                currentTrack += 1;
            }
            
            if(currentTrack >= tracks.size()) { 
                tracks.add(new LayoutTrack(tracks.size()));
            }
            
            tracks.get(currentTrack).addInterval(r);
            trackMap.put(r, tracks.get(currentTrack));
        }
    }
    
    public Iterator<Interval> iterator() { 
    	return ArrayUtils.asIterator(regions); 
    }
    
    public boolean hasTrack(Interval<X> r) { return trackMap.containsKey(r); }
    public int getNumTracks() { return tracks.size(); }
    public int getTrack(Interval<X> r) { return trackMap.get(r).getIndex(); }

	public void clear() {
		regions = null;
		clearTracks();
	}
}

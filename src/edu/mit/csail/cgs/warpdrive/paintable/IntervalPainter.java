package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.awt.font.FontRenderContext;
import java.awt.font.TextLayout;
import java.awt.font.LineMetrics;
import java.util.*;
import java.text.*;

import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.viz.DynamicAttribute;
import edu.mit.csail.cgs.datasets.chippet.ChipPetDatum;
import edu.mit.csail.cgs.datasets.chippet.ChipPetLoader;
import edu.mit.csail.cgs.datasets.expression.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.warpdrive.model.*;

public class IntervalPainter extends RegionPaintable  {
	
	private static NumberFormat nf;
	
	static {
		nf = DecimalFormat.getInstance();
		nf.setMaximumFractionDigits(2);
		nf.setMinimumFractionDigits(2);
	}

    private ScoreTrackModel model;
    private DynamicAttribute attrib;
    private static java.util.List configurationFields = null;
    private IntervalLayout layout;
    private PaintableProperties props;

    public IntervalPainter (ScoreTrackModel m) {
        super();
        model = m;
        props = new PaintableProperties();
        model.addEventListener(this);
        attrib = DynamicAttribute.getGlobalAttributes();
        layout = new IntervalLayout();
    }
    public PaintableProperties getProperties() {return props;}
    
    public void cleanup() { 
        super.cleanup();
        model.removeEventListener(this);
    }

    public java.util.List<String> configurationKeyOrder() {return configurationFields;}

    public void paintItem(Graphics2D g, 
                          int x1, int y1, 
                          int x2, int y2) {
    	int w = x2 - x1, h = y2 - y1;
        
        if(!model.isReady()) { 
            g.setColor(Color.black);
            g.drawString("Loading... Please Wait.", x1 + w / 2, y1 + h / 2);
            return;
        }
        
        Region region = model.getRegion();
        
        Stroke oldStroke = g.getStroke();
        g.setStroke(new BasicStroke((float)2));

        Iterator<ChipPetDatum> itr = model.getResults();
        LinkedList<Interval<Integer>> intvs = new LinkedList<Interval<Integer>>();
        while(itr.hasNext()) {
            ChipPetDatum intv = itr.next();
            for(int i = 0; i < intv.getCount(); i++) { 
                intvs.add(new Interval<Integer>(intv.getStart(), intv.getEnd(), i));
            }
        }
        
        layout.setRegions(intvs);
        clearLabels();
        
        int nt = layout.getNumTracks() + 1;
        int trackHeight = (int)Math.round((double)h / (double)nt);
        int ah = layout.getNumTracks() * trackHeight;


        
        g.setColor(Color.black);
        for(Interval<Integer> intv : intvs) { 
            int track = layout.getTrack(intv);
            int offset = intv.start - region.getStart();
            
            double foffset = (double)offset / (double)region.getWidth();
            double fwidth = (double)intv.getWidth() / (double)region.getWidth();
            
            int int_x = x1 + (int)Math.round(foffset * (double)w);
            int int_y = y1 + ah - (trackHeight * (track+1));
            int int_w = Math.max((int)Math.round(fwidth * (double)w), 1);
            int int_h = trackHeight;
            
            //System.out.println(String.format("%d x,y:%d,%d w,h:%d,%d", track, int_x, int_y, int_w, int_h));
            g.drawRect(int_x, int_y, int_w, int_h);
            addLabel(int_x, int_y, int_w, int_h, "Count: " + intv.data);
        }
        
        g.setStroke(oldStroke);
        //System.out.println("Interval painter painted: " + intvs.size() + " intervals.");
    }
    
    public int getMaxVertSpace() { 
        int numTracks = layout.getNumTracks();
        return Math.min(Math.max(40,numTracks * 12),120);
    }
    
    public synchronized void eventRegistered(EventObject e) {        
        if ((e.getSource() == model) &&
            model.isReady()) {
        	
            setCanPaint(true);
            setWantsPaint(true);
            notifyListeners();
        }
    }

    public void removeEventListener(Listener<EventObject> l) {
        super.removeEventListener(l);
        if (!hasListeners()) {
            model.removeEventListener(this);
        }
    }
}

class IntervalLayout<X> {
    
    private static class IntervalComparator<Y> implements Comparator<Interval<Y>> { 
        public int compare(Interval<Y> r1, Interval<Y> r2) { 
            if(r1.getWidth() > r2.getWidth()) { return -1; }
            if(r1.getWidth() < r2.getWidth()) { return 1; }
            return r1.compareTo(r2);
        }
    }
    
    private static class LayoutTrack { 
        
        private Vector<Interval> regions;
        private int index;
        
        public LayoutTrack(int ind) { 
            regions = new Vector<Interval>();
            index = ind;
        }
        
        public int getIndex() { return index; }
        public Vector<Interval> getRegions() { return regions; }
        
        public void addRegion(Interval r) { regions.add(r); }
        
        public boolean acceptsRegion(Interval r) { 
            for(Interval rt : regions) { 
                if(rt.overlaps(r)) { return false; }
            }
            return true;
        }
    }
    
    private static Comparator<Interval> comp;
    
    static { 
        comp = new IntervalComparator();
    }
    
    private Interval[] regions;
    private Vector<LayoutTrack> tracks;
    private Map<Interval,LayoutTrack> trackMap;

    public IntervalLayout() {
        regions = null;
        tracks = new Vector<LayoutTrack>();
        trackMap = new HashMap<Interval,LayoutTrack>();
    }
    
    public void setRegions(Collection<X> rs) {
        tracks.clear();
        trackMap.clear();

        regions = rs.toArray(new Interval[rs.size()]);
        Arrays.sort(regions, comp);
        
        doLayout();
    }
    
    private void doLayout() {
        for(int i = 0; i < regions.length; i++) { 
            Interval r = regions[i];
            int currentTrack = 0;
            
            while(currentTrack < tracks.size() && !tracks.get(currentTrack).acceptsRegion(r)) { 
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
    public int getTrack(Interval r) { return trackMap.get(r).getIndex(); }
}


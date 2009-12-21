package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.util.*;
import java.util.regex.*;

import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.viz.DynamicAttribute;
import edu.mit.csail.cgs.viz.colors.ColorSet;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.Named;
import edu.mit.csail.cgs.datasets.general.SpottedProbe;
import edu.mit.csail.cgs.datasets.general.Typed;
import edu.mit.csail.cgs.datasets.general.Stranded;
import edu.mit.csail.cgs.datasets.general.NamedTypedRegion;
import edu.mit.csail.cgs.warpdrive.model.Model;
import edu.mit.csail.cgs.warpdrive.model.RegionExpanderModel;
import edu.mit.csail.cgs.ewok.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.*;


/* actually works on all regions.  Has special behaviors for those that implement Stranded or Named */
public class SpottedProbePainter extends RegionPaintable {
    
    private static ColorSet colors = new ColorSet();

    private static int sMaxTracks = 12;
    private RegionExpanderModel<SpottedProbe> model;    
    private boolean dirty;
    private int trackIndex[];
    private ArrayList<SpottedProbe> probes;
    private int tracksUsed;
    private SpottedProbeProperties props;
    private double bindingThreshold;

    public SpottedProbePainter(RegionExpanderModel<SpottedProbe> model) {
        super();
        this.model = model;
        model.addEventListener(this);
        dirty = true;
        initLabels();
        props = new SpottedProbeProperties();
        bindingThreshold = 4.0;  // probes with sd >= bindingThreshold are bound
    }
    
    public SpottedProbeProperties getProperties() {
        return props;
    }

    public int getMaxVertSpace() {
        if (model.isReady()) {
            if (dirty) {
                doLayout();
            }
            return 30 * (tracksUsed > 1 ? tracksUsed : 1);
        } else {
            return 60;
        }
    }

    public void setRegion(Region r) {
        dirty = true;
        super.setRegion(r);
    }
    
    public void cleanup() { 
        super.cleanup();
        model.removeEventListener(this);
    }

    private Color getFactorColor(String factor) { 
        return colors.getColor(factor);
    }
    
    public Color getColor(Region r) {
        return colors.getColor(this.toString());
    }
    
    public String getLabel(Region r) {
        if (r instanceof Named) {
            return ((Named)r).getName();
        } else {
            return r.toString();
        }
    }

    private void doLayout() {
        Iterator<SpottedProbe> iter = model.getResults();
        probes = new ArrayList<SpottedProbe>();
        while (iter.hasNext()) {
            SpottedProbe f = iter.next();
            probes.add(f);
        }
        int lastEnd[] = new int[sMaxTracks];
        for (int i = 0; i < sMaxTracks; i++) { lastEnd[i] = -1;}
        trackIndex = new int[probes.size()];
        tracksUsed = 0;
        for (int j = 0; j < probes.size(); j++) {
            int minpos = probes.get(j).getStart();
            int maxpos = probes.get(j).getEnd();
            int minindex = -1, biggestspace = -1000000;
            for (int i = 0; i < sMaxTracks; i++) {
                if (lastEnd[i] == -1) {
                    continue;
                }
                int space = minpos - lastEnd[i];
                if (space > biggestspace) {
                    biggestspace = space;
                    minindex = i;
                }
            }
            if (minindex == -1 || biggestspace < 20) {
                for (int i = 0; i < sMaxTracks; i++) {
                    int space = minpos - lastEnd[i];
                    if (space > biggestspace) {
                        biggestspace = space;
                        minindex = i;
                    }
                }
            }
            lastEnd[minindex] = maxpos;
            trackIndex[j] = minindex;
            if (minindex >= tracksUsed) {
                tracksUsed = minindex + 1;
            }
        }        
        dirty = false;
    }
    
    public void paintItem(Graphics2D g, 
                          int ulx, int uly, 
                          int lrx, int lry) {
        if (!model.isReady()) {
            System.err.println ("  but not ready");
            return;
        }
        if (dirty) {
            doLayout();
        }
        int lineSpace = Math.abs(lry - uly);
        int bottom = lry;
        int trackY[] = new int[tracksUsed];
        int lastx[] = new int[tracksUsed];
        Font oldFont = g.getFont();
        g.setFont(DynamicAttribute.getGlobalAttributes().getRegionLabelFont(lrx-ulx,lry-uly));
        int fontsize = g.getFont().getSize();
        Stroke oldStroke = g.getStroke();
        g.setStroke(new BasicStroke(1));
        int tracksize = (lineSpace / (tracksUsed + 1));
        for (int i = 0; i < tracksUsed; i++) { 
            trackY[i] = bottom - (i + 1) * tracksize;
            lastx[i] = -1000;
        }
        boolean drawtracklabel = props.DrawTrackLabel;
        boolean drawgenenames = props.DrawGeneNames;
        boolean alwaysDrawNames = props.AlwaysDrawNames;
        boolean drewAnything = false;

        clearLabels();
        
        // this is the set of all factors which are bound in any probe in this 
        // window.  we'll use it, at the end, to draw the legend for the colored
        // circles.
        Set<String> boundFactors = new HashSet<String>();

        for (int j = 0; j < probes.size(); j++) {
            SpottedProbe probe = probes.get(j);
            Color c = getColor(probe);
            if (c == null) {
                continue;
            } else {
                g.setColor(c);
            }
            
            String label = getLabel(probe);
            
            int min = probe.getStart();
            int max = probe.getEnd();
            int x1 = getXPos(min,
                             getRegion().getStart(), getRegion().getEnd(),
                             ulx,lrx);
            int x2 = getXPos(max,
                             getRegion().getStart(), getRegion().getEnd(),
                             ulx,lrx);
            int y = trackY[trackIndex[j]];
            g.drawRect(x1, y, x2 - x1, tracksize);

            addLabel(x1,y, x2 - x1, tracksize,label);
            
            // Drawing the probe name.
            int center = (x1+x2)/2;
            int texthalf = (label.length() * fontsize)/2;
            if (((center - texthalf > lastx[trackIndex[j]]) && drawgenenames) ||
                alwaysDrawNames){
                g.setColor(Color.black);
                g.drawString(label, center - texthalf, y-2);
            }
            
            //
            // The following section draws colored circles inside the probe regions,
            // indicating which factors are bound to that spotted probe.
            //
            
            Set<String> probeBound = probe.getBoundFactors(bindingThreshold);
            int pw = x2 - x1;
            
            // first, check to make sure that we have enough space to draw these 
            // binding circles at a reasonable size.
            if(pw <= probeBound.size()*2) {
            	
            	int bfDim = tracksize;
            	int widthNeeded = bfDim * probeBound.size();
            	int rows = 1;
            	
            	while(widthNeeded >= pw && bfDim > 2) { 
            		rows += 1;
            		bfDim = Math.max(2, tracksize / rows);
            		widthNeeded = bfDim * (probeBound.size() / rows);
            	}
            	
            	int xOffset = x1, yOffset = y;
            	
            	for(String f : probeBound) { 
            		boundFactors.add(f);
            		
            		Color fc = getFactorColor(f);
            		g.setColor(fc);
            		g.fillOval(xOffset, yOffset, bfDim, bfDim);
            		
            		xOffset += bfDim; 
            		if(xOffset >= pw-bfDim) { 
            			xOffset = x1; 
            			yOffset += bfDim;
            		}
            	}
            }
            
            drewAnything = true;
            lastx[trackIndex[j]] = x2;
        }
        
        g.setStroke(oldStroke);
        if (drawtracklabel && drewAnything) {
            g.setColor(Color.BLACK);
            g.drawString(getLabel(),ulx,lry);
        }
        g.setFont(oldFont);
    }

    public synchronized void eventRegistered(EventObject e) {        
        if (e.getSource() == model &&
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

    public String toString() {
        return "NTP " + getLabel();
    }

}


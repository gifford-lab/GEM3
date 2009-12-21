package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.util.*;

import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.viz.DynamicAttribute;
import edu.mit.csail.cgs.viz.colors.ColorSet;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.Named;
import edu.mit.csail.cgs.datasets.general.Typed;
import edu.mit.csail.cgs.datasets.general.Stranded;
import edu.mit.csail.cgs.datasets.general.NamedTypedRegion;
import edu.mit.csail.cgs.warpdrive.model.Model;
import edu.mit.csail.cgs.warpdrive.model.RegionExpanderModel;
import edu.mit.csail.cgs.ewok.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.*;


/* actually works on all regions.  Has special behaviors for those that implement Stranded or Named */
public class NamedStrandedPainter extends RegionPaintable {
    private static ColorSet colors = new ColorSet();
    private static int sMaxTracks = 12;
    private RegionExpanderModel model;    
    private boolean dirty;
    private int trackIndex[];
    private ArrayList<Region> regions;
    private int tracksUsed;
    private NamedStrandedProperties props;

    public NamedStrandedPainter(RegionExpanderModel model) {
        super();
        this.model = model;
        model.addEventListener(this);
        dirty = true;
        props = new NamedStrandedProperties();
        initLabels();
    }
    public NamedStrandedProperties getProperties() {return props;}
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
    public String getLabel(Region r) {
        if (r instanceof Named) {
            return ((Named)r).getName();
        } else {
            return r.toString();
        }
    }
    public Color getColor(Region r) {
        return colors.getColor(r.toString());
    }
    private void doLayout() {
        Iterator<Region> iter = model.getResults();
        regions = new ArrayList<Region>();
        while (iter.hasNext()) {
            Region f = iter.next();
            regions.add(f);
        }
        int lastEnd[] = new int[sMaxTracks];
        for (int i = 0; i < sMaxTracks; i++) { lastEnd[i] = -1;}
        trackIndex = new int[regions.size()];
        tracksUsed = 0;
        for (int j = 0; j < regions.size(); j++) {
            int minpos = regions.get(j).getStart();
            int maxpos = regions.get(j).getEnd();
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
        boolean drawgenenames = props.DrawNames;
        boolean fillcolor = props.FillRectangles;
        boolean alwaysDrawNames = props.AlwaysDrawNames;
        boolean drewAnything = false;

        clearLabels();
        for (int j = 0; j < regions.size(); j++) {
            Region f = regions.get(j);
            Color c = colors.getColor(f.toString());
            if (c == null) {
                continue;
            } else {
                g.setColor(c);
            }
            String label = getLabel(f);
            int min = f.getStart();
            int max = f.getEnd();
            int x1 = getXPos(min,
                             getRegion().getStart(), getRegion().getEnd(),
                             ulx,lrx);
            int x2 = getXPos(max,
                             getRegion().getStart(), getRegion().getEnd(),
                             ulx,lrx);
            int y = trackY[trackIndex[j]];

            if(!fillcolor)
            	g.drawRect(x1, y, x2 - x1, tracksize);
            else {
                if (x2 == x1) {
                    x2++;
                }
            	g.fillRect(x1, y, x2 - x1, tracksize);
            }
            if (f instanceof Stranded) {
                char str = ((Stranded)f).getStrand();
                if (str == '+') {
                    g.drawLine(x1,y,x2,y + tracksize /2);
                    g.drawLine(x1,y+tracksize,x2,y + tracksize /2);
                } else if (str == '-') {
                    g.drawLine(x2,y,x1,y + tracksize /2);
                    g.drawLine(x2,y+tracksize,x1,y + tracksize /2);
                }
            }
            addLabel(x1,y, x2 - x1, tracksize,label);
            int center = (x1+x2)/2;
            int texthalf = (label.length() * fontsize)/2;
            if (((center - texthalf > lastx[trackIndex[j]]) && drawgenenames) ||
                alwaysDrawNames){
                g.setColor(Color.black);
                g.drawString(label, center - texthalf,y + tracksize);
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


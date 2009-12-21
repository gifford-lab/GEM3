package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.util.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.viz.DynamicAttribute;
import edu.mit.csail.cgs.viz.colors.ColorSet;
import edu.mit.csail.cgs.datasets.chipchip.GenericExperiment;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.warpdrive.model.ChipChipDataModel;
import edu.mit.csail.cgs.warpdrive.model.ChipChipScaleModel;
import edu.mit.csail.cgs.warpdrive.model.Model;

public class ChipChipPainter extends RegionPaintable  {

    public static ColorSet colors = new ColorSet();

    private GenericExperiment data;
    private ChipChipDataModel model;
    private ChipChipScaleModel scale;
    private DynamicAttribute attrib;
    private ChipChipProperties props;
    public ChipChipPainter (GenericExperiment data, ChipChipDataModel model) {
        super();
        this.data = data;
        this.model = model;
        props = new ChipChipProperties();
        model.addEventListener(this);
        attrib = DynamicAttribute.getGlobalAttributes();
    }
    public ChipChipProperties getProperties() {return props;}
    
    public void setScaleModel(ChipChipScaleModel s) {
        scale = s;
    }
    public ChipChipScaleModel getScaleModel() {return scale;}
    
    public void cleanup() { 
        super.cleanup();
        model.removeEventListener(this);
    }

    public void paintItem(Graphics2D g, 
                          int ulx, int uly, 
                          int lrx, int lry) {
        //        System.err.println("CCP.canpaint " + canPaint());
        if (!canPaint()) {
            return;
        }
        boolean any = true;
        int i = 0,j = 0;
        Color repcolor = colors.getColor(getLabel());
        int startindex = data.baseToIndex(getRegion().getStart());
        int endindex = data.baseToIndex(getRegion().getEnd());

        g.setColor(Color.BLACK);
        while (any || (j == 0)) {
            g.setColor(repcolor);
            repcolor = repcolor.darker();
            any = false;
            int lastx = -1;
            int lasty = -1;
            int lasti = -1;
            for (i = startindex; i < endindex; i++) {
                if (j >= data.getReplicates(i)) {lasti = -1; continue;}
                any = true;
                double maxratio;
                if (getProperties().MaxRatio > 0) {
                    maxratio = getProperties().MaxRatio;
                } else {
                    maxratio = scale.getMaxVal();
                }
                int x = getXPos(data.getPos(i),
                                getRegion().getStart(), getRegion().getEnd(),
                                ulx,lrx);
                int y = getYPos(data.getValue(i,j),
                                0, maxratio,
                                uly,lry,props.RatiosOnLogScale);
                paintDatapointAt(g,x,y,i,j);
//                  System.err.println("Painting (" + i + "=" + data.getPos(i) + "," + j +
//                                     ")=" + data.getValue(i,j) + " -> (" + x + ","+y+")    lastx=" +
//                                     lastx +"   lasty=" + lasty + "   lasti=" + lasti);
                if ((lasti != -1) && 
                    (Math.abs(data.getPos(i) - data.getPos(lasti)) < 500)) {
                    connectDatapoints(g,lastx,lasty,x,y);
                }
                lastx = x;
                lasty = y;
                lasti = i;
            }
            j++;
        }
        if (props.DrawTrackLabel) {
            g.setColor(Color.BLACK);
            g.setFont(attrib.getLargeLabelFont(lrx - ulx,lry - uly));
            g.drawString(getLabel(),ulx,uly + g.getFont().getSize() * 2);
        }
    }

    public void paintDatapointAt(Graphics2D g, int x,int y,int i,int j) {
        int ovalsize = 4;
        g.fillOval(x-ovalsize/2,y-ovalsize/2,ovalsize,ovalsize);
    }
    public void connectDatapoints(Graphics2D g, int x1, int y1, int x2, int y2) {
        g.drawLine(x1,y1,x2,y2);        
    }
    
    public synchronized void eventRegistered(EventObject e) {        
        if ((e.getSource() == model) || (e.getSource() == scale) &&
            model.isReady() &&
            ((scale == null) || scale.isReady())) {
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


package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.awt.font.LineMetrics;
import java.util.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.viz.DynamicAttribute;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.warpdrive.model.ChipChipDataModel;
import edu.mit.csail.cgs.warpdrive.model.ChipChipScaleModel;
import edu.mit.csail.cgs.warpdrive.model.Model;

public class SimplifiedChipChipScalePainter extends RegionPaintable {

    private ChipChipScaleModel model;
    private ChipChipProperties props;
    public SimplifiedChipChipScalePainter (ChipChipScaleModel m, ChipChipProperties p, Listener<EventObject> parent) {
        addEventListener(parent);
        model = m;
        props = p;
        model.addEventListener(this);
    }
    public ChipChipProperties getProperties() {
        return props;
    }
    
    public void cleanup() { 
        super.cleanup();
        model.removeEventListener(this);
    }

    public void paintItem(Graphics2D g, 
                          int ulx, int uly, 
                          int lrx, int lry) {
        int w = lrx - ulx;
        int h = lry - uly;
        if (!props.DrawRatios) {
            return;
        }

        double minval = 0;
        double maxratio = model.getMaxVal();        
        
        g.setColor(Color.BLACK);        
        Stroke oldstroke = g.getStroke();
        g.setStroke(new BasicStroke((float)2.0));
        int y = getYPos(0,0,maxratio,uly,lry,getProperties().RatiosOnLogScale);       
        //g.drawLine(ulx,y,lrx,y);

        Font oldFont = g.getFont();
        Font font = DynamicAttribute.getGlobalAttributes().getRegionLabelFont(lrx-ulx,lry - uly);
        g.setFont(font);
        
        FontMetrics fontmetrics = g.getFontMetrics();
        LineMetrics linemetrics = fontmetrics.getLineMetrics("1.0",g);
        int fontHeight= (int)(linemetrics.getHeight());
        
        for(int yh = 1; yh <= (int)maxratio; yh +=2) { 
            y = getYPos(yh, 0, maxratio, uly, lry,getProperties().RatiosOnLogScale);
            g.drawString(String.format("%d", yh),
                    ulx + 7,y + fontHeight/2);
            g.drawLine(ulx+2, y, ulx+7, y);
        }
        
        g.drawLine(ulx+2, uly, ulx+2, lry);
        g.drawLine(ulx, lry-fontHeight-5, lrx, lry-fontHeight-5);
        
        Region r = model.getRegion();
        String leftCoord = String.valueOf(r.getStart());
        String rightCoord = String.valueOf(r.getEnd());
        
        //int rightSpacing = rightCoord.length() * font.getSize();
        int rightSpacing = fontmetrics.charsWidth(rightCoord.toCharArray(), 0, rightCoord.length());
        
        g.drawString(leftCoord, ulx+2, lry-2);
        g.drawString(rightCoord, lrx-(rightSpacing), lry-2);
        
        int hashTicks = 10;
        int hashWidth = w / hashTicks;
        for(int i = 0; i < hashTicks; i++) { 
            int xh = ulx + (i+1) * hashWidth;
            g.drawLine(xh, lry-fontHeight-5, xh, lry-fontHeight);
        }

        g.setStroke(oldstroke);

        g.setFont(oldFont);
    }
    
    public boolean canPaint() {
        boolean ready = model.isReady();
        setCanPaint(ready);
        setWantsPaint(ready);
        return ready;
    }

    public synchronized void eventRegistered(EventObject e) {        
        if (canPaint()) {
            notifyListeners();
        }
    }

    public boolean wantsPaint() {
        return canPaint();
    }
    public void removeEventListener(Listener<EventObject> l) {
        super.removeEventListener(l);
        if (!hasListeners()) {
            model.removeEventListener(this);
        }
    }
    public void configurePaintable(Component c) {}

}

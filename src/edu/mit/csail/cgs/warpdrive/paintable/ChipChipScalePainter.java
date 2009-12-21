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

public class ChipChipScalePainter extends RegionPaintable {

    private ChipChipScaleModel model;
    private ChipChipProperties props;
    public ChipChipScalePainter (ChipChipScaleModel m, Listener<EventObject> parent, ChipChipProperties p) {
        addEventListener(parent);
        model = m;
        props = p;
        model.addEventListener(this);
    }
    public ChipChipProperties getProperties() {return props;}
    public void setProperties(ChipChipProperties p) {props = p;}
    
    public void cleanup() { 
        super.cleanup();
        model.removeEventListener(this);
    }

    public void paintItem(Graphics2D g, 
                          int ulx, int uly, 
                          int lrx, int lry) {
        if (props.DrawRatios.equals(Boolean.FALSE)) {
            return;
        }
        boolean logscale = props.RatiosOnLogScale;
        double minratio = 0;
        double maxratio = model.getMaxVal();        
        if (getProperties().MaxRatio > 0) {
            maxratio = getProperties().MaxRatio;
        }

        if (logscale) {
            minratio = 1.0 / maxratio;
        }

        
        
        g.setColor(Color.BLACK);        
        Stroke oldstroke = g.getStroke();
        g.setStroke(new BasicStroke((float)2.0));
        int y;
        if (!logscale) {
            y = getYPos(0,0,maxratio,uly,lry,false);       
            g.drawLine(ulx,y,lrx,y);
            g.setStroke(oldstroke);
        }

        Font oldFont = g.getFont();
        g.setFont(DynamicAttribute.getGlobalAttributes().getRegionLabelFont(lrx-ulx,lry - uly));
        FontMetrics fontmetrics = g.getFontMetrics();
        LineMetrics linemetrics = fontmetrics.getLineMetrics("1.0",g);

        boolean drawone = props.OneHorizontalLine;
        boolean drawtwo = props.TwoHorizontalLine;
        boolean drawothers = props.OtherHorizontalLines;

        if (logscale) {
            for (double ratio = 1; ratio < maxratio; ratio *= 2) {
                y = getYPos(ratio,minratio,maxratio,uly,lry, logscale);       
                String l = String.format("%.2f", Math.log(ratio));
                g.drawString(l,ulx + 5,y + (int)(linemetrics.getHeight()/2) );
                if ((Math.abs(ratio - 1) < .0001 && drawone) ||
                    (Math.abs(ratio - 2) < .0001 && drawtwo) ||
                    drawothers) {
                    g.setColor(Color.GRAY);
                    g.drawLine(ulx + 7 + fontmetrics.charsWidth("2.0".toCharArray(),0,3),y,lrx,y);
                    g.setColor(Color.BLACK);
                }
                y = getYPos(1.0 / ratio,minratio,maxratio,uly,lry, logscale);       
                l = String.format("%.2f",-1 * Math.log(ratio));
                g.drawString(l,ulx + 5,y + (int)(linemetrics.getHeight()/2) );
                if ((Math.abs(ratio - 1) < .0001 && drawone) ||
                    (Math.abs(ratio - 2) < .0001 && drawtwo) ||
                    drawothers) {
                    g.setColor(Color.GRAY);
                    g.drawLine(ulx + 7 + fontmetrics.charsWidth("2.0".toCharArray(),0,3),y,lrx,y);
                    g.setColor(Color.BLACK);
                }                
            }

        } else {
            y = getYPos(1,minratio,maxratio,uly,lry, logscale);       
            g.drawString("1.0",ulx + 5,y + (int)(linemetrics.getHeight()/2) );
            if (drawone) {
                g.drawLine(ulx + 7 + fontmetrics.charsWidth("1.0".toCharArray(),0,3),y,lrx,y);
            }
            
            linemetrics = fontmetrics.getLineMetrics("2.0",g);
            y = getYPos(2,minratio,maxratio,uly,lry, logscale);       
            g.drawString("2.0",ulx + 5,y + (int)(linemetrics.getHeight()/2) );
            if (drawtwo) {
                g.setColor(Color.GRAY);
                g.drawLine(ulx + 7 + fontmetrics.charsWidth("2.0".toCharArray(),0,3),y,lrx,y);
                g.setColor(Color.BLACK);
            }
            for (double val = 4; val < maxratio; val *= 2) {
                y = getYPos(val,minratio,maxratio,uly,lry, logscale);       
                String l = new Double(val).toString();
                g.drawString(l,ulx + 5,y + (int)(linemetrics.getHeight()/2) );
                if (drawothers) {
                    g.setColor(Color.GRAY);
                    g.drawLine(ulx + 7 + fontmetrics.charsWidth(l.toCharArray(),0,3),y,lrx,y);
                    g.setColor(Color.BLACK);
                }
            }
            g.setFont(oldFont);
        }
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

}

/*
 * Created on Nov 9, 2006
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.util.*;

import javax.swing.JMenuItem;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.warpdrive.RectangleLookup;

/**
 * @author tdanford
 */
public class MultiGenePainter extends RegionPaintable {
    
    private Vector<GenePainter> genePainters;
    private RectangleLookup<GenePainter> posLookup;
    private GeneProperties props;
    private int readyCount;
    
    public MultiGenePainter() { 
        super();
        genePainters = new Vector<GenePainter>();
        setCanPaint(false);
        setWantsPaint(false);
        readyCount = 0;
        props = new GeneProperties();
        posLookup = new RectangleLookup<GenePainter>();
    }
    public GeneProperties getProperties() {return props;}
    
    public void cleanup() { 
        super.cleanup();
        posLookup.clear();
        for(GenePainter gp : genePainters) { 
            gp.cleanup();
        }
    }
    
    public void addGenePainter(GenePainter gp) { 
        genePainters.add(gp);
        gp.addEventListener(this);
        gp.setProperties(props);
    }
    
    public void removeGenePainter(GenePainter gp) { 
        genePainters.remove(gp);
        gp.removeEventListener(this);
    }
    
    private int countCanPaint() { 
        int count = 0;
        for(GenePainter gp : genePainters) { if(gp.canPaint()) { count += 1; } }
        return count;
    }
    
    public synchronized void eventRegistered(EventObject evt) { 
        readyCount = countCanPaint();
        if(readyCount == genePainters.size()) {
            setCanPaint(true);
            setWantsPaint(true);
            
            notifyListeners();
        } else { 
            setCanPaint(false);
            setWantsPaint(false);            
        }
    }
    
    public void setRegion(Region r) { 
        super.setRegion(r);
        for(GenePainter gp : genePainters) { gp.setRegion(r); }
        readyCount = genePainters.size();
    }
    
    public int getNumTracks() { 
        int count = 0;
        for(GenePainter gp : genePainters) { count += gp.getNumTracks(); }
        return count;
    }
    
    public ArrayList<JMenuItem> mouseClickedMenu(MouseEvent e) {
        ArrayList<JMenuItem> items = super.mouseClickedMenu(e);
        if(items == null) { items = new ArrayList<JMenuItem>(); }
        
        Point p = e.getPoint();
        Collection<GenePainter> subPainters = posLookup.getAllValues(p);
        for(GenePainter gp : subPainters) {
            ArrayList<JMenuItem> subitems = gp.mouseClickedMenu(e);
            if(subitems != null) { 
                items.addAll(subitems);
            }
        }
        
        if(items.size() == 0) { return null; }
        return items;
    }
    
    public int getMaxVertSpace() {
        int numTracks = getNumTracks();
        return Math.min(Math.max(40,numTracks * 25),120);
    }

    public void paintItem(Graphics2D g, int x1, int y1, int x2, int y2) {
        int totalTracks = Math.max(1, getNumTracks());
        int w = x2 - x1, h = y2 - y1;
        int yOffset = y1;
        posLookup.clear();
        for(int i = 0; i < genePainters.size(); i++) {
            GenePainter gp = genePainters.get(i);
            int tracks = gp.getNumTracks();
            double tf = (double)tracks / (double)totalTracks;
            int height = (int)Math.floor(tf * (double)h);
            gp.paintItem(g, x1, yOffset, x2, yOffset + height);
            
            Rectangle rect = new Rectangle(x1, yOffset, x2-x1, height);
            posLookup.addValue(gp, rect);
            
            yOffset += height;
        }
    }
}

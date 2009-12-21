package edu.mit.csail.cgs.viz.paintable;

import java.awt.Graphics;
import java.util.Collection;

import javax.swing.Action;


public class FixedPaintable
extends AbstractPaintable { 
    
    private Paintable fBase;
    private Integer fX, fY, fW, fH;
    
    public FixedPaintable(Paintable p) { 
        fBase = p;
        fX = fY = fW = fH = null;
    }
    
    public void setX(int x) { fX = new Integer(x); }
    public void setY(int y) { fY = new Integer(y); }
    public void setWidth(int w) { fW = new Integer(w); }
    public void setHeight(int h) { fH = new Integer(h); }
    
    public Collection<Action> getPaintableActions() { 
        return fBase.getPaintableActions();
    }
    
    public void paintItem(Graphics g, 
            int x1, int y1, 
            int x2, int y2) { 
        
        int nx1 = x1;
        int ny1 = y1;
        int nx2 = x2;
        int ny2 = y2;
        
        if(fX != null) { nx1 = fX.intValue(); }
        if(fY != null) { ny1 = fY.intValue(); }
        if(fW != null) { nx2 = nx1 + fW.intValue(); }
        if(fH != null) { ny2 = ny1 + fH.intValue(); }
        
        fBase.paintItem(g, nx1, ny1, nx2, ny2);
    }
}

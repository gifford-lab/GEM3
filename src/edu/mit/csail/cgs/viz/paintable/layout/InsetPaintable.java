/*
 * Created on Sep 13, 2005
 */
package edu.mit.csail.cgs.viz.paintable.layout;

import java.awt.*;
import java.util.*;

import javax.swing.Action;

import edu.mit.csail.cgs.viz.paintable.AbstractPaintable;
import edu.mit.csail.cgs.viz.paintable.Paintable;
import edu.mit.csail.cgs.viz.paintable.PaintableChangedEvent;
import edu.mit.csail.cgs.viz.paintable.PaintableChangedListener;

/**
 * @author tdanford
 */
public class InsetPaintable extends AbstractPaintable implements PaintableChangedListener {
    
    private Paintable inner;

    private Integer xIn, yIn;
    private Double xFracIn, yFracIn;
    
    public InsetPaintable(int xInset, int yInset, Paintable p) { 
        super();
        xIn = xInset; yIn = yInset;
        xFracIn = yFracIn = null;
        inner = p;
        inner.addPaintableChangedListener(this);
    }
    
    public InsetPaintable(double xin, double yin, Paintable p) { 
        super();
        xIn = yIn = null;
        xFracIn = xin; yFracIn = yin;
        inner = p;
        inner.addPaintableChangedListener(this);
    }
    
    public void PaintableChanged(PaintableChangedEvent e) { 
        for(PaintableChangedListener l : fListeners) { 
            l.paintableChanged(e);
        }
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.viz.paintable.Paintable#paintItem(java.awt.Graphics, int, int, int, int)
     */
    public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
    	int w = x2 - x1, h = y2 - y1;
    	int xInset = xIn != null ? xIn : (int)Math.floor(xFracIn * (double)w) / 2;
    	int yInset = yIn != null ? yIn : (int)Math.floor(yFracIn * (double)h) / 2;
    	
    	inner.paintItem(g, x1+xInset, y1+yInset, x2-xInset, y2-yInset);
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.viz.paintable.Paintable#getPaintableActions()
     */
    public Collection<Action> getPaintableActions() {
        return super.getPaintableActions();
    }
}

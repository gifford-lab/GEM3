package edu.mit.csail.cgs.viz.paintable;

import java.awt.Graphics;
import java.util.Collection;
import java.util.LinkedList;

import javax.swing.Action;


public class BorderPaintable 
extends AbstractPaintable { 
    
    private Paintable fTop, fBottom, fLeft, fRight, fCenter;
    private double fTopFrac, fBottomFrac, fLeftFrac, fRightFrac;
    
    public BorderPaintable(Paintable center) { 
        super();
        fTopFrac = fBottomFrac = fLeftFrac = fRightFrac = 0.1;
        fTop = fBottom = fLeft = fRight = null;
        fCenter = center;
        fCenter.addPaintableChangedListener(this);
    }
    
    public double getTopFrac() { return fTopFrac; }
    public double getBottomFrac() { return fBottomFrac; }
    public double getLeftFrac() { return fLeftFrac; }
    public double getRightFrac() { return fRightFrac; }
    
    public void setTop(Paintable p) { 
        if(fTop != null) { fTop.removePaintableChangedListener(this); }
        fTop = p; dispatchChangedEvent();
        if(p != null) { p.addPaintableChangedListener(this); }
    }
    
    public void setBottom(Paintable p) { 
        if(fBottom != null) { fBottom.removePaintableChangedListener(this); }
        fBottom = p; dispatchChangedEvent();
        if(p != null) { p.addPaintableChangedListener(this); }
    }
    
    public void setLeft(Paintable p) { 
        if(fLeft != null) { fLeft.removePaintableChangedListener(this); }
        fLeft = p; dispatchChangedEvent();
        if(p != null) { p.addPaintableChangedListener(this); }
    }
    
    public void setRight(Paintable p) { 
        if(fRight != null) { fRight.removePaintableChangedListener(this); }
        fRight = p; dispatchChangedEvent();
        if(p != null) { p.addPaintableChangedListener(this); }
    }
    
    public Collection<Action> getPaintableActions() { 
        LinkedList<Action> lst = new LinkedList<Action>();
        lst.addAll(fCenter.getPaintableActions());
        if(fTop != null) { lst.addAll(fTop.getPaintableActions()); }
        if(fBottom != null) { lst.addAll(fBottom.getPaintableActions()); }
        if(fLeft != null) { lst.addAll(fLeft.getPaintableActions()); }
        if(fRight != null) { lst.addAll(fRight.getPaintableActions()); }
        return lst;
    }
    
    public void registerClick(double xf, double yf) { 
        double x1 = 0.0;
        double x2 = 1.0;
        double y1 = 0.0;
        double y2 = 1.0;
        
        if(fLeft != null) { x1 = fLeftFrac; }
        if(fRight != null) { x2 = 1.0 - fRightFrac; }
        if(fTop != null) { y1 = fTopFrac; }
        if(fBottom != null) { y2 = 1.0 - fBottomFrac; }
        
        if(x1 <= xf && xf <= x2 && y1 <= yf && yf <= y2) { 
            fCenter.registerClick((xf - x1) / (x2 - x1), (yf - y1) / (y2 - y1));
        } else { 
            if(fLeft != null && xf < x1 && 
                    y1 <= yf && yf <= y2) { 
                fLeft.registerClick(xf / x1, (yf - y1) / (y2 - y1));
            }
            
            if(fRight != null && xf > x2 && 
                    y1 <= yf && yf <= y2) { 
                fRight.registerClick((xf - x2) / (1.0 - x2), (yf - y1) / (y2 - y1));
            }
            
            if(fTop != null && yf < y1 && 
                    x1 <= xf && xf <= x2) { 
                fTop.registerClick((xf - x1) / (x2 - x1), yf / y1);
            }
            
            if(fBottom != null && yf > y2 && 
                    x1 <= xf && xf <= x2) { 
                fBottom.registerClick((xf - x1) / (x2 - x1), (yf - y2) / (1.0 - y2));
            }
        }
    }
    
    public void paintItem(Graphics g, int x1, int y1, int x2, int y2) { 
        int xi1 = x1;
        int xi2 = x2;
        int yi1 = y1;
        int yi2 = y2;
        int w = x2 - x1; 
        int h = y2 - y1;
        
        if(fLeft != null) { 
            xi1 = x1 + (int)Math.round(fLeftFrac * (double)w);
        }
        
        if(fRight != null) { 
            xi2 = x2 - (int)Math.round(fRightFrac * (double)w);
        }
        
        if(fTop != null) { 
            yi1 = y1 + (int)Math.round(fTopFrac * (double)h);
        }
        
        if(fBottom != null) { 
            yi2 = y2 - (int)Math.round(fBottomFrac * (double)h);
        }
        
        fCenter.paintItem(g, xi1, yi1, xi2, yi2);
        if(fLeft != null) { fLeft.paintItem(g, x1, yi1, xi1, yi2); }
        if(fRight != null) { fRight.paintItem(g, xi2, yi1, x2, yi2); }
        if(fTop != null) { fTop.paintItem(g, xi1, y1, xi2, yi1); }
        if(fBottom != null) { fBottom.paintItem(g, xi1, yi2, xi2, y2); }
    }
}

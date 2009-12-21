package edu.mit.csail.cgs.viz.paintable.layout;

import java.awt.Graphics;
import java.util.Collection;

import javax.swing.Action;

import edu.mit.csail.cgs.viz.paintable.AbstractPaintable;
import edu.mit.csail.cgs.viz.paintable.Paintable;


public class ScaledPaintable 
extends AbstractPaintable { 
    
    private Paintable fBase;
    private Double fXFrac, fYFrac, fWidthFrac, fHeightFrac;
    
    public ScaledPaintable(Paintable p) { 
        super();
        fBase = p;
        fXFrac = fYFrac = new Double(0.0);
        fWidthFrac = fHeightFrac = new Double(1.0);
    }
    
    public ScaledPaintable(Paintable p, 
            double x, double y,
            double w, double h) { 
        super();
        fBase = p;
        fXFrac = new Double(x);
        fYFrac = new Double(y);
        fWidthFrac = new Double(w);
        fHeightFrac = new Double(h);
    }
    
    public void setX(double x) { fXFrac = new Double(x); }
    public void setY(double y) { fYFrac = new Double(y); }
    public void setHeight(double h) { fHeightFrac = new Double(h); }
    public void setWidth(double w) { fWidthFrac = new Double(w); }
    
    public double getXFrac() { return fXFrac.doubleValue(); }
    public double getYFrac() { return fYFrac.doubleValue(); }
    public double getWidthFrac() { return fWidthFrac.doubleValue(); }
    public double getHeightFrac() { return fHeightFrac.doubleValue(); }
    
    public Collection<Action> getPaintableActions() { 
        return fBase.getPaintableActions();
    }
    
    public void paintItem(Graphics g, 
            int x1, int y1, int x2, int y2) { 
        int w = x2 - x1;
        int h = y2 - y1;
        
        int nx1 = x1;
        int ny1 = y1;
        int nx2 = y2;
        int ny2 = y2;
        
        if(fXFrac != null) { 
            nx1 = x1 + (int)Math.round(fXFrac.doubleValue() * (double)w); 
        }
        
        if(fYFrac != null) { 
            ny1 = y1 + (int)Math.round(fYFrac.doubleValue() * (double)h); 
        }
        
        if(fWidthFrac != null) { 
            nx2 = nx1 + (int)Math.round(fWidthFrac.doubleValue() * (double)w);
        }
        
        if(fHeightFrac != null) { 
            ny2 = ny1 + (int)Math.round(fHeightFrac.doubleValue() * (double)h);
        }
        
        //System.out.println("Painting Scaled Item...");
        //System.out.println(x1 + " " + y1 + " -- " + y1 + " " + y2);
        //System.out.println(nx1 + " " + ny1 + " -- " + ny1 + " " + ny2);
        
        fBase.paintItem(g, nx1, ny1, nx2, ny2);
    }
}


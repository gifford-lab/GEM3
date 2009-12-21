package edu.mit.csail.cgs.viz.paintable.layout;

import java.awt.*;
import java.util.Collection;
import java.util.LinkedList;
import java.util.Vector;

import javax.swing.Action;

import edu.mit.csail.cgs.viz.paintable.AbstractPaintable;
import edu.mit.csail.cgs.viz.paintable.Paintable;


public class StackedPaintable
extends AbstractPaintable {
    
    private boolean fVertical;
    private Vector<Paintable> fArray;
    
    public StackedPaintable() {
        super();
        fVertical = true;
        fArray = new Vector<Paintable>();
    }
    
    public StackedPaintable(Paintable... p) { 
        super();
        fVertical = true;
        fArray = new Vector<Paintable>();
        for(int i = 0; i < p.length; i++) { 
        	fArray.add(p[i]);
        	p[i].addPaintableChangedListener(this);
        }
    }
    
    public StackedPaintable(boolean vert, Paintable... p) { 
        super();
        fVertical = vert;
        fArray = new Vector<Paintable>();
        for(int i = 0; i < p.length; i++) { 
        	fArray.add(p[i]);
        	p[i].addPaintableChangedListener(this);
        }
    }
    
    public int size() { return fArray.size(); }
    public void addPaintable(Paintable p) { 
        fArray.add(p); 
        p.addPaintableChangedListener(this);
    }
    public boolean isVertical() { return fVertical; }
    
    public void toggleVertical() {
        fVertical = !fVertical;
        dispatchChangedEvent();
    }
    
    public void registerClick(double xf, double yf) { 
        if(fArray.size() > 0) { 
            double span = 1.0 / (double)fArray.size();
            if(fVertical) { 
                int index = (int)Math.floor(yf * (double)fArray.size());
                fArray.get(index).registerClick(xf, (yf - (span * index)) / span);
            } else { 
                int index = (int)Math.floor(xf * (double)fArray.size());
                fArray.get(index).registerClick((xf - (span * index)) / span, yf);
            }
        }
    }
    
    public Collection<Action> getPaintableActions() { 
        LinkedList<Action> lst = new LinkedList<Action>();
        for(Paintable p : fArray) { 
            lst.addAll(p.getPaintableActions()); 
        }
        return lst;
    }
    
    public void paintItem(Graphics g, 
            int x1, int y1, int x2, int y2) {
    	
        if(fArray.size() == 0) { return; }
        int stackedHeight = 
            (int)Math.floor((double)(y2 - y1) / (double)fArray.size());
        int stackedWidth = 
            (int)Math.floor((double)(x2 - x1) / (double)fArray.size());
        
        for(int i = 0; i < fArray.size(); i++) { 
            Paintable p = fArray.get(i);
            
            if(fVertical) { 
                p.paintItem(g, x1, y1 + (i * stackedHeight),
                        x2, y1 + ((i + 1) * stackedHeight));
            } else { 
                p.paintItem(g, x1 + (i * stackedWidth), y1, 
                        x1 + ((i + 1) * stackedWidth), y2);
            }
        }
    }
}

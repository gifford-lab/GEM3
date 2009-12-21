package edu.mit.csail.cgs.viz.paintable.layout;

import java.awt.Color;
import java.awt.Graphics;
import java.util.Collection;
import java.util.LinkedList;

import javax.swing.Action;

import edu.mit.csail.cgs.viz.paintable.AbstractPaintable;
import edu.mit.csail.cgs.viz.paintable.Paintable;


public class LayeredPaintable
extends AbstractPaintable {
    
	private Color bg;
    private LinkedList<Paintable> fLayers;
    
    public LayeredPaintable() {
    	bg = null;
        fLayers = new LinkedList<Paintable>();
    }
    
    public LayeredPaintable(Paintable... p) {
    	this();
        for(int i = 0; i < p.length; i++) { 
        	fLayers.add(p[i]);
        	p[i].addPaintableChangedListener(this);
        }
    }
    
    public void pushPaintable(Paintable p) { 
        fLayers.addLast(p); 
        p.addPaintableChangedListener(this);
    }
    
    public void setBgColor(Color c) { bg = c; }
    
    public Paintable popPaintable() { 
        Paintable p = fLayers.removeLast(); 
        p.removePaintableChangedListener(this);
        return p;
    }
    
    public int size() { return fLayers.size(); }
    
    public void paintItem(Graphics g, 
            int x1, int y1, 
            int x2, int y2) {
    	
    	int w = x2-x1, h = y2-y1;
    	if(bg != null) { 
    		g.setColor(bg);
    		g.fillRect(x1, y1, w, h);
    	}
    	
        //System.out.print("Printing LayeredPaintable ");
        for(Paintable p : fLayers) { 
            p.paintItem(g, x1, y1, x2, y2);
            //System.out.print("."); System.out.flush();
        }
        
        //System.out.print("\n");
    }
    
    /*
    public Collection<Action> getPaintableActions() { 
        LinkedList<Action> lst = new LinkedList<Action>();
        for(Paintable p : fLayers) { 
            lst.addAll(p.getPaintableActions()); 
        }
        return lst;
    }
    */
}

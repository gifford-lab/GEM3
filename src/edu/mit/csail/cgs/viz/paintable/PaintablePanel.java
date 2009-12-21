package edu.mit.csail.cgs.viz.paintable;

import java.awt.Graphics;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.Collection;
import java.util.LinkedList;

import javax.swing.JPanel;

public class PaintablePanel 
extends JPanel
implements PaintableChangedListener {
    
    private Paintable fPaintable;
    
    public PaintablePanel() {
    	super();
    	fPaintable = null;
        addMouseListener(new MouseAdapter() { 
            public void mouseClicked(MouseEvent evt) { 
                int w = getWidth();
                int h = getHeight();
                if(w > 0 && h > 0) { 
                    double xf = (double)evt.getX() / (double)w;
                    double yf = (double)evt.getY() / (double)h;
                    if(fPaintable != null) { 
                    	fPaintable.registerClick(xf, yf);
                    }
                }
            }
        });
    }
    
    public PaintablePanel(Paintable p) { 
    	this();
        fPaintable = p;
        fPaintable.addPaintableChangedListener(this);
    }
    
    public Paintable getPaintable() { return fPaintable; }
    
    public void setPaintable(Paintable p) { 
    	if(fPaintable != null) { 
    		fPaintable.removePaintableChangedListener(this);
    	}
    	fPaintable = p; 
    	if(fPaintable != null) { 
    		fPaintable.addPaintableChangedListener(this);
    	}
    	repaint(); 
    }
    
    public void paintableChanged(PaintableChangedEvent evt) { 
        //System.out.println("Repainting...");
        repaint();
    }
    
    protected void paintItem(Graphics g, int x1, int y1, int x2, int y2) { 
        if(fPaintable != null) { 
        	fPaintable.paintItem(g, x1, y1, x2, y2);
        }
    }
    
    protected void paintComponent(Graphics g) { 
        super.paintComponent(g);
        int w = getWidth();
        int h = getHeight();
        paintItem(g, 0, 0, w, h);
    }
    
    public Collection<Paintable> getPaintables() { 
        LinkedList<Paintable> lst = new LinkedList<Paintable>();
        lst.addLast(fPaintable);
        return lst;
    }
}

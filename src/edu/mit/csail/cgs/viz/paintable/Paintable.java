package edu.mit.csail.cgs.viz.paintable;

import java.awt.Graphics;
import java.util.Collection;

import javax.swing.Action;

public interface Paintable { 
    public void paintItem(Graphics g, 
            int ulx, int uly, 
            int lrx, int lry);
    public Collection<Action> getPaintableActions();
    public Action getSaveImageAction();
    public void registerClick(double xf, double yf);
    public void addPaintableChangedListener(PaintableChangedListener l);
    public void removePaintableChangedListener(PaintableChangedListener l);
}


/*
 * Created on Mar 14, 2006
 */
package edu.mit.csail.cgs.warpdrive.paintable;

import java.util.Collection;
import java.awt.Graphics2D;
import java.util.Hashtable;


/**
 * @author tdanford
 */
public interface VizPaintable {
    public void paintItem(Graphics2D g, 
            int ulx, int uly, 
            int lrx, int lry);
    /* returns true iff the VizPaintable is ready to 
       draw itself. */
    public boolean canPaint();
    /* returns true when the VizPaintable is ready to paint
       itself but has not yet done so and when the underlying data (or the
       display that will be generated) has changed.  Should be true
       only if canPaint is true */
    public boolean wantsPaint();
    // a Paintable's label is some string that it may
    // use in its display to identify itself.
    public void setLabel(String label);
    public String getLabel();
    
    public void cleanup();
}

/**
 * 
 */
package edu.mit.csail.cgs.viz.paintable;

import java.awt.*;

/**
 * @author tdanford
 *
 */
public class BackgroundPaintable extends AbstractPaintable {
	
	private Color color;
	private Paintable inner;
	
	public BackgroundPaintable() { 
		this(Color.white);
	}
	
	public BackgroundPaintable(Color c, Paintable p) { 
		color = c;
		inner = p;
		inner.addPaintableChangedListener(this);
	}
	
	public BackgroundPaintable(Color c) { 
		color = c; 
		inner = null;
	}

	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		int w = x2-x1, h = y2-y1;
		g.setColor(color);
		g.fillRect(x1, y1, w, h);
		if(inner != null) { inner.paintItem(g, x1, y1, x2, y2); }
	}
}

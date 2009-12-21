/*
 * Author: tdanford
 * Date: Nov 14, 2008
 */
package edu.mit.csail.cgs.viz.eye;

import java.awt.*;
import java.util.Collection;

public class OverlayModelPaintable extends ContainerModelPaintable {

	public OverlayModelPaintable() { 
		super();
	}
	
	public OverlayModelPaintable(ModelPaintable... pts) { 
		super(pts);
	}
	
	public OverlayModelPaintable(Collection<ModelPaintable> pts) { 
		super(pts);
	}
	
	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		g.setColor(Color.white);
		g.fillRect(x1, y1, x2-x1, y2-y1);
		for(ModelPaintable p : innerPaintables) { 
			p.paintItem(g, x1, y1, x2, y2);
		}
	}
}

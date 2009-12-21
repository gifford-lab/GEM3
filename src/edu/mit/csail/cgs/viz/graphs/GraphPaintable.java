/*
 * Author: tdanford
 * Date: May 14, 2008
 */
/**
 * 
 */
package edu.mit.csail.cgs.viz.graphs;

import edu.mit.csail.cgs.viz.paintable.AbstractPaintable;
import edu.mit.csail.cgs.viz.paintable.Paintable;

import java.awt.*;
import java.awt.geom.*;
import java.util.*;

/**
 * @author tdanford
 *
 * Wraps the GraphView object, and turns it into a Paintable implementation.
 * 
 * This is the actual class which contains all the code to turn every option stored
 * in an ObjectView of some kind, into graphical features.  
 */
public class GraphPaintable extends AbstractPaintable {
	
	private GraphView view;
	
	public GraphPaintable(GraphView gv) { 
		super();
		view = gv;
	}
	
	public void updated() { 
		dispatchChangedEvent();
	}

	public static Point rotate(Point p, double theta) { 
		double px = (double)p.x, py = (double)p.y;
		double mx = Math.cos(theta), my = Math.sin(theta);
		int x = (int)Math.round(px * mx - py * my);
		int y = (int)Math.round(px * my + py * mx);
		return new Point(x, y);
	}
	
	public void paintHightlightedNode(Graphics g, NodeView nv) { 
		int x = nv.getX(), y = nv.getY();
		int w = nv.getWidth(), h = nv.getHeight();
		g.setColor(Color.green);
		g.drawRect(x - w/2 - 2, y - h/2 - 2, w + 4, h + 4);
	}
    
	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		
		int width = 10;
		int textOffset=20;
		Graphics2D g2 = (Graphics2D)g;
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        
        g2.setColor(Color.white);
        g2.fillRect(x1, y1, x2-x1, y2-y1);
        
		Stroke oldStroke = g2.getStroke();
		g2.setStroke(new BasicStroke((float)2.0));
		
		for(EdgeView ev : view.edges()) { 
			ev.paintView(g2);
		}
		
		for(SubEdgeView ev : view.subEdges()) {
			ev.paintView(g2);
		}
					
		for(NodeView nv : view.nodes()) { 
			nv.paintView(g2);
		}
		
		g2.setStroke(oldStroke);
		g2.setColor(Color.black);

		for(NodeView nv : view.nodes()) { 
			int x = nv.getX(), y = nv.getY();
			width = (Integer)nv.getOption("width");

			nv.paintName(g2);

			if(nv.containsOption("paintable")) { 
				Paintable p = (Paintable)nv.getOption("paintable");
				int pw = (Integer)nv.getOption("paintable-width");
				int ph = (Integer)nv.getOption("paintable-height");
				
				p.paintItem(g2, x-width-pw, y, x-width, y+ph);
			}
		}
		
		g2.setColor(Color.black);
		
	} 
	
	public static void drawArrow(Graphics2D g, int side, Point loc, double rot) { 
		int n = 3;
		int[] x = new int[3];
		int[] y = new int[3];
		int side2 = side*side;
		
		int hyp = (int)Math.round(Math.sqrt((double)(side2)));
		x[2]=hyp/2; x[1]=-1*hyp/2;x[0]=-1*hyp/2;
		y[2]=0; y[1]=side/2; y[0]= -side/2;
		
		g.translate(loc.x, loc.y);
		g.rotate(rot);
		g.fillPolygon(x, y, n);
		g.rotate(-rot);
		g.translate(-loc.x, -loc.y);
	}
}
/*
 * Author: tdanford
 * Date: Jan 4, 2009
 */
package edu.mit.csail.cgs.viz.paintable;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;
import java.util.*;

public class DirectedTagPainter extends AbstractPaintable {
	
	private Color color;
	private String label;
	private boolean direction;
	
	public DirectedTagPainter() { 
		color = Color.blue;
		direction = true;
		label = null;
	}
	
	public DirectedTagPainter(Color c) {
		this();
		color = c;
	}
	
	public DirectedTagPainter(Color c, boolean dir) { 
		this(c);
		direction = dir;
	}
	
	public DirectedTagPainter(boolean dir) { 
		this(Color.blue, dir);
	}
	
	public void setFillColor(Color c) { color = c; }
	public void setLabel(String l) { label = l; }
	public void setDirection(boolean dir) { direction = dir; }
	
	/**
	 * Some experimental code for painting arbitrary arrows from point-to-point.  
	 * 
	 */
	public void paintArrow(Graphics2D g2, Point p1, Point p2, int thickness) {
		int dx = p2.x - p1.x;
		int dy = p2.y - p1.y;
		double dx2 = (double)(dx * dx); 
		double dy2 = (double)(dy * dy);
		int dist = (int)Math.round(Math.sqrt(dx2 + dy2));

		double rot = 0.0;
		if(dx == 0 && dy == 0) { 
			rot = 0.0; 
		} else { 
			if(dy < 0) { 
				rot = Math.acos((double)dx / (double)dist); 
			} else { 
				rot = 2.0 * Math.PI - Math.acos((double)dx / (double)dist);
			}
		}

		g2.translate(p1.x, p1.y);
		g2.rotate(-rot);
		
		paintItem(g2, 0, -thickness/2, dist, +thickness/2);
		
		g2.rotate(rot);
		g2.translate(-p1.x, -p1.y);
	}
	
	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) { 
		paintTag((Graphics2D)g, x1, y1, x2, y2);
	}

	public void paintTag(Graphics2D g2, int x1, int y1, int x2, int y2) { 
		Font font = g2.getFont();
		Stroke stroke = g2.getStroke();

		int thickness = 3;
		g2.setStroke(new BasicStroke((float)thickness));

		int w = x2-x1, h = y2-y1;
		int h2 = h/2, h3 = h/3;
		int w4 = w/4;
		
		int arrowDepth = Math.min(h, w);   // width of the "head" of the arrow
		
		int ya = y1+thickness;		// ya is the y of the top of the arrow base. 
		int yb = y2-thickness;      // yb is the y of the bottom of the arrow base
		int yc = y1+h2;		        // yc is the y of the middle (the point) of the arrow.
		
		int xb = direction ? x1 : x2; 	// xb is the beginning of the base.
		int xc = direction ? x2 : x1;	// xc is the x of the point of the arrow head.
		
		// xa is the x of where the arrow base meets the head
		int xa = direction ? xc - arrowDepth : xc + arrowDepth;		
		
		int bw = w-arrowDepth+1;  			  // width of the "base" of the arrow 
		int bh = h-(thickness*2); 			  // height of the "base" of the arrow
		
		int[] hx = new int[] { xa, xc, xa };   // coordinates (hx and hy) for the arrow head.
		int[] hy = new int[] { ya, yc, yb };
		
		g2.setColor(color);
		g2.fillPolygon(hx, hy, 3);     // fills in the head of the arrow
		g2.fillRect((direction ? xb : xa), ya, bw, bh);   // fills in the base of the arrow
		
		g2.setColor(Color.black);
		
		g2.drawLine(xb, ya, xa, ya);   	// draws the top line of the arrow body 
		g2.drawLine(xb, ya, xb, yb);	// draws the back line of the arrow body 
		g2.drawLine(xb, yb, xa, yb);	// draws the bottom line of the arrow body
		
		// outlines the arrow head.
		g2.drawLine(xa, ya, xa, ya);    
		g2.drawLine(xa, ya, xc, yc);
		g2.drawLine(xa, yb, xa, yb);
		g2.drawLine(xa, yb, xc, yc);
		
		g2.setStroke(stroke);
		
		if(label != null) { 
			
			int fontSize = font.getSize();
			String fontFamily = font.getFamily();
			int fontStyle = font.getStyle();
			
			FontMetrics metrics = g2.getFontMetrics();
			int lblWidth = metrics.charsWidth(label.toCharArray(), 0, label.length());
			int lblHeight = metrics.getAscent() + metrics.getDescent();
			
			while(fontSize > 7 && (lblWidth >= bw || lblHeight >= bh)) { 
				fontSize -= 1;
				Font newFont = new Font(fontFamily, fontStyle, fontSize);
				g2.setFont(newFont);
				metrics = g2.getFontMetrics();
				lblWidth = metrics.charsWidth(label.toCharArray(), 0, label.length());
				lblHeight = metrics.getAscent() + metrics.getDescent();
			}

			if(lblWidth <= bw && lblHeight <= bh) { 
				if(direction) { 
					g2.drawString(label, xb+bw-lblWidth, ya+bh-metrics.getDescent());
				} else { 
					g2.drawString(label, xa+2, ya+bh-metrics.getDescent());				
				}
			}
		}
		
		g2.setFont(font);
	}
}

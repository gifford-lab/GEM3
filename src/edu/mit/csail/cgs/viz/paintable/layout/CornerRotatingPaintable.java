/*
 * Author: tdanford
 * Date: Jul 18, 2008
 */
package edu.mit.csail.cgs.viz.paintable.layout;

import java.awt.*;
import java.awt.geom.*;

import edu.mit.csail.cgs.viz.paintable.AbstractPaintable;
import edu.mit.csail.cgs.viz.paintable.Paintable;

public class CornerRotatingPaintable extends AbstractPaintable {
	
	public static enum Corner { UL, UR, LR, LL, CENTER };
	
	private double rotation;
	private Corner corner;
	private Paintable inner;
	
	public CornerRotatingPaintable(double rot, Paintable in) {
		corner = Corner.CENTER;
		rotation = rot;
		inner = in;
		inner.addPaintableChangedListener(this);
	}
	
	public CornerRotatingPaintable(Corner c, double rot, Paintable in) {
		corner = c;
		rotation = rot;
		inner = in;
		inner.addPaintableChangedListener(this);
	}
	
	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		int w = x2-x1, h = y2-y1;
		Graphics2D g2 = (Graphics2D)g;
		AffineTransform at = AffineTransform.getRotateInstance(rotation);
		Point axis = new Point(x1+w/2, y1+h/2);
		switch(corner) { 
		case UL: 
			axis = new Point(x1, y1);
			break;
		case UR: 
			axis = new Point(x2, y1);
			break;
		case LR: 
			axis = new Point(x2, y2);
			break;
		case LL: 
			axis = new Point(x1, y2);
			break;
		case CENTER: 
			axis = new Point(x1+w/2, y1+h/2);
			break;
		}
		
		try {
			AffineTransform invAt = at.createInverse();
			g2.translate(axis.x, axis.y);
			g2.transform(at);
			inner.paintItem(g2, x1, y1, x2, y2);
			g2.transform(invAt);
			g2.translate(-axis.x, -axis.y);
		} catch (NoninvertibleTransformException e) {
			e.printStackTrace();
		}
	}
}

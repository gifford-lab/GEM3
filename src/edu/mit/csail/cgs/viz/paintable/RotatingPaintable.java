/*
 * Author: tdanford
 * Date: Jul 18, 2008
 */
package edu.mit.csail.cgs.viz.paintable;

import java.awt.*;
import java.awt.geom.*;

public class RotatingPaintable extends AbstractPaintable {
	
	private double rotation;
	private Point axis;
	private Paintable inner;
	
	public RotatingPaintable(double rot, Paintable in) {
		axis = new Point(0, 0);
		rotation = rot;
		inner = in;
		inner.addPaintableChangedListener(this);
	}
	
	public RotatingPaintable(Point ax, double rot, Paintable in) {
		axis = ax;
		rotation = rot;
		inner = in;
		inner.addPaintableChangedListener(this);
	}
	
	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		Graphics2D g2 = (Graphics2D)g;
		AffineTransform at = AffineTransform.getRotateInstance(rotation);
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

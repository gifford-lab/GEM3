package edu.mit.csail.cgs.viz.expression;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Rectangle;

import edu.mit.csail.cgs.datasets.expression.ExprMeasurement;
import edu.mit.csail.cgs.viz.paintable.AbstractPaintable;

public class ExprMeasurementSpot extends AbstractPaintable {

	private IntensityColorMap colorMap;
	private ExprMeasurement measure;
	private Rectangle rect;
	
	public ExprMeasurementSpot(ExprMeasurement em, IntensityColorMap cm) { 
		colorMap = cm;
		measure = em;
		rect = new Rectangle(10, 10);
	}
	
	public ExprMeasurement getMeasurement() { return measure; }
	public void setColorMap(IntensityColorMap cm) { colorMap = cm; }
	
	public Rectangle getBounds() { return rect; }
	public void setBounds(Rectangle r) { rect = r; }

	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		int w = x2-x1, h = y2-y1;
		
		int realizedWidth = Math.min(w, rect.width);
		int realizedHeight = Math.min(h, rect.height);
		
		Color c = colorMap.getColor(measure.getValue());
		
		g.setColor(c);
		g.fillRect(x1, y1, realizedWidth, realizedHeight);
	}
}

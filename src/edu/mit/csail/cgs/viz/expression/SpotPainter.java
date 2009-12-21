package edu.mit.csail.cgs.viz.expression;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.util.*;

import edu.mit.csail.cgs.datasets.expression.ExprMeasurement;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.viz.paintable.AbstractPaintable;

public class SpotPainter extends AbstractPaintable {
	
	protected Vector<ExprMeasurementSpot> spots;
	private Map<Integer,Pair<Double,Double>> scaledLocations;
	private double scaledSpotWidth, scaledSpotHeight;
	private IntensityColorMap colorMap;
	private boolean squarePoints;
	private Color bgColor;
	private int rows, cols;
	
	public SpotPainter(int r, int c) {
		colorMap = new CombinationColorMap(-1.0, 0.0, 1.0, 
				new IntensityColorMap.GreenMap(), 
				new IntensityColorMap.RedMap());
		scaledSpotWidth = 1.0 / (double)c;
		scaledSpotHeight = 1.0 / (double)r;
		
		squarePoints = true;
		bgColor = Color.white;
		scaledLocations = new HashMap<Integer,Pair<Double,Double>>();
		spots = new Vector<ExprMeasurementSpot>();
	}
	
	public void setIntensityColorMap(IntensityColorMap cm) { 
		colorMap = cm;
		for(ExprMeasurementSpot ems : spots) { ems.setColorMap(cm); }
	}
	
	public void setSquarePoints(boolean sp) { squarePoints = sp; }
	public void setBackgroundColor(Color c) { bgColor = c; }
	
	public void addExprMeasurement(ExprMeasurement em) { 
		ExprMeasurementSpot spot = new ExprMeasurementSpot(em, colorMap);
		spots.add(spot);
		double x = 0.5, y = 0.5;
		Pair<Double,Double> pt = new Pair<Double,Double>(x, y);
		scaledLocations.put(spots.size()-1, pt);
	}

	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		int w = x2 - x1, h = y2 - y1;
		g.setColor(bgColor);
		g.fillRect(x1, y1, w, h);
		
		int sw = Math.max(1, (int)Math.floor(scaledSpotWidth * w));
		int sh = Math.max(1, (int)Math.floor(scaledSpotHeight * h));
		
		if(squarePoints) { 
			int ptdim = Math.max(sw, sh);
			sw = sh = ptdim;
		}
		
		int xoffset = sw/2, yoffset = sh/2;
		
		for(int i = 0; i < spots.size(); i++) { 
			Pair<Double,Double> pt = scaledLocations.get(i);
			ExprMeasurementSpot spot = spots.get(i);
			
			int x = x1 + (int)Math.floor(pt.getFirst() * (double)w) - xoffset;
			int y = y1 + (int)Math.floor(pt.getLast() * (double)h) - yoffset;
			
			spot.setBounds(new Rectangle(0, 0, sw, sh));
			spot.paintItem(g, x, y, x+sw, y+sh);
		}
	}

}

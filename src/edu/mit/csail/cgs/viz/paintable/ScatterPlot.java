/*
 * Author: Timothy Danford
 * Date: Jul 29, 2008
 */
package edu.mit.csail.cgs.viz.paintable;

import java.awt.*;
import java.lang.reflect.*;
import java.util.*;
import edu.mit.csail.cgs.viz.paintable.*;
import edu.mit.csail.cgs.viz.paintable.layout.InsetPaintable;
import edu.mit.csail.cgs.viz.paintable.layout.LayeredPaintable;

/**
 * @author tdanford
 */
public class ScatterPlot extends AbstractPaintable implements PaintableScaleListener {
	
	/**
	 * Testing and debugging. 
	 */
	public static void main(String[] args) { 
		Random rand = new Random();
		Paintable p = createSquareScatterPaintable(randPoints(rand, 1000));
		PaintableFrame pf = new PaintableFrame("Scatter Plot", p);
	}
	
	/**
	 * A general-purpose method to set up a ScatterPlot as well as a vertical and horizontal scale painter.
	 * @param pts
	 * @return
	 */
	public static Paintable createScatterPaintable(Collection<double[]> pts) { 
		PaintableScale xScale = new PaintableScale(0.0, 1.0);
		PaintableScale yScale = new PaintableScale(0.0, 1.0);
		
		Paintable xPainter = new HorizontalScalePainter(xScale);
		Paintable yPainter = new VerticalScalePainter(yScale);
		
		Paintable scatter = new ScatterPlot(xScale, yScale, pts);
		Paintable p = new LayeredPaintable(xPainter, yPainter, scatter);
		p = new InsetPaintable(0.1, 0.1, p);

		return p;
	}
	
	/**
	 * Same as above, except that the ScatterPlot display is to share the same x- and y-scales.
	 * 
	 * @param pts
	 * @return
	 */
	public static Paintable createSquareScatterPaintable(Collection<double[]> pts) { 
		PaintableScale xScale = new PaintableScale(0.0, 1.0);
		PaintableScale yScale = xScale;
		
		Paintable xPainter = new HorizontalScalePainter(xScale);
		Paintable yPainter = new VerticalScalePainter(yScale);
		
		Paintable scatter = new ScatterPlot(xScale, yScale, pts);
		Paintable p = new LayeredPaintable(xPainter, yPainter, scatter);
		p = new InsetPaintable(0.1, 0.1, p);

		return p;
	}
	
	/**
	 * Generates a set of random points, linearly related to each other, for testing purposes.
	 * @param rand
	 * @param count
	 * @return
	 */
	public static Collection<double[]> randPoints(Random rand, int count) { 
		LinkedList<double[]> pts = new LinkedList<double[]>();
		
		for(int i = 0; i < count; i++) {
			double x = rand.nextDouble();
			double y = x + rand.nextDouble();
			double[] p = new double[] { x, y };
			pts.add(p);
		}
		
		return pts;
	}
	
	/**
	 * Uses the findPoints() reflective method below, but with the default parameters for the fields names of "x" and "y".
	 * @param models
	 * @return
	 */
	public static Collection<double[]> findPoints(Collection models) {
		return findPoints(models, "x", "y");
	}
	
	/**
	 * Uses reflection to pick through a list of objects, and pulls out the values under the two given field names
	 * concatenates them as a double[] array (a point), and returns a list of them.  Ignores any non-Double or inaccessible
	 * fields.
	 * 
	 * @param models
	 * @param xName
	 * @param yName
	 * @return
	 */
	public static Collection<double[]> findPoints(Collection models, String xName, String yName) { 
		LinkedList<double[]> pts = new LinkedList<double[]>();
		for(Object value : models) { 
			Class valueClass = value.getClass();
			try {
				Field xfield = valueClass.getField(xName);
				Field yfield = valueClass.getField(yName);
				Object xValue = xfield.get(value);
				Object yValue = yfield.get(value);
				if(xValue instanceof Double && yValue instanceof Double) { 
					double[] pt = new double[] { (Double)xValue, (Double)yValue };
					pts.add(pt);
				}
				
			} catch (NoSuchFieldException e) {
				// do nothing.
			} catch (IllegalAccessException e) {
				// do nothing.
			}
		}
		return pts;
	}
	
	/**
	 * Turns a pair of arrays into a collection of pairs.  Used in the constructor to ScatterPlot, below.
	 * 
	 * @param x
	 * @param y
	 * @return
	 */
	public static Collection<double[]> collectPointArrays(double[] x, double[] y) { 
		if(x.length != y.length) { 
			throw new IllegalArgumentException(String.format("%d != %d", x.length, y.length));
		}
		
		LinkedList<double[]> pts = new LinkedList<double[]>();
		for(int i = 0; i < x.length; i++) { 
			pts.add(new double[] { x[i], y[i] });
		}
		return pts;
	}
	
	private Vector<double[]> pts;
	private PaintableScale xScale, yScale;
	private int radius;
	private Color color;
	
	public ScatterPlot(PaintableScale xs, PaintableScale ys, double[] x, double[] y) { 
		this(xs, ys, collectPointArrays(x, y));
	}
	
	public ScatterPlot(PaintableScale xscale, PaintableScale yscale, Collection<double[]> points) { 
		xScale = xscale;
		yScale = yscale;
		
		xScale.addPaintableScaleListener(this);
		yScale.addPaintableScaleListener(this);
		radius = 2;
		color = Color.red;
		
		pts = new Vector<double[]>();
		for(double[] pt : points) {
			pts.add(pt.clone());
			xScale.updateScale(pt[0]);
			yScale.updateScale(pt[1]);
		}
	}
	
	public void addPoint(double[] pt) { 
		pts.add(pt.clone());
		xScale.updateScale(pt[0]);
		yScale.updateScale(pt[1]);
		dispatchChangedEvent();
	}
	
	public void setColor(Color c) { 
		color = c;
		dispatchChangedEvent();
	}

	public void setRadius(int rad) { 
		radius = rad;
		dispatchChangedEvent();
	}

	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		int w = x2-x1, h = y2-y1;
		int diam = radius*2;
		
		for(double[] pt : pts) { 
			double fx = xScale.fractionalOffset(pt[0]);
			double fy = yScale.fractionalOffset(pt[1]);
			int fxx = x1 + (int)Math.round(fx * (double)w);
			int fyy = y2 - (int)Math.round(fy * (double)h);
			g.setColor(color);
			g.fillOval(fxx-radius, fyy-radius, diam, diam);
		}
	}

	public void paintableScaleChanged(PaintableScaleChangedEvent evt) {
		dispatchChangedEvent();
	}

}

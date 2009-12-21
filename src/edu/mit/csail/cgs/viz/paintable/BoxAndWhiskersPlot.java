/*
 * Author: tdanford
 * Date: Jul 18, 2008
 */
package edu.mit.csail.cgs.viz.paintable;

import java.util.*;
import java.awt.*;

import edu.mit.csail.cgs.viz.paintable.*;
import edu.mit.csail.cgs.viz.paintable.layout.InsetPaintable;
import edu.mit.csail.cgs.viz.paintable.layout.LayeredPaintable;
import edu.mit.csail.cgs.viz.paintable.layout.StackedPaintable;

public class BoxAndWhiskersPlot extends AbstractPaintable implements PaintableScaleListener {
	
	/** 
	 * Testing whether the BoxAndWhiskersPlot works, along with the Scale painters.
	 */
	public static void main(String[] args) { 
		double min = 0.0, max = 100.0;
		Random rand = new Random();
		int count = 5;

		double[][] bounds = new double[count][];
		for(int i = 0; i < count; i++) { 
			double[] b = randomBounds(rand, min, max);
			bounds[i] = b;
		}

		PaintableScale scale = new PaintableScale(min, max);
		Paintable p = createBoxAndWhiskersPlots(scale, bounds, null);
		
		PaintableFrame pf = new PaintableFrame("Box and Whiskers", p);
	}
	
	public static Paintable createBoxAndWhiskersPlots(PaintableScale scale, double[][] bounds, double[][] pts) {
		int count = bounds.length;
		Paintable[] plots = new Paintable[count];
		
		for(int i = 0; i < count; i++) { 
			double[] b = bounds[i];
			BoxAndWhiskersPlot pp = new BoxAndWhiskersPlot(scale, b[0], b[1], b[2], b[3], b[4]);
			if(pts != null && pts[i] != null) { 
				for(int j = 0; j < pts[i].length; j++) { 
					pp.addPointValue(pts[i][j]);
				}
			}
			plots[i] = new InsetPaintable(0.25, 0.0, pp);
		}

		Paintable plot = new StackedPaintable(false, plots);
		plot = new InsetPaintable(30, 0, plot);

		VerticalScalePainter scalePainter = new VerticalScalePainter(scale);
		Paintable p = new LayeredPaintable(scalePainter, plot);
		
		return p;
	}
	
	private static double[] randomBounds(Random rand, double min, double max) {  
		double[] a = new double[5];
		
		for(int i = 0; i < a.length; i++) { a[i] = rand.nextDouble(); }
		Arrays.sort(a);
		
		double range = max-min;
		for(int i = 0; i < a.length; i++) { 
			a[i] = min + a[i] * range;
		}
		
		return a;
	}
	
	private PaintableScale scale;
	private double min, max, lowerBox, upperBox, middle;
	private LinkedList<Double> pointValues;
	
	public BoxAndWhiskersPlot(PaintableScale s, 
			double min, double lower, double middle, double upper, double max) { 
		scale = s;
		this.min = min;
		this.max = max;
		this.lowerBox = lower; 
		this.upperBox = upper;
		this.middle = middle;
		
		if(min > lower || lower > middle || middle > upper || upper > max) { 
			throw new IllegalArgumentException(String.format("%f < %f < %f < %f < %f", 
					min, lower, middle, upper, max));
		}
		
		scale.addPaintableScaleListener(this);
		pointValues = new LinkedList<Double>();
	}

	public BoxAndWhiskersPlot(PaintableScale s, double[] a) { 
		this(s, a[0], a[1], a[2], a[3], a[4]);
	}
	
	public void addPointValue(Double pt) { 
		pointValues.add(pt);
	}

	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) { 
		scale.updateScale(min);
		scale.updateScale(max);
		
		int w = x2 - x1, h = y2 - y1;
		
		double f1 = scale.fractionalOffset(min);
		double f2 = scale.fractionalOffset(lowerBox);
		double f3 = scale.fractionalOffset(middle);
		double f4 = scale.fractionalOffset(upperBox);
		double f5 = scale.fractionalOffset(max);
		
		int xline = x1 + w/2;
		
		int fy1 = y2 - (int)Math.round(f1 * (double)h);
		int fy2 = y2 - (int)Math.round(f2 * (double)h);
		int fy3 = y2 - (int)Math.round(f3 * (double)h);
		int fy4 = y2 - (int)Math.round(f4 * (double)h);
		int fy5 = y2 - (int)Math.round(f5 * (double)h);
		
		System.out.println(String.format("%.3f %.3f %.3f %.3f %.3f", f1, f2, f3, f4 ,f5));

		g.setColor(Color.black);
		g.drawLine(xline, fy1, xline, fy5); // the middle line.  
		
		g.drawLine(x1, fy1, x2, fy1);  // the top line (maximum)
		g.drawLine(x1, fy5, x2, fy5);  // the bottom line (minimum)
		
		g.setColor(Color.white);
		g.fillRect(x1, fy4, x2-x1, Math.abs(fy4-fy2));
		g.setColor(Color.black);
		g.drawRect(x1, fy4, x2-x1, Math.abs(fy4-fy2));
		
		g.drawLine(x1, fy3, x2, fy3);
		
		for(Double pt : pointValues) { 
			double fpt = scale.fractionalOffset(pt);
			int fypt = y2 - (int)Math.round(fpt * (double)h);
			int rad = 3; 
			int diam = rad*2;
			g.setColor(Color.black);
			g.fillOval(xline-rad, fypt-rad, diam, diam);
		}
	}

	public void paintableScaleChanged(PaintableScaleChangedEvent evt) {
		dispatchChangedEvent();
	}
}

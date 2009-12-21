/*
 * Author: tdanford
 * Date: Sep 16, 2008
 */
package edu.mit.csail.cgs.viz.eye;

import java.awt.*;
import java.lang.reflect.*;
import java.util.*;

import edu.mit.csail.cgs.utils.models.Model;
import edu.mit.csail.cgs.utils.models.ModelFieldAnalysis;
import edu.mit.csail.cgs.viz.paintable.*;

public class ModelHistogram extends AbstractModelPaintable {
	
	public static final String xScaleKey = "x-scale";
	public static final String yScaleKey = "y-scale";
	public static final String colorKey = "color";
	public static final String strokeKey = "stroke";
	public static final String axisColorKey = "axis-color";
	public static final String barOutlineColorKey = "bar-outline-color";
	public static final String binsKey = "bins";
	
	private String valueFieldName;
	private int maxBin;
	private Vector<Double> points;
	private Vector<Integer[]> binCounts;
	private Vector<Double[]> binBounds;
	private Vector<Model> models;
	
	public ModelHistogram() { 
		valueFieldName = "value";
		
		points = new Vector<Double>();
		models = new Vector<Model>();
		binCounts = new Vector<Integer[]>();
		binBounds = new Vector<Double[]>();
		maxBin = 0;
		
		initProperty(new PropertyValueWrapper<Integer>(binsKey, 10));
		initProperty(new PropertyValueWrapper<PaintableScale>(xScaleKey, new PaintableScale(0.0, 1.0)));
		initProperty(new PropertyValueWrapper<PaintableScale>(yScaleKey, new PaintableScale(0.0, 1.0)));
		initProperty(new PropertyValueWrapper<Color>(colorKey, Color.red));
		initProperty(new PropertyValueWrapper<Float>(strokeKey, (float)3.0));
		
		rebin();
	}
	
	public ModelHistogram(String vfield) { 
		this();
		valueFieldName = vfield;
	}
	
	public void rebin() { 
		PaintableScale scale = getPropertyValue(xScaleKey);
		Integer bins = getPropertyValue(binsKey);
		
		System.out.println(String.format("Rebinning: %d, %.3f-%.3f", bins, scale.getMin(), scale.getMax()));
		
		maxBin = 0;
		binCounts.clear();
		binBounds.clear();
		
		double lower = scale.getMin(), upper = scale.getMax();
		double binWidth = (upper-lower)/(double)bins;
		
		double start = lower, end = start+binWidth;
		for(int i = 0; i < bins; i++, start += binWidth, end += binWidth) { 
			binBounds.add(new Double[] { start, end });
			System.out.println(String.format("\t%.3f-%.3f", start, end));
			binCounts.add(new Integer[] { 0 });
		}
		
		for(Double p : points) { 
			binPoint(p);
		}
		
		PaintableScale yscale = getPropertyValue(yScaleKey);
		//yscale.updateScale((double)maxBin);
		yscale.setScale(0.0, (double)maxBin);
	}
	
	private void binPoint(double value) { 
		PaintableScale scale = getPropertyValue(xScaleKey);
		Integer bins = getPropertyValue(binsKey);
		double lower = scale.getMin(), upper = scale.getMax();
		
		scale.updateScale(value);
		if(scale.getMin() < lower || scale.getMax() > upper) { 
			rebin();
			lower = scale.getMin(); 
			upper = scale.getMax();
		}

		double binWidth = (upper-lower)/(double)bins;
		int bin = Math.min((int)Math.floor((value-lower) / binWidth), bins-1);
		int count = (binCounts.get(bin)[0] += 1);
		maxBin = Math.max(maxBin, count);
	}
	
	public void addModel(Model m) {
		Class modelClass = m.getClass();
		ModelFieldAnalysis analysis = new ModelFieldAnalysis(modelClass);

		Field vfield = analysis.findField(valueFieldName);

		if(vfield != null) { 
			try {
				Object vvalue = vfield.get(m);
				
				if(vvalue != null) { 
					Class vclass = vvalue.getClass();
					
					if(!Model.isSubclass(vclass, Number.class)) { 
						throw new IllegalArgumentException("Value must be a Number");
					}
					
					Number vnumber = (Number)vvalue;
					
					double value = vnumber.doubleValue();

					PaintableScale xScale = getPropertyValue(xScaleKey);
					xScale.updateScale(value);
					
					points.add(value);
					models.add(m);
					
					binPoint(value);
					
					dispatchChangedEvent();
					
				} else { 
					throw new IllegalArgumentException("Value was null");
				}
				
			} catch (IllegalAccessException e) {
				throw new IllegalArgumentException("Value field was inaccessible", e);
			}
		} else { 
			String msg = String.format("No Fields %s", valueFieldName);
			throw new IllegalArgumentException(msg);
		}
	}

	public void addModels(Iterator<? extends Model> itr) {
		while(itr.hasNext()) { 
			addModel(itr.next());
		}
	}
	
	public ModelPaintable setProperty(ModelPaintableProperty p) {
		super.setProperty(p);
		if(p.getKey().equals(binsKey)) { 
			rebin();
			dispatchChangedEvent();
		} else if (p.getKey().equals(xScaleKey)) { 
			rebin();
			dispatchChangedEvent();
		}
		return this; 
	}

	public void clearModels() {
		points.clear();
		models.clear();

		// this will automatically call rebin() -- see the setProperty method, above.
		setProperty(new PropertyValueWrapper<PaintableScale>(xScaleKey, new PaintableScale(0.0, 1.0)));
		setProperty(new PropertyValueWrapper<PaintableScale>(yScaleKey, new PaintableScale(0.0, 1.0)));
		
		dispatchChangedEvent();
	}

	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		PaintableScale xScale = getPropertyValue(xScaleKey);
		PaintableScale yScale = getPropertyValue(yScaleKey);
		Color barColor = getPropertyValue(colorKey, Color.red);
		Color barOutlineColor = getPropertyValue(barOutlineColorKey, Color.white);
		float stroke = getPropertyValue(strokeKey, (float)1.0);
		Color axisColor = getPropertyValue(axisColorKey, Color.black);
		
		int w = x2-x1, h = y2-y1;
		
		Graphics2D g2 = (Graphics2D)g;
		Stroke oldStroke = g2.getStroke();
		g2.setStroke(new BasicStroke(stroke));
		
		g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		
		/** Painting Code **/
		
		int strokeInt = Math.max(1, (int)Math.floor(stroke));
		HorizontalScalePainter hsp = new HorizontalScalePainter(xScale);
		
		// Axes
		g2.setColor(axisColor);
		g2.drawRect(x1, y1, w-1, h-1);
		
		for(int i = 0; i < binCounts.size(); i++) { 
			Double[] bounds = binBounds.get(i);
			int count = binCounts.get(i)[0];
			
			if(count > 0) {
				// if count > 0, then maxBin > 0
				//double countf = (double)count / (double)maxBin;
				double countf = yScale.fractionalOffset((double)count);
				int yheight = (int)Math.round(countf * (double)(h-strokeInt*2));
				
				double bx1f = xScale.fractionalOffset(bounds[0]); 
				double bx2f = xScale.fractionalOffset(bounds[1]);
				
				int bx1 = x1 + (int)Math.round(bx1f * (double)(w-strokeInt*2)) + strokeInt + 1;
				int bx2 = x1 + (int)Math.round(bx2f * (double)(w-strokeInt*2)) + strokeInt;
				
				if(barOutlineColor != null) { 
					g.setColor(barOutlineColor);
					g.fillRect(bx1, y2-yheight-strokeInt, bx2-bx1, yheight);
				}
				g.setColor(barColor);
				g.fillRect(bx1+1, y2-yheight-strokeInt+1, bx2-bx1, yheight);
			}
		}
		
		hsp.paintItem(g, x1, y2-10, x2, y2);
		
		/** End Painting Code **/
		
		g2.setStroke(oldStroke);
	}
}

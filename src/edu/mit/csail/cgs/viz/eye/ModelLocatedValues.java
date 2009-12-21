/*
 * Author: tdanford
 * Date: Sep 16, 2008
 */
package edu.mit.csail.cgs.viz.eye;

import java.awt.*;
import java.lang.reflect.*;
import java.util.*;

import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.models.Model;
import edu.mit.csail.cgs.utils.models.ModelFieldAnalysis;
import edu.mit.csail.cgs.viz.paintable.*;
import edu.mit.csail.cgs.datasets.general.ScoredPoint;

public class ModelLocatedValues extends AbstractModelPaintable {

	public static final String boundsKey = "bounds";
	public static final String scaleKey = "scale";
	public static final String radiusKey = "radius";
	public static final String colorKey = "color";
	public static final String strokeKey = "stroke";
	public static final String axisColorKey = "axis-color";
	public static final String stemKey = "stems";
	public static final String connectedKey = "connected?";

	private String xFieldName, yFieldName;
	private Vector<PointModel> points;
	private Vector<Model> models;
	
	public ModelLocatedValues() { 
		xFieldName = "location";
		yFieldName = "value";
		points = new Vector<PointModel>();
		models = new Vector<Model>();
		
		initProperty(new PropertyValueWrapper<Integer[]>(boundsKey, new Integer[] { 0, 1 }));
		initProperty(new PropertyValueWrapper<PaintableScale>(scaleKey, new PaintableScale(0.0, 1.0)));
		initProperty(new PropertyValueWrapper<Integer>(radiusKey, 3));
		initProperty(new PropertyValueWrapper<Color>(colorKey, Color.red));
		initProperty(new PropertyValueWrapper<Float>(strokeKey, (float)2.0));
		initProperty(new PropertyValueWrapper<Boolean>(stemKey, (Boolean)true));
		initProperty(new PropertyValueWrapper<Boolean>(connectedKey, (Boolean)false));
		
		startDrawingPoints();
	}
	
	public ModelLocatedValues(String xfield, String yfield) { 
		this();
		xFieldName = xfield;
		yFieldName = yfield;
	}
	
	public void setBounds(int start, int end) { 
		PropertyValueWrapper<Integer[]> wrapper = getProperty(boundsKey);
		wrapper.setValue(new Integer[] { start, end });
	}
	
	public void addValue(Object value) { 
		if (value instanceof ScoredPoint) {  
			ScoredPoint p = (ScoredPoint)value;
			int location = p.getLocation();
			double val = p.getScore();
			
			Color c = getPropertyValue(colorKey);
			models.add(new PointModel(location, val, null, c));
			addLocatedValue(location, val, null, c);
			
		} else { 
			super.addValue(value);
		}
	}
	
	public void addModel(Model m) {
		Class modelClass = m.getClass();
		ModelFieldAnalysis analysis = new ModelFieldAnalysis(modelClass);

		Field xfield = analysis.findField(xFieldName);
		Field yfield = analysis.findField(yFieldName);

		if(xfield != null && yfield != null) { 
			try {
				Object xvalue = xfield.get(m);
				Object yvalue = yfield.get(m);
				
				if(xvalue != null && yvalue != null) { 
					Class xclass = xvalue.getClass();
					Class yclass = yvalue.getClass();
					
					if(!Model.isSubclass(xclass, Integer.class)) { 
						throw new IllegalArgumentException("Location value must be an Integer");
					}
					if(!Model.isSubclass(yclass, Number.class)) { 
						throw new IllegalArgumentException("Value must be a Number");
					}
					
					Integer xnumber = (Integer)xvalue;
					Number ynumber = (Number)yvalue;
					
					int x = xnumber.intValue();
					double y = ynumber.doubleValue();
					
					Color color = getPropertyValue(colorKey);
					Field colorfield = analysis.findField("color");
					if(colorfield != null && Model.isSubclass(colorfield.getType(), Color.class)) { 
						Color c = (Color)colorfield.get(m);
						if(c != null) { 
							color = c;
						}
					}
					
					models.add(m);
					addLocatedValue(x, y, m, color);
					
				} else { 
					throw new IllegalArgumentException("location or value was null");
				}
				
			} catch (IllegalAccessException e) {
				throw new IllegalArgumentException("location or value field was inaccessible", e);
			}
		} else { 
			String msg = "No Fields:";
			if(xfield == null) { 
				msg += String.format(" %s", xFieldName);
			}
			if(yfield == null) { 
				msg += String.format(" %s", yFieldName);
			}
			throw new IllegalArgumentException(msg);
		}
	}
	
	private void addLocatedValue(int x, double y, Model m, Color c) { 
		PaintableScale yScale = getPropertyValue(scaleKey);
		ModelPaintableProperty boundsProp = getProperty(boundsKey);
		Integer[] bounds = (Integer[])boundsProp.getValue();
		
		PointModel pm = new PointModel(x, y, m, c);
		points.add(pm);

		if(x < bounds[0] || x > bounds[1]) {
			bounds[0] = Math.min(x, bounds[0]);
			bounds[1] = Math.max(x, bounds[1]);
			setProperty(new PropertyValueWrapper<Integer[]>(boundsKey, bounds));
		}
		yScale.updateScale(y);
		
		dispatchChangedEvent();
	}

	public void addModels(Iterator<? extends Model> itr) {
		while(itr.hasNext()) { 
			addModel(itr.next());
		}
	}

	public void clearModels() {
		points.clear();
		models.clear();

		dispatchChangedEvent();
	}

	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		Integer[] bounds = getPropertyValue(boundsKey);
		PaintableScale yScale = getPropertyValue(scaleKey);
		int radius = getPropertyValue(radiusKey, 3);
		Color color = getPropertyValue(colorKey, Color.red);
		float stroke = getPropertyValue(strokeKey, (float)1.0);
		int strokeWidth = Math.max(1, (int)Math.floor(stroke));
		Color axisColor = getPropertyValue(axisColorKey, Color.black);
		Boolean stems = getPropertyValue(stemKey, true);
		Boolean connected = getPropertyValue(connectedKey, true);
		
		int diam = radius*2;
		int length = Math.max(1, bounds[1] - bounds[0]);
		
		int w = x2-x1, h = y2-y1;
		
		Graphics2D g2 = (Graphics2D)g;
		Stroke oldStroke = g2.getStroke();
		g2.setStroke(new BasicStroke(stroke));
		
		g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		
		clearDrawnPoints();
		
		/** Painting Code **/
		
		// Axes
		//g2.setColor(axisColor);
		//g2.drawRect(x1, y1, w-1, h-1);
		
		int wArea = w - strokeWidth*2;
		
		if(connected) {
			int ppx = -1, ppy = -1;
			
			g2.setColor(color);

			for(PointModel p : points) { 
				int location = p.location;
				double value = p.value;
				
				double xf = (double)(location-bounds[0]) / (double)length; 
				double yf = yScale.fractionalOffset(value);
				
				int px = x1 + (int)Math.round(xf * (double)wArea) + strokeWidth;
				int py = y2 - (int)Math.round(yf * (double)h);
				
				if(ppx != -1 && ppy != -1) {
					g.setColor(p.color);
					g.drawLine(ppx, ppy, px, py);
				}
			
				ppx = px; ppy = py;
			}
		}
		
		// Points
		for(PointModel p : points) { 
			int location = p.location;
			double value = p.value;
			//System.out.println(String.format("=> %d, %.2f", location, value));
			
			double xf = (double)(location-bounds[0]) / (double)length; 
			double yf = yScale.fractionalOffset(value);
			boolean outOfBounds = value < yScale.getMin() || value > yScale.getMax();
			
			int px = x1 + (int)Math.round(xf * (double)wArea) + strokeWidth;
			int py = y2 - (int)Math.round(yf * (double)h);

			if(!outOfBounds) { 
				g2.setColor(p.color);
				if(stems) { 
					g2.drawLine(px, py, px, y2);
				}
				//g2.setColor(Color.white);
				//g2.fillOval(px-radius, py-radius, diam, diam);
				g2.setColor(p.color);
				g2.drawOval(px-radius, py-radius, diam, diam);

				if(p.model != null) {
					drawPoint(new Point(px, py), p.model);
				}
			}
		}
		
		g2.setStroke(oldStroke);
	}
	
	public static class PointModel extends Model { 
		
		public Integer location;
		public Double value;
		public Model model;
		public Color color;
		
		public PointModel(int l, double v, Model m, Color c) { 
			location = l; value = v;
			model = m;
			color = c;
		}
	}
}



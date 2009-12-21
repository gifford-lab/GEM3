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

public class ModelLocatedDifferenceValues extends AbstractModelPaintable {

	public static final String boundsKey = "bounds";
	public static final String scaleKey = "scale";
	public static final String radiusKey = "radius";
	public static final String y1ColorKey = "y1color";
	public static final String y2ColorKey = "y2color";
	public static final String strokeKey = "stroke";
	public static final String axisColorKey = "axis-color";

	private String xFieldName, y1FieldName, y2FieldName;
	private Vector<DiffModel> points;
	private Vector<Model> models;
	
	public ModelLocatedDifferenceValues() { 
		xFieldName = "location";
		y1FieldName = "value1";
		y2FieldName = "value2";
		points = new Vector<DiffModel>();
		models = new Vector<Model>();
		
		initProperty(new PropertyValueWrapper<Integer[]>(boundsKey, new Integer[] { 0, 1 }));
		initProperty(new PropertyValueWrapper<PaintableScale>(scaleKey, new PaintableScale(0.0, 1.0)));
		initProperty(new PropertyValueWrapper<Integer>(radiusKey, 3));
		initProperty(new PropertyValueWrapper<Color>(y1ColorKey, Color.red));
		initProperty(new PropertyValueWrapper<Color>(y2ColorKey, Color.green));
		initProperty(new PropertyValueWrapper<Float>(strokeKey, (float)3.0));
		
		startDrawingPoints();
	}
	
	public ModelLocatedDifferenceValues(String xfield, String y1, String y2) { 
		this();
		xFieldName = xfield;
		y1FieldName = y1;
		y2FieldName = y2;
	}
	
	public void setBounds(int start, int end) { 
		PropertyValueWrapper<Integer[]> wrapper = getProperty(boundsKey);
		wrapper.setValue(new Integer[] { start, end });
	}
	
	public void addModel(Model m) {
		Class modelClass = m.getClass();
		ModelFieldAnalysis analysis = new ModelFieldAnalysis(modelClass);

		Field xfield = analysis.findField(xFieldName);
		Field y1field = analysis.findField(y1FieldName);
		Field y2field = analysis.findField(y2FieldName);

		if(xfield != null && y1field != null && y2field != null) { 
			try {
				Object xvalue = xfield.get(m);
				Object y1 = y1field.get(m);
				Object y2 = y2field.get(m);
				
				if(xvalue != null && y1 != null && y2 != null) {  
					Class xclass = xvalue.getClass();
					Class y1class = y1.getClass();
					Class y2class = y2.getClass();
					
					if(!Model.isSubclass(xclass, Integer.class)) { 
						throw new IllegalArgumentException("Location value must be an Integer");
					}
					if(!Model.isSubclass(y1class, Number.class)) { 
						throw new IllegalArgumentException("Value must be a Number");
					}
					if(!Model.isSubclass(y2class, Number.class)) { 
						throw new IllegalArgumentException("Value must be a Number");
					}
					
					Integer xnumber = (Integer)xvalue;
					Number y1num = (Number)y1;
					Number y2num = (Number)y2;
					
					int x = xnumber.intValue();
					double v1 = y1num.doubleValue();
					double v2 = y2num.doubleValue();
					
					models.add(m);
					addLocatedValue(x, v1, v2, m);
					
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
			if(y1field == null) { 
				msg += String.format(" %s", y1FieldName);
			}
			if(y2field == null) { 
				msg += String.format(" %s", y2FieldName);
			}
			throw new IllegalArgumentException(msg);
		}
	}
	
	private void addLocatedValue(int x, double y1, double y2, Model m) { 
		PaintableScale yScale = getPropertyValue(scaleKey);
		ModelPaintableProperty boundsProp = getProperty(boundsKey);
		Integer[] bounds = (Integer[])boundsProp.getValue();
		
		points.add(new DiffModel(x, y1, y2, m));

		if(x < bounds[0] || x > bounds[1]) {
			bounds[0] = Math.min(x, bounds[0]);
			bounds[1] = Math.max(x, bounds[1]);
			setProperty(new PropertyValueWrapper<Integer[]>(boundsKey, bounds));
		}
		yScale.updateScale(y1);
		yScale.updateScale(y2);
		
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
		Color y1color = getPropertyValue(y1ColorKey, Color.red);
		Color y2color = getPropertyValue(y2ColorKey, Color.green);
		float stroke = getPropertyValue(strokeKey, (float)1.0);
		int strokeWidth = Math.max(1, (int)Math.floor(stroke));
		Color axisColor = getPropertyValue(axisColorKey, Color.black);
		
		
		int length = Math.max(1, bounds[1] - bounds[0]);
		
		int w = x2-x1, h = y2-y1;
		
		int numPoints = 0;
		for(DiffModel p : points) { 
			if(p.location >= bounds[0] && p.location <= bounds[1]) { 
				numPoints += 1;
			}
		}

		// Basically, we don't want to draw the lines so thick that they overlap too much.  
		int maxRadius = Math.max(1, (int)Math.floor((double)w / (double)Math.max(1, numPoints/2)));
		radius = Math.min(radius, maxRadius);
		
		int diam = radius*2;
		
		Graphics2D g2 = (Graphics2D)g;
		Stroke oldStroke = g2.getStroke();
		g2.setStroke(new BasicStroke(stroke));
		
		g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		
		clearDrawnPoints();
		
		/** Painting Code **/
		
		// Axes
		g2.setColor(axisColor);
		g2.drawRect(x1, y1, w-1, h-1);
		
		int wArea = w - strokeWidth*2;
				
		// Points
		for(DiffModel p : points) { 
			int location = p.location;
			//System.out.println(String.format("=> %d, %.2f", location, value));

			double v1 = p.value1, v2 = p.value2;
			Color color = p.diff >= 0.0 ? y1color : y2color;
			
			double xf = (double)(location-bounds[0]) / (double)length; 
			double y1f = yScale.fractionalOffset(v1);
			double y2f = yScale.fractionalOffset(v2);
			
			int px = x1 + (int)Math.round(xf * (double)wArea) + strokeWidth;
			int py1 = y2 - (int)Math.round(y1f * (double)h);
			int py2 = y2 - (int)Math.round(y2f * (double)h);

			g2.setColor(color);
			//g2.drawLine(px, py1, px, py2);
			g2.fillRoundRect(px-radius, Math.min(py1, py2)-radius, diam, Math.abs(py2-py1)+diam, diam, diam);
			
			//g2.setColor(Color.lightGray);
			//g2.drawOval(px-radius, py1-radius, diam, diam);
			//g2.drawOval(px-radius, py2-radius, diam, diam);
			
			/*
			g2.setColor(Color.white);
			g2.fillOval(px-radius, py-radius, diam, diam);
			g2.setColor(color);
			g2.drawOval(px-radius, py-radius, diam, diam);
			*/
			
			if(p.model != null) {
				drawPoint(new Point(px, Math.min(py1, py2)), p.model);
			}
			
		}
		
		g2.setStroke(oldStroke);
	}
	
	public static class DiffModel extends Model { 
		
		public Integer location;
		public Double value1, value2, diff;
		public Model model;
		
		public DiffModel(int l, double v1, double v2, Model m) { 
			location = l; 
			value1 = v1; value2 = v2; 
			diff = value1 - value2; 
			model = m;
		}
	}
}



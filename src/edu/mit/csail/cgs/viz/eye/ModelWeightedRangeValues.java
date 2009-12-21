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
import edu.mit.csail.cgs.viz.colors.Coloring;
import edu.mit.csail.cgs.viz.paintable.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.ScoredPoint;
import edu.mit.csail.cgs.datasets.general.ScoredRegion;

public class ModelWeightedRangeValues extends AbstractModelPaintable {

	public static final String boundsKey = "bounds";
	public static final String scaleKey = "scale";
	public static final String colorKey = "color";
	public static final String strokeKey = "stroke";
	public static final String axisColorKey = "axis-color";
	public static final String drawWeightKey = "draw-weight";

	private String xFieldName, yFieldName, scoreFieldName;
	private Vector<Pair<Integer,Integer>> ranges;
	private Vector<Double> scores;
	private Vector<Model> models;
	
	public ModelWeightedRangeValues() { 
		xFieldName = "start";
		yFieldName = "end";
		scoreFieldName = "score";
		ranges = new Vector<Pair<Integer,Integer>>();
		scores = new Vector<Double>();
		models = new Vector<Model>();
		
		Color transblue = Color.blue;
		transblue = Coloring.clearer(Coloring.clearer(transblue));
		
		initProperty(new PropertyValueWrapper<Integer[]>(boundsKey, new Integer[] { 0, 1 }));
		initProperty(new PropertyValueWrapper<PaintableScale>(scaleKey, new PaintableScale(0.0, 1.0)));
		initProperty(new PropertyValueWrapper<Color>(colorKey, transblue));
		initProperty(new PropertyValueWrapper<Float>(strokeKey, (float)3.0));
		initProperty(new PropertyValueWrapper<Boolean>(drawWeightKey, Boolean.FALSE));
	}
	
	public ModelWeightedRangeValues(String xfield, String yfield, String scoreField) { 
		this();
		xFieldName = xfield;
		yFieldName = yfield;
		scoreFieldName = scoreField;
	}
	
	public void addValue(Object value) { 
		if (value instanceof ScoredRegion) {  
			ScoredRegion r = (ScoredRegion)value;
			int start = r.getStart(), end = r.getEnd();
			double s = r.getScore();
			
			models.add(new ScoredRangeModel(start, end, s));
			addScoredRangeValue(start, end, s);
			
		} else { 
			super.addValue(value);
		}
	}
	
	public void addModel(Model m) {
		Class modelClass = m.getClass();
		ModelFieldAnalysis analysis = new ModelFieldAnalysis(modelClass);

		Field xfield = analysis.findField(xFieldName);
		Field yfield = analysis.findField(yFieldName);
		Field scorefield = analysis.findField(scoreFieldName);

		if(xfield != null && yfield != null && scorefield != null) { 
			try {
				Object xvalue = xfield.get(m);
				Object yvalue = yfield.get(m);
				Object svalue = scorefield.get(m);
				
				if(xvalue != null && yvalue != null && scorefield != null) { 
					Class xclass = xvalue.getClass();
					Class yclass = yvalue.getClass();
					Class sclass = svalue.getClass();
					
					if(!Model.isSubclass(xclass, Integer.class)) { 
						throw new IllegalArgumentException("Start value must be an Integer");
					}
					if(!Model.isSubclass(yclass, Integer.class)) { 
						throw new IllegalArgumentException("End value must be an Integer");
					}
					if(!Model.isSubclass(sclass, Number.class)) { 
						throw new IllegalArgumentException("Score value must be an Number");
					}
					
					Integer xnumber = (Integer)xvalue;
					Integer ynumber = (Integer)yvalue;
					Number snumber = (Number)svalue;
					
					int x = xnumber.intValue();
					int y = ynumber.intValue();
					double sc = snumber.doubleValue();
					
					models.add(m);
					addScoredRangeValue(x, y, sc);
					
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
			if(scorefield == null) { 
				msg += String.format(" %s", scoreFieldName);
			}
			throw new IllegalArgumentException(msg);
		}
	}
	
	private void addScoredRangeValue(int x, int y, double s) {
		if(x > y) { throw new IllegalArgumentException(); }
		
		PaintableScale scale = getPropertyValue(scaleKey);
		
		ModelPaintableProperty boundsProp = getProperty(boundsKey);
		Integer[] bounds = (Integer[])boundsProp.getValue();
		
		ranges.add(new Pair<Integer,Integer>(x, y));
		scores.add(s);

		if(x < bounds[0] || y > bounds[1]) {
			bounds[0] = Math.min(x, bounds[0]);
			bounds[1] = Math.max(y, bounds[1]);
			setProperty(new PropertyValueWrapper<Integer[]>(boundsKey, bounds));
		}

		scale.updateScale(s);
		
		dispatchChangedEvent();
	}

	public void addModels(Iterator<? extends Model> itr) {
		while(itr.hasNext()) { 
			addModel(itr.next());
		}
	}

	public void clearModels() {
		ranges.clear();
		models.clear();
		scores.clear();
		
		dispatchChangedEvent();
	}

	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		Integer[] bounds = getPropertyValue(boundsKey);
		Color color = getPropertyValue(colorKey, Color.red);
		float stroke = getPropertyValue(strokeKey, (float)1.0);
		int strokeWidth = Math.max(1, (int)Math.floor(stroke));
		Color axisColor = getPropertyValue(axisColorKey, Color.black);
		PaintableScale yScale = getPropertyValue(scaleKey);
		Boolean drawWeights = getPropertyValue(drawWeightKey);
		
		int length = Math.max(1, bounds[1] - bounds[0] + 1);
		
		int w = x2-x1, h = y2-y1;
		
		Graphics2D g2 = (Graphics2D)g;
		Stroke oldStroke = g2.getStroke();
		g2.setStroke(new BasicStroke(stroke));
		
		g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		
		/** Painting Code **/
		
		// Axes
		g2.setColor(axisColor);
		g2.drawRect(x1, y1, w-1, h-1);
		
		int wArea = w - strokeWidth*2;
		
		// Points
		for(int i = 0; i < ranges.size(); i++) { 
			Pair<Integer,Integer> r = ranges.get(i);
			Double score = scores.get(i);
			
			int start = r.getFirst();
			int end = r.getLast();
			//System.out.println(String.format("=> %d, %d : %.3f", start, end, score));
			
			double xf = (double)(start-bounds[0]) / (double)length; 
			double yf = (double)(end-bounds[0]+1) / (double)length;
			double zf = yScale.fractionalOffset(score);
			
			int px = x1 + (int)Math.round(xf * (double)wArea);
			int py = x1 + (int)Math.round(yf * (double)wArea);
			int pz = y2 - (int)Math.round(zf * (double)h);

			g2.setColor(color);
			g2.fillRect(px, pz, py-px, y2-pz);
			
			g2.setColor(Color.black);
			g2.drawLine(px, pz, py, pz);
			
			if(drawWeights) { 
				g2.drawString(String.format("%.2", score), px+(py-px)/2, y2-pz-1);
			}
		}
		
		g2.setStroke(oldStroke);
	}
	
	public static class ScoredRangeModel extends Model { 
		public Integer start, end;
		public Double score;
		
		public ScoredRangeModel(int s, int e, double ss) { 
			start = s; end = e;
			score = ss;
		}
	}
}



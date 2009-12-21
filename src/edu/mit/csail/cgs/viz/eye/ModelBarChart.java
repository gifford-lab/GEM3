/*
 * Author: tdanford
 * Date: Jun 3, 2009
 */
package edu.mit.csail.cgs.viz.eye;

import java.util.*;
import java.lang.reflect.*;
import java.awt.*;
import java.io.*;

import edu.mit.csail.cgs.utils.models.*;
import edu.mit.csail.cgs.utils.models.data.DataFrame;
import edu.mit.csail.cgs.viz.paintable.PaintableFrame;
import edu.mit.csail.cgs.viz.paintable.PaintableScale;
import edu.mit.csail.cgs.viz.paintable.VerticalScalePainter;

public class ModelBarChart extends AbstractModelPaintable {
	
	public static void main(String[] args) { 
		File f = new File(args[0]);
		
		try {
			DataFrame<BarValue> frame = new DataFrame<BarValue>(BarValue.class, f, false, "x", "y");
			System.out.println(String.format("Loaded %d values.", frame.size()));
			
			ModelBarChart chart = new ModelBarChart();
			chart.addModels(frame.iterator());
			
			//chart.setProperty(new PropertyValueWrapper(labelKey, Boolean.FALSE));
			chart.setProperty(new PropertyValueWrapper(sortedKey, Boolean.TRUE));
			
			PaintableFrame pf = new PaintableFrame(f.getName(), chart);
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static final String scaleKey = "scale";
	public static final String colorKey = "color";
	public static final String labelKey = "labels";
	public static final String outlineKey = "outlines";
	public static final String sortedKey = "sorted";

	private LinkedList<Model> models;
	private Map<String,Double> values;
	
	private String xField, yField;
	
	public ModelBarChart() { 
		this("x", "y");
	}
	
	public ModelBarChart(String x, String y) { 
		xField = x; 
		yField = y;
		models = new LinkedList<Model>();
		values = new LinkedHashMap<String,Double>();
		initProperty(new PropertyValueWrapper(scaleKey, new PaintableScale(0.0, 1.0)));
		initProperty(new PropertyValueWrapper(colorKey, Color.blue));
		initProperty(new PropertyValueWrapper(labelKey, Boolean.TRUE));
		initProperty(new PropertyValueWrapper(outlineKey, Boolean.TRUE));
		initProperty(new PropertyValueWrapper(sortedKey, Boolean.FALSE));
	}

	public void addModel(Model m) {
		ModelFieldAnalysis analysis = new ModelFieldAnalysis(m.getClass());
		Field xfld = analysis.findTypedField(xField, String.class);
		Field yfld = analysis.findTypedField(yField, Number.class);
		
		if(xfld != null && yfld != null) { 
			try {
				String xvalue = (String) xfld.get(m);
				Number yvalue = ((Number)yfld.get(m));
				
				if(xvalue != null && yvalue != null) { 
					models.add(m);
					double y = yvalue.doubleValue();
					PaintableScale scale = getPropertyValue(scaleKey);
					scale.updateScale(y);
					
					values.put(xvalue, y);
					dispatchChangedEvent();
				}
			} catch (IllegalAccessException e) {
				e.printStackTrace();
			}
		}
	}

	public void addModels(Iterator<? extends Model> itr) {
		while(itr.hasNext()) { 
			addModel(itr.next());
		}
	}
	
	public BarValue[] values() { 
		BarValue[] varray = new BarValue[values.size()];
		int i = 0;
		for(String key : values.keySet()) { 
			varray[i++] = new BarValue(key, values.get(key));
		}
		Arrays.sort(varray);
		return varray;
	}
	
	private Collection<String> sortedKeys() { 
		BarValue[] varray = values();
		LinkedList<String> keys = new LinkedList<String>();
		for(int i = 0; i < varray.length; i++) { 
			keys.add(varray[i].x);
		}
		return keys;
	}

	public void clearModels() {
		models.clear();
		values.clear();
	}

	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		PaintableScale scale = getPropertyValue(scaleKey);
		Color color = getPropertyValue(colorKey);
		Boolean labels = getPropertyValue(labelKey);
		Boolean outlines = getPropertyValue(outlineKey);
		Boolean sorted = getPropertyValue(sortedKey);
		
		int w = x2-x1, h = y2-y1;
		
		Graphics2D g2 = (Graphics2D)g;
		Font oldFont = g2.getFont();
		int fontSize = 28;
		Font newFont = new Font("Times New Roman", Font.PLAIN, fontSize);
		
		g2.setFont(newFont);
		g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, 
				RenderingHints.VALUE_ANTIALIAS_ON);
		
		FontMetrics fm = g2.getFontMetrics();
		int textHeight = fm.getAscent() + fm.getDescent();
		
		int numValues = Math.max(1, values.size());
		int strips = (values.size()*2) + 1;
		int stripWidth = (int)Math.floor((double)w / (double)strips);
		
		int labelHeight = labels ? (int)Math.floor((double)h / 5.0) : 0;
		int dataHeight = h - labelHeight;
		
		while(textHeight > stripWidth*2 && fontSize > 1) { 
			fontSize -= 1;
			newFont = new Font(newFont.getFamily(), Font.PLAIN, fontSize);
			g2.setFont(newFont);
			fm = g2.getFontMetrics();
			textHeight = fm.getAscent() + fm.getDescent();
		}
		
		int i = 0;
		
		Collection<String> keys = sorted ? sortedKeys() : values.keySet();
		
		for(String key : keys) { 
			int x = stripWidth + (i*2*stripWidth);
			double value = values.get(key);
			double valuef = scale.fractionalOffset(value);
			
			int barHeight = (int)Math.round(valuef * (double)dataHeight);
			int y = y1 + dataHeight - barHeight;

			g.setColor(color);
			g.fillRect(x, y, stripWidth, barHeight);

			if(outlines) { 
				g.setColor(Color.black);
				Stroke oldStroke = g2.getStroke();
				Stroke newStroke = new BasicStroke((float)2.0);
				g2.drawRect(x, y, stripWidth, barHeight);
				g2.setStroke(oldStroke);
			}

			if(labels) { 
				g.setColor(Color.black);
				
				int textWidth = fm.charsWidth(key.toCharArray(), 0, key.length());
				FontMetrics cfm = fm;
				Font currentFont = g2.getFont();
				int currentFontSize = fontSize;
				
				while(textWidth > (h - dataHeight) && 
						currentFontSize > 1) {
					
					currentFontSize -= 1;
					currentFont = new Font(currentFont.getFamily(), Font.PLAIN, currentFontSize);
					
					g2.setFont(currentFont);
					cfm = g2.getFontMetrics();
					
					textWidth = cfm.charsWidth(key.toCharArray(), 0, key.length());
				}
				
				//int tx = x + stripWidth/2 - textWidth/2;
				//int ty = y + barHeight + fm.getAscent() + 2;
				
				int tx = x + stripWidth/2 - cfm.getDescent();
				int ty = y2 - labelHeight + 2;
				
				g2.translate(tx, ty);
				g2.rotate(Math.PI/2.0);
				
				//g2.drawString(key, tx, ty);
				g2.drawString(key, 0, 0);
				
				g2.rotate(-Math.PI / 2.0);
				g2.translate(-tx, -ty);
				
				g2.setFont(newFont);
			}
			
			i += 1;
		}
		
		VerticalScalePainter vsp = new VerticalScalePainter(scale);
		g.setColor(Color.black);
		g.drawLine(x1, y1+dataHeight, x2, y1+dataHeight);
		vsp.paintItem(g, x1, y1, x2, y1+dataHeight);
		
		g2.setFont(oldFont);
	}
	
	public static class BarValue extends Model implements Comparable<BarValue> {
		
		public String x; 
		public Double y; 
		
		public BarValue() {}
		
		public BarValue(String x, Double y) { 
			this.x = x; 
			this.y = y;
		}

		public int hashCode() { return x.hashCode(); }
		
		public boolean equals(Object o) { 
			if(!(o instanceof BarValue)) { return false; }
			BarValue v = (BarValue)o;
			if(!x.equals(v.x)) { return false; }
			return true;
		}
		
		public int compareTo(BarValue v) { 
			if(y > v.y) { return -1; }
			if(y < v.y) { return 1; }
			return x.compareTo(v.x);
		}
	}
}



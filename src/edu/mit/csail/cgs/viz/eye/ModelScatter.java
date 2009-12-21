/*
 * Author: tdanford
 * Date: Sep 16, 2008
 */
package edu.mit.csail.cgs.viz.eye;

import java.awt.*;
import java.lang.reflect.*;
import java.util.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.models.Model;
import edu.mit.csail.cgs.utils.models.ModelFieldAnalysis;
import edu.mit.csail.cgs.viz.paintable.*;

public class ModelScatter extends AbstractModelPaintable {
	
	public static final String xScaleKey = "x-scale";
	public static final String yScaleKey = "y-scale";
	public static final String radiusKey = "radius";
	public static final String colorKey = "color";
	public static final String strokeKey = "stroke";
	public static final String axisColorKey = "axis-color";
	public static final String showScaleKey = "show-scales";

	private String xFieldName, yFieldName;
	private Vector<ScatterPoint> points;
	private Vector<Model> models;
	
	public ModelScatter() { 
		xFieldName = "x"; 
		yFieldName = "y";
		points = new Vector<ScatterPoint>();
		models = new Vector<Model>();
		
		initProperty(new PropertyValueWrapper<PaintableScale>(xScaleKey, new PaintableScale(0.0, 1.0)));
		initProperty(new PropertyValueWrapper<PaintableScale>(yScaleKey, new PaintableScale(0.0, 1.0)));
		initProperty(new PropertyValueWrapper<Integer>(radiusKey, 3));
		initProperty(new PropertyValueWrapper<Color>(colorKey, Color.red));
		initProperty(new PropertyValueWrapper<Float>(strokeKey, (float)3.0));
		initProperty(new PropertyValueWrapper<Boolean>(showScaleKey, Boolean.TRUE));
		
		startDrawingPoints();
	}
	
	public ModelScatter(String xfield, String yfield) { 
		this();
		xFieldName = xfield;
		yFieldName = yfield;
	}
	
	public void addModel(Model m) {
		Class modelClass = m.getClass();
		ModelFieldAnalysis analysis = new ModelFieldAnalysis(modelClass);

		Field xfield = analysis.findField(xFieldName);
		Field yfield = analysis.findField(yFieldName);
		Field vectorField = analysis.findTypedField("vector", Boolean.class);
		Field colorField = analysis.findTypedField("color", Color.class);
		
		if(xfield != null && yfield != null) { 
			try {
				Object xvalue = xfield.get(m);
				Object yvalue = yfield.get(m);
				
				if(xvalue != null && yvalue != null) { 
					Class xclass = xvalue.getClass();
					Class yclass = yvalue.getClass();
					
					if(!Model.isSubclass(xclass, Number.class)) { 
						throw new IllegalArgumentException("X value must be a Number");
					}
					if(!Model.isSubclass(yclass, Number.class)) { 
						throw new IllegalArgumentException("Y value must be a Number");
					}
					
					Number xnumber = (Number)xvalue;
					Number ynumber = (Number)yvalue;
					
					double x = xnumber.doubleValue();
					double y = ynumber.doubleValue();
					
					PaintableScale xScale = getPropertyValue(xScaleKey);
					PaintableScale yScale = getPropertyValue(yScaleKey);
				
					ScatterPoint sp = new ScatterPoint(x,y);
					if(vectorField != null) { 
						Boolean isVector = (Boolean)vectorField.get(m);
						if(isVector != null && isVector) { 
							sp.vector = true;
						}
					}
					
					if(colorField != null) { 
						sp.color = (Color)colorField.get(m);
					}
					
					points.add(sp);
					models.add(m);
					
					xScale.updateScale(x);
					yScale.updateScale(y);
					
					if(sp.vector) { 
						xScale.updateScale(0.0);
						yScale.updateScale(0.0);
					}
					
					dispatchChangedEvent();
					
				} else { 
					throw new IllegalArgumentException("x or y value was null");
				}
				
			} catch (IllegalAccessException e) {
				throw new IllegalArgumentException("x or y field was inaccessible", e);
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

	public void addModels(Iterator<? extends Model> itr) {
		while(itr.hasNext()) { 
			addModel(itr.next());
		}
	}
	
	public <T> Set<T> findTags(double x1, double y1, double x2, double y2, String tagName) {
		System.out.println(String.format("findTags(%.2f,%.2f,%.2f,%.2f,%s)",
				x1,y1, x2,y2, tagName));
		
		Set<T> tags = new HashSet<T>();

		Collection<Model> contained = findModels(x1, y1, x2, y2);
		for(Model m : contained) { 
			try {
				Field f = m.getClass().getField(tagName);
				T value = (T)f.get(m);
				tags.add(value);
				
			} catch (NoSuchFieldException e) {
				// silently fail. 
			} catch (IllegalAccessException e) {
				e.printStackTrace();
			}
		}
		
		return tags;
	}
	
	public Collection<Model> findModels(double x1, double y1, double x2, double y2) {
		ArrayList<Model> ms = new ArrayList<Model>();
		for(int i = 0; i < points.size(); i++) { 
			ScatterPoint pt = points.get(i);
			Double x = pt.x, y = pt.y;
			if(x >= x1 && y >= y1 && x < x2 && y < y2) { 
				ms.add(models.get(i));
			}
		}
		return ms;
	}

	public void clearModels() {
		points.clear();
		models.clear();

		setProperty(new PropertyValueWrapper<PaintableScale>(xScaleKey, new PaintableScale(0.0, 1.0)));
		setProperty(new PropertyValueWrapper<PaintableScale>(yScaleKey, new PaintableScale(0.0, 1.0)));
		
		dispatchChangedEvent();
	}

	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		PaintableScale xScale = getPropertyValue(xScaleKey);
		PaintableScale yScale = getPropertyValue(yScaleKey);
		
		int radius = getPropertyValue(radiusKey, 3);
		Color color = getPropertyValue(colorKey, Color.red);
		float stroke = getPropertyValue(strokeKey, (float)1.0);
		Color axisColor = getPropertyValue(axisColorKey, Color.black);
		Boolean showScales = getPropertyValue(showScaleKey, Boolean.TRUE);
		
		int diam = radius*2;
		
		int w = x2-x1, h = y2-y1;
		
		Graphics2D g2 = (Graphics2D)g;
		Stroke oldStroke = g2.getStroke();
		g2.setStroke(new BasicStroke(stroke));
		
		g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		
		/** Painting Code **/
		
		VerticalScalePainter vsp = new VerticalScalePainter(yScale);
		HorizontalScalePainter hsp = new HorizontalScalePainter(xScale);
		
		clearDrawnPoints();
		
		// We may need the coordinates for the origin, so let's go ahead and 
		// get their scale and coordinates now. 
		double zxf = xScale.fractionalOffset(0.0);
		double zyf = yScale.fractionalOffset(0.0);
		int zx = x1 + (int)Math.round(zxf * (double)w);
		int zy = y2 - (int)Math.round(zyf * (double)h);
		
		if(showScales) { 
			g.setColor(Color.black);
			g.drawLine(x1, zy, x2, zy);
			hsp.paintItem(g, x1, zy, x2, y2);

			g.setColor(Color.black);
			g.drawLine(zx, y1, zx, y2);
			vsp.paintItem(g, zx, y1, x2, y2);
		}
		
		// Points
		for(int i = 0; i < points.size(); i++) { 
			ScatterPoint p = points.get(i);
			Model m = models.get(i);
			
			double xf = xScale.fractionalOffset(p.x);
			double yf = yScale.fractionalOffset(p.y);
			
			int px = x1 + (int)Math.round(xf * (double)w);
			int py = y2 - (int)Math.round(yf * (double)h);

			g2.setColor(p.color != null ? p.color : color);
			
			if(p.vector) { 
				// Some points are "vectors", and get an arrow drawn for them.
				int dx = px-zx, dy = py-zy;
				
				double len = Math.sqrt((double)(dx*dx + dy*dy));
				double rot = p.y >= 0.0 ? 
						Math.acos((double)dx/len) : Math.PI + Math.acos((double)dx/len);
				
				g2.translate(zx, zy);
				g2.rotate(rot);

				int off = (int)Math.round(len);
				int arrow = Math.max(5, Math.min(15, off/20));
				
				g2.drawLine(0, 0, off, 0);
				g2.drawLine(off, 0, off-arrow, -arrow);
				g2.drawLine(off, 0, off-arrow, arrow);
				
				g2.rotate(-rot);
				g2.translate(-zx, -zy);
				
			} else {
				// Usually, we just draw a normal point.
				g2.drawOval(px-radius, py-radius, diam, diam);
			}
			
			if(m != null) { 
				drawPoint(new Point(px, py), m);
			}
		}
		
		g2.setStroke(oldStroke);
	}
	
	public static class InteractiveFrame extends JFrame {
		
		private InteractivePanel panel;
		private ModelScatter scatter;
		
		public InteractiveFrame(ModelScatter sc, String tagName) { 
			super(tagName);
			
			Container c = (Container)getContentPane();
			c.setLayout(new BorderLayout());
		
			scatter = sc;
			c.add(panel = new InteractivePanel(scatter, tagName), BorderLayout.CENTER);
			
			panel.setPreferredSize(new Dimension(300, 300));
			setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			
			SwingUtilities.invokeLater(new Runnable() {
				public void run() { 
					setVisible(true);
					pack();
				}
			});
		}
	}
	
	public static class InteractivePanel extends PaintablePanel {
		
		private ModelScatter scatter;
		private int lastWidth, lastHeight;
		private String tagName;
		private Point p1, p2;
		
		private Point mouse, models;
		private Set<Model> mouseModels;
		
		public InteractivePanel(ModelScatter sc, String tagn) { 
			//super(new DoubleBufferedPaintable(sc));
			super(new DoubleBufferedPaintable(new OverlayModelPaintable(sc)));
			scatter = sc;
			lastWidth = lastHeight = 1;
			tagName = tagn;
			p1 = p2 = null;

			addMouseListener(new MouseAdapter() {
				public void mousePressed(MouseEvent e) {
					if(tagName != null) { 
						p1 = e.getPoint();
						p2 = p1;
						repaint();
					}
				}
				
				public void mouseReleased(MouseEvent e) {
					if(p1 != null && p2 != null && tagName !=null) {
						System.out.println(String.format("Searching %s, %s", p1.toString(), p2.toString()));
						Set tags = findTags(p1, p2, tagName);
						System.out.println("Found: " + tags.size());
						for(Object t : tags) { 
							System.out.println(t.toString());
						}
					}
				}
				
				public void mouseClicked(MouseEvent e) { 
					p1 = p2 = null;
					repaint();
				}
			});
			
			addMouseMotionListener(new MouseMotionAdapter()  { 
				public void mouseDragged(MouseEvent e) { 
					if(p1 != null) { 
						p2 = e.getPoint();
						repaint();
					}
				}
				
				public void mouseMoved(MouseEvent e) { 
					if(e.getButton() == MouseEvent.NOBUTTON) { 
						mouse = e.getPoint();
						Pair<Point,Set<Model>> nearby = scatter.findNearestDrawnPoint(mouse);
						models = nearby.getFirst();
						mouseModels = nearby.getLast();
						repaint();
						
					} else { 
						mouse = models = null;
						mouseModels = null;
					}
				}
			});
		}
		
		public <T> Set<T> findTags(Point p1, Point p2, String tagName) { 
			int minX = Math.min(p1.x, p2.x), maxX = Math.max(p1.x, p2.x);
			int minY = Math.min(p1.y, p2.y), maxY = Math.max(p1.y, p2.y);
			double x1 = xPixToCoord(minX), x2 = xPixToCoord(maxX);
			double y1 = yPixToCoord(lastHeight-maxY), y2 = yPixToCoord(lastHeight-minY);
			return scatter.findTags(x1, y1, x2, y2, tagName);
		}
		
		private double xPixToCoord(int x) {
			PaintableScale scale = scatter.getPropertyValue(ModelScatter.xScaleKey);
			double f = (double)x / (double)lastWidth;
			double range = scale.getRange();
			return scale.getMin() + f*range;
		}
		
		private double yPixToCoord(int y) { 
			PaintableScale scale = scatter.getPropertyValue(ModelScatter.yScaleKey);
			double f = (double)y / (double)lastHeight;
			double range = scale.getRange();
			return scale.getMin() + f*range;
		}
		
		protected void paintComponent(Graphics g) {
			super.paintComponent(g);
			
			lastWidth = getWidth();
			lastHeight = getHeight();
			
			if(p1 != null && p2 != null) { 
				int x1 = Math.min(p1.x, p2.x), x2 = Math.max(p1.x, p2.x);
				int y1 = Math.min(p1.y, p2.y), y2 = Math.max(p1.y, p2.y);
				g.setColor(Color.blue);
				g.drawRect(x1, y1, x2-x1, y2-y1);
			}
			
			if(mouse != null && models != null) { 
				g.setColor(Color.black);
				g.drawLine(mouse.x, mouse.y, models.x, models.y);
				g.drawString(mouseModels.toString(), mouse.x, mouse.y);
			}
		}
	}
	
	public static class ScatterPoint extends Model { 
		public Double x, y;
		public Boolean vector;
		public Color color;
		public String name;
		
		public ScatterPoint() {}
		
		public ScatterPoint(Double _x, Double _y) { 
			x = _x; y = _y; vector = false;
			color = null;
			name = super.toString();
		}
		
		public ScatterPoint(Double x, Double y, Color c, String n) { 
			this(x, y);
			color = c;
			name = n;
		}
		
		public String toString() { return name; }
	}
}

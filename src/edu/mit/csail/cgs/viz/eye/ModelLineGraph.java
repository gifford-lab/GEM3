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

public class ModelLineGraph extends AbstractModelPaintable {
	
	public static final String xScaleKey = "x-scale";
	public static final String yScaleKey = "y-scale";
	public static final String radiusKey = "radius";
	public static final String colorKey = "color";
	public static final String strokeKey = "stroke";
	public static final String axisColorKey = "axis-color";

	private String xFieldName, yFieldName;
	private Map<Double,LinePoint> points;
	private Map<Double,Model> models;
	
	public ModelLineGraph() { 
		xFieldName = "x"; 
		yFieldName = "y";
		points = new TreeMap<Double,LinePoint>();
		models = new TreeMap<Double,Model>();
		
		initProperty(new PropertyValueWrapper<PaintableScale>(xScaleKey, new PaintableScale(0.0, 1.0)));
		initProperty(new PropertyValueWrapper<PaintableScale>(yScaleKey, new PaintableScale(0.0, 1.0)));
		initProperty(new PropertyValueWrapper<Integer>(radiusKey, 3));
		initProperty(new PropertyValueWrapper<Color>(colorKey, Color.red));
		initProperty(new PropertyValueWrapper<Float>(strokeKey, (float)3.0));
		
		startDrawingPoints();
	}
	
	public ModelLineGraph(String xfield, String yfield) { 
		this();
		xFieldName = xfield;
		yFieldName = yfield;
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
				
					LinePoint sp = new LinePoint(x,y);
					
					points.put(x, sp);
					models.put(x, m);
					
					xScale.updateScale(x);
					yScale.updateScale(y);
					
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
			LinePoint pt = points.get(i);
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
		
		int diam = radius*2;
		
		int w = x2-x1, h = y2-y1;
		
		Graphics2D g2 = (Graphics2D)g;
		Stroke oldStroke = g2.getStroke();
		g2.setStroke(new BasicStroke(stroke));
		
		g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		
		/** Painting Code **/
		
		g2.setColor(Color.white);
		g2.fillRect(x1, y1, w, h);
		
		// Axes
		g2.setColor(axisColor);
		g2.drawRect(x1, y1, w-1, h-1);
		
		VerticalScalePainter vsp = new VerticalScalePainter(yScale);
		HorizontalScalePainter hsp = new HorizontalScalePainter(xScale);
		
		vsp.paintItem(g, x1, y1, x2, y2);
		hsp.paintItem(g, x1, y1, x2, y2);
		
		clearDrawnPoints();
		
		Vector<Point> pts = new Vector<Point>();
		Vector<Double> keys = new Vector<Double>();
		Point previous = null;
		
		// Points
		for(Double key : points.keySet()) { 
			LinePoint p = points.get(key);
			keys.add(key);
			
			double xf = xScale.fractionalOffset(p.x);
			double yf = yScale.fractionalOffset(p.y);
			
			int px = x1 + (int)Math.round(xf * (double)w);
			int py = y2 - (int)Math.round(yf * (double)h);
			
			Point pp = new Point(px, py);
			pts.add(pp);
			
			if(previous != null) { 
				g2.setColor(color);
				g2.drawLine(previous.x, previous.y, pp.x, pp.y);
			}
			
			previous = pp;
		}

		for(int i = 0; i < pts.size(); i++) {
			Point pt = pts.get(i);
			Model m = models.get(keys.get(i));

			int px = pt.x, py = pt.y;
			g2.setColor(color);
			g2.drawOval(px-radius, py-radius, diam, diam);
			
			if(m != null) { 
				drawPoint(new Point(px, py), m);
			}
		}
		
		g2.setStroke(oldStroke);
	}
	
	public static class InteractiveFrame extends JFrame {
		
		private InteractivePanel panel;
		private ModelLineGraph scatter;
		
		public InteractiveFrame(ModelLineGraph sc, String tagName) { 
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
		
		private ModelLineGraph scatter;
		private int lastWidth, lastHeight;
		private String tagName;
		private Point p1, p2;
		
		private Point mouse, models;
		private Set<Model> mouseModels;
		
		public InteractivePanel(ModelLineGraph sc, String tagn) { 
			super(new DoubleBufferedPaintable(sc));
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
			PaintableScale scale = scatter.getPropertyValue(ModelLineGraph.xScaleKey);
			double f = (double)x / (double)lastWidth;
			double range = scale.getRange();
			return scale.getMin() + f*range;
		}
		
		private double yPixToCoord(int y) { 
			PaintableScale scale = scatter.getPropertyValue(ModelLineGraph.yScaleKey);
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
	
	public static class LinePoint extends Model { 
		public Double x, y;
		
		public LinePoint(Double _x, Double _y) { 
			x = _x; y = _y;  
		}
	}
}

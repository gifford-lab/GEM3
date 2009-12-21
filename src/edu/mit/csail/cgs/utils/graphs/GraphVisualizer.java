package edu.mit.csail.cgs.utils.graphs;

import java.util.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.font.FontRenderContext;
import java.awt.font.TextLayout;

import javax.swing.*;
import javax.swing.border.*;

import edu.mit.csail.cgs.utils.Pair;

public class GraphVisualizer extends JPanel {
	
	private Graph graph;
	private Map<String,Object> parameters;
	private LinkedList<Pair<Rectangle,String>> nodeBounds;
	
	private String movingVertex;

	public GraphVisualizer(Graph g) { 
		graph = g;
		parameters = new HashMap<String,Object>();
		parameters.put("showDirected", true);
		nodeBounds = new LinkedList<Pair<Rectangle,String>>();
		movingVertex = null;
		
		rebuildLocations();
		
		addMouseListener(new MouseAdapter() { 
			public void mousePressed(MouseEvent e) { 
				movingVertex = retrieveVertexAtPoint(e.getPoint());
			}
			
			public void mouseReleased(MouseEvent e) { 
				if(movingVertex != null) {
					Point p = (Point)parameters.get("location:" + movingVertex);
					int radius = (Integer)parameters.get("radius:" + movingVertex);
					Rectangle rect = new Rectangle(p.x-radius, p.y-radius, 2*radius, 2*radius);
					nodeBounds.addFirst(new Pair<Rectangle,String>(rect, movingVertex));
					movingVertex = null;
					repaint();
				}
			}
		});
		
		addMouseMotionListener(new MouseMotionAdapter() { 
			public void mouseDragged(MouseEvent e) { 
				if(movingVertex != null) { 
					parameters.put("location:" + movingVertex, e.getPoint());
					repaint();
				}
			}
		});
	}
	
	public JFrame putInFrame() { 
		JFrame f = new JFrame();
		
		Container c = (Container)f.getContentPane();
		c.setLayout(new BorderLayout());
		c.add(this, BorderLayout.CENTER);
		
		f.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		f.setVisible(true);
		f.pack();
		return f;
	}
	
	private String retrieveVertexAtPoint(Point p) { 
		Iterator<Pair<Rectangle,String>> itr = nodeBounds.iterator();
		while(itr.hasNext()) { 
			Pair<Rectangle,String> pair = itr.next();
			if(pair.getFirst().contains(p)) { 
				itr.remove();
				return pair.getLast();
			}
		}
		return null;
	}
	
	public String getVertexAtPoint(Point p) { 
		for(Pair<Rectangle,String> pr : nodeBounds) { 
			if(pr.getFirst().contains(p)) { 
				return pr.getLast();
			}
		}
		return null;
	}
	
	public Dimension getPreferredSize() { 
		return new Dimension(400, 200);
	}
	
	public void rebuildLocations() {
		int w = getWidth(), h = getHeight();
		w = Math.max(w, 100);
		h = Math.max(h, 100);
		
		Random r = new Random();		
		nodeBounds.clear();
		for(String vertex : graph.getVertices()) { 
			String key = "location:" + vertex;
			if(!parameters.containsKey(key)) {
				Point p = new Point(r.nextInt(w), r.nextInt(h));
				parameters.put(key, p);
				parameters.put("radius:" + vertex, 10);
			}
			
			Point p = (Point)parameters.get("location:" + vertex);
			int radius = (Integer)parameters.get("radius:" + vertex);
			Rectangle rect = new Rectangle(p.x-radius, p.y-radius, 2*radius, 2*radius);
			nodeBounds.addLast(new Pair<Rectangle,String>(rect, vertex));
		}
		
	}
	
	protected void paintComponent(Graphics g) { 
		super.paintComponent(g);
		
        Graphics2D g2 = (Graphics2D)g;
        Map oldHints = g2.getRenderingHints();
        
        Map newHints = new HashMap(oldHints);
        newHints.put(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2.setRenderingHints(newHints);
		
		int w = getWidth(), h = getHeight();
		
		g.setColor(Color.white);
		g.fillRect(0, 0, w, h);
		
		boolean showDirected = (Boolean)parameters.get("showDirected");
		for(String head : graph.getVertices()) { 
			for(String tail : graph.getNeighbors(head)) { 
				if(head.compareTo(tail) < 1 || showDirected) { 
					paintEdge((Graphics2D)g, head, tail);
				}
			}
		}
		
		Stroke oldStroke = g2.getStroke();
		g2.setStroke(new BasicStroke((float)2.0));
		
		for(String v : graph.getVertices()) { 
			paintNode(g2, v);
		}
		
		g2.setStroke(oldStroke);
		g2.setRenderingHints(oldHints);
	}
	
	private void paintNode(Graphics2D g, String v) { 
		Point p = (Point)parameters.get("location:" + v);
		int radius = (Integer)parameters.get("radius:" + v);
		
		g.setColor(Color.white);
		g.fillOval(p.x-radius, p.y-radius, 2*radius, 2*radius);
		g.setColor(Color.black);
		g.drawOval(p.x-radius, p.y-radius, 2*radius, 2*radius);
		
		g.drawString(v, p.x-radius+2, p.y);
	}
	
	private static final double twopi = 2.0 * Math.PI;
	
	private void paintEdge(Graphics2D g, String h, String t) { 
		Point p1 = (Point)parameters.get("location:" + h);
		Point p2 = (Point)parameters.get("location:" + t);
		int rad = (Integer)parameters.get("radius:" + h); 

		int dx = p2.x - p1.x, dy = p1.y - p2.y;
		double theta = 0.0;

		if((dx != 0 || dy != 0)) { 
			int cx = (p1.x + p2.x)/2, cy = (p1.y + p2.y)/2;

			double r = Math.sqrt(dx*dx + dy * dy);
			int r2 = (int)Math.round(r / 2.0);

			if(dx != 0) {
				if(dx > 0) {  
					theta = Math.atan((double)dy / (double)dx);
					theta = twopi - theta;
				} else { 
					theta = - Math.atan((double)dy / (double)dx) - Math.PI;
				}
			} else { 
				theta = dy > 0 ? Math.PI / 2.0 : 3.0 * Math.PI / 2.0;
				theta += Math.PI;
			}

			g.translate(cx, cy);
			g.rotate(theta);

			int depth = 3;
			int sep = depth*2;
			int offset = r2/2;

			Stroke oldStroke = g.getStroke();
			g.setStroke(new BasicStroke((float)3.0));
	
			depth = rad;

			g.setColor(Color.cyan);
			int ptx = r2-(int)(rad*1.5);
			g.drawLine(-r2, 0, ptx, 0);
			g.drawLine(ptx, 0, ptx-depth, -depth);
			g.drawLine(ptx, 0, ptx-depth, depth);

			g.setStroke(oldStroke);

			g.rotate(-theta);
			g.translate(-cx, -cy);

		}
	}
}

package edu.mit.csail.cgs.utils.graphs.ui;

import java.util.*;
import java.awt.Graphics;
import java.io.*;
import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

import edu.mit.csail.cgs.utils.graphs.*;

public class Visualizer {
	
	public static void main(String[] args) { 
		DirectedGraph dg = new DirectedGraph();
		Frame f = new Frame(dg);
	}
	
	public static class Frame extends JFrame { 
		
		private JButton in,out,left,right,up,down,addVert,addEdge,recenter;
		private Graph graph;
		private Visualizer viz;
		
		public Frame(Graph g) { 
			super("Graph Frame");
			graph = g;
			viz = new Visualizer(graph, (graph instanceof DirectedGraph));
			in = new JButton("+");
			out = new JButton("-");
			left = new JButton("<<");
			right = new JButton(">>");
			up = new JButton("^^");
			down = new JButton("vv");
			addVert = new JButton("VERT+");
			addEdge = new JButton("EDGE+");
			recenter = new JButton("Arr");
			
			JPanel zoomPanel = new JPanel(); zoomPanel.setLayout(new GridLayout(1, 2));
			zoomPanel.add(out); zoomPanel.add(in);

			JPanel movePanel = new JPanel();
			movePanel.setLayout(new BorderLayout());
			
			movePanel.add(left, BorderLayout.WEST);
			movePanel.add(right, BorderLayout.EAST);
			movePanel.add(up, BorderLayout.NORTH);
			movePanel.add(down,BorderLayout.SOUTH);
			movePanel.add(zoomPanel, BorderLayout.CENTER);
			
			JPanel buttons = new JPanel(); buttons.setLayout(new GridLayout(1, 4));
			buttons.add(movePanel);
			buttons.add(addVert);
			buttons.add(addEdge);
			buttons.add(recenter);
			
			JPanel mainPanel = new JPanel();
			mainPanel.setLayout(new BorderLayout());
			mainPanel.add(buttons, BorderLayout.SOUTH);
			mainPanel.add(new Panel(viz), BorderLayout.CENTER);
			
			Container c = (Container)getContentPane();
			c.setLayout(new BorderLayout());
			c.add(mainPanel, BorderLayout.CENTER);
			
			in.addActionListener(new ActionListener() { 
				public void actionPerformed(ActionEvent e) { 
					viz.addAction(new ScaleAction(ScaleAction.ZOOM_IN));
					repaint();
				}
			});
			out.addActionListener(new ActionListener() { 
				public void actionPerformed(ActionEvent e) { 
					viz.addAction(new ScaleAction(ScaleAction.ZOOM_OUT));
					repaint();					
				}
			});
			up.addActionListener(new ActionListener() { 
				public void actionPerformed(ActionEvent e) { 
					viz.addAction(new ScaleAction(ScaleAction.SHIFT_UP));
					repaint();					
				}
			});
			down.addActionListener(new ActionListener() { 
				public void actionPerformed(ActionEvent e) { 
					viz.addAction(new ScaleAction(ScaleAction.SHIFT_DOWN));
					repaint();					
				}
			});
			left.addActionListener(new ActionListener() { 
				public void actionPerformed(ActionEvent e) { 
					viz.addAction(new ScaleAction(ScaleAction.SHIFT_LEFT));
					repaint();					
				}
			});
			right.addActionListener(new ActionListener() { 
				public void actionPerformed(ActionEvent e) { 
					viz.addAction(new ScaleAction(ScaleAction.SHIFT_RIGHT));					
					repaint();
				}
			});
			addVert.addActionListener(new ActionListener() { 
				public void actionPerformed(ActionEvent e) { 
					int nodeNum = graph.getVertices().size() + 1;
					String name = "Node#" + nodeNum;
					viz.addVertex(name);
					repaint();
				}
			});
			addEdge.addActionListener(new ActionListener() { 
				public void actionPerformed(ActionEvent e) { 
					Random rand = new Random();
					Vector<String> names = new Vector<String>(graph.getVertices());
					String n1 = names.get(rand.nextInt(names.size()));
					String n2 = names.get(rand.nextInt(names.size()));
					viz.addEdge(n1, n2);
					repaint();
				}
			});
			recenter.addActionListener(new ActionListener() { 
				public void actionPerformed(ActionEvent e) { 
					viz.recenterVertices();
					repaint();
				}
			});
			
			setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
			setVisible(true);
			pack();
		}
	}
	
	public static class Panel extends JPanel {
		
		private Visualizer viz;
		private String selected;
		
		public Panel(Visualizer v) { 
			super();
			viz = v;
			selected = null;
			addMouseListener(new MouseAdapter() { 
				public void mousePressed(MouseEvent e) { 
					selected = viz.getVertexAtPoint(e.getPoint());
				}
				
				public void mouseReleased(MouseEvent e) { 
					viz.saveOverrides();
					viz.clearOverrides();
					repaint();
				}
			});
			
			this.addMouseMotionListener(new MouseMotionAdapter() { 
				public void mouseDragged(MouseEvent e) { 
					if(selected != null) { 
						viz.setOverride(selected, e.getPoint());
						repaint();
					}
				}
			});
		}
		
		public Dimension getPreferredSize() { 
			Dimension pd = super.getPreferredSize();
			return new Dimension(Math.max(100, pd.width), Math.max(100, pd.height));
		}
		
		protected void paintComponent(Graphics g) { 
			super.paintComponent(g);
			g.setColor(Color.white);
			int w = getWidth(), h = getHeight();
			g.fillRect(0, 0, w, h);
			viz.paintItem(g, 0, 0, w, h);
		}
	}
	
	private static Random rand;
	
	static { 
		rand = new Random();
	}
	
	private boolean directed;
	private Graph graph;
	private Map<String,VizVertex> vizVertices;
	private Map<String,Point> overridePoints;
	private LinkedList<ScaleAction> actions;
	
	private Point ulPoint;
	private double scale;
	private int lastWidth, lastHeight;
	private ScaleWindow lastWin;
	
	public Visualizer(Graph g, boolean dir) {
		graph = g;
		directed = dir;
		vizVertices = new HashMap<String,VizVertex>();
		overridePoints = new HashMap<String,Point>();
		ulPoint = new Point(0, 0);
		lastWidth = lastHeight = 100;
		lastWin = null;
		scale = 1.0;
		actions = new LinkedList<ScaleAction>();
		int diam = 6;
		for(String vertex : graph.getVertices()) { 
			int x = (int)Math.floor(rand.nextDouble() * (double)lastWidth); 
			int y = (int)Math.floor(rand.nextDouble() * (double)lastHeight); 
			vizVertices.put(vertex, new VizVertex(vertex, x, y, diam));
		}
	}
	
	public void addAction(ScaleAction sa) { actions.addLast(sa); }
	public void removeAction() { if(!actions.isEmpty()) { actions.removeLast(); } }
	
	public void setOverride(String n, Point pt) { overridePoints.put(n, pt); }
	public void clearOverrides() { overridePoints.clear(); }
	
	public void saveOverrides() {
		if(lastWin != null && ulPoint != null) {
			int lx = ulPoint.x, uy = ulPoint.y;
			int rx = lx + lastWidth, ly = uy + lastHeight;
			for(String n : overridePoints.keySet()) { 
				VizVertex vv = vizVertices.get(n);
				Point unpt = lastWin.getUnscaledPoint(overridePoints.get(n), lx, uy, rx, ly);
				vv.setPoint(unpt);
			}
		}
	}
	
	public String getVertexAtPoint(Point pt) { 
		for(String name : vizVertices.keySet()) { 
			VizVertex vv = vizVertices.get(name);
			int rad = vv.getScaledDiam() / 2;
			Point vpt = vv.getLastPoint();
			if(Math.abs(pt.x - vpt.x) <= rad && Math.abs(pt.y - vpt.y) <= rad) { 
				return name;
			}
		}
		return null;
	}
	
	public void recenterVertices() { 
		if(lastWin != null && ulPoint != null) {
			Point pt = lastWin.getUpperLeftPoint();
			int w = lastWin.getWidth(), h = lastWin.getHeight();

			for(String n : vizVertices.keySet()) {
				Point p = new Point(vizVertices.get(n).getX(), vizVertices.get(n).getY());
				if(!lastWin.containsPoint(p)) { 
					int x = rand.nextInt(w) + pt.x;
					int y = rand.nextInt(h) + pt.y;
					vizVertices.get(n).setPoint(new Point(x, y));
				}
			}
		}
	}
	
	public void addVertex(String name) {
		if(directed) { 
			DirectedGraph dg = (DirectedGraph)graph;
			dg.addVertex(name);
		} else { 
			UndirectedGraph ug = (UndirectedGraph)graph;
			ug.addVertex(name);
		}

		int diam = 6;
		int w = lastWin != null ? lastWin.getWidth() : 100;
		int h = lastWin != null ? lastWin.getHeight() : 100;
		int wx = lastWin != null ? lastWin.getUpperLeftPoint().x : 0;
		int wy = lastWin != null ? lastWin.getUpperLeftPoint().y : 0;
		
		int x = wx + (int)Math.floor(rand.nextDouble() * (double)w); 
		int y = wy + (int)Math.floor(rand.nextDouble() * (double)h); 
		vizVertices.put(name, new VizVertex(name, x, y, diam));
	}
	
	public void addEdge(String n1, String n2) { 
		if(directed) { 
			DirectedGraph dg = (DirectedGraph)graph;
			dg.addEdge(n1, n2);
		} else { 
			UndirectedGraph ug = (UndirectedGraph)graph;
			ug.addEdge(n1, n2);
		}		
	}

	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		Graphics2D g2 = (Graphics2D)g;
		int w = x2 - x1, h = y2 - y1;
		lastWidth = w; lastHeight = h;
		
		ScaleWindow win = new ScaleWindow(ulPoint.x, ulPoint.y, ulPoint.x + w, ulPoint.y + h, scale);
		for(ScaleAction action : actions) { 
			win.handleAction(action);
		}
		lastWin = win;
		
		int minDiam = -1;
		
		// This locates each vertex, according to the scale of the window.
		for(String vn : vizVertices.keySet()) { 
			VizVertex vv = vizVertices.get(vn);
			vv.locateVertex(x1, y1, x2, y2, win);
			if(minDiam == -1 || vv.getScaledDiam() < minDiam) { 
				minDiam = vv.getScaledDiam();
			}
		}

		Stroke oldStroke = g2.getStroke();
		float thickness = (float)(Math.max(minDiam / 8, 1));
		if(minDiam > 8) { g2.setStroke(new BasicStroke(thickness)); }
		

		// This draws the edges.
		g.setColor(Color.black);
		for(String vn : vizVertices.keySet()) {
			VizVertex vvn = vizVertices.get(vn);
			Point s = vvn.getLastPoint();
			if(overridePoints.containsKey(vn)) { 
				s = overridePoints.get(vn);
			}
			
			for(String vt : graph.getNeighbors(vn)) { 
				VizVertex vvt = vizVertices.get(vt);
				Point e = vvt.getLastPoint();
				if(overridePoints.containsKey(vt)) { 
					e = overridePoints.get(vt); 
				}
				
				g.drawLine(s.x, s.y, e.x, e.y);
				if(directed) {
					int dx = e.x - s.x, dy = e.y - s.y;
					double dist = Math.sqrt((double)((dx * dx) + (dy * dy)));
					if(dist > 0.0) { 
						double theta = Math.asin((double)dy / dist);
						if(dx < 0) { theta = Math.PI - theta; }

						g2.translate(s.x, s.y);
						g2.rotate(theta);

						int nx = (int)Math.round(dist);
						int targetRad = vvt.getScaledDiam() / 2;
						int ah = targetRad/2;
						int ax = nx - targetRad;
						g.drawLine(ax, 0, ax-ah, ah);
						g.drawLine(ax, 0, ax-ah, -ah);
						//g.drawLine(ax-ah, ah, ax-ah, -ah);

						g2.rotate(-theta);
						g2.translate(-s.x, -s.y);
					}
				}
			}
		}
		
		g2.setStroke(oldStroke);
		
		// This draws the vertices themselves.
		for(String vn : vizVertices.keySet()) { 
			VizVertex vv = vizVertices.get(vn);
			if(overridePoints.containsKey(vn)) { 
				vv.paintAt(g, overridePoints.get(vn));
			} else { 
				vv.paint(g);
			}
		}
	}
}

class VizVertex { 
	
	private int cx, cy;
	private int diam;
	private String name;
	
	private Point lastPoint;
	private int scaledDiam;
	
	public VizVertex(String n, int x, int y, int d) { 
		name = n;
		cx = x; cy = y;
		diam = d;
		lastPoint = null;
		scaledDiam = diam;
	}
	
	public int getX() { return cx; }
	public int getY() { return cy; }
	public int getDiam() { return diam; }
	public String getName() { return name; }
	public void setPoint(Point p) { cx = p.x; cy = p.y; }
	
	public Point getLastPoint() { return lastPoint; }
	public int getScaledDiam() { return scaledDiam; }
	
	public boolean containsSpot(int x, int y) { 
		if(lastPoint != null) { 
			int rad = Math.max(1, scaledDiam / 2);
			return Math.abs(x-lastPoint.x) <= rad && Math.abs(y-lastPoint.y) <= rad;
		}
		return false;
	}
	
	public void locateVertex(int x1, int y1, int x2, int y2, ScaleWindow win) { 
		lastPoint = win.getScaledPoint(new Point(cx, cy), x1, y1, x2, y2);
		scaledDiam = win.getScaledDiameter(diam);		
	}
	
	public void paint(Graphics g) { 
		paintAt(g, lastPoint);
	}
	
	public void paintAt(Graphics g, Point pt) {
		Graphics2D g2 = (Graphics2D)g;
		Stroke oldStroke = g2.getStroke();
		float thickness = (float)Math.max(1, scaledDiam / 10);
		g2.setStroke(new BasicStroke(thickness));
		
		int rad = Math.max(1, scaledDiam/2);
		g.setColor(Color.pink);
		g.fillOval(pt.x-rad, pt.y-rad, scaledDiam, scaledDiam);
		g.setColor(Color.black);
		g.drawOval(pt.x-rad, pt.y-rad, scaledDiam, scaledDiam);
		g2.setStroke(oldStroke);
	}
}

class ScaleAction {
	
	public static final int ZOOM_IN = 0;
	public static final int ZOOM_OUT = 1;
	public static final int SHIFT_LEFT = 2;
	public static final int SHIFT_RIGHT = 3;
	public static final int SHIFT_UP = 4;
	public static final int SHIFT_DOWN = 5;

	public int type;
	public ScaleAction(int t) { 
		type = t;
	}
}

class ScaleWindow { 
	
	private int x1, y1, x2, y2;
	private double scale;
	
	public ScaleWindow(int x1, int y1, int x2, int y2, double sc) { 
		this.x1 = x1;
		this.y1 = y1;
		this.x2 = x2;
		this.y2 = y2;
		scale = sc;
	}
	
	public Point getUpperLeftPoint() { return new Point(x1, y1); }
	public int getWidth() { return x2 - x1; }
	public int getHeight() { return y2 - y1; }
	public boolean containsPoint(Point p) { 
		return p.x >= x1 && p.x < x2 && p.y >= y1 && p.y < y2;
	}
	
	public Point getUnscaledPoint(Point p, int lx, int uy, int rx, int ly) { 
		int ww = x2 - x1, wh = y2 - y1;
		int fw = rx - lx, fh = ly - uy;
		double xf = (double)(p.x-lx) / (double)fw;
		double yf = (double)(p.y-uy) / (double)fh;
		int px = x1 + (int)Math.round(xf * (double)ww);
		int py = y1 + (int)Math.round(yf * (double)wh);
		return new Point(px, py);		
	}
	
	public Point getScaledPoint(Point p, int lx, int uy, int rx, int ly) { 
		int ww = x2 - x1, wh = y2 - y1;
		int fw = rx - lx, fh = ly - uy;
		double xf = (double)(p.x-x1) / (double)ww;
		double yf = (double)(p.y-y1) / (double)wh;
		int px = lx + (int)Math.round(xf * (double)fw);
		int py = uy + (int)Math.round(yf * (double)fh);
		return new Point(px, py);
	}
	
	public void handleAction(ScaleAction sa) { 
		switch(sa.type) { 
		case ScaleAction.ZOOM_IN:
			zoomIn();
			break;
		case ScaleAction.ZOOM_OUT:
			zoomOut();
			break;
		case ScaleAction.SHIFT_LEFT:
			shiftLeft();
			break;
		case ScaleAction.SHIFT_RIGHT:
			shiftRight();
			break;
		case ScaleAction.SHIFT_UP:
			shiftUp();
			break;
		case ScaleAction.SHIFT_DOWN:
			shiftDown();
			break;
		}
	}
	
	public void translate(int dx, int dy) { 
		x1 += dx; x2 += dx;
		y1 += dy; y2 += dy;
	}
	
	public void shiftLeft() { translate(-(x2 - x1) / 4, 0); }
	public void shiftRight() { translate((x2 - x1) / 4, 0); }
	public void shiftUp() { translate(0, -(y2 - y1) / 4); }
	public void shiftDown() { translate(0, (y2 - y1) / 4); }
	
	public void resize(int dx, int dy) { 
		x1 -= dx;
		y1 -= dy;
	}
	
	public void zoomIn() { 
		int w = x2 - x1, h = y2 - y1;
		int w4 = Math.max(1, (int)Math.floor((double)w / 4.0));
		int h4 = Math.max(1, (int)Math.floor((double)h / 4.0));

		if(w4 >= 1 && h4 >= 1) { 
			x1 += w4; x2 -= w4;
			y1 += h4; y2 -= h4;
			int nw = x2 - x1, nh = y2 - y1;
			if(nw >= nh) { 
				scale *= (double)w / (double)nw;
			} else { 				
				scale *= (double)h / (double)nh;
			}
		}
	}
	
	public void zoomOut() { 
		int w = x2 - x1, h = y2 - y1;
		int w2 = Math.max(1, (int)Math.floor((double)w / 2.0));
		int h2 = Math.max(1, (int)Math.floor((double)h / 2.0));
		
		x1 -= w2; x2 += w2;
		y1 -= h2; y2 += h2;
		
		int nw = x2 - x1, nh = y2 - y1;
		if(nw >= nh) { 
			scale *= (double)w / (double)nw;
		} else { 				
			scale *= (double)h / (double)nh;
		}
	}
	
	public int getScaledDiameter(int diam) {
		return (int)Math.round(diam * scale);
	}	
}

/**
 * @author tdanford
 */
package edu.mit.csail.cgs.viz.graphs;

import java.util.*;
import java.io.*;

import edu.mit.csail.cgs.utils.graphs.DirectedGraph;
import edu.mit.csail.cgs.viz.paintable.*;
import edu.mit.csail.cgs.viz.utils.FileChooser;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;

public class GraphView extends ObjectView {
	
	public static void main(String[] args) {
		GraphView gv = new GraphView();
				
		NodeView n1 = gv.createNode();
		NodeView n2 = gv.createNode();
		n1.setOption("width", 50);
		n1.setName("node1");
		n1.setX(50);n1.setY(50);
		n2.setOption("width", 20);
		n2.setName("node2");
		n2.setX(200);n2.setY(200);
		EdgeView e1 = gv.createEdge(n1, n2);
		//e1.setDirected(true);
		e1.setOption("edgeWidth", new Integer(40));
		//e1.setOption("arrowSize", new Integer(15));
		e1.setOption("curved", new Boolean(true));
		e1.setOption("name", "hello");
		
		SubEdgeView s1 = gv.createSubEdge(n1, n2);
		s1.setTotalWidth(40);
		s1.addSubEdge(10, Color.green);
		s1.setOption("curved", new Boolean(true));
		s1.addSubEdge(15, Color.yellow);
		
		
		NodeView n = gv.createNode();
		n.setOption("width", 20);
		n.setX(100);n.setY(300);
		
		EdgeView e = gv.createEdge(n,n);
		e.setDirected(true);
		double angle = (2*Math.PI)/(8);
		e.setOption("selfAngle", new Double(angle-Math.toRadians(90)));
		e.setOption("edgeWidth", new Integer(40));
		e.setOption("arrowSize", new Integer(20));
		
		SubEdgeView s2 = gv.createSubEdge(n, n);
		s2.setTotalWidth(40);
		s2.addSubEdge(10, Color.green);
		s2.setOption("selfAngle", new Double(angle-Math.toRadians(90)));
		s2.addSubEdge(15, Color.yellow);

		TestGraphFrame f = new TestGraphFrame(gv.createInteractive());
	}
	
	private LinkedList<NodeView> nodes;
	private LinkedList<EdgeView> edges;
	private LinkedList<SubEdgeView> subedges;
	
	private ObjectView defaultEdgeOptions, defaultNodeOptions;

	public GraphView() {
		super();
		nodes = new LinkedList<NodeView>();
		edges = new LinkedList<EdgeView>();
		subedges = new LinkedList<SubEdgeView>();
		defaultEdgeOptions = new ObjectView();
		defaultNodeOptions = new ObjectView();
		
		defaultNodeOptions.setOption("width", 20);
		defaultEdgeOptions.setOption("color", Color.red);
        defaultEdgeOptions.setOption("arrowSize", 10);
	}
	
	public GraphView(DirectedGraph g) { 
		this();
		
		Map<String,NodeView> views = new HashMap<String,NodeView>();
		for(String v : g.getVertices()) { 
			views.put(v, createNode());
			views.get(v).setName(v);
		}
		
		for(String v1 : g.getVertices()) { 
			for(String v2 : g.getNeighbors(v1)) { 
				EdgeView ev = createEdge(views.get(v1), views.get(v2));
				ev.setDirected(true);
			}
		}
	}

	public void clear() {
		nodes.clear();
		edges.clear();
		subedges.clear();
	}

	public void removeNode(NodeView n) { 
		nodes.remove(n);
	}
	
	public void removeEdge(EdgeView v) { 
		edges.remove(v);
	}
	
	public void removeSubEdge(SubEdgeView v) { 
		subedges.remove(v);
	}
	
	public void removeTextView(TextView v) { 
		removeNode(v);
	}
	
	public Collection<NodeView> nodes() { return nodes; }
	public Collection<EdgeView> edges() { return edges; }
	public Collection<SubEdgeView> subEdges() { return subedges; }
	
	public NodeView findTopNode(int x, int y) {
		Point p = new Point(x, y);
		for(NodeView nv : nodes) { 
			if(nv.containsPoint(p)) { 
				return nv;
			}
		}
		return null;
	}
	
	public void setDefaultNodeOption(String k, Object v) { defaultNodeOptions.setOption(k, v); }
	public void setDefaultEdgeOption(String k, Object v) { defaultEdgeOptions.setOption(k, v); }

	public NodeView createNode() { 
		NodeView nv = new NodeView(defaultNodeOptions, this);
		nodes.addLast(nv);
		return nv;
	}
	
	public TextView createText() { 
		TextView tv = new TextView(defaultNodeOptions, this);
		nodes.addLast(tv);
		return tv;
	}
	
	public TextView createText(String para) { 
		TextView tv = new TextView(defaultNodeOptions, this, para);
		nodes.addLast(tv);
		return tv;
	}
	
	public EdgeView createEdge(NodeView start, NodeView end) {
		if(start.getGraph() != this || end.getGraph() != this) { 
			throw new IllegalArgumentException();
		}
		EdgeView ev = new EdgeView(defaultEdgeOptions, this, start, end);
		start.addEdge(ev);
		edges.addLast(ev);
		return ev;
	}
	
	public SubEdgeView createSubEdge(NodeView start, NodeView end) {
		if(start.getGraph() != this || end.getGraph() != this) { 
			throw new IllegalArgumentException();
		}
		SubEdgeView sev = new SubEdgeView(defaultEdgeOptions, this, start, end);
		start.addEdge(sev);
		subedges.addLast(sev);
		return sev;
	}
	
	public Paintable createPaintable() { return new GraphPaintable(this); }
	
	public InteractiveGraphPanel createInteractive() { 
		return new InteractiveGraphPanel(this); 
	}
	
	public static class TestGraphFrame extends JFrame {
		
		private InteractiveGraphPanel p;
		private JButton addNode, addEdges, twistSelf;
		
		public TestGraphFrame(InteractiveGraphPanel pn) {
			super("Interactive Graph");
			p = pn;
            
            setJMenuBar(createMenuBar());
            
			Container c = (Container)getContentPane();
			c.setLayout(new BorderLayout());
			c.add(p, BorderLayout.CENTER);
			JPanel buttonPanel = new JPanel();
			
			buttonPanel.setLayout(new GridLayout(1, 2));
			buttonPanel.add(addNode = new JButton("+ NODE"));
			buttonPanel.add(addEdges = new JButton("+ EDGES"));
			buttonPanel.add(twistSelf = new JButton("TWIST-SELF"));
			c.add(buttonPanel, BorderLayout.SOUTH);
			
			addNode.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent arg0) {
					p.createNode();
				} 
			});
			
			addEdges.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent arg0) {
					p.createHighlightedEdges();
				} 
			});
			twistSelf.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent arg0) {
					p.twistSelfEdge();
				} 
			});
			setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
			setVisible(true);
			pack();
			setLocation(getX() + 100, getY() + 100);
		}
        
        private JMenuBar createMenuBar() { 
            JMenuBar bar = new JMenuBar();
            JMenu menu;
            JMenuItem item;
            
            bar.add(menu = new JMenu("Image"));
            menu.add(item = new JMenuItem("Save Snapshot..."));
            
            item.addActionListener(new ActionListener() { 
                public void actionPerformed(ActionEvent e) { 
                    p.saveImage();
                }
            });
            
            return bar;
        }
	}

}

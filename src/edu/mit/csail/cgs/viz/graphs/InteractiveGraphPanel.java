/*
 * Author: tdanford
 * Date: May 27, 2008
 */
package edu.mit.csail.cgs.viz.graphs;

import java.util.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

import edu.mit.csail.cgs.viz.utils.FileChooser;

import java.io.*;

public class InteractiveGraphPanel extends JPanel {
	
	private GraphView view;
	private GraphPaintable gp;
	private NodeView selectedNode;
	private Vector<NodeView> highlightedNodes;
	private Boolean displayPath=false;
	
	public InteractiveGraphPanel(GraphView v) { 
		super();
		view = v;
		gp = new GraphPaintable(view);
		selectedNode = null;
		highlightedNodes = new Vector<NodeView>();
		
		addMouseMotionListener(new MouseMotionAdapter() {
			public void mouseDragged(MouseEvent arg0) {
				moved(arg0);
			}
		});
		
		addMouseListener(new MouseAdapter() {
			public void mouseClicked(MouseEvent e) { 
				clicked(e);
			}
			
			public void mousePressed(MouseEvent arg0) {
				pressed(arg0);
			}

			public void mouseReleased(MouseEvent arg0) {
				released(arg0);
			}
		});
	}
    
    public void saveImage() { 
        FileChooser chooser = new FileChooser(null);
        int w = getWidth(), h = getHeight();
        File f = chooser.choose();
        try {
            gp.saveImage(f, w, h, true);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
	
	public void createNode() {
		Random r = new Random();
		NodeView nv = view.createNode();
		nv.setX(r.nextInt(getWidth()));
		nv.setY(r.nextInt(getHeight()));
		repaint();
	}
	
	public void createText() {
		Random r = new Random();
		TextView tv = view.createText();
		tv.setX(r.nextInt(getWidth()));
		tv.setY(r.nextInt(getHeight()));
		repaint();
	}
	
	public void createHighlightedEdges() { 
		for(int i = 0; i < highlightedNodes.size(); i++) {
			NodeView nv1 = highlightedNodes.get(i);
			for(int j = i + 1; j < highlightedNodes.size(); j++) { 
				NodeView nv2 = highlightedNodes.get(j);
				EdgeView ev = view.createEdge(nv1, nv2);
                ev.setDirected(true);
			}
		}
		if(highlightedNodes.size() > 0) { 
			repaint();
		}
	}
	
	public void twistSelfEdge() { 
		for(int i = 0; i < highlightedNodes.size(); i++) {
			NodeView nv = highlightedNodes.get(i);
			for(int e=0; e<nv.getNumEdges(); e++){
				EdgeView ev =nv.getEdge(e);
				if(ev.containsOption("self") && (Boolean)ev.getOption("self")==true){
					ev.twistSelfEdge(Math.PI/4);						
				}
			}
		}
		if(highlightedNodes.size() > 0) { 
			repaint();
		}
	}
	
	private void clicked(MouseEvent e) { 
		if(e.getButton() == MouseEvent.BUTTON3) { 
			if(e.getClickCount() == 1) { 
				NodeView nv = view.findTopNode(e.getX(), e.getY());
				//if(nv != null && !highlightedNodes.contains(nv)) {
				if(nv != null) {
					highlightedNodes.add(nv);
					repaint();
				}
			} else if (e.getClickCount() == 2) { 
				highlightedNodes.clear();
				repaint();
			}
		}
	}
	
	private void pressed(MouseEvent e) { 
		if(e.getButton() == MouseEvent.BUTTON1) { 
			int x = e.getX(), y = e.getY();
			selectedNode = view.findTopNode(x, y);
		}
	}
	
	private void released(MouseEvent e) { 
		if(e.getButton() == MouseEvent.BUTTON1) { 
			selectedNode = null;
		}
	}
	
	private void moved(MouseEvent e) { 
		if(selectedNode != null) { 
			selectedNode.setX(e.getX());
			selectedNode.setY(e.getY());
            
            for(int i = 0; i < selectedNode.getNumEdges(); i++) { 
                selectedNode.getEdge(i).clearDynamicAttributes();
            }
            
            for(int i = 0; i < selectedNode.getNumIncomingEdges(); i++) { 
                selectedNode.getIncomingEdge(i).clearDynamicAttributes();
            }
            
			repaint();
		}
	}
	
	public void paintComponent(Graphics g) { 
		super.paintComponent(g);
		gp.paintItem(g, 0, 0, getWidth(), getHeight());
		for(NodeView nv : highlightedNodes) { 
			gp.paintHightlightedNode(g, nv);
		}
	}
	
	public Vector<NodeView> getHighlightedNodes(){return highlightedNodes;}
	public void unHighlight(){highlightedNodes.clear();}
	public void setDisplayPath(Boolean dp){displayPath=dp;}
	
	public NodeView getTopNode(int x, int y) {
		return view.findTopNode(x, y);
	}
		
}
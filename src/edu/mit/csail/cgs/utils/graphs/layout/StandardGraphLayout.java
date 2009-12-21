package edu.mit.csail.cgs.utils.graphs.layout;

import java.util.*;
import java.awt.*;

import edu.mit.csail.cgs.utils.graphs.*;

public class StandardGraphLayout implements GraphLayout<DirectedGraph> { 

	private DirectedGraph graph;
	private Rectangle bounds;
	private Map<String,Point> nodeLocations;
	private Map<String,Object> parameters;

	public StandardGraphLayout(DirectedGraph dg, Map<String,Point> locs) { 
		graph = dg;
		bounds = new Rectangle(0, 0, 0, 0);
		nodeLocations = new HashMap<String,Point>(locs);
		parameters = new HashMap<String,Object>();
	}

	public void setScaling(boolean s) { 
		parameters.put("scale?", s);
	}

	public boolean isScaling() { 
		return parameters.containsKey("scale?") ? (Boolean)parameters.get("scale?") : true; 
	}

	public DirectedGraph getGraph() { return graph; }

	public void displayGraph(Graphics2D g2, Rectangle bounds) { 
		for(String node : graph.getVertices()) { 
			for(String target : graph.getNeighbors(node)) { 
				displayEdge(node, target, g2, bounds);
			}
		}

		for(String node : graph.getVertices()) { 
			displayNode(node, g2, bounds);
		}
	}

	public void displayNode(String node, Graphics2D g2, Rectangle bounds) { 
	}

	public void displayEdge(String head, String tail, Graphics2D g2, Rectangle bounds) { 
	}
}


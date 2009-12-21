package edu.mit.csail.cgs.utils.graphs.layout;

import java.util.*;
import java.awt.*;

import edu.mit.csail.cgs.utils.graphs.*;

public interface GraphLayout<G extends Graph> { 
	public G getGraph();

	public void displayGraph(Graphics2D g2, Rectangle bounds);
	public void displayNode(String node, Graphics2D g2, Rectangle bounds);
	public void displayEdge(String head, String tail, Graphics2D g2, Rectangle bounds);
}


/*
 * Created on Feb 14, 2008
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.viz.graphs;

import java.util.*;

import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.graphs.*;

public class UndirectedGraphAdapter {

    private UndirectedGraph graph;
    private GraphView view;
    
    private Map<String,NodeView> nodeViews;
    private Map<Pair<String,String>,EdgeView> edgeViews;
    
    public UndirectedGraphAdapter(UndirectedGraph dg) { 
        graph = dg;
        view = new GraphView();
        
        nodeViews = new TreeMap<String,NodeView>();
        edgeViews = new HashMap<Pair<String,String>,EdgeView>();
        
        buildViews();
    }
    
    private void buildViews() { 
        for(String node : graph.getVertices()) { 
            NodeView nview = view.createNode();
            nview.setOption("name", node);
            nodeViews.put(node, nview);
        }
        
        for(String from : graph.getVertices()) {
            NodeView fromView = nodeViews.get(from);
            
            for(String to : graph.getNeighbors(from)) {
            	if(from.compareTo(to) <= 0) { 
            		NodeView toView = nodeViews.get(to);
                
                	EdgeView edge = view.createEdge(fromView, toView);
                	edge.setDirected(false);
                	edgeViews.put(new Pair<String,String>(from, to), edge);
            	}
            }
        }
    }
    
    public UndirectedGraph getGraph() { 
        return graph; 
    }

    public GraphView getView() { 
        return view;
    }
}

/*
 * Author: tdanford
 * Date: Nov 15, 2008
 */
package edu.mit.csail.cgs.utils.graphs;

import java.util.*;

public class DirectedWeightedGraph extends DirectedGraph implements WeightedGraph {
	
	private Double defaultWeight;
	private Map<String,Double> nodeWeights;
	private Map<String,Map<String,Double>> edgeWeights;
	
	public DirectedWeightedGraph() { 
		this(1.0);
	}

	public DirectedWeightedGraph(Double defaultWeight) { 
		super();
		this.defaultWeight = defaultWeight;
		nodeWeights = new HashMap<String,Double>();
		edgeWeights = new HashMap<String,Map<String,Double>>();
	}
	
	public DirectedWeightedGraph(DirectedGraph g, Double dweight) { 
		super(g);
		defaultWeight = dweight;
		nodeWeights = new HashMap<String,Double>();
		edgeWeights = new HashMap<String,Map<String,Double>>();
	}
	
	public void removeVertex(String node) { 
		super.removeVertex(node);
		if(nodeWeights.containsKey(node)) { nodeWeights.remove(node); }
		if(edgeWeights.containsKey(node)) { edgeWeights.remove(node); }
		for(String n : edgeWeights.keySet()) { 
			if(edgeWeights.get(n).containsKey(node)) { 
				edgeWeights.get(n).remove(node);
			}
		}
	}
	
	public void removeEdge(String n1, String n2) { 
		super.removeEdge(n1, n2);
		if(edgeWeights.containsKey(n1) && edgeWeights.get(n1).containsKey(n2)) { 
			edgeWeights.get(n1).remove(n2);
			if(edgeWeights.get(n1).isEmpty()) { 
				edgeWeights.remove(n1);
			}
		}
	}
	
	public void removeEdges(String n1) { 
		if(edgeWeights.containsKey(n1)) { 
			edgeWeights.remove(n1);
		}
		
		for(String n2 : getNeighbors(n1)) { 
			removeEdge(n1, n2);
		}
	}
	
	public void removeAllEdges() { 
		super.removeAllEdges();
		edgeWeights.clear();
	}
	
	public void addWeightedEdge(String n1, String n2, Double w) { 
		addEdge(n1, n2);
		setWeight(n1, n2, w);
	}
	
	public void setWeight(String node, Double w) {
		if(!getVertices().contains(node)) { throw new IllegalArgumentException(); }
		nodeWeights.put(node, w);
	}
	
	public void setWeight(String n1, String n2, Double w) { 
		if(!isNeighbor(n1, n2)) { throw new IllegalArgumentException(); }
		if(!edgeWeights.containsKey(n1)) { edgeWeights.put(n1, new HashMap<String,Double>()); }
		edgeWeights.get(n1).put(n2, w);
	}

	public Double weight(String node) {
		return nodeWeights.containsKey(node) ? nodeWeights.get(node) : defaultWeight;
	}

	public Double weight(String n1, String n2) {
		if(edgeWeights.containsKey(n1) && edgeWeights.get(n1).containsKey(n2)) { 
			return edgeWeights.get(n1).get(n2);
		} else { 
			return defaultWeight;
		}
	}
}

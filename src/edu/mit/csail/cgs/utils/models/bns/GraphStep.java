/*
 * Author: tdanford
 * Date: Dec 4, 2008
 */
package edu.mit.csail.cgs.utils.models.bns;

import edu.mit.csail.cgs.utils.graphs.DirectedGraph;

/**
 * A GraphStep is a reversible transformation to a DirectedGraph.  
 * 
 * Finding a graph step for a graph is a necessary part of doing the random graph search, 
 * which will let us determine structure.
 * 
 * @author tdanford
 */
public interface GraphStep {
	public void forward(DirectedGraph g);
	public void reverse(DirectedGraph g);
}

class FlipEdgeStep implements GraphStep {
	
	private String v1, v2;
	
	public FlipEdgeStep(String vf, String vs) { 
		v1 = vf; v2 = vs;
	}

	public void forward(DirectedGraph g) {
		if(!g.containsEdge(v1, v2)) { throw new IllegalArgumentException(); }
		g.removeEdge(v1, v2);
		g.addEdge(v2, v1);
	}
	
	public void reverse(DirectedGraph g) {
		if(!g.containsEdge(v2, v1)) { throw new IllegalArgumentException(); }
		g.removeEdge(v2, v1);
		g.addEdge(v1, v2);
	}
	
	public String toString() { 
		return String.format("FLIP %s->%s", v1, v2);
	}
}

class AddEdgeStep implements GraphStep {
	
	private String v1, v2;
	
	public AddEdgeStep(String vf, String vs) { 
		v1 = vf; v2 = vs;
	}

	public void forward(DirectedGraph g) {
		if(g.containsEdge(v1, v2)) { throw new IllegalArgumentException(); }
		g.addEdge(v1, v2);
	}

	public void reverse(DirectedGraph g) {
		if(!g.containsEdge(v1, v2)) { throw new IllegalArgumentException(); }
		g.removeEdge(v1, v2);
	}
	public String toString() { 
		return String.format("ADD %s->%s", v1, v2);
	}
}

class RemoveEdgeStep implements GraphStep {
	
	private String v1, v2;
	
	public RemoveEdgeStep(String vf, String vs) { 
		v1 = vf; v2 = vs;
	}

	public void reverse(DirectedGraph g) {
		if(g.containsEdge(v1, v2)) { throw new IllegalArgumentException(); }
		g.addEdge(v1, v2);
	}

	public void forward(DirectedGraph g) {
		if(!g.containsEdge(v1, v2)) { throw new IllegalArgumentException(); }
		g.removeEdge(v1, v2);
	}
	public String toString() { 
		return String.format("REMOVE %s->%s", v1, v2);
	}
}

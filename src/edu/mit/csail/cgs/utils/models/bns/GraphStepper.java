/*
 * Author: tdanford
 * Date: Dec 4, 2008
 */
package edu.mit.csail.cgs.utils.models.bns;

import java.util.*;
import edu.mit.csail.cgs.utils.graphs.*;

public interface GraphStepper {
	public GraphStep step(DirectedGraph g);
}

/**
 * Tries to add, delete, or flip an edge.
 * 
 * @author tdanford
 */
class SimpleGraphStepper implements GraphStepper {
	
	private DirectedAlgorithms algos;
	private Random rand;
	
	public SimpleGraphStepper() { 
		algos = null;
		rand = new Random();
	}
	
	public GraphStep step(DirectedGraph g) {
		algos = new DirectedAlgorithms(g);
		GraphStep step = null;
		
		while((step = makeStep(g)) == null) { 
			// repeat.
		}
		
		return step;
	}
	
	private GraphStep makeStep(DirectedGraph g) {
		Vector<String> topo = algos.getTopologicalOrdering();
		String v1 = topo.get(rand.nextInt(topo.size()));
		String v2 = topo.get(rand.nextInt(topo.size()));
		
		if(!v1.equals(v2)) { 
			if(g.containsEdge(v1, v2)) {
				int action = rand.nextInt(2);
				if(action == 0) { 
					// we're going to delete the edge.
					return new RemoveEdgeStep(v1, v2);
					
				} else { 
					// we'll try to flip the edge.
					g.removeEdge(v1, v2);
					if(!g.getAncestors(v2).contains(v1)) {
						 g.addEdge(v1, v2);
						 return new FlipEdgeStep(v1, v2);
						 
					} else { 
						g.addEdge(v1, v2);
						return null;
					}
				}

			} else { 
				// try to add v1 to v2
				if(!g.getAncestors(v1).contains(v2)) {
					return new AddEdgeStep(v1, v2);
				}
			}
		}
		
		return null;
	}
	
	
}

/**
 * 
 */
package edu.mit.csail.cgs.utils.graphs;

import java.util.*;
import java.io.*;

/**
 * @author Timothy Danford
 *
 */
public class DirectedAlgorithms extends Algorithms {

	private DirectedGraph dgraph;
	
	public DirectedAlgorithms(DirectedGraph dg) { super(dg); dgraph = dg; }
	
	public Set<String> getAncestors(Set<String> vs) { 
		HashSet<String> ancs = new HashSet<String>();
		for(String v : vs) { 
			ancs.addAll(dgraph.getAncestors(v));
		}
		return ancs;
	}
	
	public boolean hasCycle() { 
		for(String v : dgraph.getVertices()) { 
			for(String a : dgraph.getAncestors(v)) { 
				if(dgraph.getAncestors(a).contains(v)) { 
					return true;
				}
			}
		}
		return false;
	}
	
	public Vector<String> getTopologicalOrdering() { 
		Vector<String> order = new Vector<String>();
		HashSet<String> seen = new HashSet<String>();
		Set<String> roots = dgraph.getRoots();
		order.addAll(roots);
		seen.addAll(roots);
		LinkedList<String> remaining = new LinkedList<String>(tools.subtract(dgraph.getVertices(), seen));
		
		while(!remaining.isEmpty()) { 
			String first = remaining.removeFirst();
			if(tools.subtract(dgraph.getParents(first), seen).isEmpty()) {  
				seen.add(first);
				order.add(first);
			} else { 
				remaining.addLast(first);
			}
		}
		
		return order;
	}

	

}

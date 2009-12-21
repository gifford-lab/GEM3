/**
 * 
 */
package edu.mit.csail.cgs.utils.graphs;

import java.util.*;
import java.io.*;

/**
 * @author Timothy Danford
 */
public class DirectedCycleChecker implements CycleChecker {
	
	public static void main(String[] args) {
		File inputFile = new File(args[0]);
		Parser p = new Parser(inputFile);
		try { 
			DirectedGraph g = p.parseDirectedGraph();
			CycleChecker checker = new DirectedCycleChecker(g);
			System.out.println(args[0] + " " + 
					(checker.containsCycle() ? "is cyclic" : "isn't cyclic"));
		} catch(IOException ie) { 
			ie.printStackTrace(System.err);
		}
	}

	private DirectedGraph dgraph;
	private Map<String,Set<String>> ancestors;
	
	public DirectedCycleChecker(DirectedGraph dg) { 
		dgraph = dg;
		ancestors = new HashMap<String,Set<String>>();
		rebuild();
	}
	
	public void rebuild() { 
		ancestors.clear();
		for(String vertex : dgraph.getVertices()) {
			ancestors.put(vertex, dgraph.getAncestors(vertex));
		}
	}
	
	/* (non-Javadoc)
	 * @see edu.mit.csail.cgs.utils.graphs.CycleChecker#containsCycle()
	 */
	public boolean containsCycle() { 
		for(String vertex : ancestors.keySet()) { 
			if(ancestors.get(vertex).contains(vertex)) { 
				return true;
			}
		}
		return false;
	}
	
	/* (non-Javadoc)
	 * @see edu.mit.csail.cgs.utils.graphs.CycleChecker#checkCycleWithEdge(java.lang.String, java.lang.String)
	 */
	public boolean checkCycleWithEdge(String v1, String v2) { 
		return containsCycle() || ancestors.get(v1).contains(v2);
	}	
}

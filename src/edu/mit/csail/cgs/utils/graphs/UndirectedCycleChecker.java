package edu.mit.csail.cgs.utils.graphs;

import java.util.*;

public class UndirectedCycleChecker implements CycleChecker {
	
	private UndirectedAlgorithms algs;
	private UndirectedGraph graph;
	private Map<String,Set<String>> connected;
	private GraphSearch searcher;
	private boolean cyclic;
	
	public UndirectedCycleChecker(UndirectedGraph ug) { 
		graph = ug;
		searcher = new GraphSearch(graph);
		algs = new UndirectedAlgorithms(graph);
		connected = new HashMap<String,Set<String>>();
		rebuild();
	}

	public boolean checkCycleWithEdge(String v1, String v2) {
		return connected.get(v2).contains(v1) || containsCycle();
	}

	public boolean containsCycle() {
		return cyclic;
	}

	public void rebuild() { 
		connected.clear();	
		cyclic = false;
		for(String v : graph.getVertices()) { 
			TotalNeighborSearch tns = new TotalNeighborSearch(v);
			searcher.ConnectedBFS(v, tns);
			connected.put(v, tns.getTotal());
			cyclic = cyclic || connected.get(v).contains(v);
		}
	}

	private static class TotalNeighborSearch implements SearchInterface { 
		private String start;
		private Set<String> total;

		public TotalNeighborSearch(String s) { 
			start = s;
			total = new HashSet<String>();
		}

		public String getStartVertex() { return start; }
		public boolean isCyclic() { return total.contains(start); }
		public Set<String> getTotal() { return total; }

		public boolean searchNode(Graph g, String node) { 
			total.add(node);
			return true;
		}
	}
}

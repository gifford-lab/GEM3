package edu.mit.csail.cgs.utils.graphs;

import java.util.*;
import java.io.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.models.Model;

public class UndirectedGraph implements Graph { 

	private HashSet<String> nodes;
	private Map<String,Set<String>> edges;
	private static SetTools<String> tools;
	
	static { 
		tools = new SetTools<String>();
	}
	
	public UndirectedGraph(UndirectedGraph ug) { 
		nodes = new HashSet<String>(ug.nodes);
		edges = new HashMap<String,Set<String>>();
		for(String k : ug.edges.keySet()) { 
			edges.put(k, new HashSet<String>(ug.edges.get(k))); 
		}
	}

	public UndirectedGraph() { 
		nodes = new HashSet<String>();
		edges = new HashMap<String,Set<String>>();
	}
	
	public UndirectedGraph(GraphModel m) { 
		this();
		for(int i = 0; i < m.nodes.length; i++) { 
			addVertex(m.nodes[i]);
		}
		
		for(int i = 0; i < m.edges.length; i++) { 
			addEdge(m.edges[i][0], m.edges[i][1]);
		}
	}
	
	public GraphModel asModel() { return new GraphModel(this); }
	
	public UndirectedGraph(Collection<String> n, Map<String,Set<String>> e) { 
		nodes = new HashSet<String>(n);
		edges = new HashMap<String,Set<String>>();

		for(String node : nodes) { 
			edges.put(node, new HashSet<String>());
		}

		for(String v1 : e.keySet()) { 
			for(String v2 : e.get(v1)) { 
				if(nodes.contains(v1) && nodes.contains(v2)) { 
					addEdge(v1, v2);
				}
			}
		}
	}

	public void printGraph(PrintStream ps) { 
		TreeSet<String> nodesOrdered = new TreeSet<String>(nodes);
		for(String n : nodesOrdered) { 
			ps.print(n + " : ( ");
			for(String nn : edges.get(n)) { 
				ps.print(nn + " ");
			}
			ps.println(")");
		}
	}
	
	public boolean isNeighbor(String v, String n) { 
		return edges.get(v).contains(n);
	}

	public void addVertex(String v) { 
		if(!nodes.contains(v)) { 
			nodes.add(v);
			edges.put(v, new HashSet<String>());
		}
	}

	public void addEdge(String v1, String v2) { 
		if(nodes.contains(v1) && nodes.contains(v2)) { 
			edges.get(v1).add(v2);
			edges.get(v2).add(v1);
		}
	}
	
	public boolean containsEdge(String v1, String v2) { 
		return nodes.contains(v1) && nodes.contains(v2) && edges.get(v1).contains(v2);
	}
	
	public void removeEdge(String v1, String v2) { 
		if(nodes.contains(v1) && nodes.contains(v2)) { 
			edges.get(v1).remove(v2);
			edges.get(v2).remove(v1);
		}
	}
	
	public void removeVertex(String n) { 
		nodes.remove(n);
		edges.remove(n);
		for(String k : edges.keySet()) { edges.get(k).remove(n); }
	}

	public Set<String> getVertices() { 
		return new HashSet<String>(nodes);
	}

	public Set<String> getNeighbors(String vertex) { 
		return edges.get(vertex);
	}

	public Graph getSubgraph(Set<String> vs) { 
		return new UndirectedGraph(vs, edges);
	}
}


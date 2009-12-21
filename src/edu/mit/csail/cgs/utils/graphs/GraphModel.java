/*
 * Author: tdanford
 * Date: Apr 19, 2009
 */
package edu.mit.csail.cgs.utils.graphs;

import java.util.*;

import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.models.Model;

public class GraphModel extends Model {

	public String[] nodes; 
	public String[][] edges;
	
	private Map<String,Set<String>> neighbors; 
	
	public GraphModel() { 
		nodes = null;
		edges = null;
		neighbors = null;
	}
	
	public GraphModel(UndirectedGraph g) { 
		nodes = g.getVertices().toArray(new String[0]);
		Arrays.sort(nodes);
		int ne = 0;
		for(int i = 0; i < nodes.length; i++) { 
			String n1 = nodes[i];
			Set<String> nbs = g.getNeighbors(n1);
			for(String n2 : nbs) { 
				if(n1.compareTo(n2) <= 0) { 
					ne += 1;
				}
			}
		}
		
		edges = new String[ne][];
		int k = 0;
		for(int i = 0; i < nodes.length; i++) { 
			String n1 = nodes[i];
			Set<String> nbs = g.getNeighbors(n1);
			for(String n2 : nbs) { 
				if(n1.compareTo(n2) <= 0) {
					edges[k++] = new String[] { n1, n2 };
				}
			}
		}
	}
	
	public GraphModel(DirectedGraph g) { 
		nodes = g.getVertices().toArray(new String[0]);
		Arrays.sort(nodes);
		int ne = 0;
		for(int i = 0; i < nodes.length; i++) { 
			String n1 = nodes[i];
			Set<String> nbs = g.getNeighbors(n1);
			ne += nbs.size();
		}
		
		edges = new String[ne][];
		int k = 0;
		for(int i = 0; i < nodes.length; i++) { 
			String n1 = nodes[i];
			Set<String> nbs = g.getNeighbors(n1);
			for(String n2 : nbs) { 
				edges[k++] = new String[] { n1, n2 };
			}
		}
	}
	
	public GraphModel(String[] es) { 
		this();
		neighbors = new TreeMap<String,Set<String>>();
		Set<String> ns = new TreeSet<String>();
		int ne = 0;
		for(String e : es) {
			String n1, n2;
			String[] array = e.split(",");
			n1 = array[0];
			n2 = array[1];
			
			ns.add(n1); 
			ns.add(n2);
			if(!neighbors.containsKey(n1)) { neighbors.put(n1, new TreeSet<String>()); }
			if(!neighbors.containsKey(n2)) { neighbors.put(n2, new TreeSet<String>()); }
			neighbors.get(n1).add(n2);
			ne += 1;
		}
		
		nodes = ns.toArray(new String[0]);
		edges = new String[ne][2];
		
		int i = 0;
		for(String n1 : neighbors.keySet()) { 
			for(String n2 : neighbors.get(n1)) { 
				edges[i][0] = n1; edges[i][1] = n2;
				i++;
			}
		}
	}
	
	public GraphModel(Collection es) {
		this();
		neighbors = new TreeMap<String,Set<String>>();
		Set<String> ns = new TreeSet<String>();
		int ne = 0;
		for(Object e : es) {
			String n1, n2;
			if(e instanceof Pair) { 
				Pair p = (Pair)e;
				if(!(p.getFirst() instanceof String) ||
					!(p.getLast() instanceof String)) { 
					throw new IllegalArgumentException();
				}
				n1 = (String)p.getFirst(); 
				n2 = (String)p.getLast();
			} else if (e instanceof String[]) {
				String[] earray = (String[])e;
				n1 = earray[0]; 
				n2 = earray[1];
			} else { 
				throw new IllegalArgumentException(e.toString());
			}
			
			ns.add(n1); 
			ns.add(n2);
			if(!neighbors.containsKey(n1)) { neighbors.put(n1, new TreeSet<String>()); }
			if(!neighbors.containsKey(n2)) { neighbors.put(n2, new TreeSet<String>()); }
			neighbors.get(n1).add(n2);
			ne += 1;
		}
		
		nodes = ns.toArray(new String[0]);
		edges = new String[ne][2];
		
		int i = 0;
		for(String n1 : neighbors.keySet()) { 
			for(String n2 : neighbors.get(n1)) { 
				edges[i][0] = n1; edges[i][1] = n2;
				i++;
			}
		}
	}
	
	private void rebuildNeighbors() { 
		neighbors = new TreeMap<String,Set<String>>();
		for(String n : nodes) { 
			neighbors.put(n, new TreeSet<String>());
		}
		for(String[] e : edges) { 
			if(neighbors.containsKey(e[0])) { 
				neighbors.get(e[0]).add(e[1]);
			}
		}
	}
	
	public DirectedGraph asDirectedGraph() { 
		DirectedGraph g = new DirectedGraph();
		if(neighbors==null) { rebuildNeighbors(); }
		for(String n : nodes) { g.addVertex(n); }
		for(String n1 : neighbors.keySet()) { 
			for(String n2 : neighbors.get(n1)) { 
				g.addEdge(n1, n2);
			}
		}
		return g;
	}
	
	public int numNodes() { return nodes.length; }
	public String node(int i) { return nodes[i]; }
	public int numEdges() { return edges.length; }
	public String[] edge(int i) { return edges[i]; }

	public boolean hasNode(String n) {
		if(neighbors==null) { rebuildNeighbors(); }
		return neighbors.containsKey(n);
	}

	public boolean hasEdge(String n1, String n2) {
		if(neighbors==null) { rebuildNeighbors(); }
		return neighbors.containsKey(n1) ? neighbors.get(n1).contains(n2) : false;
	}
	
	public Set<String> forward(String n) { 
		if(neighbors==null) { rebuildNeighbors(); }
		return neighbors.containsKey(n) ? neighbors.get(n) : new TreeSet<String>(); 
	}
	
	public Set<String> reverse(String n) { 
		if(neighbors==null) { rebuildNeighbors(); }
		TreeSet<String> nbs = new TreeSet<String>();
		for(String n1 : neighbors.keySet()) { 
			if(neighbors.get(n1).contains(n)) { 
				nbs.add(n1);
			}
		}
		return nbs;
	}

}

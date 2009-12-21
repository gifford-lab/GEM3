/*
 * Author: tdanford
 * Date: Oct 31, 2008
 */
package edu.mit.csail.cgs.utils.graphs;

import java.io.*;
import java.util.*;

import edu.mit.csail.cgs.utils.SetTools;

public class DirectedMultiGraph implements MultiGraph { 

	private static SetTools<String> tools;
	
	static { 
		tools = new SetTools<String>();
	}
	
	private HashSet<String> nodes;
	private Set<String> tags;
	
	// node -> outgoing/incoming -> tags
	private Map<String,Map<String,Set<String>>> edges, parents;

	public DirectedMultiGraph() { 
		nodes = new HashSet<String>();
		tags = new TreeSet<String>();
		edges = new HashMap<String,Map<String,Set<String>>>();
		parents = new HashMap<String,Map<String,Set<String>>>();
	}
	
	public DirectedMultiGraph(DirectedGraph dg, String t) { 
		tags = new TreeSet<String>();
		tags.add(t);
		nodes = new HashSet<String>(dg.getVertices());
		edges = new HashMap<String,Map<String,Set<String>>>();
		parents = new HashMap<String,Map<String,Set<String>>>();
		for(String v : nodes) { 
			edges.put(v, new TreeMap<String,Set<String>>());
			parents.put(v, new TreeMap<String,Set<String>>());
			for(String p : dg.getParents(v)) { 
				parents.get(v).put(p, new HashSet<String>());
				parents.get(v).get(p).add(t);
			}
			for(String n : dg.getNeighbors(v)) { 
				edges.get(v).put(n, new HashSet<String>());
				edges.get(v).get(n).add(t);
			}
		}
	}

	public void printGraph(PrintStream ps) { 
		TreeSet<String> nodesOrdered = new TreeSet<String>(nodes);
		for(String n : nodesOrdered) { 
			ps.print(n + " : ( ");
			for(String nn : edges.get(n).keySet()) {
				for(String tag : edges.get(n).get(nn)) { 
					ps.print(String.format("%s:%s ", nn, tag));
				}
			}
			ps.println(")");
		}
	}

	public void addVertex(String v) { 
		if(!nodes.contains(v)) { 
			nodes.add(v);
			edges.put(v, new TreeMap<String,Set<String>>());
			parents.put(v, new TreeMap<String,Set<String>>());
		}
	}

	public void addEdge(String v1, String v2, String t) { 
		if(nodes.contains(v1) && nodes.contains(v2)) {
			if(!edges.get(v1).containsKey(v2)) { 
				edges.get(v1).put(v2, new HashSet<String>());
			}
			if(!parents.get(v2).containsKey(v1)) { 
				parents.get(v2).put(v1, new HashSet<String>());
			}
			
			edges.get(v1).get(v2).add(t);
			parents.get(v2).get(v1).add(t);
			tags.add(t);
		}
	}
	
	public void removeEdge(String v1, String v2, String t) { 
		if(containsEdge(v1, v2)) {
			if(edges.get(v1).containsKey(v2)) { 
				edges.get(v1).get(v2).remove(t);
				if(parents.get(v1).get(v2).isEmpty()) { 
					parents.get(v1).remove(v2);
				}
			}
			if(parents.get(v2).containsKey(v1)) { 
				parents.get(v2).get(v1).remove(t);
				if(parents.get(v2).get(v1).isEmpty()) { 
					parents.get(v2).remove(v1);
				}
			}
			
			rebuildTags();
		}
	}
	
	public void removeVertex(String n) { 
		if(nodes.contains(n)) { 
			nodes.remove(n);
			parents.remove(n);
			edges.remove(n);

			for(String k : edges.keySet()) {
				if(parents.get(k).containsKey(n)) { 
					parents.get(k).remove(n);
				}
				if(edges.get(k).containsKey(n)) { 
					edges.get(k).remove(n);
				}
			}
			
			rebuildTags();
		}
	}
	
	private void rebuildTags() { 
		tags.clear();
		for(String v1 : nodes) { 
			for(String v2 : edges.get(v1).keySet()) { 
				tags.addAll(edges.get(v1).get(v2));
			}
		}
	}

	public boolean isNeighbor(String v, String n) { return containsEdge(v, n); }
	
	public boolean containsEdge(String v1, String v2) {
		return edges.get(v1).containsKey(v2);
	}
	
	public boolean containsEdge(String v1, String v2, String t) { 
		return tags.contains(t) && 
			edges.get(v1).containsKey(v2) && edges.get(v1).get(v2).contains(t);
	}

	public Set<String> getAncestors(String vertex) { 
		Set<String> searched = new HashSet<String>(); 
		LinkedList<String> pending = new LinkedList<String>(); 
		pending.addAll(parents.get(vertex).keySet());
		
		while(!pending.isEmpty()) { 
			String cv = pending.removeFirst();
			searched.add(cv);
			Set<String> ns = tools.subtract(parents.get(cv).keySet(), searched);
			ns.removeAll(pending);
			pending.addAll(ns);
		}
		return searched;
	}

	public Set<String> getDescendants(String vertex) { 
		Set<String> searched = new HashSet<String>(); 
		LinkedList<String> pending = new LinkedList<String>(); 
		pending.addAll(edges.get(vertex).keySet());
		
		while(!pending.isEmpty()) { 
			String cv = pending.removeFirst();
			searched.add(cv);
			Set<String> ns = tools.subtract(edges.get(cv).keySet(), searched);
			ns.removeAll(pending);
			pending.addAll(ns);
		}
		return searched;
	}
	
	public Set<String> getRoots() { 
		HashSet<String> roots = new HashSet<String>();
		for(String node : nodes) { 
			if(parents.get(node).size() == 0) { 
				roots.add(node);
			}
		}
		return roots;
	}

	public Set<String> getVertices() { 
		return new HashSet<String>(nodes);
	}

	public Set<String> getParents(String vertex) { 
		return parents.get(vertex).keySet();
	}

	public Set<String> getNeighbors(String vertex) { 
		return edges.get(vertex).keySet();
	}

	public Set<String> getNeighbors(String vertex, String tag) {
		TreeSet<String> ns = new TreeSet<String>();
		for(String v2 : edges.get(vertex).keySet()) { 
			if(edges.get(vertex).get(v2).contains(tag)) { 
				ns.add(v2);
			}
		}
		return ns;
	}

	public Set<String> getTags() {
		return tags;
	}

	public boolean isNeighbor(String v, String n, String tag) {
		return edges.containsKey(v) && edges.get(v).containsKey(n) && edges.get(v).get(n).contains(tag);
	}

	public Set<String> findNeighborTags(String v1, String v2) {
		if(edges.containsKey(v1) && edges.get(v1).containsKey(v2)) {  
			return edges.get(v1).get(v2);
		} else { 
			return new TreeSet<String>();
		}
	}
}



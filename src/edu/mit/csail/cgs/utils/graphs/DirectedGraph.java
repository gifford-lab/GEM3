package edu.mit.csail.cgs.utils.graphs;

import java.util.*;
import java.io.*;
import edu.mit.csail.cgs.utils.*;

public class DirectedGraph implements Graph { 

	private static SetTools<String> tools;
	
	static { 
		tools = new SetTools<String>();
	}
	
	private HashSet<String> nodes;
	private Map<String,Set<String>> edges, parents;
	
	public DirectedGraph() { 
		nodes = new HashSet<String>();
		edges = new HashMap<String,Set<String>>();
		parents = new HashMap<String,Set<String>>();
	}
	
	public DirectedGraph(Collection<String> ns) { 
		nodes = new HashSet<String>(ns);
		edges = new HashMap<String,Set<String>>();
		parents = new HashMap<String,Set<String>>();
	}
	
	public DirectedGraph(GraphModel m) { 
		this();
		for(int i = 0; i < m.nodes.length; i++) { 
			addVertex(m.nodes[i]);
		}
		
		for(int i = 0; i < m.edges.length; i++) { 
			addEdge(m.edges[i][0], m.edges[i][1]);
		}
	}
	
	public GraphModel asModel() { return new GraphModel(this); }

	public DirectedGraph(Collection<String> n, Map<String,Set<String>> e) { 
		nodes = new HashSet<String>(n);
		edges = new HashMap<String,Set<String>>();
		parents = new HashMap<String,Set<String>>();

		for(String node : nodes) { 
			edges.put(node, new HashSet<String>());
			parents.put(node, new HashSet<String>());
		}

		for(String v1 : e.keySet()) { 
			for(String v2 : e.get(v1)) { 
				if(nodes.contains(v1) && nodes.contains(v2)) { 
					addEdge(v1, v2);
				}
			}
		}
	}
	
	public DirectedGraph(DirectedGraph dg) { 
		nodes = new HashSet<String>(dg.nodes);
		edges = new HashMap<String,Set<String>>();
		parents = new HashMap<String,Set<String>>();
		for(String v : nodes) { 
			edges.put(v, new HashSet<String>(dg.edges.get(v)));
			parents.put(v, new HashSet<String>(dg.parents.get(v)));
		}
	}

	public DirectedGraph reverse() { 
		DirectedGraph dg = new DirectedGraph();
		dg.nodes.addAll(nodes);
		for(String n : edges.keySet()) { 
			dg.parents.put(n, new HashSet<String>(edges.get(n)));
		}
		for(String n : parents.keySet()) { 
			dg.edges.put(n, new HashSet<String>(parents.get(n)));
		}
		return dg;
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

	public void addVertex(String v) { 
		if(!nodes.contains(v)) { 
			nodes.add(v);
			edges.put(v, new HashSet<String>());
			parents.put(v, new HashSet<String>());
		}
	}

	public void addEdge(String v1, String v2) { 
		if(nodes.contains(v1) && nodes.contains(v2)) { 
			edges.get(v1).add(v2);
			parents.get(v2).add(v1);
		} else if (!nodes.contains(v1)) { 
			throw new IllegalArgumentException(v1);
		} else { 
			throw new IllegalArgumentException(v2);			
		}
	}
	
	public void removeEdge(String v1, String v2) { 
		if(containsEdge(v1, v2)) { 
			edges.get(v1).remove(v2);
			parents.get(v2).remove(v1);
		}
	}
	
	public void removeVertex(String n) { 
		if(nodes.contains(n)) { 
			nodes.remove(n);
			parents.remove(n);
			edges.remove(n);
			for(String k : edges.keySet()) { 
				if(parents.get(k).contains(n)) { parents.get(k).remove(n); }
				if(edges.get(k).contains(n)) { edges.get(k).remove(n); }
			}
		}
	}
	
	public void removeAllEdges() { 
		edges.clear();
		parents.clear();
	}
	
	public String chooseNeighbor(String n) { 
		if(edges.containsKey(n) && !edges.get(n).isEmpty()) { 
			Iterator<String> nitr = edges.get(n).iterator();
			return nitr.next();
		} else { 
			return null;
		}
	}

	public boolean isNeighbor(String v, String n) { return containsEdge(v, n); }
	
	public boolean containsEdge(String v1, String v2) { 
		return edges.get(v1).contains(v2);
	}

	public boolean containsVertex(String searchOrf) {
		return nodes.contains(searchOrf);
	}

	public Set<String> getAncestors(String vertex) { 
		Set<String> searched = new HashSet<String>(); 
		LinkedList<String> pending = new LinkedList<String>(); 
		pending.addAll(parents.get(vertex));
		
		while(!pending.isEmpty()) { 
			String cv = pending.removeFirst();
			searched.add(cv);
			Set<String> ns = tools.subtract(parents.get(cv), searched);
			ns.removeAll(pending);
			pending.addAll(ns);
		}
		return searched;
	}

	public Set<String> getDescendants(String vertex) { 
		Set<String> searched = new HashSet<String>(); 
		LinkedList<String> pending = new LinkedList<String>(); 
		pending.addAll(edges.get(vertex));
		
		while(!pending.isEmpty()) { 
			String cv = pending.removeFirst();
			searched.add(cv);
			Set<String> ns = tools.subtract(edges.get(cv), searched);
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
		return parents.get(vertex);
	}

	public Set<String> getNeighbors(String vertex) { 
		return edges.get(vertex);
	}

	public UndirectedGraph createUndirected() { 
		return new UndirectedGraph(nodes, edges);
	}

	public UndirectedGraph moralize() { 
		UndirectedGraph analog = createUndirected();
		for(String v : nodes) { 
			Vector<String> pv = new Vector<String>(parents.get(v));
			for(int i = 0; i < pv.size(); i++) { 
				String p1 = pv.get(i);
				for(int j = i + 1; j < pv.size(); j++) { 
					String p2 = pv.get(j);
					analog.addEdge(p1, p2);
				}
			}
		}
		return analog;
	}

	public Graph getSubgraph(Set<String> vs) { 
		return new DirectedGraph(vs, edges);
	}

	public int size() {
		return nodes.size();
	}
}


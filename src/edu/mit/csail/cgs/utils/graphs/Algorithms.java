package edu.mit.csail.cgs.utils.graphs;

import java.util.*;
import java.io.*;
import edu.mit.csail.cgs.utils.*;

public class Algorithms { 

	public static void main(String[] args) { 
		try { 
			Parser p = new Parser(new File(args[0]));
			DirectedGraph dg = p.parseDirectedGraph();
			Algorithms algs = new Algorithms(dg);

			String v1 = args[1], v2 = args[2];
			System.out.println("Paths (" + v1 + " to " + v2 + "):");
			LinkedList<LinkedList<String>> paths = algs.findAllPaths(v1, v2);

			for(LinkedList<String> path : paths) { 
				for(String v : path) { 
					System.out.print(v + " ");
				}
				System.out.println();
			}

		} catch(IOException ie) { 
			ie.printStackTrace(System.err);
		}
	}

	protected Graph graph;
    protected GraphSearch searcher;
	protected static SetTools<String> tools;
	
	static { 
		tools = new SetTools<String>();
	}

	public Algorithms(Graph g) { 
		graph = g;
        searcher = new GraphSearch(g);
	}
	
	public void printPath(LinkedList<String> path, PrintStream ps) { 
		boolean first = true;
		for(String n : path) { 
			ps.print((first ? "" : " --> ") + n);
		}
	}
	
	public boolean isPath(LinkedList<String> path) {
		String prev = null;
		for(String node : path) { 
			if(!graph.getVertices().contains(node)) { return false; }
			if(prev != null) { 
				if(!graph.getNeighbors(prev).contains(node)) { 
					return false;
				}
			}
			prev = node;
		}
		return true;
	}

	public LinkedList<LinkedList<String>> findAllPaths(String v1, String v2) { 
		LinkedList<LinkedList<String>> allPaths = new LinkedList<LinkedList<String>>();
		searchAllPaths(v1, v2, new HashSet<String>(), allPaths);
		return allPaths;
	}

	private void searchAllPaths(String v1, String v2, Set<String> searched, 
			LinkedList<LinkedList<String>> allPaths) { 

		searched.add(v1);
		if(v1.equals(v2)) { 
			LinkedList<String> path = new LinkedList<String>();
			path.addLast(v1);
			allPaths.addLast(path);
			return;
		}

		Set<String> next = tools.subtract(graph.getNeighbors(v1), searched);
		for(String n : next) { 
			LinkedList<LinkedList<String>> nextPaths = new LinkedList<LinkedList<String>>();
			searchAllPaths(n, v2, new HashSet<String>(), nextPaths);
			for(LinkedList<String> path : nextPaths) { path.addFirst(v1); }
			allPaths.addAll(nextPaths);
		}
	}

	public LinkedList<String> findPath(String v1, String v2) { 
		LinkedList<String> path = new LinkedList<String>();
		return searchPath(v1, v2, new HashSet<String>(), path) ? path : null;
	}

	private boolean searchPath(String v1, String v2, Set<String> searched, LinkedList<String> path) { 
		if(v1.equals(v2)) { path.addLast(v1); return true; }
		searched.add(v1);
		Set<String> ns = tools.subtract(graph.getNeighbors(v1), searched);
		for(String n : ns) { 
			if(searchPath(n, v2, searched, path)) { 
				path.addFirst(v1);
				return true;
			}
		}
		return false;
	}
	
	public static class ConnectedSetSearcher implements SearchInterface { 
		private Set<String> total;
		public ConnectedSetSearcher() { total = new HashSet<String>(); }
		public Set<String> getConnectedComponent() { return total; }
		
		public boolean searchNode(Graph g, String n) {
			total.add(n);
			return true;
		}
	}
    
    public static class ConnectedSearcher implements SearchInterface {
        private String target;
        private boolean found;
        public ConnectedSearcher(String t) { target = t; found = false; }
        
        public boolean isFound() { return found; }
        
        public boolean searchNode(Graph g, String node) {
            found = node.equals(target);
            return !found;
        } 
    }
    
    public Set<String> getConnected(Set<String> vs) {
		ConnectedSetSearcher css = new ConnectedSetSearcher();
    	for(String v : vs) { 
    		searcher.BFS(v, css);
    	}
    	return css.getConnectedComponent();
    }

	public boolean isConnected(String v1, String v2) {
        ConnectedSearcher cs = new ConnectedSearcher(v2);
        searcher.DFS(v1, cs);
        return cs.isFound();
	}
}

/*
 * Created on Oct 6, 2006
 */
package edu.mit.csail.cgs.utils.graphs;

import java.util.*;

import edu.mit.csail.cgs.utils.SetTools;

/**
 * @author tdanford
 */
public class GraphSearch {
    
    private static SetTools<String> tools;
    
    static { 
        tools = new SetTools<String>();
    }
    
    private Graph graph;

    public GraphSearch(Graph g) {
        graph = g;
    }

    public void ConnectedBFS(String start, SearchInterface searcher) { 
        HashSet<String> seen = new HashSet<String>();
        LinkedList<String> pending = new LinkedList<String>();
        if(!graph.getVertices().contains(start)) { throw new IllegalArgumentException(); }

		pending.addAll(graph.getNeighbors(start));
        boolean keepSearching = true;

        while(keepSearching && !pending.isEmpty()) { 
            String current = pending.removeFirst();
            seen.add(current);
            keepSearching = searcher.searchNode(graph, current);

            if(keepSearching) { 
                Set<String> neighbors = tools.subtract(graph.getNeighbors(current), seen);
                neighbors.removeAll(pending);
				pending.addAll(neighbors);
            }
        }
    }

    public void BFS(String start, SearchInterface searcher) { 
        HashSet<String> seen = new HashSet<String>();
        LinkedList<String> pending = new LinkedList<String>();
        pending.addLast(start);
        
        if(!graph.getVertices().contains(start)) { throw new IllegalArgumentException(); }
        boolean keepSearching = true;
        while(keepSearching && !pending.isEmpty()) { 
            String current = pending.removeFirst();
            seen.add(current);
            keepSearching = searcher.searchNode(graph, current);

            if(keepSearching) { 
                Set<String> neighbors = tools.subtract(graph.getNeighbors(current), seen);
                neighbors.removeAll(pending);
                for(String n : neighbors) { 
                    pending.addLast(n);
                }
            }
        }
    }

    public void DFS(String start, SearchInterface searcher) { 
        HashSet<String> seen = new HashSet<String>();
        LinkedList<String> pending = new LinkedList<String>();
        pending.addLast(start);
        
        if(!graph.getVertices().contains(start)) { throw new IllegalArgumentException(); }
        boolean keepSearching = true;
        while(keepSearching && !pending.isEmpty()) { 
            String current = pending.removeFirst();
            seen.add(current);
            keepSearching = searcher.searchNode(graph, current);

            if(keepSearching) { 
                Set<String> neighbors = tools.subtract(graph.getNeighbors(current), seen);
                neighbors.removeAll(pending);
                for(String n : neighbors) { 
                    pending.addFirst(n);
                }
            }
        }
    }
}


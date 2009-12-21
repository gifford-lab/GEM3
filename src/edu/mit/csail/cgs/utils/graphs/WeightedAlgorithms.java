/*
 * Author: tdanford
 * Date: Nov 15, 2008
 */
package edu.mit.csail.cgs.utils.graphs;

import java.io.*;
import java.util.*;

import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.datastructures.Heap;

public class WeightedAlgorithms extends Algorithms {
	
	public static void main(String[] args) { 
		//dijkstra_test(args);
		//bellman_ford_test(args);
		floyd_warshall_test(args);
	}
	
	public static void floyd_warshall_test(String[] args) { 
		WeightedGraph wg = figure_26_1();
		
		WeightedAlgorithms algs = new WeightedAlgorithms(wg);
		WeightedGraph[] paths = algs.floydWarshall();
		
		String[] vertices = algs.vertices;
		int source = 1, target = 0;
		
		Pair<List<String>,Double> pathpair = algs.parsePathTree(paths[source], vertices[target]);
		System.out.println(String.format("%s -> %s", vertices[source], vertices[target]));
		System.out.println(String.format("Weight: %f", pathpair.getLast()));
		System.out.println(String.format("Path: %s", pathpair.getFirst()));
	}
	
	/**
	 * Tests the Bellman-Ford implementation using the example from CLRS page 589 (Figure 24.4).  
	 * 
	 * @param args
	 */
	public static void bellman_ford_test(String[] args) { 
		WeightedGraph wg = figure_24_4();
		String source = "s", target = "z";
		
		WeightedAlgorithms algs = new WeightedAlgorithms(wg);
		WeightedGraph paths = algs.bellmanFord(source);
		
		Pair<List<String>,Double> pathpair = algs.parsePathTree(paths, target);
		System.out.println(String.format("%s -> %s", source, target));
		System.out.println(String.format("Weight: %f", pathpair.getLast()));
		System.out.println(String.format("Path: %s", pathpair.getFirst()));
	}
	
	/**
	 * Runs, as a test of the dijkstra implementation, the example from CLRS on page 596 of 
	 * the 2nd edition.
	 * 
	 * @param args
	 */
	public static void dijkstra_test(String[] args) { 
		DirectedWeightedGraph wg = figure_25_5();
		String source = "s", target = "x";
		
		WeightedAlgorithms algs = new WeightedAlgorithms(wg);
		WeightedGraph paths = algs.dijkstra(source);
		
		Pair<List<String>,Double> pathpair = algs.parsePathTree(paths, target);
		System.out.println(String.format("%s -> %s", source, target));
		System.out.println(String.format("Weight: %f", pathpair.getLast()));
		System.out.println(String.format("Path: %s", pathpair.getFirst()));
	}
	
	private WeightedGraph graph;
	private String[] vertices;

	public WeightedAlgorithms(WeightedGraph g) { 
		super(g);
		graph = g;
		vertices = graph.getVertices().toArray(new String[0]);
		Arrays.sort(vertices);
	}
	
	/**
	 * A utility method, for manipulating weighted path-tree objects (which are 
	 * represented as WeightedGraphs in this implementation).  
	 * 
	 * For a given tree and a given "target" node, this method will follow links in the 
	 * path tree back to the root -- this gives the shortest path (in reverse order, of course)
	 * of the path from root --> target in the original graph.
	 * 
	 * This method also returns the ultimate weight of that (shortest) path.  
	 * 
	 * @param g The path-tree, which must have been calculated with respect to a given source.  
	 * @param target The ultimate target in the graph.  
	 * @return
	 */
	public Pair<List<String>, Double> parsePathTree(WeightedGraph g, String target) { 
		LinkedList<String> path = new LinkedList<String>();
		Double weight = g.weight(target);
		
		while(target != null) { 
			path.addFirst(target);
			Set<String> nbs = g.getNeighbors(target);
			Iterator<String> itr = nbs.iterator();
			target = itr.hasNext() ? itr.next() : null;
		}
		
		return new Pair<List<String>,Double>(path, weight);
	}
	
	/**
	 * Performs the Floyd-Warshall all-pairs shortest-paths algorithm.  Returns an array of path 
	 * tree matrices.  
	 * @return A path-tree matrix, where matrix[i] is the path-tree for paths originating from vertices[i].
	 */
	public WeightedGraph[] floydWarshall() { 

		DirectedWeightedGraph[][] karray = new DirectedWeightedGraph[][] { 
				createFloydWarshallMatrix(), 
				createFloydWarshallMatrix()
		};
		
		System.out.println("Init:");
		printFloydWarshallMatrix(karray[0], System.out);
		
		int thisk = 0, lastk = 0;
		for(int k = 1; k <= vertices.length; k++) { 
			thisk = k % 2;
			lastk = 1-thisk;
			
			floydWarshallStep(k-1, karray[thisk], karray[lastk]);
			
			System.out.println(String.format("%d:", k));
			printFloydWarshallMatrix(karray[thisk], System.out);
		}

		WeightedGraph[] finalMatrix = new WeightedGraph[karray[0].length];
		for(int i = 0; i < karray[0].length; i++) { 
			finalMatrix[i] = karray[thisk][i];
		}
		
		return finalMatrix;
	}
	
	public void printFloydWarshallMatrix(DirectedWeightedGraph[] karray, PrintStream ps) { 
		for(int j = 0; j < vertices.length; j++) { 
			ps.print(String.format("\t%s", vertices[j]));
		}
		ps.println();
		
		for(int i = 0; i < karray.length; i++) { 
			ps.print(vertices[i]);
			for(int j = 0; j < vertices.length; j++) { 
				Double w = karray[i].weight(vertices[j]);
				String wstr = w != null ? String.format("%.1f", w) : "*";
				ps.print(String.format("\t%s", wstr));
			}
			
			ps.println();
		}
	}
	
	private void floydWarshallStep(int k, DirectedWeightedGraph[] thisK, DirectedWeightedGraph[] lastK) {
		for(int i = 0; i < vertices.length; i++) { 
			for(int j = 0; j < vertices.length; j++) { 

				String lastPath = lastK[i].chooseNeighbor(vertices[j]);
				Double lastWeight = lastPath != null ? graph.weight(lastPath, vertices[j]) : null;
				
				thisK[i].removeEdges(vertices[j]);

				Double dij = lastK[i].weight(vertices[j]);

				Double dik = lastK[i].weight(vertices[k]);
				Double dkj = lastK[k].weight(vertices[j]);
				Double ikjWeight = add(dik, dkj);
	
				Boolean isBetter = lessThan(ikjWeight, dij);

				if(isBetter) { 
					thisK[i].addWeightedEdge(vertices[j], vertices[k], graph.weight(vertices[k], vertices[j]));
					thisK[i].setWeight(vertices[j], ikjWeight);
				} else if (lastPath != null) { 
					thisK[i].addWeightedEdge(vertices[j], lastPath, lastWeight);
					thisK[i].setWeight(vertices[j], dij);
				}
			}
		}
	}
	
	private Double add(Double v1, Double v2) { 
		return v1 == null || v2 == null ? null : v1 + v2;
	}
	
	private Boolean lessThan(Double v1, Double v2) { 
		if(v1 == null) { 
			return false;
		} else if (v2 == null) { 
			return true;
		} else { 
			return v1 < v2;
		}
	}
	
	private DirectedWeightedGraph[] createFloydWarshallMatrix() { 
		DirectedWeightedGraph[] matrix = new DirectedWeightedGraph[vertices.length];
		for(int i = 0; i < vertices.length; i++) { 
			matrix[i] = new DirectedWeightedGraph(null);
			for(int j = 0; j < vertices.length; j++) { 
				matrix[i].addVertex(vertices[j]);
			} 
			
			for(int j = 0; j < vertices.length; j++) { 
				if(i == j) { 
					matrix[i].setWeight(vertices[j], 0.0);
				} else if (graph.isNeighbor(vertices[i], vertices[j])) { 
					Double w = graph.weight(vertices[i], vertices[j]);

					matrix[i].setWeight(vertices[j], w);
					matrix[i].addWeightedEdge(vertices[j], vertices[i], w);
				}
			}
		}
		return matrix;
	}
	
	public WeightedGraph bellmanFord(String source) { 

		if(!graph.getVertices().contains(source)) { 
			throw new IllegalArgumentException(source);
		}
		DirectedWeightedGraph pathtree = initializePathTree();
		pathtree.setWeight(source, 0.0);
		boolean hasNegativeLoop = false;
		
		int V = graph.getVertices().size();

		for(int i = 0; i < V-1; i++) { 
			for(String n1 : graph.getVertices()) { 
				for(String n2 : graph.getNeighbors(n1)) { 
					double w = graph.weight(n1, n2);
					relax(pathtree, n1, n2, w);
				}
			}
		}

		outerLoop: for(String u : graph.getVertices()) { 
			double du = pathtree.weight(u);
			for(String v : graph.getNeighbors(u)) { 
				double dv = pathtree.weight(v);
				if(dv > du + graph.weight(u, v)) { 
					hasNegativeLoop = true; 
					break outerLoop;  // ends the outer-most for loop.
				}
			}
		}

		return hasNegativeLoop ? null : pathtree;
	}
	
	/**
	 * Returns a weighted graph of shortest-path links.  
	 * 
	 * For a given v2 \in G = dijkstra(v1), the 
	 * weight(v2) in G is the weight of the shortest path v1 -> v2. 
	 * 
	 * Furthermore, v2 has exactly one neighbor in G, v3, which is the previous
	 * node on that shortest path from v2 back to v2.  The weight(v2, v3) is the 
	 * weight of the edge v3->v2 in the *original* graph which is on that shortest path.
	 * 
	 * dijkstra() throws an IllegalArgumentException if it encounters a negative-weight edge.
	 * 
	 * @param source The source node in the original graph.
	 * @return The weighted path-tree.  
	 */
	public WeightedGraph dijkstra(String source) {
		
		if(!graph.getVertices().contains(source)) { 
			throw new IllegalArgumentException(source);
		}
		
		DirectedWeightedGraph pathtree = initializePathTree();
		pathtree.setWeight(source, 0.0);
				
		Heap<WeightedNode> Q = new Heap<WeightedNode>(-1);
		for(String vertex : pathtree.getVertices()) { 
			Q.insert(new WeightedNode(vertex, pathtree.weight(vertex)));
		}
		
		while(Q.size() > 0) { 
			//System.out.println(String.format("\nQ: %s", Q.asList().toString()));
			WeightedNode u = Q.removeFirst();
			//System.out.println(String.format("u: %s", u.toString()));
			if(u.weight == null) {
				// The rest of the graph is unconnected -- so we just return the pathtree.
				return pathtree;
			} else { 
				for(String v : graph.getNeighbors(u.node)) { 
					Double w = graph.weight(u.node, v);
					if(w < 0.0) { 
						// Dijkstra's algorithm doesn't work, when there are negative weight edges.
						throw new IllegalArgumentException(String.format("%s->%s has negative weight %f",
								u.node, v, w));
					}
					relax(Q, pathtree, u, v, w);
				}
			}
		}
		
		return pathtree;
	}
	
	/**
	 * This is totally ugly.  
	 * 
	 * The RELAX() method described in CLRS is actually tricky -- 
	 * it requires an (implicit) update of the Q heap, and that requires knowing where each 
	 * old key was in the heap so we can call increaseKey() on it ... so I do that here with a
	 * linear search through the heap (implemented in the increase() method, in Heap.java), but 
	 * this is probably sub-optimal... sigh.    
	 * 
	 * @param Q See the note above.
	 * @param pathtree  Stores the 'd' and 'pi' variables from CLRS, implicitly.
	 * @param u 
	 * @param v
	 * @param w
	 */
	private void relax(Heap<WeightedNode> Q, DirectedWeightedGraph pathtree, 
			WeightedNode u, String v, Double w) { 
		
		Double uweight = pathtree.weight(u.node);
		Double oldvweight = pathtree.weight(v);
		Double newvweight = uweight != null ? uweight + w : null;

		WeightedNode oldvn = new WeightedNode(v, oldvweight);
		WeightedNode newvn = new WeightedNode(v, newvweight);
		
		if(newvn.compareTo(oldvn) == -1) {
			for(String nn : pathtree.getNeighbors(v)) { 
				pathtree.removeEdge(v, nn);
			}
			
			pathtree.addEdge(v, u.node);
			pathtree.setWeight(v, newvweight);
			pathtree.setWeight(v, u.node, w);
			
			//System.out.println(String.format("%s => %s", oldvn.toString(), newvn.toString()));
			Q.increase(oldvn, newvn);
		}
	}
	
	private void relax(DirectedWeightedGraph pathtree, String u, String v, Double w) { 
		Double uweight = pathtree.weight(u);
		Double oldvweight = pathtree.weight(v);
		Double newvweight = uweight != null ? uweight + w : null;
		
		if(newvweight == null) { return; }
		
		if(oldvweight == null || newvweight < oldvweight) { 
			for(String nn : pathtree.getNeighbors(v)) { 
				pathtree.removeEdge(v, nn);
			}
			
			pathtree.addEdge(v, u);
			pathtree.setWeight(v, newvweight);
			pathtree.setWeight(v, u, w);
		}		
	}
	
	private DirectedWeightedGraph initializePathTree() { 
		DirectedWeightedGraph tree = new DirectedWeightedGraph(null);
		for(String node : graph.getVertices()) { 
			tree.addVertex(node);
		}
		return tree;
	}
	
	private class WeightedNode implements Comparable<WeightedNode> {
		
		public String node; 
		public Double weight;  // a null weight will represent 'infinity.'
		
		public WeightedNode(String n, Double w) { 
			node = n; weight = w;
		}
		
		public WeightedNode(String n) { 
			node = n;
			weight = null;
		}
		
		public String toString() { 
			if(weight != null) { 
				return String.format("%s (%.2f)", node, weight);
			} else { 
				return String.format("%s (inf)", node);
			}
		}
		
		public int compareTo(WeightedNode n) {
			if(weight != null || n.weight != null) { 
				if(weight != null && n.weight != null) { 
					if(weight < n.weight) { return -1; }
					if(weight > n.weight) { return 1; }
				} else { 
					if(weight == null) { return 1; }
					if(n.weight == null) { return -1; }
				}
			}
			//return node.compareTo(n.node);
			return 0;
		}
		
		public int hashCode() { 
			return node.hashCode();
		}
		
		public boolean equals(Object o) { 
			if(!(o instanceof WeightedNode)) { return false; }
			WeightedNode w = (WeightedNode)o;
			return w.node.equals(node);
		}
	}

	/**
	 * The example from CLRS, used to illustrate Dijkstra's Algorithm.
	 * @return
	 */
	public static DirectedWeightedGraph figure_25_5() { 
		DirectedWeightedGraph wg = new DirectedWeightedGraph(1.0);
		wg.addVertex("s");
		wg.addVertex("t");
		wg.addVertex("x");
		wg.addVertex("y");
		wg.addVertex("z");
		wg.addWeightedEdge("s", "t", 10.0);
		wg.addWeightedEdge("s", "y", 5.0);
		wg.addWeightedEdge("t", "y", 2.0);
		wg.addWeightedEdge("y", "t", 3.0);
		wg.addWeightedEdge("y", "x", 9.0);
		wg.addWeightedEdge("t", "x", 1.0);
		wg.addWeightedEdge("y", "z", 2.0);
		wg.addWeightedEdge("x", "z", 4.0);
		wg.addWeightedEdge("z", "x", 6.0);
		wg.addWeightedEdge("z", "s", 7.0);
		
		return wg;
	}
	
	public static WeightedGraph figure_26_1() { 
		DirectedWeightedGraph wg = new DirectedWeightedGraph(null);
		
		wg.addVertex("1");
		wg.addVertex("2");
		wg.addVertex("3");
		wg.addVertex("4");
		wg.addVertex("5");
		
		wg.addWeightedEdge("1", "5", -4.0);
		wg.addWeightedEdge("1", "2", 3.0);
		wg.addWeightedEdge("1", "3", 8.0);
		wg.addWeightedEdge("2", "5",  7.0);
		wg.addWeightedEdge("2", "4", 1.0);
		wg.addWeightedEdge("3", "2", 4.0);
		wg.addWeightedEdge("4", "3", -5.0);
		wg.addWeightedEdge("4", "1", 2.0);
		wg.addWeightedEdge("5", "4", 6.0);
		
		return wg;
	}
	
	public static WeightedGraph figure_24_4() { 
		DirectedWeightedGraph wg = new DirectedWeightedGraph(null);
		wg.addVertex("s");
		wg.addVertex("t");
		wg.addVertex("x");
		wg.addVertex("y");
		wg.addVertex("z");
		wg.addWeightedEdge("s", "t", 6.0);
		wg.addWeightedEdge("s", "y", 7.0);
		wg.addWeightedEdge("t", "y", 8.0);
		wg.addWeightedEdge("t", "x", 5.0);
		wg.addWeightedEdge("t", "z", -4.0);
		wg.addWeightedEdge("y", "x", -3.0);
		wg.addWeightedEdge("y", "z", 9.0);
		wg.addWeightedEdge("x", "t", -2.0);
		wg.addWeightedEdge("z", "s", 2.0);
		wg.addWeightedEdge("z", "x", 7.0);

		return wg;
	}	
}

/*
 * Author: tdanford
 * Date: Nov 9, 2008
 */
package edu.mit.csail.cgs.utils.graphs;

import java.util.*;
import Jama.*;

public class AdjacencyMatrix {

	private Matrix adj;
	private String[] vertices;
	
	public AdjacencyMatrix(Graph g) { 
		Set<String> verts = g.getVertices();
		vertices = verts.toArray(new String[0]);
		Arrays.sort(vertices);
		
		adj = new Matrix(vertices.length, vertices.length, 0.0);
		
		for(int i = 0; i < vertices.length; i++) { 
			for(int j = 0; j < vertices.length; j++) { 
				if(i == j) { 
					adj.set(i, j, 1.0);
				} else { 
					if(g.getNeighbors(vertices[i]).contains(vertices[j])) { 
						adj.set(i, j, 1.0);
					}
				}
			}
		}
	}
	
	public void path() { 
		adj = adj.times(adj);
	}
	
	public void path(int n) { 
		for(int i = 0; i < n; i++) { 
			path();
		}
	}
	
	public int size() { return vertices.length; } 
	public String vertex(int i) { return vertices[i]; }
	
	public boolean connected(int i, int j) { 
		return adj.get(i, j) >= 0.5;
	}
	
	public int value(int i, int j) { 
		return (int)Math.round(adj.get(i, j));
	}
	
	public Matrix matrix() { 
		return adj;
	}
}

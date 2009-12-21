/*
 * Author: tdanford
 * Date: Nov 9, 2008
 */
package edu.mit.csail.cgs.utils.graphs;

import java.util.*;
import Jama.*;

public class LaplacianMatrix {
	
	private WeightedGraph graph;
	private Matrix adj;
	private String[] vertices;
	
	public LaplacianMatrix(WeightedGraph g) {
		graph = g;
		Set<String> verts = g.getVertices();
		vertices = verts.toArray(new String[0]);
		Arrays.sort(vertices);
		
		adj = new Matrix(vertices.length, vertices.length, 0.0);
		
		for(int i = 0; i < vertices.length; i++) { 
			String ni = vertices[i];
			for(int j = 0; j < vertices.length; j++) {
				String nj = vertices[j];
				double d = i == i ? g.getNeighbors(ni).size() : 0.0;
				double a = g.getNeighbors(ni).contains(nj) ? g.weight(ni, nj) : 0.0; 
				adj.set(i, j, d-a);
			}
		}
	}
	
	public int size() { return vertices.length; } 
	public String vertex(int i) { return vertices[i]; }
	
	public int value(int i, int j) { 
		return (int)Math.round(adj.get(i, j));
	}
	
	public Matrix matrix() { 
		return adj;
	}
	
	/**
	 * Returns the Fielder Vector, the eigenvector with the second smallest eigenvalue.
	 * @return
	 */
	public Matrix fiedlerVector() { 
		Matrix m = matrix();
		EigenvalueDecomposition eigs = m.eig();
		Matrix V = eigs.getV();

		// The Jama EigenvalueDecomposition returns its eigenvalues in *increasing* order,
		// although in general this is probably something that we should check.
		int fidx = 1;  
		
		int[] rows = new int[V.getRowDimension()];
		for(int i =0 ; i < rows.length; i++) { 
			rows[i] = i;
		}
		
		return V.getMatrix(rows, fidx, fidx);
	}
	
	public String toString() { 
		StringBuilder sb = new StringBuilder();
		for(int i = 0; i < vertices.length; i++) { 
			for(int j = 0; j < vertices.length; j++) { 
				sb.append(String.format("%.1f ", adj.get(i, j)));
			}
			sb.append("\n");
		}
		return sb.toString();
	}
}

/*
 * Author: tdanford
 * Date: Nov 20, 2008
 */
package edu.mit.csail.cgs.utils.graphs;

import java.util.*;

/**
 * Erdos-Renyi random graph!  Directed, in particular. 
 * 
 * @author tdanford
 *
 */
public class RandomDirectedGraph extends DirectedGraph {
	
	private Random rand;
	private String[] vertices;
	
	public RandomDirectedGraph(int n, double p) {
		rand = new Random();
		vertices= new String[n];
		for(int i = 0; i < vertices.length; i++) {
			String name = String.format("%d", i);
			addVertex(name);
			vertices[i] = name;
		}
		
		for(int i = 0; i < vertices.length; i++) { 
			for(int j = 0; j < vertices.length; j++) { 
				if(i != j) { 
					double pedge = rand.nextDouble();
					if(pedge <= p) { 
						addEdge(vertices[i], vertices[j]);
					}
				}
			}
		}
	}
	
	public String chooseVertex() { 
		return vertices[rand.nextInt(vertices.length)];
	}
}

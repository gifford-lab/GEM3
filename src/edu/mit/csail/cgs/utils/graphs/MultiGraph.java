/*
 * Author: tdanford
 * Date: Oct 31, 2008
 */
package edu.mit.csail.cgs.utils.graphs;

import java.util.Set;

public interface MultiGraph extends Graph {
	
	/*
	 * From Graph:
	public Set<String> getVertices();
	public Set<String> getNeighbors(String vertex);
	public boolean isNeighbor(String v, String n);
	*/

	// "Tags" are the labels we use on the multi-edges.  
	public Set<String> getTags();
	public Set<String> getNeighbors(String vertex, String tag);
	public boolean isNeighbor(String v, String n, String tag);
	public Set<String> findNeighborTags(String v1, String v2);
}


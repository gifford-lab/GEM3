package edu.mit.csail.cgs.clustering;

import java.util.Collection;

/**
 * ClusteringMethod is a means of turning a collection of objects
 * into a collection of Cluster objects.  Both the kmeans and 
 * hierarchical sub-packages provide implementations of this interface.
 * 
 * @author Timothy Danford
 */
public interface ClusteringMethod<X> {
	public Collection<Cluster<X>> clusterElements(Collection<X> elmts);
}

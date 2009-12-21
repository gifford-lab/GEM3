package edu.mit.csail.cgs.clustering;

import java.util.Set;

/**
 * Cluster is one of the key interfaces for the clustering package.  
 * An implementation of Cluster is a collection of ClusterElements, 
 * and is returned by a ClusteringMethod.
 * 
 * @author Timothy Danford
 */
public interface Cluster <X extends Object> {
    public Set<X> getElements();
    public int size();	
}

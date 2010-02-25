package edu.mit.csail.cgs.clustering;

import java.util.Set;

//-2rT4v  -DDuplex:true  --header='$n|$%/$=|%W' --columns=1  --portrait  --line-numbers  --font=Courier9  --highlight=java  
//--color  --word-wrap  --mark-wrapped-lines=arrow  --verbose  -dxerox5 ${container_loc}/${resource_name}

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

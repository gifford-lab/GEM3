package edu.mit.csail.cgs.clustering;

/**
 * ClusterRepresentative is an interface which is implemented by any class
 * which is used to turn a Cluster into a "representative" ClusterElement 
 * used for further clustering -- this doesn't necessarily need to be an 
 * element of the Cluster itself.
 * 
 * @author Timothy Danford
 *
 */
public interface ClusterRepresentative<X> {	
    public X getRepresentative(Cluster<X> c);       
}

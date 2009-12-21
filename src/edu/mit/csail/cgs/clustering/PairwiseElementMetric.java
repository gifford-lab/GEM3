package edu.mit.csail.cgs.clustering;

/**
 * PairwiseElementMetric is the interface whose implementations are 
 * metrics on pairs of ClusterElement -- implementations of this are 
 * used in several ClusteringMethod implementations.
 * 
 * @author Timothy Danford
 */
public interface PairwiseElementMetric<X> {
	
    public double evaluate(X e1, X e2);
	
}

package edu.mit.csail.cgs.clustering.vectorcluster;

import edu.mit.csail.cgs.clustering.PairwiseElementMetric;

public class ChebyshevDistance<X extends VectorClusterElement> implements
		PairwiseElementMetric<X> {

	public double evaluate(X e1, X e2) {
		if(e1.dimension() != e2.dimension()) { throw new IllegalArgumentException(); }
		double maxvalue = 0;
		for(int i = 0; i < e1.dimension(); i++) {
			if(!e1.isMissingValue(i) && !e2.isMissingValue(i)) {
				double s = Math.abs(e1.getValue(i) - e2.getValue(i));
				if (s>maxvalue) {
					maxvalue = s;
				}
			}
		}
		return maxvalue;
	}

}

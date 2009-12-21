package edu.mit.csail.cgs.tools.motifs;

import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.clustering.PairwiseElementMetric;

/**
 * compares two weight matrices and returns the euclidian distance between them 
 * normalized over the length of the query matrix
 */
public class WMNormalizedDistanceComparator implements WMComparator, PairwiseElementMetric<WeightMatrix> {
    public double compare(WeightMatrix query, WeightMatrix target) {        
        if (target.matrix.length < query.matrix.length) {
            return Double.MAX_VALUE;
        }
        double minscore = Double.MAX_VALUE;
        for (int offset = 0; offset <= target.matrix.length - query.matrix.length; offset++) {
            double score = 0;
            for (int i = 0; i < query.matrix.length; i++) {
                score += Math.sqrt(Math.pow(query.matrix[i]['A'] - target.matrix[offset+i]['A'],2) +
                                   Math.pow(query.matrix[i]['C'] - target.matrix[offset+i]['C'],2) +
                                   Math.pow(query.matrix[i]['G'] - target.matrix[offset+i]['G'],2) +
                                   Math.pow(query.matrix[i]['T'] - target.matrix[offset+i]['T'],2));
            }
            if (score < minscore) {
                minscore = score;
            }
            score = 0;
            for (int i = 0; i < query.matrix.length; i++) {
                score += Math.sqrt(Math.pow(query.matrix[query.matrix.length - i -1]['T'] - target.matrix[offset+i]['A'],2) +
                                   Math.pow(query.matrix[query.matrix.length - i -1]['G'] - target.matrix[offset+i]['C'],2) +
                                   Math.pow(query.matrix[query.matrix.length - i -1]['C'] - target.matrix[offset+i]['G'],2) +
                                   Math.pow(query.matrix[query.matrix.length - i -1]['A'] - target.matrix[offset+i]['T'],2));
            }
            if (score < minscore) {
                minscore = score;
            }
        }
        return minscore / query.length();
    }

    public double evaluate(WeightMatrix query, WeightMatrix target) {
        if (query.length() > target.length()) {
            return compare(target,query);
        } else {
            return compare(query,target);
        }
    }

}

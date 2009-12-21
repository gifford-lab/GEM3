package edu.mit.csail.cgs.tools.motifs;

import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.clustering.PairwiseElementMetric;

/**
 * compares two weight matrices and returns the euclidian distance between them 
 */
public class WMDistanceComparator implements WMComparator, PairwiseElementMetric<WeightMatrix> {
    private boolean normalize;
    private int compareLength;
    public WMDistanceComparator(boolean n, int c) {
        normalize = n;
        compareLength = -1;
    }

    public double compare(WeightMatrix query, WeightMatrix target) {        
        double minscore = Double.MAX_VALUE;
        int complen = compareLength;
        if (compareLength < 0) {
            complen = query.matrix.length;
        }
        if (target.matrix.length < complen) {
            return Double.MAX_VALUE;
        }
        for (int queryoffset = 0; queryoffset <= query.matrix.length - complen ; queryoffset++) {
            for (int targetoffset = 0; targetoffset <= target.matrix.length - complen; targetoffset++) {
                double score = 0;
                for (int i = 0; i < complen; i++) {
                    score += Math.sqrt(Math.pow(query.matrix[queryoffset + i]['A'] - target.matrix[targetoffset+i]['A'],2) +
                                       Math.pow(query.matrix[queryoffset + i]['C'] - target.matrix[targetoffset+i]['C'],2) +
                                       Math.pow(query.matrix[queryoffset + i]['G'] - target.matrix[targetoffset+i]['G'],2) +
                                       Math.pow(query.matrix[queryoffset + i]['T'] - target.matrix[targetoffset+i]['T'],2));
                }
                if (score < minscore) {
                    minscore = score;
                }
                score = 0;
                for (int i = 0; i < complen; i++) {
                    score += Math.sqrt(Math.pow(query.matrix[queryoffset + complen - i -1]['T'] - target.matrix[targetoffset+i]['A'],2) +
                                       Math.pow(query.matrix[queryoffset + complen - i -1]['G'] - target.matrix[targetoffset+i]['C'],2) +
                                       Math.pow(query.matrix[queryoffset + complen - i -1]['C'] - target.matrix[targetoffset+i]['G'],2) +
                                       Math.pow(query.matrix[queryoffset + complen - i -1]['A'] - target.matrix[targetoffset+i]['T'],2));
                }
                if (score < minscore) {
                    minscore = score;
                }
            
            }
        }
        if (normalize) {
            return minscore / complen;
        } else {
            return minscore;
        }
    }

    public double evaluate(WeightMatrix query, WeightMatrix target) {
        if (query.length() > target.length()) {
            return compare(target,query);
        } else {
            return compare(query,target);
        }
    }

}

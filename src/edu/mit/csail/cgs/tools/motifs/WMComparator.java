package edu.mit.csail.cgs.tools.motifs;

import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.clustering.PairwiseElementMetric;

/** compares two weight matrices and returns a score to describe how
 * close they are.  Zero is the minumum score and represents a perfect match.
 *   Worse matches should receive larger values 
 */
public interface WMComparator extends PairwiseElementMetric<WeightMatrix> {

    public double compare(WeightMatrix query, WeightMatrix target);
}

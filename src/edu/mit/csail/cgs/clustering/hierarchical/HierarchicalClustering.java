package edu.mit.csail.cgs.clustering.hierarchical;

import java.util.Collection;
import java.util.Vector;

import edu.mit.csail.cgs.clustering.Cluster;
import edu.mit.csail.cgs.clustering.ClusterRepresentative;
import edu.mit.csail.cgs.clustering.ClusteringMethod;
import edu.mit.csail.cgs.clustering.PairwiseElementMetric;
import edu.mit.csail.cgs.clustering.SingletonCluster;

/**
 * @author Timothy Danford
 *
 */
public class HierarchicalClustering<X> implements ClusteringMethod<X> {
    
    private ClusterRepresentative<X> repr;
    private PairwiseElementMetric<X> metric;
    private double maxDistanceToAccept;
	
    public HierarchicalClustering(ClusterRepresentative<X> rep, PairwiseElementMetric<X> m) { 
        repr = rep;
        metric = m;
        maxDistanceToAccept = Double.MAX_VALUE;
    }

    public void setMaxDistanceToAccept(double d) {
        maxDistanceToAccept = d;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.clustering.ClusteringMethod#clusterElements(java.util.Collection)
     */
    public Collection<Cluster<X>> clusterElements(Collection<X> elmts) {
        Vector<Cluster<X>> clusters = new Vector<Cluster<X>>();
        Vector<X> reps = new Vector<X>();
        Double distances[][] = new Double[elmts.size()][elmts.size()];
        for(X ce : elmts) { 
            Cluster c= new SingletonCluster<X>(ce);
            clusters.add(c);
            reps.add(repr.getRepresentative(c));
        }
        for (int i = 0; i < elmts.size(); i++) {
            for (int j = 0; j < elmts.size(); j++) {
                X r1 = reps.get(i);
                X r2 = reps.get(j);
                distances[i][j] = metric.evaluate(r1,r2);
            }
        }
        
        int nclusters = clusters.size();
        while(nclusters > 1) { 
            int mini = -1, minj = -1;
            double mindist = Double.MAX_VALUE;

            for(int i = 0; i < clusters.size() - 1; i++) { 
                if (clusters.get(i) == null) {continue;}
                for(int j = i + 1; j < clusters.size(); j++) {
                    if (clusters.get(j) == null) {continue;}
                    double d;
                    if (Double.isNaN(distances[i][j])) {
                        X r1 = reps.get(i);
                        X r2 = reps.get(j);
                        d = metric.evaluate(r1, r2);
                        distances[i][j] = d;
                    } else {
                        d = distances[i][j];
                    }
                    if(!Double.isNaN(d) && d < mindist) { 
                        mindist = d;
                        mini = i; minj = j;
                    }
                }
            }
		
            if (mini == -1) {
                break;
            }
            if (mindist > maxDistanceToAccept) {
                break;
            }
            Cluster<X> left = clusters.get(mini), right = clusters.get(minj);
            if (left == null || right == null) {
                throw new NullPointerException("left is " + left + " + from " + mini +".  right is " + right
                                               + " from " + minj);
            }
            ClusterNode<X> node = new ClusterNode<X>(left, right);
            clusters.set(minj,null);
            reps.set(minj,null);
            for (int i = 0; i < elmts.size(); i++) {
                distances[i][minj] = Double.NaN;
                distances[minj][i] = Double.NaN;
                distances[i][mini] = Double.NaN;
                distances[mini][i] = Double.NaN;
            }
            clusters.set(mini,node);
            reps.set(mini,repr.getRepresentative(node));
            nclusters--;
            //System.out.println("# Clusters: " + nclusters + "(" + mindist + ")");
        }
        Vector<Cluster<X>> output = new Vector<Cluster<X>>();
        for (int i = 0; i < clusters.size(); i++) {
            if (clusters.get(i) != null) {
                output.add(clusters.get(i));
            }
        }
        return output;
    }

}

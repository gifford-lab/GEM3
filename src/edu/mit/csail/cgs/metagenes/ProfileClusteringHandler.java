package edu.mit.csail.cgs.metagenes;

import java.util.*;

import edu.mit.csail.cgs.clustering.*;
import edu.mit.csail.cgs.clustering.vectorcluster.*;
import edu.mit.csail.cgs.clustering.hierarchical.*;

public class ProfileClusteringHandler {

	private ProfileClusterRepresentative repr;
	private PairwiseElementMetric<ProfileClusterable> metric;
	private HierarchicalClustering<ProfileClusterable> clustering;
	
	public ProfileClusteringHandler(BinningParameters bps) { 
		repr = new ProfileClusterRepresentative(bps);
		metric = new EuclideanDistance<ProfileClusterable>();
		clustering = new HierarchicalClustering<ProfileClusterable>(repr, metric);
	}
	
	public Vector<Integer> runClustering(Vector<Profile> profs) { 
		Vector<Integer> indices = new Vector<Integer>();
		
		Vector<ProfileClusterable> clusterables = new Vector<ProfileClusterable>();
		for(int i = 0; i < profs.size(); i++) { 
			clusterables.add(new ProfileClusterable(i, profs.get(i)));
		}
		
		Collection<Cluster<ProfileClusterable>> clusterTree = clustering.clusterElements(clusterables);
		for(Cluster<ProfileClusterable> tree : clusterTree) { 
			appendIndices(indices, tree);
		}
		
		return indices;
	}
	
	private void appendIndices(Vector<Integer> indices, Cluster<ProfileClusterable> tree) { 
		if(tree instanceof ClusterNode) { 
			ClusterNode<ProfileClusterable> treeNode = (ClusterNode<ProfileClusterable>)tree;
			appendIndices(indices, treeNode.getLeft());
			appendIndices(indices, treeNode.getRight());
		} else { 
			for(ProfileClusterable pc : tree.getElements()) { 
				Integer idx = pc.getIndex();
				if(idx != null) { 
					indices.add(idx);
				} else { 
					System.err.println("Null index encountered...");
				}
			}
		}
	}
}

class ProfileClusterRepresentative implements ClusterRepresentative<ProfileClusterable> {
	
	private int index;
	private BinningParameters params;
	
	public ProfileClusterRepresentative(BinningParameters bps) { 
		index = 0;
		params = bps;
	}

	public ProfileClusterable getRepresentative(Cluster<ProfileClusterable> c) {
		String name = String.format("cluster-rep-%d", index++);
		MetaProfile mp = new MetaProfile(name, params);
		
		for(ProfileClusterable pc : c.getElements()) { 
			Profile p = pc.getProfile();
			mp.addProfile(p);
		}
		
		return new ProfileClusterable(mp);
	} 
	
}

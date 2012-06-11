package edu.mit.csail.cgs.clustering.pairedhitcluster;

import edu.mit.csail.cgs.clustering.Cluster;
import edu.mit.csail.cgs.clustering.ClusterRepresentative;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.projects.readdb.PairedHit;

public class PairedHitClusterRepresentative implements
		ClusterRepresentative<PairedHitClusterable> {

	public PairedHitClusterable getRepresentative(
			Cluster<PairedHitClusterable> c) {
		int leftchrom = 0;
		long leftpos = 0;
		short leftlen = 0;
		int rightchrom = 0;
		long rightpos = 0;
		short rightlen = 0;
		float weight = 0f;
		Genome g = null;
		for (PairedHitClusterable phc : c.getElements()) {
			PairedHit hit = phc.getHit();
			leftchrom = hit.leftChrom;
			leftpos += hit.leftPos;
			leftlen = hit.leftLength;
			rightchrom = hit.rightChrom;
			rightpos += hit.rightPos;
			rightlen = hit.rightLength;
			weight = hit.weight;
			g = phc.getGenome();
		}
		int avgleftpos = (int)(leftpos / c.size());
		int avgrightpos = (int)(rightpos / c.size());
		//leftpos /= c.size();
		//rightpos /= c.size();
		return new PairedHitClusterable(new PairedHit(leftchrom, avgleftpos, true, leftlen, rightchrom, avgrightpos, true, rightlen, c.getElements().size()),g);
	}

}

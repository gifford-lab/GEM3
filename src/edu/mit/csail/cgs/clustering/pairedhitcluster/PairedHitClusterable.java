package edu.mit.csail.cgs.clustering.pairedhitcluster;

import edu.mit.csail.cgs.clustering.vectorcluster.VectorClusterElement;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.projects.readdb.PairedHit;

public class PairedHitClusterable implements VectorClusterElement {
	
	private PairedHit hit;
	private Genome g;
	
	public PairedHitClusterable(PairedHit hit, Genome g) {
		if (hit.leftStrandedPoint(g).compareTo(hit.rightStrandedPoint(g))<0) {
			this.hit = hit.flippedCopy();
		} else {
			this.hit = hit;
		}
		this.g = g;
	}
	
	public PairedHit getHit() {
		return hit;
	}
	
	public Genome getGenome() {
		return g;
	}

	public int dimension() {
		return 2;
	}

	public double getValue(int i) {
		return i==0 ? hit.leftPos : hit.rightPos;
	}

	public boolean isMissingValue(int i) {
		return false;
	}

	public int numMissingValues() {
		return 0;
	}

	public String getTag(String k) {
		return null;
	}

	public boolean hasTag(String k) {
		return false;
	}

}

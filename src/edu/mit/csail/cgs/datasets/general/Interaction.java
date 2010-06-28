package edu.mit.csail.cgs.datasets.general;

import edu.mit.csail.cgs.datasets.species.Genome;

public class Interaction {

	private Region left, right;
	private int count;
	private double pval;
	
	public Interaction(Region left, Region right) {
		this.left = left;
		this.right = right;
	}
	
	public Interaction(Region left, Region right, int count, double pval) {
		this.left = left;
		this.right = right;
		this.count = count;
		this.pval = pval;
	}

	public Region getLeft() {
		return left;
	}

	public Region getRight() {
		return right;
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof Interaction)) {
			return false;
		} else {
			Interaction tmp = (Interaction)obj;
			return (this.left.equals(tmp.left) && this.right.equals(tmp.right)) || (this.left.equals(tmp.right) && this.right.equals(tmp.left));
		}
	}

	@Override
	public int hashCode() {
		return left.hashCode() + right.hashCode();
	}
	
	public String toString() {
		return left.toString()+"\t"+right.toString()+"\t"+count+"\t"+pval;
	}
	
}

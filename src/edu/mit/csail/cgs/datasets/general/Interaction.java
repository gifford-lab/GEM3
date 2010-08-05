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
	
	public boolean similarTo(Interaction other, int distance) {
		if (((this.left.getChrom().equals(other.left.getChrom()) && this.left.distance(other.left) <= distance) && 
				(this.right.getChrom().equals(other.right.getChrom()) && this.right.distance(other.right) <= distance)) || 
				((this.left.getChrom().equals(other.right.getChrom()) && this.left.distance(other.right) <= distance) && 
						(this.right.getChrom().equals(other.left.getChrom()) && this.right.distance(other.left) <= distance))) {
			return true;
		} else {
			return false;
		}
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

package edu.mit.csail.cgs.utils.probability;

public class Sample implements Comparable {

	private Double value;
	private int group;
	
	public Sample(double value, int group) {
		this.value = value;
		this.group = group;
	}
	
	public double getValue() {
		return value;
	}
	
	public int getGroup() {
		return group;
	}
	
	public int compareTo(Object o) {
		return this.value.compareTo(((Sample)o).value);
	}

}

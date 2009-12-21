package edu.mit.csail.cgs.utils.stats;

public class DoubleSort implements Comparable {
	public int pos = 0;
	public double value = 0.0;
	public int compareTo(Object o) {
		DoubleSort d = (DoubleSort) o;
		
		if (value > d.value)
			return 1;
		
		if (value == d.value)
			return 0;
		
		return -1;
	}
}

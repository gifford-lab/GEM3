package edu.mit.csail.cgs.clustering;

public class SimpleClusterable implements Clusterable {

	private String name;
	
	public SimpleClusterable(String name) {
		this.name = name;
	}
	
	public String name() {
		// TODO Auto-generated method stub
		return name;
	}
	
	public boolean equals(Object o) {
		if (o instanceof SimpleClusterable) {
			return this.name.equals(((SimpleClusterable)o).name);
		} else {
			return false;
		}
	}

}

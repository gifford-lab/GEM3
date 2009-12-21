package edu.mit.csail.cgs.metagenes;

import edu.mit.csail.cgs.clustering.vectorcluster.VectorClusterElement;

public class ProfileClusterable implements VectorClusterElement {
	
	private Integer index;
	private Profile profile;
	
	public ProfileClusterable(int idx, Profile p) { 
		index = idx;
		profile = p;
	}
	
	public ProfileClusterable(Profile p) { 
		profile = p;
	}
	
	public Integer getIndex() { return index; }
	
	public Profile getProfile() { 
		return profile;
	}

	public int dimension() {
		return profile.length();
	}

	public String getTag(String k) {
		return null;
	}

	public double getValue(int i) {
		return profile.value(i);
	}

	public boolean hasTag(String k) {
		return false;
	}

	public boolean isMissingValue(int i) {
		return false;
	}

	public int numMissingValues() {
		return 0;
	}

}

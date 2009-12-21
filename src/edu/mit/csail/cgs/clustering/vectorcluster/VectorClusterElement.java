package edu.mit.csail.cgs.clustering.vectorcluster;

/**
 * interface for a clusterable object that acts as a vector of doubles.
 * 
 * @author Timothy Danford
 */
public interface VectorClusterElement {
	public int dimension();
	
	/**
	 * getValue() should return a double for every integer index in the
	 * range [0, dimension()).  The only exception to this is if 
	 * isMissingValue(i) returns true.
	 * 
	 * @param i the index of the value
	 * @return the value
	 */
	public double getValue(int i);

	/**
	 * Checks whether the value for a given index is "missing."
	 * 
	 * @param i the index
 	 * @return true, if the value is missing; false otherwise.
	 */
	public boolean isMissingValue(int i);
	
	/**
	 * Counts the number of values which are missing.
	 * @return An integer-value of the number of missing values.
	 */
	public int numMissingValues();
	
	/**
	 * Gets the "tag" associated with a particular string key.
	 * @param k The tag's key.
	 * @return The value of the tag.
	 */
	public String getTag(String k);
	
	/**
	 * Checks whether there is a tag associated with the given key.
	 * @param k The tag-key.
	 * @return True if the key has a tag, false otherwise.
	 */
	public boolean hasTag(String k);
}

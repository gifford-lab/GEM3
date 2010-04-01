package edu.mit.csail.cgs.datasets.motifs;


public interface BackgroundModelFrequencySupport {

	/**
	 * Get the probability for the specified kmer from the distribution of all
   * kmers of that length
	 * @param kmer the kmer to lookup
	 * @return the probability for the kmer, or 0 if it's not in the map
	 */
  public abstract double getFrequency(String kmer);  
  
  /**
   * Get the probability for the kmer corresponding to the specified intVal and
   * kmerLen from the distribution of all kmers of that length
   * @param intVal the integer representation of the kmer to lookup
   * @param kmerLen the length of the kmer to lookup
   * @return
   */
  public abstract double getFrequency(int intVal, int kmerLen);
	
  
  /**
   * Remove strandedness from the model by setting reverse-complements to have
   * equal probabilities/counts
   */
  public abstract void degenerateStrands();
  
}

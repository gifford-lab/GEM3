package edu.mit.csail.cgs.datasets.motifs;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.stats.Fmath;

/**
 * @author mahony
 * @author rca This class represents a background model where the values are
 *         frequencies. For example the probabilities of all 16 possible
 *         dinucleotides sum to 1.
 */
public class FrequencyBackgroundModel extends BackgroundModel implements BackgroundModelFrequencySupport {
  
  
  public FrequencyBackgroundModel(BackgroundModelMetadata md) throws NotFoundException {
    super(md);
    if (!BackgroundModelImport.FREQUENCY_TYPE_STRING.equals(md.getDBModelType())) {
      throw new IllegalArgumentException("Metadata model type must be FREQUENCY");
    }
  }
  
  public FrequencyBackgroundModel(String name, Genome gen) {
    super(name, gen);
  }

  
  public FrequencyBackgroundModel(String name, Genome gen, int modelLength) {
    super(name, gen, modelLength);
  }
  
  
  public FrequencyBackgroundModel(CountsBackgroundModel cbg) {
  	super(cbg);
  	for (int i = 1; i <= cbg.getMaxKmerLen(); i++) {
  		cbg.computeFrequencies(i);
  		modelProbs[i].putAll(cbg.modelProbs[i]);
  	}
  }
  
  
  protected void init() {
    this.setDBModelType(BackgroundModelImport.FREQUENCY_TYPE_STRING);
  }
  
  /**
   * @see BackgroundModel
   */
  public Set<String> getKmers(int kmerLen) {
  	return modelProbs[kmerLen].keySet();
  }

  
  /**
   * @see BackgroundModel
   * On the fly compute the markov probability from the known frequencies
   */
	public double getMarkovProb(String kmer) {
		String prevBases = kmer.substring(0, kmer.length() - 1);
		double total = 0.0;
		for (char currBase : BackgroundModel.BASE_ORDER) {
			String currKmer = prevBases + currBase;
			total = total + this.getFrequency(currKmer);
		}
		if (total == 0.0) {
			return 0.0;
		}
		else {
			return (this.getFrequency(kmer) / total);
		}
	}
	
	
  /**
   * @see BackgroundModelFrequencySupport
   */
  public double getFrequency(String kmer) {
  	if (modelProbs[kmer.length()].containsKey(kmer)) {
  		return (modelProbs[kmer.length()].get(kmer));
  	}
  	else {
  		return 0.0;
  	}			
  }
  
  
  /**
   * @see BackgroundModelFrequencySupport
   */
  public double getFrequency(int intVal, int kmerLen) {
    return this.getFrequency(BackgroundModel.int2seq(intVal, kmerLen));
  }
	
  
  /**
   * Set the model to have the frequency probabilities specified in the map 
   * @param probs A map of kmers of identical length and their probabilities,
   * which must sum to 1.
   */
  public void setKmerFrequencies(Map<String, Double> probs) {
  	//check that the specified probabilties sum to 1 and are for kmers of the
  	//correct length
  	double total = 0.0;
  	int kmerLen = probs.keySet().iterator().next().length();
  	if (kmerLen < 1) {
  		throw new IllegalArgumentException("Kmer length must be >= 1.");
  	}
  	
  	for (String kmer : probs.keySet()) {
  		if (kmer.length() == kmerLen) {
  			total += probs.get(kmer);
  		}
  		else {
  			throw new IllegalArgumentException("Kmers in map must all have the same length");
  		}
  	}
  	
  	//set the model to use the specified probabilities
  	if (Fmath.isEqualWithinLimits(total, 1.0, BackgroundModel.EPSILON)) {
  		modelProbs[kmerLen] = new HashMap<String, Double>(probs);
  	}
		else {
			throw new IllegalArgumentException("Probabilities must sum to 1, but instead sum to " + total);
		}
  	
		//reset the isStranded variable to null to indicate unknown strandedness
  	isStranded = null;
  }
  
  
  /**
   * @see BackgroundModel
   * Check whether the reverse complement kmers have equal frequencies for all
   * the kmer lengths 
   */
  public boolean checkAndSetIsStranded() {
    int currKmerLen = 1;
    while (currKmerLen <= modelProbs.length) {
      List<Pair<Integer, Integer>> revCompPairs = BackgroundModel.computeDistinctRevCompPairs(currKmerLen);

      for (Pair<Integer, Integer> rcPair : revCompPairs) {
        if (!Fmath.isEqualWithinLimits(this.getFrequency(currKmerLen, rcPair.car()), 
            this.getFrequency(currKmerLen,rcPair.cdr()), BackgroundModel.EPSILON)) {
          isStranded = true;
          return isStranded;
        }
      }
    }
    
    /**
     * If all kmerlens have been checked and the method hasn't returned then the
     * model is not stranded
     */
    isStranded = false;
    return isStranded;
  }
  

  /**
   * @see BackgroundModelFrequencySupport
   * For a frequency based model this is accomplished by setting reverse
   * complements to have the average of their frequencies, so there is no harm
   * in running this method even if the model is already destranded.
   */
  public void degenerateStrands() {
    for (int currKmerLen = 1; currKmerLen <= this.getMaxKmerLen(); currKmerLen++) {
    	List<Pair<Integer, Integer>> revCompPairs = BackgroundModel.computeDistinctRevCompPairs(currKmerLen);

    	for (Pair<Integer, Integer> rcPair : revCompPairs) {
   			double freq = this.getFrequency(rcPair.car(), currKmerLen);
    		double revCompFreq = this.getFrequency(rcPair.cdr(), currKmerLen);
    		double newFreq = (freq + revCompFreq) / 2.0;
    		modelProbs[currKmerLen].put(BackgroundModel.int2seq(rcPair.car(), currKmerLen), newFreq);
    		modelProbs[currKmerLen].put(BackgroundModel.int2seq(rcPair.cdr(), currKmerLen), newFreq);
    	}    	
    }
    isStranded = false;
  }
  
  
  /**
   * Check if this model is normalized properly.
   * @return the shortest kmerlen for which the model is not normalized, or -1 
   * if normalized properly.
   */
  public int verifyNormalization() {
    //iterate over each order level of the model
    for (int i = 1; i <= this.getMaxKmerLen(); i++) {
      Double total = 0.0;
      //iterate over all the n-mers of the order summing up the values
      for (int k = 0; k < (int) Math.pow(4, i); k++) {
        String currMer = int2seq(k, i);
        total += modelProbs[i].get(currMer);
      }
      if (!Fmath.isEqualWithinLimits(total, 1.0, 1E-6)) {
        return i;
      }
    }
    return -1;
  }
}

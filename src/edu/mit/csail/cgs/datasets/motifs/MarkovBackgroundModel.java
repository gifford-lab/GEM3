package edu.mit.csail.cgs.datasets.motifs;

import java.util.List;
import java.util.Set;

import cern.colt.matrix.*;
import cern.colt.matrix.impl.*;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.stats.Fmath;

/**
 * @author mahony
 * @author rca This class represents a background model where the values are
 *         conditional probabilities. For dinucleotides, for example, p(A|A) +
 *         p(C|A) + p(G|A) + p(T|A) = 1.
 */
public class MarkovBackgroundModel extends BackgroundModel {
  
  public MarkovBackgroundModel(String name, Genome gen) {
    super(name, gen);
  }


  public MarkovBackgroundModel(String name, Genome gen, int modelLength) {
    super(name, gen, modelLength);
  }
  
  
  /**
   * Construct a Markov Background Model from an existing Frequency Background
   * Model. 
   * @param fbg
   */
  public MarkovBackgroundModel(FrequencyBackgroundModel fbg) {
  	super(fbg);
    
    //iterate over each order level of the model
    for (int i = 1; i <= this.getMaxKmerLen(); i++) {
      //iterate over all sets of conditions for that order
      for (int k = 0; k < (int) Math.pow(4, i); k += 4) {
        Double total = 0.0;
        //iterate over the 4 outcomes for the conditional probability 
        //summing up the values
        for (int b = 0; b < 4; b++) {
          String currKmer = int2seq(k + b, i);
          total += fbg.getFrequency(currKmer);
        }
        //if the total is 0 skip these kmers and keep the model sparse
        if (total > 0) {
        	//iterate over the 4 outcomes normalizing 
        	for (int b = 0; b < 4; b++) {
        		String currKmer = int2seq(k + b, i);
        		modelProbs[i].put(currKmer, fbg.getFrequency(currKmer) / total);
        	}
        }
      }
    }    
  }
  
  
  protected void init() {
    this.setDBModelType(BackgroundModelImport.MARKOV_TYPE_STRING);
  }
  
  
  /**
   * @see BackgroundModel
   */
  public Set<String> getKmers(int kmerLen) {
  	return modelProbs[kmerLen].keySet();
  }

  
  /**
   * Construct a Markov Background Model from an existing Counts Background
   * Model. 
   * @param cbg
   */
  public MarkovBackgroundModel(CountsBackgroundModel cbg) {
  	this(new FrequencyBackgroundModel(cbg));
  }

  
  /**
   * Return the markov probability for the last base of the specified kmer
   * conditioned upon the preceding bases.
   */
	public double getMarkovProb(String kmer) {
		if (modelProbs[kmer.length()].containsKey(kmer)) {
			return modelProbs[kmer.length()].get(kmer);
		}
		else {
			return 0.0;
		}
	}

	
	/**
	 * Sets the 4 markov probabilities that are conditioned on the specified 
	 * string of previous bases.   
	 * @param prevBases The bases that the probabilities are conditionally
	 * dependent upon
	 * @param aProb P( A | prevBases)
	 * @param cProb P( C | prevBases)
	 * @param gProb P( G | prevBases)
	 * @param tProb P( T | prevBases)
	 */
	public void setMarkovProb(String prevBases, double aProb, double cProb, double gProb, double tProb) {
		double total = aProb + cProb + gProb + tProb;
		if (Fmath.isEqualWithinLimits(total, 1.0, BackgroundModel.EPSILON)) {
			int kmerLen = prevBases.length() + 1;
			modelProbs[kmerLen].put(prevBases + "A", aProb);
			modelProbs[kmerLen].put(prevBases + "C", cProb);
			modelProbs[kmerLen].put(prevBases + "G", gProb);
			modelProbs[kmerLen].put(prevBases + "T", tProb);
		}
		else {
			throw new IllegalArgumentException("Probabilities must sum to 1, but instead sum to " + total);
		}
		
		//reset the isStranded variable to null to indicate unknown strandedness
		isStranded = null;
	}
	
	
  /**
   * Check if this model is normalized properly.
   * @return an array of string containing the kmers for which the model isn't
   * normalized, or null if it is normalized correctly
   */
  public String[] verifyNormalization() {
    //iterate over each order level of the model
    for (int i = 1; i <= this.getMaxKmerLen(); i++) {
      //iterate over all sets of conditions for that order
      for (int k = 0; k < (int) Math.pow(4, i); k += 4) {
        Double total = 0.0;
        //iterate over the 4 outcomes for the conditional probability 
        //summing up the values
        String[] currMers = new String[4];
        for (int b = 0; b < 4; b++) {
          String currMer = int2seq(k + b, i);
          currMers[b] = currMer;
          total += this.getMarkovProb(currMer);
        }
        if (!Fmath.isEqualWithinLimits(total, 1.0, 1E-6)) {
          return currMers;
        }
      }
    }
    return null;
  }

  
  /**
   * Checks whether this model is based on a single strand of sequence or
   * double-stranded sequence. For the Markov model a linear system can be 
   * created from the log-space version of the equations for normalizing the 
   * model and equations indicating the equality of reverse complement kmers.
   * The system will be underconstrained, so there will always be a solution,
   * but once converted back to linear-space additional constraints can be 
   * verified to determine whether the solution is valid. 
   */
  public boolean checkAndSetIsStranded() {
    int currKmerLen = 1;
    while (currKmerLen <= modelProbs.length) {      
      //set up and solve the log-space linear system for the  current kmers
      List<Pair<Integer, Integer>> revCompPairs = BackgroundModel.computeDistinctRevCompPairs(currKmerLen);
      int numCurrKmers = (int)Math.pow(4, currKmerLen);
      DoubleMatrix2D linSys = MarkovBackgroundModel.createStrandednessLinearSystem(numCurrKmers, revCompPairs);
      
      DoubleMatrix1D logMarkov = new DenseDoubleMatrix1D(numCurrKmers);
      for (int i = 0; i < logMarkov.size(); i++) {
        Double prob = this.getMarkovProb(currKmerLen, i);
        if (prob != null) {
          logMarkov.setQuick(i, Math.exp(prob));
        }
        else {
          logMarkov.setQuick(i, -Double.MAX_VALUE);
        }
      }
            
      Algebra algebra = new Algebra();
      DoubleMatrix1D sol = algebra.mult(algebra.inverse(linSys), logMarkov).assign(Functions.exp);
      
      //now check that the rev. comp components of the solution are matched       
      for (Pair<Integer, Integer> rcPair : revCompPairs) {
        if (!Fmath.isEqualWithinLimits(sol.getQuick(rcPair.car()), sol.getQuick(rcPair.cdr()), BackgroundModel.EPSILON)) {
          isStranded = true;
          return isStranded;
        }
      }
      
      //finally, check that the sums used to normalize match correctly
      for (int i = 0; i < numCurrKmers; i+= 4) {
        double sum = sol.viewPart(i, 4).zSum();
        int solSumIndex = numCurrKmers + (i / 4);
        if (!Fmath.isEqualWithinLimits(sum, sol.getQuick(solSumIndex), BackgroundModel.EPSILON)) {
          isStranded = true;
          return isStranded;
        }
      }
      
      currKmerLen++;
    }
    
    /**
     * If all kmerlens have been checked and the method hasn't returned then the
     * model is not stranded
     */
    isStranded = false;
    return isStranded;
  }
  
  
  /**
   * Create a linear system of the log-space versions of the equations used to
   * normalize a frequency/count model as a markov model.
   * @param numKmers
   * @param revCompPairs
   * @return
   */
  private static DoubleMatrix2D createStrandednessLinearSystem(int numKmers, List<Pair<Integer, Integer>> revCompPairs) {       
    int numSeqGroups = numKmers / 4;
    DoubleMatrix2D linSys = new SparseDoubleMatrix2D(numKmers + revCompPairs.size(), (numKmers + numSeqGroups));
    for (int i = 0; i < numKmers; i++) {
      //set the kmer's entry for itself
      linSys.setQuick(i, i, 1);
      //set the kmer's entry for the kmer group sum
      int groupIndex = numKmers + (i / 4);
      linSys.setQuick(i, groupIndex, -1);      
    }
    
    for (int i = 0; i < revCompPairs.size(); i++) {
      Pair<Integer, Integer> revComps = revCompPairs.get(i);
      linSys.setQuick(numKmers + i, revComps.car(), 1);
      linSys.setQuick(numKmers + i, revComps.cdr(), -1);
    }
    
    return linSys;
  }     
}

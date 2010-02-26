package edu.mit.csail.cgs.datasets.motifs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import cern.colt.matrix.*;
import cern.colt.matrix.impl.*;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.utils.stats.Fmath;

/**
 * @author mahony
 * @author rca This class represents a background model where the values are
 *         conditional probabilities. For dinucleotides, for example, p(A|A) +
 *         p(C|A) + p(G|A) + p(T|A) = 1.
 */
public class MarkovBackgroundModel extends BackgroundModel {
  
  public MarkovBackgroundModel() {
    super();
  }


  public MarkovBackgroundModel(int modelLength) {
    super(modelLength);
  }
  
  
  public int getMarkovOrder() {
    return model.length - 2;
  }
  
  
  public boolean checkAndSetIsStranded() {
    int currKmerLen = 1;
    while (currKmerLen <= model.length) {      
      //set up and solve the log-space linear system for the  current kmers
      List<Pair<Integer, Integer>> revCompPairs = BackgroundModel.computeRevCompPairs(currKmerLen);
      int numCurrKmers = (int)Math.pow(4, currKmerLen);
      DoubleMatrix2D linSys = MarkovBackgroundModel.createStrandednessLinearSystem(numCurrKmers, revCompPairs);
      
      DoubleMatrix1D logMarkov = new DenseDoubleMatrix1D(numCurrKmers);
      for (int i = 0; i < logMarkov.size(); i++) {
        Double prob = this.getModelProb(currKmerLen, i);
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
      
      //finally, check that the sums match correctly
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
    
  
  /**
   * Check with the specified Markov Background Model is normalized properly.
   * @param mbg
   * @return an array of string containing the kmers for which the model isn't
   * normalized, or null if it is normalized correctly
   */
  public static String[] isModelNormalized(MarkovBackgroundModel mbg) {
    //iterate over each order level of the model
    for (int i = 1; i <= mbg.getMaxKmerLen(); i++) {
      //iterate over all sets of conditions for that order
      for (int k = 0; k < (int) Math.pow(4, i); k += 4) {
        Double total = 0.0;
        //iterate over the 4 outcomes for the conditional probability 
        //summing up the values
        String[] currMers = new String[4];
        for (int b = 0; b < 4; b++) {
          String currMer = int2seq(k + b, i);
          currMers[b] = currMer;
          if (mbg.model[i].containsKey(currMer)) {
            total += mbg.model[i].get(currMer).car();
          }
        }
        if (Fmath.isEqualWithinLimits(total, 1.0, 1E-6)) {
          return currMers;
        }
      }
    }
    return null;
  }
}

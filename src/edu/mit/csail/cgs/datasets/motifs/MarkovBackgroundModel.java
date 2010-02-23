package edu.mit.csail.cgs.datasets.motifs;

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
  
  
  public int getModelOrder() {
    return model.length - 2;
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
          String currMer = intToSeq(k + b, i);
          currMers[b] = currMer;
          if (mbg.model[i].containsKey(currMer)) {
            total += mbg.model[i].get(currMer);
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

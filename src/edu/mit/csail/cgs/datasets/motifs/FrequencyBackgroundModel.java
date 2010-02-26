package edu.mit.csail.cgs.datasets.motifs;

import java.util.Arrays;
import java.util.List;

import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.utils.stats.Fmath;

/**
 * @author mahony
 * @author rca This class represents a background model where the values are
 *         frequencies. For example the probabilities of all 16 possible
 *         dinucleotides sum to 1.
 */
public class FrequencyBackgroundModel extends BackgroundModel {
  
  public FrequencyBackgroundModel() {
    super();
  }

  
  public FrequencyBackgroundModel(int modelLength) {
    super(modelLength);
  }
  
  
  public boolean checkAndSetIsStranded() {
    int currKmerLen = 1;
    while (currKmerLen <= model.length) {
      List<Pair<Integer, Integer>> revCompPairs = BackgroundModel.computeRevCompPairs(currKmerLen);
      
      for (Pair<Integer, Integer> rcPair : revCompPairs) {
        if (!Fmath.isEqualWithinLimits(this.getModelProb(currKmerLen, rcPair.car()), this.getModelProb(currKmerLen, rcPair.cdr()), BackgroundModel.EPSILON)) {
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
   * Remove strandedness from the model by setting reverse-complements to have
   * indentical probabilities equal to the average of both
   */
//  public void degenerateStrands() {
//    for (int i = 1; i <= this.getMaxKmerLen(); i++) {
//      boolean[] check = new boolean[(int) Math.pow(4, i)];
//      Arrays.fill(check, false);
//      for (int k = 0; k < Math.pow(4, i); k++) {
//        String currMer = intToSeq(k, i);
//        if (!check[k]) {
//          String revMer = SequenceUtils.reverseComplement(currMer);
//          int revID = seqToInt(revMer);
//          check[k] = true;
//          check[revID] = true;
//
//          double newVal = (model[i].get(currMer) + model[i].get(revMer)) / 2.0;
//          model[i].put(currMer, newVal);
//          model[i].put(revMer, newVal);
//        }
//      }
//    }
//  }
//  
  
  /**
   * Check with the specified frequency background model is normalized properly.
   * @param fbg
   * @return the kmerlen for which the model is not normalized, or -1 if 
   * normalized properly.
   */
  public static int isModelNormalized(FrequencyBackgroundModel fbg) {
    //iterate over each order level of the model
    for (int i = 1; i <= fbg.getMaxKmerLen(); i++) {
      Double total = 0.0;
      //iterate over all the n-mers of the order summing up the values
      for (int k = 0; k < (int) Math.pow(4, i); k++) {
        String currMer = int2seq(k, i);
        if (fbg.model[i].containsKey(currMer)) {
          total += fbg.model[i].get(currMer).car();
        }
      }
      if (!Fmath.isEqualWithinLimits(total, 1.0, 1E-6)) {
        return i;
      }
    }
    return -1;
  }
}

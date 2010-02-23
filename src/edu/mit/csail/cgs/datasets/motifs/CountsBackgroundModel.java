/**
 * 
 */
package edu.mit.csail.cgs.datasets.motifs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.verbs.ChromRegionIterator;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.io.parsing.FASTAStream;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;

/**
 * @author rca
 *
 */
public class CountsBackgroundModel extends BackgroundModel {
  

  public CountsBackgroundModel() {
    super();
  }


  public CountsBackgroundModel(int modelLength) {
    super(modelLength);
  }

  
  /**
   * Examine all the appropriately sized kmers from the specified sequence and
   * add them to this model
   * @param sequence
   */
  protected void addKmerCountsFromSequence(String sequence) {
    
    //handle all the positions with complete length kmers
    String currSub = "";
    int maxLen = this.getMaxKmerLen();
    for (int i = 0; i <= sequence.length() - maxLen; i++) {
      currSub = sequence.substring(i, i + maxLen);
      for (int k = 1; k <= maxLen; k++) {
        String currKmer = currSub.substring(0, k);
        if (model[k].containsKey(currKmer)) {
          model[k].put(currKmer, model[k].get(currKmer));
        }
        else {
          model[k].put(currKmer, 1.0);
        }
      }
    }

    //Put in the last few kmers
    int start = sequence.length() - maxLen + 1;
    for (int i = start; i < sequence.length(); i++) {
      currSub = sequence.substring(i);
      for (int k = 1; k <= currSub.length(); k++) {
        String currKmer = currSub.substring(0, k);
        if (model[k].containsKey(currKmer)) {
          model[k].put(currKmer, model[k].get(currKmer));
        }
        else {
          model[k].put(currKmer, 1.0);
        }
      }
    }
  }
  
    

  /**
   * Remove strandedness from the model by setting reverse-complements to have
   * equal probabilities
   */
  public void degenerateStrands() {
    for (int i = 1; i <= this.getMaxKmerLen(); i++) {
      boolean[] check = new boolean[(int) Math.pow(4, i)];
      Arrays.fill(check, false);
      for (int k = 0; k < Math.pow(4, i); k++) {
        String currMer = intToSeq(k, i);
        if (!check[k]) {
          String revMer = SequenceUtils.reverseComplement(currMer);
          int revID = seqToInt(revMer);
          check[k] = true;
          check[revID] = true;

          double newVal = (model[i].get(currMer) + model[i].get(revMer)) / 2.0;
          model[i].put(currMer, newVal);
          model[i].put(revMer, newVal);
        }
      }
    }
  }
  
  
  public MarkovBackgroundModel normalizeToMarkovModel() {
    MarkovBackgroundModel mbg = new MarkovBackgroundModel(this.getMaxKmerLen());
    mbg.gen = this.gen;
    
    //iterate over each order level of the model
    for (int i = 1; i <= this.getMaxKmerLen(); i++) {
      //iterate over all sets of conditions for that order
      for (int k = 0; k < (int) Math.pow(4, i); k += 4) {
        Double total = 0.0;
        //iterate over the 4 outcomes for the conditional probability 
        //summing up the values
        for (int b = 0; b < 4; b++) {
          String currMer = intToSeq(k + b, i);
          if (model[i].containsKey(currMer)) {
            total += model[i].get(currMer);
          }
        }
        //iterate over the 4 outcomes normalizing 
        for (int b = 0; b < 4; b++) {
          String currMer = intToSeq(k + b, i);
          if (model[i].containsKey(currMer)) {
            Double prev = model[i].get(currMer);
            mbg.model[i].put(currMer, prev / total);
          }
          else {
            mbg.model[i].put(currMer, 0.0);
          }
        }
      }
    }
    return mbg;
  }
  
  
  public FrequencyBackgroundModel normalizeToFrequencyModel() {
    FrequencyBackgroundModel fbg = new FrequencyBackgroundModel(this.getMaxKmerLen());
    fbg.gen = this.gen;
    
    //iterate over each order level of the model
    for (int i = 1; i <= this.getMaxKmerLen(); i++) {
      Double total = 0.0;
      //iterate over all the n-mers of the order summing up the values
      for (int k = 0; k < (int) Math.pow(4, i); k++) {
        String currMer = intToSeq(k, i);
        if (model[i].containsKey(currMer)) {
          total += model[i].get(currMer);
        }
      }
      //iterate over all the n-mers of the order normalizing
      for (int k = 0; k < (int) Math.pow(4, i); k++) {
        String currMer = intToSeq(k, i);
        if (model[i].containsKey(currMer)) {
          Double prev = model[i].get(currMer);
          fbg.model[i].put(currMer, prev / total);
        }
        else {
          fbg.model[i].put(currMer, 0.0);
        }
      }
    }
    return fbg;
  }
  
  
  
  /**
   * 
   * @param gen
   * @return
   */
  public static CountsBackgroundModel modelFromRegionList(Genome gen){
    ArrayList <Region>chromList = new ArrayList<Region>();
    Iterator<NamedRegion> chroms = new ChromRegionIterator(gen);
    while (chroms.hasNext()) {
      chromList.add(chroms.next());
    }
    return CountsBackgroundModel.modelFromRegionList(gen, chromList);
  } 
  
 
  /**
   * 
   * @param gen
   * @param regionList
   * @return
   */
  public static CountsBackgroundModel modelFromRegionList(Genome gen, List<Region> regionList) {
    SequenceGenerator<Region> seqgen = new SequenceGenerator<Region>();
    seqgen.useCache(false);

    CountsBackgroundModel cbg = new CountsBackgroundModel();
    cbg.gen = gen;

    for (Region currR : regionList) {
      String tmpSeq = seqgen.execute(currR);
      String regionSeq = tmpSeq.toUpperCase();

      cbg.addKmerCountsFromSequence(regionSeq);
    }
    return cbg;
  }

  
  /**
   * 
   * @param stream
   * @return
   */
  public static CountsBackgroundModel modelFromFASTAStream(FASTAStream stream) {
    CountsBackgroundModel cbg = new CountsBackgroundModel();
    while (stream.hasNext()) {
      Pair<String, String> currSeq = stream.next();

      String tmpSeq = currSeq.cdr();
      String regionSeq = tmpSeq.toUpperCase();

      cbg.addKmerCountsFromSequence(regionSeq);
    }
    stream.close();
    return cbg;
  }
}

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
   * Overrides the superclass method so that only the model variable is used to
   * hold data
   */
  public Double getModelCount(String mer) {
    if (mer.length() > 0) {
      return (model[mer.length()].get(mer));
    }
    else {
      throw new IllegalArgumentException("Zero length kmer.");
    }
  }


  /**
   * Overrides the superclass method so that only the model variable is used to
   * hold data
   */
  public Double getModelCount(int kmerLen, int intVal) {
    if (kmerLen > 0) {
      return (model[kmerLen].get(BackgroundModel.intToSeq(intVal, kmerLen)));
    }
    else {
      throw new IllegalArgumentException("kmerLen must be greater than zero.");
    }
  }


  /**
   * Overrides the superclass method so that only the model variable is used to
   * hold data
   */
  public void setModelCount(String mer, double val) {
    if (mer.length() <= this.getMaxKmerLen() && mer.length() > 0) {
      model[mer.length()].put(mer, val);
    }
    else if (mer.length() < 1) {     
      throw new IllegalArgumentException("Zero length kmer.");      
    }
    else {
      throw new IllegalArgumentException("Kmer " + mer + " must have length less than model length (" + this.getMaxKmerLen() + ").");
    }
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
   * counts equal to the average of both. Strictly this is incorrect, because
   * if the reverse strand had been counted the total counts would be the sum 
   * of both, but setting to the average prevents the count from increasing if
   * the method is called repeatedly, and generally it won't cause problems.
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
  
  
  /**
   * Create a Markov Background Model by normalizing this model. 
   * 
   * Note: It will often be appropriate to call degenerateStrands before calling
   * this method. The Markov background model does not have this method.
   * @return
   */
  public MarkovBackgroundModel normalizeToMarkovModel() {
    MarkovBackgroundModel mbg = new MarkovBackgroundModel(this.getMaxKmerLen());
    mbg.gen = this.gen;
    mbg.counts = this.model.clone();
    
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
  
  
  /**
   * Create a Frequency Background Model by normalizing this model.
   * @return
   */
  public FrequencyBackgroundModel normalizeToFrequencyModel() {
    FrequencyBackgroundModel fbg = new FrequencyBackgroundModel(this.getMaxKmerLen());
    fbg.gen = this.gen;
    fbg.counts = this.model.clone();
    
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
   * Create a model from the entirety of the specified genome
   * @param gen
   * @return
   */
  public static CountsBackgroundModel modelFromWholeGenome(Genome gen){
    ArrayList <Region>chromList = new ArrayList<Region>();
    Iterator<NamedRegion> chroms = new ChromRegionIterator(gen);
    while (chroms.hasNext()) {
      chromList.add(chroms.next());
    }
    return CountsBackgroundModel.modelFromRegionList(gen, chromList);
  } 
  
 
  /**
   * Create a model from a list of regions from the specified genome
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
   * Create a model from a FASTAStream object
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

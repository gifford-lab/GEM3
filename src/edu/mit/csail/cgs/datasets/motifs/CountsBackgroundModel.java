/**
 * 
 */
package edu.mit.csail.cgs.datasets.motifs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.verbs.ChromRegionIterator;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.io.parsing.FASTAStream;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.utils.stats.Fmath;

/**
 * @author rca
 *
 */
public class CountsBackgroundModel extends BackgroundModel implements BackgroundModelFrequencySupport {
  
	 /**
   * Map for holding the kmer counts for the model. Each 
   * element of the array holds kmers whose length is the index of that element.
   * i.e. model[2] holds 2-mers. Accordingly, model[0] should always be empty.
   * 
   * Note: This may be null if count information is unavailable (as may be the
   * case for a frequency model or markov model)
   */
  private Map<String, Long>[] kmerCounts;

  
  public CountsBackgroundModel(String name, Genome gen) {
    this(name, gen, BackgroundModel.DEFAULT_MAX_KMER_LEN);
  }


  public CountsBackgroundModel(String name, Genome gen, int maxKmerLen) {
    super(name, gen, maxKmerLen);
    kmerCounts = new HashMap[maxKmerLen + 1];
    for (int i = 1; i <= maxKmerLen; i++) {
      kmerCounts[i] = new HashMap<String, Long>();
    }
    
    //set the model map array to be empty so that it can be used to lazily
    //shadow the counts with frequencies
    Arrays.fill(modelProbs, null);
  }
  
  
  protected void init() {
    //do nothing (for now)
  }
  
  /**
   * @see BackgroundModel
   */
  public Set<String> getKmers(int kmerLen) {
  	return kmerCounts[kmerLen].keySet();
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
		//build the cache if it's not ready
		if (modelProbs[kmer.length()] == null) { 
			this.computeFrequencies(kmer.length());
		}
		
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
   * 
   * @param mer
   * @return
   */
  public long getKmerCount(String kmer) {
		if (kmerCounts[kmer.length()].containsKey(kmer)) {
 			return (kmerCounts[kmer.length()].get(kmer));
 		}
 		else {
 			return 0;
 		}
  }


  /**
   * 
   * @param kmerLen
   * @param intVal
   * @return
   */
  public long getKmerCount(int kmerLen, int intVal) {
    return (this.getKmerCount(BackgroundModel.int2seq(intVal, kmerLen)));
  }

  
  public void setKmerCount(String kmer, long count) {  	
  	kmerCounts[kmer.length()].put(kmer, count);

  	//clear the frequency map for the kmer's length
  	modelProbs[kmer.length()] = null;
  	
  	//reset the isStranded variable to null to indicate unknown strandedness
  	isStranded = null;
  }
  
  
  public void addToKmerCount(String kmer, long addlCount) {
  	kmerCounts[kmer.length()].put(kmer, this.getKmerCount(kmer) + addlCount);
  	
  	//clear the frequency map for the kmer's length
  	modelProbs[kmer.length()] = null;

  	//reset the isStranded variable to null to indicate unknown strandedness
  	isStranded = null;
  }

  
  /**
   * @see addKmerCountsFromSequence(String sequence, boolean addRevComp)
   * @param sequence
   */
  public void addKmerCountsFromSequence(String sequence) {
  	this.addKmerCountsFromSequence(sequence, false);
  }

  
  /**
   * Examine all the appropriately sized kmers from the specified sequence and
   * add them to this model
   * @param sequence
   * @param addRevComp if true also add counts for the reverse complement
   */
  public void addKmerCountsFromSequence(String sequence, boolean addRevComp) {
    
    //handle all the positions with complete length kmers
    String currSub = "";
    int maxLen = this.getMaxKmerLen();
    for (int i = 0; i <= sequence.length() - maxLen; i++) {
      currSub = sequence.substring(i, i + maxLen);
      for (int k = 1; k <= maxLen; k++) {
        String currKmer = currSub.substring(0, k);
        this.addToKmerCount(currKmer, 1);
        this.addToKmerCount(SequenceUtils.reverseComplement(currKmer), 1);
      }
    }

    //Put in the last few kmers
    int start = sequence.length() - maxLen + 1;
    for (int i = start; i < sequence.length(); i++) {
      currSub = sequence.substring(i);
      for (int k = 1; k <= currSub.length(); k++) {
        String currKmer = currSub.substring(0, k);
        this.addToKmerCount(currKmer, 1);
        this.addToKmerCount(SequenceUtils.reverseComplement(currKmer), 1);
      }
    }
  	isStranded = !addRevComp;
  	
    //set the model map array to be empty so that it can be used to lazily
    //shadow the counts with frequencies
    Arrays.fill(modelProbs, null);

  }
  
  
  /**
   * @see BackgroundModel
   * Check whether the reverse complement kmers have equal counts for all
   * the kmer lengths 

   */
  public boolean checkAndSetIsStranded() {
    int currKmerLen = 1;
    while (currKmerLen <= modelProbs.length) {
      List<Pair<Integer, Integer>> revCompPairs = BackgroundModel.computeDistinctRevCompPairs(currKmerLen);
      
      for (Pair<Integer, Integer> rcPair : revCompPairs) {
        if (this.getKmerCount(currKmerLen, rcPair.car()) == this.getKmerCount(currKmerLen, rcPair.cdr())) {
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
   * For a count based model this is accomplished by setting reverse
   * complements to have the sum of their counts.
   */
  public void degenerateStrands() {
  	if (isStranded == null) {
  		this.checkAndSetIsStranded();
  	}
  	//only destrand the model if its stranded, otherwise the counts will blow up
  	if (isStranded) {
  		for (int currKmerLen = 1; currKmerLen <= this.getMaxKmerLen(); currKmerLen++) {
  			List<Pair<Integer, Integer>> revCompPairs = BackgroundModel.computeDistinctRevCompPairs(currKmerLen);

  			for (Pair<Integer, Integer> rcPair : revCompPairs) {
  				long count = this.getKmerCount(rcPair.car(), currKmerLen);;
  				long revCompCount = this.getKmerCount(rcPair.cdr(), currKmerLen);;
  				long newCount = count + revCompCount;
  				this.setKmerCount(BackgroundModel.int2seq(rcPair.car(), currKmerLen), newCount);
  				this.setKmerCount(BackgroundModel.int2seq(rcPair.cdr(), currKmerLen), newCount);
  			}    	
  		}
  	}
  	
    //set the model map array to be empty so that it can be used to lazily
    //shadow the counts with frequencies
    Arrays.fill(modelProbs, null);
  }
  
  
  /**
   * Compute the frequencies corresponding to the counts for kmers of the
   * specified length
   * @param kmerLen
   */
  protected void computeFrequencies(int kmerLen) {
  	double total = 0;
  	for (long count : kmerCounts[kmerLen].values()) {
  		total += count;
  	}
  	
  	modelProbs[kmerLen] = new HashMap<String, Double>();
  	if (total > 0) {
  		for (String kmer : kmerCounts[kmerLen].keySet()) {
  			modelProbs[kmerLen].put(kmer, ((double)this.getKmerCount(kmer) / total));
  		}
  	}
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

    CountsBackgroundModel cbg = new CountsBackgroundModel(null, gen);
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
    CountsBackgroundModel cbg = new CountsBackgroundModel(null, null);
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

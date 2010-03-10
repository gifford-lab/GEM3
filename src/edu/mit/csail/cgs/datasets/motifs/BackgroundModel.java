package edu.mit.csail.cgs.datasets.motifs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;

public abstract class BackgroundModel extends BackgroundModelMetadata {

  public static final double EPSILON = 1E-6;
  
  public static final char[] BASE_ORDER = new char[] {'A', 'C', 'G', 'T'};

  public static final int DEFAULT_MAX_KMER_LEN = 3;

  //a reference to a genome object in addition to the metadata genomeID
  protected Genome gen;

  /**
   * For keeping track of whether or not both strands are accounted for in 
   * the model.
   * True if the model was based on a single strand
   * False if the model was based on both strands
   */
  protected Boolean isStranded = null;

  /**
   * Map for holding the kmer probability values for the model. Each 
   * element of the array holds kmers whose length is the index of that element.
   * i.e. model[2] holds 2-mers. Accordingly, model[0] should always be empty.
   * 
   * Note: this should only be used for probabilties, and not for counts. The
   * CountBasedModel uses a separate Map to hold the counts and this map to
   * keep the normalized counts (frequencies).
   */
  protected Map<String, Double>[] modelProbs;


  
  /**
   * 
   */
  public BackgroundModel(String name, Genome gen) {
    this(name ,gen, DEFAULT_MAX_KMER_LEN);
    this.init();
  }


  /**
   * 
   * @param maxKmerLen
   */
  public BackgroundModel(String name, Genome gen, int maxKmerLen) {
    super(name, maxKmerLen, gen.getDBID());
  	this.gen = gen;
    modelProbs = new HashMap[maxKmerLen + 1];
    for (int i = 1; i <= maxKmerLen; i++) {
      modelProbs[i] = new HashMap<String, Double>();
    }
    this.init();
  }
  

  /**
   * Construct a Background Model from an existing Background Model. Set the
   * instance variables of this class to match the source model. Subclasses
   * should only set the dbid if the type of the subclass model matches the type
   * of the source model.
   *  
   * @param source the Background Model on which to base the one being constructed
   */
  public BackgroundModel(BackgroundModel source) {
    //FIXME
  	this(source.name, source.gen, source.getMaxKmerLen());
  	this.isStranded = source.isStranded;
  	this.init();
  }
  
  //TODO
  protected abstract void init();
  
  
  /**
   * Returns this model's genome
   * @return the genome, which is null if it's not set
   */
  public Genome getGenome() {
  	return gen;
  }
  
  
  /**
   * Set this model to have the specified genome
   * @param gen
   */
  public void setGenome(Genome gen) {
  	this.gen = gen;
  }
    

  /**
   * Get the Markov-order of this model
   * @return
   */
  public int getMarkovOrder() {
    return modelProbs.length - 2;
  }
  
  
  public boolean isStranded() {
  	if (isStranded == null) {
  		this.checkAndSetIsStranded();  		
  	}
  	return isStranded;
  }
  
  
  /**
   * Check whether this model is making use of both strands or just a single 
   * strand (i.e. check if all probabilities/counts match their reverse
   * complement
   * @return
   */
  public abstract boolean checkAndSetIsStranded();
    
  
  /**
   * Return the set of kmers of the specified length for which this model has
   * probabilities or counts
   * @param kmerLen
   * @return
   */
  public abstract Set<String> getKmers(int kmerLen);

  
  /**
   * Return the markov probability for the last base of the specified kmer
   * conditioned upon the preceding bases. 
   * @param kmer
   * @return
   */
  public abstract double getMarkovProb(String kmer);
  
  
  /**
   * Return the markov probability for the last base of the kmer represented by
   * the specified int, conditioned upon the preceding bases.
   * @param intVal
   * @param kmerLen
   * @return
   */
  public double getMarkovProb(int intVal, int kmerLen) {
    return (this.getMarkovProb(BackgroundModel.int2seq(intVal, kmerLen)));
  }
  
    
  /**
   * Convert a base to an int value
   * @param base
   * @return
   */
	public static int base2int(char base) {
		int intVal = -1;
		switch (base) {
		case 'A':
			intVal = 0;
			break;
		case 'C':
			intVal = 1;
			break;
		case 'G':
			intVal = 2;
			break;
		case 'T':
			intVal = 3;
			break;
		default:
			throw new IllegalArgumentException("Invalid character: " + base);
		}
		return intVal;
	}


	/**
	 * Return a base for the specified integer
	 * @param x
	 * @return
	 */
	public static char int2base(int x) {
		char base;
		switch (x) {
		case 0:
			base = 'A';
			break;
		case 1:
			base = 'C';
			break;
		case 2:
			base = 'G';
			break;
		case 3:
			base = 'T';
			break;
		default:
			throw new IllegalArgumentException("Invalid int: " + x);
		}
		return (base);
	}


	/**
	 * Convert a nucleotide sequence to an integer value
	 * @param seq
	 * @return
	 */
	public static int seq2int(String seq) {
		int intVal = 0;
		int len = seq.length();

		for (int i = 0; i < len; i++) {
			long currInt = BackgroundModel.base2int(seq.charAt(i));
			if (currInt == -1) {
				return -1;
			}
			intVal = intVal << 2;
			intVal += currInt;
		}
		return intVal;
	}


	/**
	 * 
	 * @param x 
	 * @return
	 */
	public static String int2seq(long x, int kmerLen) {	
		/**
		 * check that the x is valid for the specified maxKmerLen. 
		 * Note: 4 << (2 * (kmerLen - 1)) = 4^kmerLen
		 */
	  if (x > ((4 << (2 * (kmerLen - 1))) - 1)) {
		  throw new IllegalArgumentException("Invalid int value, " + x + ", for kmerLen " + kmerLen);
		}
    StringBuffer seq = new StringBuffer(kmerLen);
    for (int i = 0; i < kmerLen; i++) {
      int baseVal = (int)(x % 4);
      seq.append(BackgroundModel.int2base(baseVal));
      x = x >> 2;             
    }
		return seq.reverse().toString();
	}


  /**
   * Assemble a list of distinct pairs of kmers and non-palindrome reverse complements.
   * @param kmerLen
   * @return
   */
  protected static List<Pair<Integer, Integer>> computeDistinctRevCompPairs(int kmerLen) {
    int numKmers = (int)Math.pow(4, kmerLen);
    int[] revCompMap = new int[numKmers];
    Arrays.fill(revCompMap, -1);
    ArrayList<Pair<Integer, Integer>> revCompPairs = new ArrayList<Pair<Integer, Integer>>();
    
    for (int i = 0; i < revCompMap.length; i++) {
      if (revCompMap[i] == -1) {
        int revComp = BackgroundModel.seq2int(SequenceUtils.reverseComplement(BackgroundModel.int2seq(i, kmerLen)));
        revCompMap[i] = revComp;
        if (i != revComp) {
          revCompMap[revComp] = i;
          revCompPairs.add(new Pair<Integer, Integer>(i, revComp));
        }        
      }
    }
    
    return revCompPairs;
  }
  
  
  /**
   ***********************************************************************
   */
	
	public static void main(String[] args) {
	  for (int i = 1; i <= 3; i++) {
	    for (int j = 0; j < Math.pow(4,i); j++) {
	      String seq = BackgroundModel.int2seq(j, i);
	      System.out.println(j + "\t" + seq + "\t" + BackgroundModel.seq2int(seq)); 
	      assert (j == BackgroundModel.seq2int(seq));
	    }	    
	  }
	}
}


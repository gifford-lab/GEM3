package edu.mit.csail.cgs.datasets.motifs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;

public abstract class BackgroundModel {

  public static final double EPSILON = 1E-6;
  
  public static final char[] BASE_ORDER = new char[] {'A', 'C', 'G', 'T'};

  public static final int DEFAULT_MODEL_LENGTH = 3;

  public String name;
  public int dbid;
  public boolean hasdbid;
  
  protected Genome gen;

  /**
   * Map for holding the kmer probability (and count) values for the model. Each 
   * element of the array holds kmers whose length is the index of that element.
   * i.e. model[2] holds 2-mers. Accordingly, model[0] should always be empty. 
   */
  protected Map<String, Pair<Double, Integer>>[] model;


  /**
   * For keeping track of whether or not both strands are accounted for in 
   * the model.
   * True if the model was based on a single strand
   * False if the model was based on both strands
   */
  protected boolean isStranded = true;
  
  /**
   * 
   */
  public BackgroundModel() {
    this(DEFAULT_MODEL_LENGTH);    
  }


  /**
   * 
   * @param modelLength
   */
  public BackgroundModel(int modelLength) {
    model = new HashMap[DEFAULT_MODEL_LENGTH + 1];
    for (int i = 0; i <= DEFAULT_MODEL_LENGTH; i++) {
      model[i] = new HashMap<String, Pair<Double, Integer>>();
    }
  }
  
  
  /**
   * Return the order of this model.
   * @return
   */
  public int getMaxKmerLen() {
    return model.length;
  }
  
  
  /**
   * 
   * @param mer
   * @return
   */
  public Pair<Double, Integer> getModelValuePair(String mer) {
    if (mer.length() > 0) {      
      return (model[mer.length()].get(mer));
    }
    else {
      throw new IllegalArgumentException("Zero length kmer.");
    }
  }


  /**
   * 
   * @param kmerLen
   * @param intVal
   * @return
   */
  public Pair<Double, Integer> getModelValuePair(int kmerLen, int intVal) {
    if (kmerLen > 0) {
      return (this.getModelValuePair(BackgroundModel.int2seq(intVal, kmerLen)));
    }
    else {
      throw new IllegalArgumentException("kmerLen must be greater than zero.");
    }
  }
  
  
  /**
   * 
   * @param mer
   * @param val
   */
  public void setModelValuePair(String mer, Double prob, Integer count) {
    if (mer.length() > 0) {     
      Pair<Double, Integer> values = new Pair<Double, Integer>(prob, count);
      model[mer.length()].put(mer, values);
    }
    else {     
      throw new IllegalArgumentException("Zero length kmer.");      
    }
  }
  
  
  /**
   * 
   * @param mer
   * @return
   */
  public Integer getModelCount(String mer) {
    Pair<Double, Integer> values = this.getModelValuePair(mer);
    if (values != null) {
      return values.cdr(); 
    }
    else {
      return null;
    }
  }


  /**
   * 
   * @param kmerLen
   * @param intVal
   * @return
   */
  public Integer getModelCount(int kmerLen, int intVal) {
    Pair<Double, Integer> values = this.getModelValuePair(kmerLen, intVal);
    if (values != null) {
      return values.cdr(); 
    }
    else {
      return null;
    }
  }


  /**
   * 
   * @param mer
   * @param val
   */
  public void setModelCount(String mer, int count) {
    this.setModelValuePair(mer, this.getModelProb(mer), count);
  }
  
  
  /**
   * 
   * @param mer
   * @return
   */
  public Double getModelProb(String mer) {
    Pair<Double, Integer> values = this.getModelValuePair(mer);
    if (values != null) {
      return values.car(); 
    }
    else {
      return null;
    }
  }


  /**
   * 
   * @param kmerLen
   * @param intVal
   * @return
   */
  public Double getModelProb(int kmerLen, int intVal) {
    Pair<Double, Integer> values = this.getModelValuePair(kmerLen, intVal);
    if (values != null) {
      return values.car(); 
    }
    else {
      return null;
    }
  }


  /**
   * 
   * @param mer
   * @param val
   */
  public void setModelProb(String mer, double prob) {
    this.setModelValuePair(mer, prob, this.getModelCount(mer));
  }
    
  
  /**
   * Check whether this model is making use of both strands or just a single 
   * strand (i.e. check if all probabilities/counts match their reverse
   * complement
   * @return
   */
  public abstract boolean checkAndSetIsStranded();
  
  
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
		 * check that the x is valid for the specified kmerlen. 
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
  protected static List<Pair<Integer, Integer>> computeRevCompPairs(int kmerLen) {
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


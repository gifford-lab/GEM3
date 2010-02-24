package edu.mit.csail.cgs.datasets.motifs;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;

public abstract class BackgroundModel {

  public static final char[] BASE_ORDER = new char[] {'A', 'C', 'G', 'T'};

  public static final int DEFAULT_MODEL_LENGTH = 3;

  public String name;
  public int dbid;
  public boolean hasdbid;
  
  protected Genome gen;

  /**
   * Map for holding the kmer count/probability values for the model. Each 
   * element of the array holds kmers whose length is the index of that element.
   * i.e. model[2] holds 2-mers. So, model[0] should always be empty. 
   */
  protected Map<String, Double>[] model;

  /**
   * This map is meant to store the counts from which a non-count model was
   * derived. It may be null. For a count based model this should be ignored and
   * the methods for accessing and modifying counts should be overridden so that
   * only the model map is used. 
   */
  protected Map<String, Double>[] counts = null;
  
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
      model[i] = new HashMap<String, Double>();
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
  public Double getModelCount(String mer) {
    if (mer.length() > 0) {
      return (counts[mer.length()].get(mer));
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
  public Double getModelCount(int kmerLen, int intVal) {
    if (kmerLen > 0) {
      return (counts[kmerLen].get(BackgroundModel.intToSeq(intVal, kmerLen)));
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
  public void setModelCount(String mer, double val) {
    if (mer.length() <= this.getMaxKmerLen() && mer.length() > 0) {
      counts[mer.length()].put(mer, val);
    }
    else if (mer.length() < 1) {     
      throw new IllegalArgumentException("Zero length kmer.");      
    }
    else {
      throw new IllegalArgumentException("Kmer " + mer + " must have length less than model length (" + this.getMaxKmerLen() + ").");
    }
  }
  
  
  
  /**
   * 
   * @param mer
   * @return
   */
  public Double getModelVal(String mer) {
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
  public Double getModelVal(int kmerLen, int intVal) {
    if (kmerLen > 0) {
      return (model[kmerLen].get(BackgroundModel.intToSeq(intVal, kmerLen)));
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
  public void setModelVal(String mer, double val) {
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
	public static int seqToInt(String seq) {
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
	public static String intToSeq(long x, int kmerLen) {	
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
   ***********************************************************************
   */
	
	public static void main(String[] args) {
	  for (int i = 1; i <= 3; i++) {
	    for (int j = 0; j < Math.pow(4,i); j++) {
	      String seq = BackgroundModel.intToSeq(j, i);
	      System.out.println(j + "\t" + seq + "\t" + BackgroundModel.seqToInt(seq)); 
	      assert (j == BackgroundModel.seqToInt(seq));
	    }	    
	  }
	}
}


/**
 * 
 */
package edu.mit.csail.cgs.deepseq.discovery.hmm;

import java.util.*;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.DoubleMatrix3D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix3D;
import cern.jet.math.Arithmetic;
import cern.jet.math.Functions;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.deepseq.BindingModel;
import edu.mit.csail.cgs.deepseq.ReadHit;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.utils.numeric.Numerical;

/**
 * @author rca
 *
 */
public class HMMAnalysisRegion {

  //the genomic region which this analysis region spans
  private Region genomicRegion;
  
  //the sequence for the '+' strand of that region
  private String sequence;
  public int[] sequenceAsInts;
  
  //the IP channel reads for the region
  private ArrayList<ReadHit> signalReads;
  
  //boolean to keep track of whether a set of WCE extract reads for this region exist
  boolean hasWCE;
  
  //the control channel reads for the region
  private ArrayList<ReadHit> controlReads;
  
  //the binding frequency values calculated from the controlReads
  double[] wceBindingFreq;
  
  //the binding model to use for determining read probabilities
  private BindingModel bindingModel;

  //some hashes for keeping the reads organized by assignment
  private HashMap<ReadHit, Integer> readAssignments;
  private HashMap<Integer, List<ReadHit>> readsByLoc;
  private TreeMap<Integer, Integer> numReadsByLoc;
  
  
  private HashMap<ReadHit, Double> logReadProbCache;
  private HashMap<Integer, Double> locLogReadProbCache;
  
  //expectations
  private DoubleMatrix2D logAlpha;
  private DoubleMatrix2D logBeta;
  private DoubleMatrix3D logXiNumerator;

  private DoubleMatrix2D logGamma;
  private DoubleMatrix3D logXi;
  
  //ML state sequence matrix
  private DoubleMatrix2D viterbiLogProb;
  private int[][] viterbiTraceback;
  
  //log likelihood for this region only, not for the whole data set
  private double logLikelihood;
  
  public HMMAnalysisRegion(Region genomicRegion, List<ReadHit> signalReads, Double minFreq, BindingModel model) {
    this(genomicRegion, signalReads, null, null, null, model);
  }
  
  
  public HMMAnalysisRegion(Region genomicRegion, List<ReadHit> signalReads, List<ReadHit> controlReads, Double wceScaleFactor, Double minFreq, BindingModel model) {
    this.genomicRegion = genomicRegion;
    this.signalReads = new ArrayList<ReadHit>(signalReads);
    this.bindingModel = model;

    if (controlReads == null) {
      hasWCE = false;
      wceBindingFreq = null;
    }
    else {
      hasWCE = true;
      this.controlReads = new ArrayList<ReadHit>(controlReads);
      this.computeWCEBindingFreqs(controlReads, wceScaleFactor.doubleValue(), minFreq.doubleValue());
    }
    
    SequenceGenerator<Region> seqGen = new SequenceGenerator<Region>();
    sequence = seqGen.execute(genomicRegion).toUpperCase();
    sequenceAsInts = this.convertSequenceToInts(sequence);
    
    readAssignments = new HashMap<ReadHit, Integer>(signalReads.size());
    readsByLoc = new HashMap<Integer, List<ReadHit>>();
    numReadsByLoc = new TreeMap<Integer, Integer>();
    
    logReadProbCache = new HashMap<ReadHit, Double>(signalReads.size());
    locLogReadProbCache = new HashMap<Integer, Double>();
  }
  

  /**
   * 
   * @return
   */
  public Region getGenomicRegion() {
    return genomicRegion;
  }


  /**
   * Accessor for the width of the analysis region
   * @return
   */
  public int getWidth() {
    return genomicRegion.getWidth();
  }
  
  
  /**
   * Accessor for the genomic sequence of the region
   * @return
   */
  public String getSequence() {
    return sequence;
  }


  /**
   * Accessor for the list of reads from the signal channel
   * @return
   */
  public List<ReadHit> getSignalReads() {
    return signalReads;
  }
  
  
  /**
   * Gets a character at the specified 0-based offset from the start of the region
   * @param pos
   * @return
   */
  public char baseAt(int pos) {
    return sequence.charAt(pos);
  }
  
  
  /**
   * 
   * @param pos
   * @return
   */
  public double getWCEBindingFreq(int pos) {
    if (hasWCE) {
      return wceBindingFreq[pos];
    }
    else {
      throw new NullPointerException("No WCE Binding Freq is available!");
    }
  }
  
  
  /**
   * getter for alpha
   * @return
   */
  public DoubleMatrix2D getLogAlpha() {
    return logAlpha;
  }


  /**
   * getter for beta
   * @return
   */
  public DoubleMatrix2D getLogBeta() {
    return logBeta;
  }


  /**
   * getter for xi numerator
   * @return
   */
  public DoubleMatrix3D getLogXiNumerator() {
    return logXiNumerator;
  }
  
  
  /**
   * getter for gamma
   * @return
   */
  public DoubleMatrix2D getLogGamma() {
    return logGamma;
  }

  
  /**
   * getter for xi
   * @return
   */
  public DoubleMatrix3D getLogXi() {
    return logXi;
  }

  
  /**
   * getter for logLikelihood
   * @return
   */
  public double getLogLikelihood() {
    return logLikelihood;
  }
  
  
  public DoubleMatrix2D getViterbiLogProb() {
    return viterbiLogProb;
  }


  public int[][] getViterbiTraceback() {
    return viterbiTraceback;
  }
  
  
  /**
   * setter for alpha
   * @param logAlpha
   */
  public void setLogAlphaBeta(DoubleMatrix2D logAlpha, DoubleMatrix2D logBeta, double logLikelihood) {
    this.logAlpha = logAlpha;
    this.logBeta = logBeta;
    this.logLikelihood = logLikelihood;
  }

  
  /**
   * setter for the numerator of formula for xi
   * @param logXiNumerator
   */
  public void setLogXiNumerator(DoubleMatrix3D logXiNumerator) {
    this.logXiNumerator = logXiNumerator;
  }


  /**
   * setter for gamma
   * @param logGamma
   */
  public void setLogGamma(DoubleMatrix2D logGamma) {
    this.logGamma = logGamma;
  }


  /**
   * setter for xi
   * @param logXi
   */
  public void setLogXi(DoubleMatrix3D logXi) {
    this.logXi = logXi;
  }
  
  
  public void setViterbiLogProb(DoubleMatrix2D viterbiLogProb) {
    this.viterbiLogProb = viterbiLogProb;
  }

  
  public void setViterbiTraceback(int[][] viterbiTraceback) {
    this.viterbiTraceback = viterbiTraceback;
  }

  /**
   * Get a list of reads assigned to the specified location
   * @param loc
   * @return
   */
  public List<ReadHit> readsAssignedTo(int loc) {
    List<ReadHit> reads = readsByLoc.get(loc);
    if (reads == null) {
      reads = new ArrayList<ReadHit>(); 
    }
    return reads;
  }
  
  
  /**
   * Gets the number of reads assigned to the specified location
   * @param loc
   * @return
   */
  public int numReadsAssignedTo(int loc) {
    Integer numReads = numReadsByLoc.get(loc);
    if (numReads != null) {
    	return numReads.intValue();
    }
    else {
    	return 0;
    }
  }
  
  
  public TreeMap<Integer, Integer> getNumReadsByLoc() {
    return numReadsByLoc;
  }


  /**
   * Updates the binding event location to which a hit is assigned
   * @param hit
   * @param newLoc location where hit is being assigned, as a 0-based offset from
   * the start of this region
   */
  public void assignRead(ReadHit hit, int newLoc) {
    //check if this hit has been assigned before
    if (readAssignments.containsKey(hit)) {
      //check that the new assignment is different
      int oldLoc = readAssignments.get(hit);
      if (oldLoc != newLoc) {
        readAssignments.put(hit, newLoc);
        
        //clear the cache
        logReadProbCache.put(hit, null);
        locLogReadProbCache.put(oldLoc, null);
        
        //update the hash of reads by location
        List<ReadHit> oldLocHits = readsByLoc.get(oldLoc);
        if (oldLocHits != null) {
          oldLocHits.remove(hit);
          numReadsByLoc.put(oldLoc, oldLocHits.size());
        }
        List<ReadHit> newLocHits = readsByLoc.get(newLoc);
        if (newLocHits == null) {
         newLocHits = new ArrayList<ReadHit>();
         readsByLoc.put(newLoc, newLocHits);
        }
        newLocHits.add(hit);
        numReadsByLoc.put(newLoc, newLocHits.size());
      }
    }
    else {
      //the hit hasn't been assigned before, so add it
      readAssignments.put(hit, newLoc);
      List<ReadHit> newLocHits = readsByLoc.get(newLoc);
      if (newLocHits == null) {
       newLocHits = new ArrayList<ReadHit>();
       readsByLoc.put(newLoc, newLocHits);
      }
      newLocHits.add(hit);
      numReadsByLoc.put(newLoc, newLocHits.size());
    }
  }
  
  
  /**
   * Returns the location to which the specified hit is assigned, or -1 if
   * unassigned
   * @param hit
   * @return
   */
  public int getReadAssignment(ReadHit hit) {
    Integer assignedLoc = readAssignments.get(hit);
    if (assignedLoc != null) {
      return assignedLoc.intValue();     
    }
    else {
      return -1;
    }
  }
  
  
  /**
   * compute the probability of all of the hits assigned to the specified location
   * coming from a binding event at that location - does _not_ include poisson 
   * frequency. Use caching.
   * @param loc
   * @return
   */
  public double computeLogReadProb(int loc) {
    Double cachedLogProb = locLogReadProbCache.get(loc);
    if (cachedLogProb != null) {
      return cachedLogProb.doubleValue();
    }
    else {
      List<ReadHit> hits = readsByLoc.get(loc);
      double logProb = this.computeLogReadProb(hits, loc);
      locLogReadProbCache.put(loc, logProb);
      return logProb;
    }
  }

  
  /**
   * compute the probability of all of the specified hits coming from a binding
   * event at the specified location - does _not_ include the poisson probability
   * @param hits
   * @param loc
   * @return
   */
  public double computeLogReadProb(List<ReadHit> hits, int loc) {
    double totalProb = 0;
    for (ReadHit hit : hits) {
      totalProb = totalProb + this.computeLogReadProb(hit, loc);      
    }
    return totalProb;
  }
  
  
  /**
   * Compute the log probability of the specified location producing the 
   * specified hit. Cache results.
   * @param hit
   * @param loc
   * @return 
   */
  public double computeLogReadProb(ReadHit hit, int loc) {
    Double cachedLogProb = logReadProbCache.get(hit);
    Integer assignedLoc = readAssignments.get(hit);
    if ((cachedLogProb != null) && (assignedLoc != null) && (assignedLoc.intValue() == loc)) {      
        return cachedLogProb.doubleValue();
    } 
    else {
      int distance = genomicRegion.getStart() + loc;
      if (hit.getStrand() == '+') {
        distance = distance - hit.getFivePrime();
      }
      else {
        distance = hit.getFivePrime() - distance;
      }
      double logProb = bindingModel.logProbability(distance);

      //if this is for the assigned location then cache the value
      if ((assignedLoc != null) && (assignedLoc.intValue() == loc)) {
        logReadProbCache.put(hit, logProb);
      }
      return logProb;
    }
  }  
  
  
  /**
   * Compute the poisson binding frequencies for the WCE. 
   * 
   * This implementation computes the binding freq. as the count of reads within
   * the window of the binding model around the location. 
   * @param controlReads
   */
  private void computeWCEBindingFreqs(List<ReadHit> controlReads, double scaleFactor, double minFreq) {
    int regStart = genomicRegion.getStart();
    int regionWidth = genomicRegion.getWidth();
    int modelMin = bindingModel.getMin();
    int modelMax = bindingModel.getMax();
    
    wceBindingFreq = new double[regionWidth];
    
    for (ReadHit currentHit : controlReads) {
      //calculate the range of locations the read could be assigned
      int mStart, mEnd;
      if (currentHit.getStrand() == '+') {
        mStart = currentHit.getFivePrime() - regStart + modelMin;
        mEnd = currentHit.getFivePrime() -regStart + modelMax;
      }
      else {
        //subtract the min and max in order to keep mStart as the lesser value
        mStart = currentHit.getFivePrime() - regStart - modelMax;
        mEnd = currentHit.getFivePrime() - regStart - modelMin;
      }
      //ensure the values are valid
      mStart = Math.max(0, mStart);
      mEnd = Math.min(regionWidth - 1, mEnd);
      
      //at each position add the probability from the binding model to the 
      //count for that position
      int m = mStart;
      int distance = regStart + m;
      int delta;
      if (currentHit.getStrand() == '+') {
        distance = distance - currentHit.getFivePrime();
        delta = 1;
      }
      else {
        distance = currentHit.getFivePrime() - distance;
        delta = -1;
      }
      while (m <= mEnd) {
        double prob = bindingModel.probability(distance);
        wceBindingFreq[m] = wceBindingFreq[m] + prob;
        distance += delta;
        m++;
      }
    }
    for (int m = 0; m < regionWidth; m++) {
      wceBindingFreq[m] = Math.max(wceBindingFreq[m], minFreq);
      wceBindingFreq[m] = wceBindingFreq[m] * scaleFactor;
    }
  }
  
  
  public void printNumReadsByLoc() {
    StringBuffer buf = new StringBuffer();
    for (int loc : numReadsByLoc.keySet()) {
      buf.append(loc + ": " + numReadsByLoc.get(loc) + "\t");
    }
    System.out.println(buf.toString());
  }
  
  
  private int[] convertSequenceToInts(String sequence) {
  	int[] intSeq = new int[sequence.length()]; 
  	for (int i = 0; i < sequence.length(); i++) {
  		intSeq[i] = BindingMotifHMM.baseToIndex(sequence.charAt(i));
  	}
  	return intSeq;
  }
}

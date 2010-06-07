/**
 * 
 */
package edu.mit.csail.cgs.deepseq.discovery.hmm;

import java.util.*;

import org.apache.log4j.Logger;

import cern.colt.matrix.*;
import cern.colt.matrix.impl.*;
import cern.jet.math.Functions;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.deepseq.BindingModel;
import edu.mit.csail.cgs.deepseq.DeepSeqExpt;
import edu.mit.csail.cgs.deepseq.ReadHit;
import edu.mit.csail.cgs.deepseq.discovery.SingleConditionFeatureFinder;
import edu.mit.csail.cgs.deepseq.features.EnrichedFeature;
import edu.mit.csail.cgs.deepseq.features.Feature;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.numeric.Numerical;
import edu.mit.csail.cgs.utils.stats.StatUtil;

/**
 * @author rca
 * 
 */
public class BindingMotifHMM extends SingleConditionFeatureFinder {

  private static Logger logger = Logger.getLogger(BindingMotifHMM.class);

  private static final String UPDATE_MOTIF_KEY = "update_motif";
  private static final String BG_BINDING_FREQ_KEY = "bg_binding_freq";
  private static final String MOTIF_BINDING_FREQ_KEY = "motif_binding_freq";
  private static final String MOTIF_INIT_STATE_PROB_KEY = "motif_init_state_prob";
  private static final String MOTIF_END_START_TRANS_PROB_KEY = "motif_end_start_prob";
  private static final String MAX_ITER_KEY = "max_iter";
  private static final String CONV_THRESH_KEY = "conv_thresh";
  
  private static final String UPDATE_MOTIF_DFLT = "TRUE";
  //private static final double BG_BINDING_FREQ_DFLT = 1.0;
  private static final double MOTIF_BINDING_FREQ_DFLT = 45.0; //FIXME make a CL param
  private static final double MOTIF_INIT_STATE_PROB_DFLT = 0.0;
  private static final double MOTIF_END_START_TRANS_PROB_DFLT = 1.0/10000.0;
  private static final int MAX_ITER_DFLT = 10000;
  private static final int MIN_ITER_DFLT = 5000;
  private static final double CONV_THRESH_DFLT = 1.0e-5; //default is run to max iter

  //pseudocount for sequence emission probabilities
  private static final double PSEUDOCOUNT = 0.01;
  
//  //to filter out towers and needles.
//  protected final int MAX_DUPLICATES = 3; 
//
//  //number of reads that need to be in a region for it to possibly be a tower 
//  protected final int TOWER_THRESHOLD = 1000; 

  /*************************************************************************
   * Data
   *************************************************************************/
  List<HMMAnalysisRegion> analysisRegions;
  int regionSizeTotal = 0;

  
  /*************************************************************************
   * HMM Parameters
   *************************************************************************/
  /**
   * Conventions for HMM state numbering:
   * The HMM States will be indexed 0 to K-1. 
   * The background states will be 0 to numBG-1.
   * 
   * For each motif the states will be in their sequential order followed by
   * the reverse complement states in reverse order 
   * (e.g. m1, m2, m3, \bar{m3}, \bar{m2}, \bar{m1})
   * so looking just at the '+' strand the states would be iterated over in
   * ascending order whether the occurence is the motif or its reverse complement
   */
  int hmm_numStates; //aka "K"
  
  //the number of states that would exist if this were a regular HMM not a GHMM
  int hmm_numEffectiveStates; 
  
  int hmm_numBGStates; //the number of background states, starting at 0
  int hmm_numMotifs;
  int[] hmm_motifFwdStates;
  int[] hmm_motifLengths;
  boolean[] isBGState;
  boolean[] isMotifState;
  boolean[] isMotifFwdState;
  boolean[] isStateFixedDuration;
  int[] stateDurations;
  int maxStateDuration;
  
  Hashtable<Integer, Integer> stateToMotifMap;
  Hashtable<HMMAnalysisRegion, int[][]> regionToMKReadCounts = new Hashtable<HMMAnalysisRegion, int[][]>();
   
    
  //initial state log-probabilities
  double motifInitialStateProb;
  DoubleMatrix1D hmm_LogPi;
  
  //state transition probabilities
  double motifEndStartInitialTransProb;
  DoubleMatrix2D hmm_LogA;
  
  double[][] hmm_transitionPriorParams; //alpha params for dirichlet priors 
  double hmm_transitionPriorWeight = 1.0;
  
  TreeMap<Integer, DoubleMatrix2D> hmm_LogAdjustedTransProb;
  
  //sequence emission probabilities
  DenseDoubleMatrix2D[] hmm_phi;
  
  
  //state binding frequencies
  double minBGFreq;
  double minSignalFreq;
  double bgBindingFreq;  
  double motifBindingFreq;
  DoubleMatrix1D hmm_binding;
  boolean hmm_useWCEforBGBinding;
  
  /**
   * this scale factor relates the number of reads that will be at a particular
   * spot in the signal channel to the number of reads that are present in the
   * same spot in the WCE channel. This depends on the ratio of the number of
   * reads in each channel, and also the percentage of signal channel reads that
   * are "noise".
   */
  double hmm_WCEBindingFreqScaleFactor;
  
  //whether or not to update the motif PWM
  boolean hmm_updateMotif;
  
  //whether or not to update read assignments stochastically
  boolean isReadUpdateStochastic = true;
  //TODO make this parameter adjustable from the CL
  
  /**
   * Binding model representing the empirical distribution of a binding event
   */
  private BindingModel model;

  private ArrayList<Feature> insignificantFeatures;

  private int maxIterations;
  private double convergenceThreshold;
  
  /**
   * Log-likelihood from current iteration
   */
  private double totalLogLikelihood;
  private double prevTotalLogLikelihood;
   
  
  private Random randomSrc = new Random(137);
  
  /**
   * 
   * @param signal
   * @param control
   * @param selectedRegions
   * @param args
   */
  public BindingMotifHMM(DeepSeqExpt signal, DeepSeqExpt control, Double wceScaleFactor, List<Region> selectedRegions, BindingModel model, String[] args) {
    super(signal, control);
    String updateMotif = Args.parseString(args, UPDATE_MOTIF_KEY, UPDATE_MOTIF_DFLT).toUpperCase();
    if (updateMotif.equals("FALSE")) {
      hmm_updateMotif = false;
    }
    else {
      hmm_updateMotif = true;
    }
    
    /**
     * Compute a minimum read frequency because we can't have 0.
     * Use the number of reads in the control divided by the length of the genome
     */    
    minSignalFreq = (double)signal.getHitCount() / (double)signal.getGenome().getGenomeLength();
    Double minFreqDbl = null;
    if (control != null) {
      minBGFreq = (double)control.getHitCount() / (double)control.getGenome().getGenomeLength();
      minFreqDbl = new Double(minBGFreq);
    }
    else {
      minBGFreq = minSignalFreq;
    }
    
    
    bgBindingFreq = Args.parseDouble(args, BG_BINDING_FREQ_KEY, minBGFreq);
    motifBindingFreq = Args.parseDouble(args, MOTIF_BINDING_FREQ_KEY, MOTIF_BINDING_FREQ_DFLT);
    
    //In the GHMM version of this model we force the model to start and end in the BG state
    motifInitialStateProb = 0;//Args.parseDouble(args, MOTIF_INIT_STATE_PROB_KEY, MOTIF_INIT_STATE_PROB_DFLT);
    motifEndStartInitialTransProb = Args.parseDouble(args, MOTIF_END_START_TRANS_PROB_KEY, MOTIF_END_START_TRANS_PROB_DFLT);
    maxIterations = Args.parseInteger(args, MAX_ITER_KEY, MAX_ITER_DFLT);
    convergenceThreshold = Args.parseDouble(args, CONV_THRESH_KEY, CONV_THRESH_DFLT);
        
    this.model = model; 
    
    
    
    analysisRegions = new ArrayList<HMMAnalysisRegion>(selectedRegions.size());
    for (Region region : selectedRegions) {
      List<ReadHit> signalReads = signal.loadHits(region);
      List<ReadHit> controlReads = null;
      if (control != null) {
        controlReads = control.loadHits(region);
      }
      
      HMMAnalysisRegion hmmRegion = new HMMAnalysisRegion(region, signalReads, controlReads, wceScaleFactor, minFreqDbl, model);
      analysisRegions.add(hmmRegion);
      regionSizeTotal += region.getWidth();
    }
    logger.debug(selectedRegions.size() + " regions loaded for analysis.");

    if (control != null) {
      hmm_useWCEforBGBinding = true;
      hmm_WCEBindingFreqScaleFactor = wceScaleFactor.doubleValue();
    }
    
    logger.debug("Binding Mixture Constructed.");
  }
  
  
  /**
   * Constructor for running without wce channel
   * @param signal
   * @param selectedRegions
   * @param args
   */
  public BindingMotifHMM(DeepSeqExpt signal, List<Region> selectedRegions, BindingModel model, String[] args) {
    this(signal, null, null, selectedRegions, model, args);
  }
  
  
  public void setHMMBinding(DoubleMatrix1D hmm_binding) {
    assert (hmm_binding.size() == hmm_numStates) : 
      "Matrix of binding values is length " + hmm_binding.size() + " instead of " + hmm_numStates;  
    this.hmm_binding = hmm_binding;
  }

  
  public void setHMMInitialStateLogProb(DoubleMatrix1D hmm_LogPi) {
    assert (hmm_LogPi.size() == hmm_numStates) :
      "Matrix of initial probabilities is length " + hmm_LogPi.size() + " instead of " + hmm_numStates;
    this.hmm_LogPi = hmm_LogPi;
  }
  
  
  public void setHMMStateTransitionLogProb(DoubleMatrix2D hmm_LogA) {
    assert ((hmm_LogA.rows() == hmm_numStates) && (hmm_LogA.columns() == hmm_numStates)):
      "Matrix of transition probabilities is " + hmm_LogA.rows() + "x" + hmm_LogA.columns() + " instead of " + hmm_numStates + "x" + hmm_numStates;
    this.hmm_LogA = hmm_LogA;
  }

  
  public void setMaxIterations(int maxIterations) {
    this.maxIterations = maxIterations;
  }
  
  
  public void setConvergenceThreshold(double convThresh) {
    this.convergenceThreshold = convThresh;
  }
  
  public void initializeHMM(List<WeightMatrix> motifPWMs, WeightMatrix bgPWM) {
    hmm_numBGStates = 1;
    hmm_numMotifs = motifPWMs.size();
    hmm_motifLengths = new int[hmm_numMotifs];
    for (int i = 0; i < motifPWMs.size(); i++) {
      hmm_motifLengths[i] = motifPWMs.get(i).length();
    }

    //set up the hmm state structure    
    this.createHMMStates();
      
    //state binding frequencies
    this.initializeStateBindingFrequencies();

    /**
     * initialize initial state probabilities and 
     * state transition probabilities
     * Note: binding frequencies must be initialized first!
     */
    this.initializeInitialStateProbabilities();
    this.initializeTransitionProbabilities();        
    
    //initialize sequence emission probabilities
    List<WeightMatrix> bgPWMs = new ArrayList<WeightMatrix>();
    bgPWMs.add(bgPWM);
    this.initializeEmissionProbabilities(motifPWMs, bgPWMs);
    
    //initialize read assignments
    this.initializeReadAssignments();
  }
  
  
  public void initializeHMM(List<WeightMatrix> motifPWMs, List<WeightMatrix> bgPWMs,
      DoubleMatrix1D hmm_binding, DoubleMatrix1D hmm_LogPi, 
      DoubleMatrix2D hmm_LogA) {
    
    hmm_numBGStates = bgPWMs.size();
    hmm_numMotifs = motifPWMs.size();
    hmm_motifLengths = new int[hmm_numMotifs];
    for (int i = 0; i < motifPWMs.size(); i++) {
      hmm_motifLengths[i] = motifPWMs.get(i).length();
    }
    
    //set up the hmm state structure    
    this.createHMMStates();
    
    //set the binding frequencies and transition probabilities as specified
    this.hmm_binding = hmm_binding;
    this.hmm_LogPi = hmm_LogPi;
    this.hmm_LogA = hmm_LogA; 
    
    //initialize sequence emission probabilities
    this.initializeEmissionProbabilities(motifPWMs, bgPWMs);
    
    //initialize read assignments
    this.initializeReadAssignments();
  }
  

  /**
   * Create an HMM with a single background state
   */
  private void createHMMStates() {
    hmm_motifFwdStates = new int[hmm_numMotifs];
    
    //this works because the first bg state is #0
    int currFirstState = hmm_numBGStates;
    int numMotifStates = 0;
    //in the GHMM version each motif is just a single state for the fwd. motif
    //and a single state for the rev. motif
    for (int i = 0; i < hmm_numMotifs; i++) {
      hmm_motifFwdStates[i] = currFirstState + (2*i);
      numMotifStates += 2;
    }      
    
    hmm_numStates = hmm_numBGStates + numMotifStates;
    
    isBGState = new boolean[hmm_numStates];
    isMotifState = new boolean[hmm_numStates];
    isMotifFwdState = new boolean[hmm_numStates];
    stateToMotifMap = new Hashtable<Integer, Integer>();

    isStateFixedDuration = new boolean[hmm_numStates];
    stateDurations = new int[hmm_numStates];
    maxStateDuration = 0;

    Arrays.fill(isBGState, false);
    Arrays.fill(isMotifState, false);
    Arrays.fill(isMotifFwdState, false);
    Arrays.fill(isStateFixedDuration, false);
    Arrays.fill(stateDurations, 1);
    hmm_numEffectiveStates = 0;
    
    for (int i = 0; i < hmm_numBGStates; i++) {
      isBGState[i] = true;
      hmm_numEffectiveStates += 1;
    }
    
    for (int i = 0; i < hmm_motifFwdStates.length; i++) {
      int firstState = hmm_motifFwdStates[i];
      isMotifState[firstState] = true;
      isMotifState[firstState + 1] = true;
      
      isMotifFwdState[firstState] = true;
      
      stateToMotifMap.put(firstState, i);
      stateToMotifMap.put(firstState + 1, i);
      
      isStateFixedDuration[firstState] = true;
      isStateFixedDuration[firstState + 1] = true;
      stateDurations[firstState] = hmm_motifLengths[i];
      stateDurations[firstState + 1] = hmm_motifLengths[i];
      maxStateDuration = Math.max(maxStateDuration, hmm_motifLengths[i]);
      
      hmm_numEffectiveStates += 2 * hmm_motifLengths[i];
    }
  }
  
  
  /**
   * Initialize binding probabilities for an HMM with a single background state
   */
  private void initializeStateBindingFrequencies() {
    hmm_binding = new DenseDoubleMatrix1D(hmm_numStates);
    
    //first make sure everything is at least the minimum frequency
    int motifBGBindingFreq = 2;
//    hmm_binding.assign(minSignalFreq);
    hmm_binding.assign(motifBGBindingFreq);
    
    //set binding freq for the bg states
    hmm_binding.setQuick(0, bgBindingFreq);
    
    //set binding freq for the motif states   
    for (int j = 0; j < hmm_numMotifs; j++) {
      int fwdBindingState = hmm_motifFwdStates[j];
      int revBindingState = hmm_motifFwdStates[j] + 1;
      hmm_binding.setQuick(fwdBindingState, motifBindingFreq);
      hmm_binding.setQuick(revBindingState, motifBindingFreq);    
    }
  }
  
  /**
   * Initialize initial state probabilities for an HMM with a single BG state
   */
  private void initializeInitialStateProbabilities() {
    hmm_LogPi = new DenseDoubleMatrix1D(hmm_numStates);
        
    //set initial probabilities for the binding states
    hmm_LogPi.setQuick(0, Math.log(1.0));
//    if (hmm_numBGStates == 1) {
//      hmm_LogPi.setQuick(0, Math.log(1.0 - motifInitialStateProb)); 
//    }
//    else {
//      /**
//       * probably should do something more complicated if there are 
//       * multiple BG states, which presumably will have distinct binding freqs.
//       * 
//       * For now, weight the probability of transitioning into the background
//       * states as inversely proportional to their initial binding frequency
//       */
//      double bgTotalProb = 1.0 - motifInitialStateProb;
//      double[] bgBindingLogProb = new double[hmm_numBGStates];
//      double invBindingFreqSum = 0.0;
//      for (int k = 0; k < hmm_numBGStates; k++) {
//        bgBindingLogProb[k] = 1.0 / hmm_binding.getQuick(k);
//        invBindingFreqSum += bgBindingLogProb[k];
//      }
//      for (int k = 0; k < hmm_numBGStates; k++) {
//        //take the log here so it doesn't need to be done repeatedly later
//        hmm_LogPi.setQuick(k, Math.log(bgTotalProb * bgBindingLogProb[k] / invBindingFreqSum));
//      }
//    }
    
    //set initial probabilities for motif states
    /**
     * set them to a portion of the motifInitialStateProb weighted by the number
     * of motifs and the lengths of the motifs
     */    
    for (int j = 0; j < hmm_numMotifs; j++) {
      hmm_LogPi.setQuick(hmm_motifFwdStates[j], Double.NEGATIVE_INFINITY);
      hmm_LogPi.setQuick(hmm_motifFwdStates[j] + 1, Double.NEGATIVE_INFINITY);
    }    
  }
  
  
  /**
   * Initialize the transition probabilities and for the states that don't have
   * a fixed duration the continuation probabilities
   */
  private void initializeTransitionProbabilities() {
    hmm_LogA = new DenseDoubleMatrix2D(hmm_numStates, hmm_numStates);
    hmm_LogA.assign(Double.NEGATIVE_INFINITY);
    hmm_transitionPriorParams = new double[hmm_numStates][hmm_numStates];
    
    /**
     * Start by initializing the transitions from BG states...
     * 
     * If we expect approximately one motif/binding event occurence per sequence
     * then because the transition probabilities are multinomials we can 
     * calculate that the probability of transitioning to the motif from BG is 
     * numRegions / (totalNumBases * numMotifs * 2.0)
     * the 2.0 is because the probability needs to be split among the motif and
     * its reverse complement. 
     * 
     * Note: there may be multiple close matches to the motif per sequence, but
     * motif matches that aren't responsible for reads will be accounted for
     * by the background state. 
     */
    double bgMotifTotalTransProb = (double)analysisRegions.size() / (double)regionSizeTotal;
    double bgBGTotalTransProb = 1.0 - bgMotifTotalTransProb;
    double bgMotifIndivTransLogProb = Math.log(bgMotifTotalTransProb / (2.0 * (double)hmm_numMotifs));
    double[][] bgBGTransLogProb = new double[hmm_numBGStates][hmm_numBGStates];
    bgBGTransLogProb[0][0] = Math.log(bgBGTotalTransProb);
//    if (hmm_numBGStates == 1) {
//     bgBGTransLogProb[0][0] = Math.log(bgBGTotalTransProb);
//    }
//    else {
//      /**
//       * probably should do something more complicated if there are 
//       * multiple BG states, which presumably will have distinct binding freqs.
//       * 
//       * For now, weight the probability of transitioning into the background
//       * states as inversely proportional to their initial binding frequency
//       */
//      double[] bgBindingLogProb = new double[hmm_numBGStates];
//      double invBindingFreqSum = 0.0;
//      for (int k = 0; k < hmm_numBGStates; k++) {
//        bgBindingLogProb[k] = 1.0 / hmm_binding.getQuick(k);
//        invBindingFreqSum += bgBindingLogProb[k];
//      }
//      for (int k = 0; k < hmm_numBGStates; k++) {
//        //take the log here so it doesn't need to be done repeatedly later
//        bgBindingLogProb[k] = Math.log(bgBGTotalTransProb * bgBindingLogProb[k] / invBindingFreqSum);
//      }
//      for (int j = 0; j < hmm_numBGStates; j++) {
//        for (int k = 0; k < hmm_numBGStates; k++) {
//          bgBGTransLogProb[j][k] = bgBindingLogProb[k];
//        }
//      }
//    }

    /** 
     * set the transition from the background states to the calculated values
     */
    for (int j = 0; j < hmm_numBGStates; j++) {
      for (int k = 0; k < hmm_numStates; k++) {
        if (isBGState[k]) {
          hmm_LogA.setQuick(j, k, bgBGTransLogProb[j][k]);
          hmm_transitionPriorParams[j][k] = Math.exp(bgBGTransLogProb[j][k]) * regionSizeTotal; 
        }
        else if (isMotifState[k]) {
          hmm_LogA.setQuick(j, k, bgMotifIndivTransLogProb);
          hmm_transitionPriorParams[j][k] = Math.exp(bgMotifIndivTransLogProb) * regionSizeTotal;
        }
      }
    }
        
    
    /**
     * Initialize the transtions from motif states 
     */
    for (int i = 0; i < hmm_numMotifs; i++) {
      int fwdMotifState = hmm_motifFwdStates[i];
      int revMotifState = hmm_motifFwdStates[i] + 1;
      
      hmm_LogA.setQuick(fwdMotifState, fwdMotifState, Math.log(motifEndStartInitialTransProb / 2.0));
      hmm_LogA.setQuick(fwdMotifState, revMotifState, Math.log(motifEndStartInitialTransProb / 2.0));
      hmm_LogA.setQuick(revMotifState, fwdMotifState, Math.log(motifEndStartInitialTransProb / 2.0));
      hmm_LogA.setQuick(revMotifState, revMotifState, Math.log(motifEndStartInitialTransProb / 2.0));
      hmm_transitionPriorParams[fwdMotifState][fwdMotifState] = motifEndStartInitialTransProb / 2.0 * regionSizeTotal; 
      hmm_transitionPriorParams[fwdMotifState][revMotifState] = motifEndStartInitialTransProb / 2.0 * regionSizeTotal;
      hmm_transitionPriorParams[revMotifState][fwdMotifState] = motifEndStartInitialTransProb / 2.0 * regionSizeTotal; 
      hmm_transitionPriorParams[revMotifState][revMotifState] = motifEndStartInitialTransProb / 2.0 * regionSizeTotal;
      
      //total amount of probability left for motif to BG transitions
      double motifBGTransTotalProb = 1.0 - motifEndStartInitialTransProb;
      
      hmm_LogA.setQuick(fwdMotifState, 0, Math.log(motifBGTransTotalProb));
      hmm_LogA.setQuick(revMotifState, 0, Math.log(motifBGTransTotalProb));
      hmm_transitionPriorParams[fwdMotifState][0] = motifBGTransTotalProb * regionSizeTotal; 
      hmm_transitionPriorParams[revMotifState][0] = motifBGTransTotalProb * regionSizeTotal;
                                               
//      if (hmm_numBGStates == 1) {
//        hmm_LogA.setQuick(endState, 0, Math.log(motifBGTransTotalProb));
//        hmm_LogA.setQuick(revEndState, 0, Math.log(motifBGTransTotalProb));
//      }
//      else {
//        /**
//         * probably should do something more complicated if there are 
//         * multiple BG states, which presumably will have distinct binding freqs.
//         * 
//         * For now, weight the probability of transitioning into the background
//         * states as inversely proportional to their initial binding frequency
//         */
//        double[] motifBGBindingLogProb = new double[hmm_numBGStates];
//        double invBindingFreqSum = 0.0;
//        for (int k = 0; k < hmm_numBGStates; k++) {
//          motifBGBindingLogProb[k] = 1.0 / hmm_binding.getQuick(k);
//          invBindingFreqSum += motifBGBindingLogProb[k];
//        }
//        for (int k = 0; k < hmm_numBGStates; k++) {
//          //take the log here so it doesn't need to be done repeatedly later
//          motifBGBindingLogProb[k] = Math.log(motifBGTransTotalProb * motifBGBindingLogProb[k] / invBindingFreqSum);
//        }
//
//        for (int k = 0; k < hmm_numBGStates; k++) {
//          hmm_LogA.setQuick(endState, revStartState, motifBGBindingLogProb[k]);
//          hmm_LogA.setQuick(revEndState, startState, motifBGBindingLogProb[k]);
//        }
//      }      
    }
    
    this.computeAdjustedTransitionProbabilities();
  }
    
  
  private void computeAdjustedTransitionProbabilities() {
    /**
     * compute adjusted transition probabilities for when various states are
     * inaccessible on account of the state duration and the position within
     * a region
     * 
     * These probabilities go into effect when the available duration is
     * _less than_ the duration used to compute the adjusted probabilities, 
     * but the CeilingEntry method of TreeMap finds entries with a key >= the 
     * specified key, so to compensate we set the keys to be the duration - 1.
     * For example, suppose the longest state has duration 20 and another 
     * state has duration 10. Suppose the current available duration is 20, 
     * Then ceilingEntry will return null, i.e. all states are available options. 
     * 
     * Suppose instead that the available duration is 10 to 19, then 
     * ceilingEntry will return the adjusted probabilities for when the 
     * longest state (only) is unavailable.
     */
    hmm_LogAdjustedTransProb = new TreeMap<Integer, DoubleMatrix2D>();
    TreeSet<Integer> durations = new TreeSet<Integer>();
    for (int j = 0; j < hmm_numStates; j++) {
      if (isStateFixedDuration[j] && (stateDurations[j] > 1)) {
        durations.add(stateDurations[j]);
      }
    }
    
    for (int duration : durations) {
      DoubleMatrix2D logAdjustedTransProb = hmm_LogA.copy();
      for (int j = 0; j < hmm_numStates; j++) {
        double totalLogProb = Double.NEGATIVE_INFINITY;
        for (int k = 0; k < hmm_numStates; k++) {          
          if (isStateFixedDuration[k] && (stateDurations[k] >= duration)) {
              logAdjustedTransProb.setQuick(j, k, Double.NEGATIVE_INFINITY);
          }
          else {
            totalLogProb = Numerical.log_add(totalLogProb, logAdjustedTransProb.getQuick(j, k));
          }          
        }
        /**
         * adjustment to renormalize probabilities is: orig. prob / totalProb
         * which in log space becomes: orig. log prob - totalLogProb  
         */        
        logAdjustedTransProb.viewRow(j).assign(Functions.minus(totalLogProb));
      }
      hmm_LogAdjustedTransProb.put((duration - 1), logAdjustedTransProb);
    }
  }
  
  /**
   * 
   * @param motifPWMs
   * @param bgPWMs
   */
  private void initializeEmissionProbabilities(List<WeightMatrix> motifPWMs, List<WeightMatrix> bgPWMs) {
    hmm_phi = new DenseDoubleMatrix2D[hmm_numStates];
    //set up the emission probabilities for all the bg states to be the bgPWM
    for (int i = 0; i < hmm_numBGStates; i++) {
      hmm_phi[i] = new DenseDoubleMatrix2D(1, 4);
    	WeightMatrix currentPWM = bgPWMs.get(i);
      for (int j = 0; j < WeightMatrix.letters.length; j++) {
        hmm_phi[i].setQuick(i, this.baseToIndex(WeightMatrix.letters[j]), currentPWM.matrix[0][WeightMatrix.letters[j]]);
      }
    }    
    
    //set up the emission probabilities for the motifs
    for (int i = 0; i < hmm_numMotifs; i++) {
      WeightMatrix currentPWM = motifPWMs.get(i);
      int fwdState = hmm_motifFwdStates[i];
      int revState = fwdState + 1;
      hmm_phi[fwdState] = new DenseDoubleMatrix2D(currentPWM.length(), 4);
      hmm_phi[revState] = new DenseDoubleMatrix2D(currentPWM.length(), 4);
      
      for (int j = 0; j < currentPWM.length(); j++) {
        for (int k = 0; k < WeightMatrix.letters.length; k++) {
          hmm_phi[fwdState].setQuick(j, this.baseToIndex(WeightMatrix.letters[k]), currentPWM.matrix[j][WeightMatrix.letters[k]]);
          
          
          int revCompPosIndex = currentPWM.length() - 1 - j;         
          hmm_phi[revState].setQuick(j, this.baseToIndex(WeightMatrix.letters[k]), currentPWM.matrix[revCompPosIndex][WeightMatrix.revCompLetters[k]]);
        }
      }
    }
  }
  
  
  private void initializeReadAssignments() {
    for (HMMAnalysisRegion region : analysisRegions) {
      int currentRegionStart = region.getGenomicRegion().getStart();
      List<ReadHit> hits = region.getSignalReads();
      for (int n = 0; n < hits.size(); n++) {
        ReadHit currentHit = hits.get(n);
        int hitStart;
        int newLoc;
        //FIXME 49 is the maximum of the shear distribution, but this should
        //be fixed to avoid using a constant
        hitStart = currentHit.getFivePrime() - currentRegionStart;
        if (currentHit.getStrand() == '+') {
          newLoc = Math.min(hitStart + 50, region.getWidth()-1);
        }
        else {
          //subtract the min and max in order to keep mStart as the lesser value
          newLoc = Math.max(hitStart - 50, 0);
        }
        region.assignRead(currentHit, newLoc);
      }
      region.printNumReadsByLoc();
      
      
      //Also initialize the regionToMKReadCount hash
      regionToMKReadCounts.put(region, new int[region.getWidth()][hmm_numStates]);
    }    
  }

  public List<Feature> execute() {        
    logger.debug("Starting EM");

    logger.debug("======Initial Parameters======");
    this.printPi();
    this.printTransitionProbabilities();
    this.printSequenceEmission();
    this.printBindingFrequencies();
    this.printReadAssignments();              
    this.printMLStateSequence();        
    this.printFeatures();

    
    int numIterations = 0;
    boolean converged = false;
    while ((numIterations < maxIterations) && !converged) {

      boolean computeViterbi = true;//TODO (numIterations % 50 == 0);
      
      //TODO come up with a strategy for alternating between these two updates
      if ((numIterations > 0) && (numIterations % 2 == 1/*0*/)) { //TODO set 1 back to 0
        this.iterateEM(true, false, computeViterbi);
      }
      else {
        this.iterateEM(false, true, computeViterbi);
      }
      
      signalFeatures = this.getFeaturesFromViterbi();
      if (numIterations % 50 == 0) {
        
        
        this.printPi();
        this.printTransitionProbabilities();
        this.printSequenceEmission();
        this.printBindingFrequencies();
        //this.printReadAssignments();              
        this.printMLStateSequence();        
        this.printFeatures();
      }
      
      numIterations++;
      if ((Math.abs(totalLogLikelihood - prevTotalLogLikelihood) < convergenceThreshold) && (numIterations >= MIN_ITER_DFLT)) {
        converged = true;
      }
      logger.debug("Finished iteration " + numIterations + " LL: " + totalLogLikelihood);
    }
    
    /**
     * Run Forward-Backward one last time to compute the likelihood and the
     * Viterbi decoding for the final parameter estimates
     */
    prevTotalLogLikelihood = totalLogLikelihood;
    totalLogLikelihood = 0;
    for (HMMAnalysisRegion currentRegion : analysisRegions) {
      this.runForwardBackward(currentRegion, true);
      totalLogLikelihood = totalLogLikelihood + currentRegion.getLogLikelihood();
    }
    
    signalFeatures = this.getFeaturesFromViterbi();        

    this.printPi();
    this.printTransitionProbabilities();
    this.printSequenceEmission();
    this.printBindingFrequencies();
    this.printReadAssignments();
    this.printMLStateSequence();
    

    return signalFeatures;
  }
  
  
  private void iterateEM(boolean updateBindingFreq, boolean updateReadAssignments, boolean computeViterbi) {
    //1. compute expectations
    //run forward-backward on each region to compute alpha and beta
  	prevTotalLogLikelihood = totalLogLikelihood;
    totalLogLikelihood = 0;
    for (HMMAnalysisRegion currentRegion : analysisRegions) {
      this.runForwardBackward(currentRegion, computeViterbi);
      totalLogLikelihood = totalLogLikelihood + currentRegion.getLogLikelihood();
    }
    
    //next compute gamma and xi using the totalLogLikelihood
    for (HMMAnalysisRegion currentRegion : analysisRegions) {
    	//FIXME This should probably be the region's likelihood rather than the total likelihood
      //compute gamma = alpha*beta/Likelihood
      currentRegion.setLogGamma(currentRegion.getLogAlpha().copy().assign(currentRegion.getLogBeta(), Functions.plus).assign(Functions.minus(totalLogLikelihood)));
      
      //compute xi(m-1,m) = alpha(m-1) * emission * transition * beta(m)/Likelihood      
      currentRegion.setLogXi(currentRegion.getLogXiNumerator().copy().assign(Functions.minus(totalLogLikelihood)));
    }
    
    
    //2. Update initial state (pi) parameters
    this.updatePi();
    
    
    //3. Update transition parameters
    this.updateTransitionProbabilities();
    
    
    //4. Update sequence emission parameters
    this.updateSequenceEmission();
    
    
    //5. Alternate between updating hit assignments and updating binding frequencies        
    if (updateBindingFreq) {
      if (!hmm_useWCEforBGBinding) {
        this.updateBGStateBindingFrequencies();
      }      
      this.updateMotifBindingFrequencies(true);
    }
    if (updateReadAssignments) {
      //this.batchUpdateReadAssignments();
      this.nonBatchUpdateReadAssignments(isReadUpdateStochastic);
    }
    
  }
  
  
  
  /**
   * run the forward backward algorithm on the specified region
   * @param region
   */
  private void runForwardBackward(HMMAnalysisRegion region, boolean computeViterbi) {
    int regionWidth = region.getWidth();
    DenseDoubleMatrix3D logXiNumerator = new DenseDoubleMatrix3D(regionWidth - 1, hmm_numStates, hmm_numStates);
    
    /**
     * Run Forward algorithm
     */
    DoubleMatrix2D logAlpha = new DenseDoubleMatrix2D(regionWidth, hmm_numStates);
    DoubleMatrix2D logEmissionProbs = new DenseDoubleMatrix2D(regionWidth, hmm_numStates);
    DoubleMatrix2D viterbiLogProb = new DenseDoubleMatrix2D(regionWidth, hmm_numStates);
    
    logXiNumerator.assign(Double.NEGATIVE_INFINITY);
    logAlpha.assign(Double.NEGATIVE_INFINITY);
    logEmissionProbs.assign(Double.NEGATIVE_INFINITY);
    viterbiLogProb.assign(Double.NEGATIVE_INFINITY);
    
    
    int[][] viterbiTraceback = new int[regionWidth][hmm_numStates];
    for (int m = 0; m < regionWidth; m++) {
      Arrays.fill(viterbiTraceback[m], -1);
    }

    int[][] mkReadCounts = regionToMKReadCounts.get(region);
    
    /**
     * variables for holding values computed during forward backward alg.
     */
    double logReadProb = 0;
    double logSequenceProb = 0;    
    double logBindingProb = 0;
    double motifLogReadProb = 0;
    double motifLogBindingProb = 0;
    double logEmissionProb = 0;
    double logAlpha_z_1k = 0;
    

    //compute alpha(z_1)
    List<ReadHit> readsAtLoc = region.readsAssignedTo(0);
    logReadProb = region.computeLogReadProb(readsAtLoc, 0);
    for (int k = 0; k < hmm_numStates; k++) {
      //if the initial transition probability is 0 then don't bother computing
      //the sequence or read emission probabilities
      double initLogProb = hmm_LogPi.getQuick(k); 
      if (initLogProb != Double.NEGATIVE_INFINITY) {
        logSequenceProb = this.computeLogSequenceProb(hmm_phi[k], region.sequenceAsInts, 0);
        
      	//FIXME need to adjust the read count when the state duration is > 1
        mkReadCounts[0][k] = readsAtLoc.size();
        logBindingProb = this.computeLogBindingProb(region, 0, k, mkReadCounts[0][k]);
        
        
        logEmissionProbs.setQuick(0, k, logSequenceProb + logReadProb + logBindingProb);
        logAlpha_z_1k = initLogProb + logSequenceProb + logReadProb + logBindingProb; 
        logAlpha.setQuick(0, k, logAlpha_z_1k);
        
        //initialize the viterbi matrix the same way
        viterbiLogProb.setQuick(0, k, logAlpha_z_1k);
      }
      else {
        logAlpha.setQuick(0, k, Double.NEGATIVE_INFINITY);
        viterbiLogProb.setQuick(0, k, Double.NEGATIVE_INFINITY);
        mkReadCounts[0][k] = -1;
      }
    }
    
    
    //matrix for holding recurrence values computed during algorithm
//    DoubleMatrix1D logRecurrenceProbsMat = new DenseDoubleMatrix1D(hmm_numStates);
//    DoubleMatrix1D logViterbiRecurrenceProbsMat = new DenseDoubleMatrix1D(hmm_numStates);
    
    double maxLogRecurrenceProb;
    int maxLogRecurrenceProbIdx;    
    double maxViterbiLogRecurrenceProb;
    int maxViterbiLogRecurrenceProbIdx;
    
    double logTransProb;
    double logAlpha_z_prev;
    double viterbiLogPrev = 0;
    double viterbiLogRecurrenceProb;
    double logAlpha_z_mk;
    
    DoubleMatrix2D currentTransProb;
    
    /**
     * compute alpha values for the rest of the locations
     */
    for (int m = 1; m < regionWidth; m++) {
      currentTransProb = this.getAdjustedTransProbs(m, regionWidth);
      
      readsAtLoc = region.readsAssignedTo(m);      
      logReadProb = region.computeLogReadProb(readsAtLoc, m);
      
			for (int k = 0; k < hmm_numStates; k++) {
				/**
				 * Check that state k isn't too long to start at this location
				 * The condition should be <= regionWidth for BG states because the 
				 * indexing is 0-based and position m is the first position of the state.
				 * For non-BG states the condition should be < regionWidth, because we
				 * want to force the HMM to always end in a BG state   
				 */
				if ((m + stateDurations[k]) <= regionWidth) { 				    
					//the probability of emitting the upcoming sequence from state k
					logSequenceProb = this.computeLogSequenceProb(hmm_phi[k], region.sequenceAsInts, m);

					if (stateDurations[k] == 1) {
						mkReadCounts[m][k] = readsAtLoc.size();
						logBindingProb = this.computeLogBindingProb(region, m, k, mkReadCounts[m][k]);

						// emission prob is sequence emission, read emission, and binding prob
						logEmissionProb = logSequenceProb + logReadProb + logBindingProb;
					}
					else {
						/**
						 * if this is a reverse motif state reuse the calculations from the
						 * previous iteration
						 */					
						if ((k > 1) && isMotifFwdState[k-1]) {
						  mkReadCounts[m][k] = mkReadCounts[m][k-1];
							logEmissionProb = logSequenceProb + motifLogReadProb + motifLogBindingProb;
						}
						else {
							//compute them new and store them for the next iteration
							mkReadCounts[m][k] = readsAtLoc.size();
							motifLogReadProb = logReadProb;
							for (int n = m + 1; n < m + stateDurations[k]; n++) {
								List<ReadHit> reads = region.readsAssignedTo(n);
								motifLogReadProb += region.computeLogReadProb(reads, n);
								mkReadCounts[m][k] += reads.size();
							}
							motifLogBindingProb = this.computeLogBindingProb(region, m, k, mkReadCounts[m][k]);
							logEmissionProb = logSequenceProb + motifLogReadProb + motifLogBindingProb;
						}
					}
					logEmissionProbs.setQuick(m, k, logEmissionProb);
					
					/**
					 * need to store the probabilities of getting to this state from the
					 * possible previous states and then add them up using the logsumexp
					 * method...
					 */
					double[] logRecurrenceProbs = new double[hmm_numStates];

					/**
					 * keep track of both the max. previous alpha probability (recurrence
					 * probability) and the index of it to speed the log-sum-exp calculation
					 */
					maxLogRecurrenceProb = Double.NEGATIVE_INFINITY;
					maxLogRecurrenceProbIdx = -1;

					/**
					 * and similarly for max viterbi prob
					 */
					maxViterbiLogRecurrenceProb = Double.NEGATIVE_INFINITY;
					maxViterbiLogRecurrenceProbIdx = -1;

					for (int j = 0; j < hmm_numStates; j++) {
						//figure out how many positions back to look for the start of state j
						int jStart = m - 1;
						if (isStateFixedDuration[j]) {
							jStart = m - stateDurations[j];
						}

						if ((jStart) >= 0) {
							logTransProb = currentTransProb.getQuick(j, k);

							logAlpha_z_prev = logAlpha.getQuick(jStart, j);
							logRecurrenceProbs[j] = logAlpha_z_prev + logTransProb;
							if (logRecurrenceProbs[j] > maxLogRecurrenceProb) {
								maxLogRecurrenceProb = logRecurrenceProbs[j];
								maxLogRecurrenceProbIdx = j;
							}          

							viterbiLogPrev = viterbiLogProb.getQuick(jStart, j);
							viterbiLogRecurrenceProb = viterbiLogPrev + logTransProb;
							if (viterbiLogRecurrenceProb > maxViterbiLogRecurrenceProb) {
								maxViterbiLogRecurrenceProb = viterbiLogRecurrenceProb;
								maxViterbiLogRecurrenceProbIdx = j;
							}

							/**
							 * do the first part of the calculation for Xi - everything but beta term
							 *
							 * The indexing for Xi is offset, so logXi[0,j,k] should hold
							 * the value for xi(z_{0,j}, z_{1,k} and in general
							 * logXi[m,j,k] =  xi(z_{m,j}, z_{m+1,k} 
							 * or to keep consistent with the notation jStart, m
							 * logXi[jStart,j,k] = xi(z_{jStart,j}, z_{m,k}) 
							 * = alpha(z_{jStart,j}) * emissionProb * transitionProb * beta(z_{m,k})
							 */
							logXiNumerator.setQuick(jStart, j, k, logRecurrenceProbs[j] + logEmissionProb);
						}
						else {
							logRecurrenceProbs[j] = Double.NEGATIVE_INFINITY;
						}
					}

					if (maxLogRecurrenceProb != Double.NEGATIVE_INFINITY) {
						logAlpha_z_mk = logEmissionProb + Numerical.log_sum_exp(logRecurrenceProbs, maxLogRecurrenceProb); 
					}
					else {
						logAlpha_z_mk = Double.NEGATIVE_INFINITY;
					}
					//          double matTest = logEmissionProb + Numerical.log_sum_exp(logRecurrenceProbsMat);
					logAlpha.setQuick(m, k, logAlpha_z_mk);

					//update the viterbi probabilities
					viterbiLogProb.setQuick(m, k, logEmissionProb + maxViterbiLogRecurrenceProb);
					//          double vitTest = logViterbiRecurrenceProbsMat.aggregate(Functions.max, Functions.identity);
					viterbiTraceback[m][k] = maxViterbiLogRecurrenceProbIdx;
				}
				else {
					//state k is too long to start at position m, so fill in accordingly 
					logAlpha.setQuick(m, k, Double.NEGATIVE_INFINITY);
					viterbiLogProb.setQuick(m, k, Double.NEGATIVE_INFINITY);
					viterbiTraceback[m][k] = -1;
					mkReadCounts[m][k] = -1;
				}
			} 
    }     
    region.setViterbiLogProb(viterbiLogProb);
    region.setViterbiTraceback(viterbiTraceback);
    
    /**
     * Run Backward Algorithm
     */
    double[] logSumComponents = new double[hmm_numStates];
    double maxComponent = Double.NEGATIVE_INFINITY;
    double logBeta_z_next;
    
    DoubleMatrix2D logBeta = new DenseDoubleMatrix2D(regionWidth, hmm_numStates);
    logBeta.assign(Double.NEGATIVE_INFINITY);
    //initialize log beta(z_m), set to 0 ( 0 = ln(1) ) if it's possible to finish in that state 
    for (int k = 0; k < hmm_numStates; k++) {
    	if (isStateFixedDuration[k] && (stateDurations[k] > 1)) {
    		logBeta.setQuick(regionWidth - 1, k, Double.NEGATIVE_INFINITY);
    	}
    	else { 
    		logBeta.setQuick(regionWidth - 1, k, 0.0);
    	}
    }
    
    //compute the beta values for the rest of the locations
    for (int m = regionWidth - 2; m >= 0; m--) {      
      for (int j = 0; j < hmm_numStates; j++) {         
        int betaNextOffset = 1;
				if (isStateFixedDuration[j]) {
					betaNextOffset = stateDurations[j];
				}
				/**
				 * If the duration is > 1 we need to check that it isn't too long to
				 * have started at position m. Also, normally we would use Beta(m+1), 
				 * for the probability of the emission coming after this state, so 
				 * instead we use Beta(m + duration)
				 */
				if((m + betaNextOffset) < regionWidth) {
          currentTransProb = this.getAdjustedTransProbs(m + betaNextOffset, regionWidth);
	      	maxComponent = Double.NEGATIVE_INFINITY;
	      	
					for (int k = 0; k < hmm_numStates; k++) {
	          logTransProb = currentTransProb.getQuick(j, k);
	          logEmissionProb = logEmissionProbs.getQuick(m + betaNextOffset, k);          
	          logBeta_z_next = logBeta.getQuick(m + betaNextOffset, k);
	          logSumComponents[k] = logTransProb + logEmissionProb + logBeta_z_next;
	          maxComponent = Math.max(maxComponent, logSumComponents[k]);
	          
	          /**
	           * finish the calculation of xi here...
	           */
	           logXiNumerator.setQuick(m, j, k, logXiNumerator.getQuick(m,j,k) + logBeta_z_next);          
	        }
					
	        double logBeta_z_mj;
	        if (maxComponent != Double.NEGATIVE_INFINITY) {
	          logBeta_z_mj = Numerical.log_sum_exp(logSumComponents, maxComponent);
	        }
	        else {
	          logBeta_z_mj = Double.NEGATIVE_INFINITY;
	        }
	        
	        logBeta.setQuick(m, j, logBeta_z_mj);        					
				}
				else if ((m + betaNextOffset) == regionWidth) {
					logBeta.setQuick(m, j, 0.0);
				}
				else {
				  logBeta.setQuick(m, j, Double.NEGATIVE_INFINITY);
				}
      }
    }
    
    
    /**
     * Check that the likelihoods match
     */    
    this.checkLikelihood(logAlpha, logBeta);
    double regionLogLikelihood = Numerical.log_sum_exp(logAlpha.viewRow(0).copy().assign(logBeta.viewRow(0), Functions.plus));
    
    region.setLogAlphaBeta(logAlpha, logBeta, regionLogLikelihood);
    region.setLogXiNumerator(logXiNumerator);
  }
  
  
  /**
   * Update the initial state probabilities. 
   * Calculations done in log-space.
   */
  private void updatePi() {
    //set up a matrix whose columns will hold the gamma(z_0) vector of each
    //analysis region
    DoubleMatrix2D logGamma_zeros = new DenseDoubleMatrix2D(hmm_numStates, analysisRegions.size());
    for (int i = 0; i < analysisRegions.size(); i++) {
      HMMAnalysisRegion region = analysisRegions.get(i);
      logGamma_zeros.viewColumn(i).assign(region.getLogGamma().viewRow(0));      
    }
    
    //the sums along each of the rows will be the numerators for the update eqs  
    DoubleMatrix1D hmm_LogPi_new = new DenseDoubleMatrix1D(hmm_numStates);
    for (int k = 0; k < hmm_numStates; k++) {
      hmm_LogPi_new.setQuick(k, Numerical.log_sum_exp(logGamma_zeros.viewRow(k)));
    }
    
    //summing all of those numerators will give the denominator for the update eq
    double logGammaSum = Numerical.log_sum_exp(hmm_LogPi_new);
    hmm_LogPi_new.assign(Functions.minus(logGammaSum));
    
    hmm_LogPi = hmm_LogPi_new;
  }
  
  
  /**
   * Calculate updated values for the transition probabilities. 
   * Calculations done in log-space.
   * 
   */
  private void updateTransitionProbabilities() {
    //iterate over each transition start state and update the probabilities for
    //each state it could transition to
    for (int j = 0; j < hmm_numStates; j++) {
       //create an array to hold all the log-space values that need to be summed    
      double[][] logXiExps = new double[hmm_numStates][regionSizeTotal - analysisRegions.size()];
      
      int regionOffset = 0;
      for (int i = 0; i < analysisRegions.size(); i++) {
        HMMAnalysisRegion region = analysisRegions.get(i);
        int regionWidth = region.getWidth();
        DoubleMatrix2D logXiRow = region.getLogXi().viewRow(j);
        /**
         * copy all of the relevant
         * xi terms to the appropriate space in the logXiExps array
         */
        for (int k = 0; k < hmm_numStates; k++) {
          System.arraycopy(logXiRow.viewColumn(k).toArray(), 0, 
              logXiExps[k], regionOffset, regionWidth - 1);
        }
        regionOffset = regionOffset + regionWidth - 1;
      }

      //perform the summations using log_sum_exp
      double[] pseudocounts = new double[hmm_numStates];
      StatUtil.dirichlet_rnd(pseudocounts, hmm_transitionPriorParams[j], hmm_numStates);
      double[] logXiSums = new double[hmm_numStates];
      double maxLogXiSum = Double.NEGATIVE_INFINITY;
      for (int k = 0; k < hmm_numStates; k++) {
        logXiSums[k] = Numerical.log_sum_exp(logXiExps[k]);
        double pseudocount = Math.log(pseudocounts[k] * hmm_transitionPriorWeight * regionSizeTotal);
        logXiSums[k] = Numerical.log_add(logXiSums[k], pseudocount);
        maxLogXiSum = Math.max(maxLogXiSum, logXiSums[k]);
      }
      double logXiSumsTotal = Numerical.log_sum_exp(logXiSums, maxLogXiSum);
      
      
      //update the transition probabilities
      for (int k = 0; k < hmm_numStates; k++) {
        hmm_LogA.setQuick(j, k, logXiSums[k] - logXiSumsTotal);
      }
      if (Math.abs(Numerical.log_sum_exp(hmm_LogA.viewRow(j))) > 1e-6) {
        logger.error("Transition probabilities don't sum to 1");
      }
    }
    
    this.computeAdjustedTransitionProbabilities();
  }
  
  
  /**
   * 
   */
  private void updateSequenceEmission() {
    /**
     * these are handled in separate subfunctions because if the motif is already
     * known it shouldn't be updated, and when it is updated the update for 
     * the emission probability for both the motif state and the reverse
     * complement state need to be updated together
     */
    this.updateBGStateSequenceEmission();  
    if (hmm_updateMotif) {
      this.updateMotifStateSequenceEmission();
    }
  }
  
  
  /**
   * 
   */
  private void updateBGStateSequenceEmission() {
    for (int k = 0; k < hmm_numBGStates; k++) {
      //initialize an array to hold the sums of the gamma values
      double[] logGammaSums = new double[4];
      Arrays.fill(logGammaSums, Double.NEGATIVE_INFINITY);
      //iterate over all the analysis regions summing the gamma values
      for (int i = 0; i < analysisRegions.size(); i++) {
        HMMAnalysisRegion currentRegion = analysisRegions.get(i);
        
        //iterate over all the gamma values adding each to the sum for the
        //base at the corresponding location
        DoubleMatrix1D logGammas = currentRegion.getLogGamma().viewColumn(k);
        String sequence = currentRegion.getSequence();
        for (int m = 0; m < logGammas.size(); m++) {
          double logGamma = logGammas.getQuick(m);
          switch(sequence.charAt(m)) {
            case 'A': logGammaSums[0] = Numerical.log_add(logGammaSums[0], logGamma);
              break;
            case 'C': logGammaSums[1] = Numerical.log_add(logGammaSums[1], logGamma);     
              break;
            case 'G': logGammaSums[2] = Numerical.log_add(logGammaSums[2], logGamma);
              break;
            case 'T': logGammaSums[3] = Numerical.log_add(logGammaSums[3], logGamma);
              break;
            default: logger.error("Invalid character, " + sequence.charAt(m) 
                + ", in region " + currentRegion.getGenomicRegion().regionString()
                + " at offset " + m + ".");
          }
        }
      }

      //normalize and update the probabilities
      double totalLogGammaSum = Numerical.log_sum_exp(logGammaSums);
      for (int c = 0; c < logGammaSums.length; c++) {
        double prob = (Math.exp(logGammaSums[c] - totalLogGammaSum) + PSEUDOCOUNT) / (1.0 + 4 * PSEUDOCOUNT);
        hmm_phi[k].setQuick(k, c, prob);
      }
    }
  }
  
  
  /**
   * 
   */
  private void updateMotifStateSequenceEmission() {
    for (int j = 0 ; j < hmm_motifFwdStates.length; j++) {
      int currentMotifLength = hmm_motifLengths[j];
      int motifFwdState = hmm_motifFwdStates[j];
      int motifRevState = motifFwdState + 1;
      
      //initialize an array to hold the sums of the gamma values
      double[][] logGammaSums = new double[currentMotifLength][4];
      for (int l = 0; l < logGammaSums.length; l++) {
       	Arrays.fill(logGammaSums[l], Double.NEGATIVE_INFINITY);
      }    

      //iterate over all the analysis regions summing the gamma values
      for (int i = 0; i < analysisRegions.size(); i++) {
        HMMAnalysisRegion currentRegion = analysisRegions.get(i);
          
        //iterate over all the gamma values adding each to the sum for the
        //base at the corresponding location
        DoubleMatrix1D logGammasMotif = currentRegion.getLogGamma().viewColumn(motifFwdState);
        DoubleMatrix1D logGammasRevComp = currentRegion.getLogGamma().viewColumn(motifRevState);
          
        /**
         * note: logGammasMotif.size() should be the same as currentRegion.size()
         * it is _not_ the motif length
         */
        for (int m = 0; m < logGammasMotif.size(); m++) {
          double logGammaMotif = logGammasMotif.getQuick(m);
          double logGammaRevComp = logGammasRevComp.getQuick(m);
          for (int l = 0; l < currentMotifLength; l++) {
          	int revCompPosIndex = currentMotifLength - 1 - l; 
          	int baseIndex = currentRegion.sequenceAsInts[m+l];
          	int revCompBaseIndex = BindingMotifHMM.baseIndexToRevCompIndex(baseIndex);
          	logGammaSums[l][baseIndex] = Numerical.log_add(logGammaSums[l][baseIndex], logGammaMotif);
          	logGammaSums[revCompPosIndex][revCompBaseIndex] =  Numerical.log_add(logGammaSums[revCompPosIndex][revCompBaseIndex], logGammaRevComp);
          }
        }

        //normalize and update the probabilities
        double[] totalLogGammaSum = new double[currentMotifLength];
        for (int l = 0; l < currentMotifLength; l++) {
        	totalLogGammaSum[l] = Numerical.log_sum_exp(logGammaSums[l]);
        }
        for (int l = 0; l < currentMotifLength; l++) {
        	int revCompPosIndex = currentMotifLength - 1 - l;        
        	for (int c = 0; c < logGammaSums[l].length; c++) {
        		int revCompBaseIndex = BindingMotifHMM.baseIndexToRevCompIndex(c);      
        		double prob = (Math.exp(logGammaSums[l][c] - totalLogGammaSum[l]) + PSEUDOCOUNT) / (1.0 + 4.0 * PSEUDOCOUNT);
        		hmm_phi[motifFwdState].setQuick(l, c, prob);
        		/**
        		 * for the state representing the reverse complement of the same motif
        		 * position the same probability should be used for the reverse complement
        		 * of the current base
        		 */
        		hmm_phi[motifRevState].setQuick(revCompPosIndex, revCompBaseIndex, prob);
        	}
        }        
      }
    }
  }
  
  
  /**
   * 
   */
  private void updateBGStateBindingFrequencies() {
    /**
     * Update the binding frequencies for the background states and for the
     * one motif state (and its reverse complement) that have a non-zero
     * binding frequency
     */
    //TODO move this snippet of code to initialization so it only gets run once
    List<Integer> bgFrequencyStates = new ArrayList<Integer>();
    //add the bg states if the WCE isn't being used for the binding freqs
    if (!hmm_useWCEforBGBinding) {
      for (int k = 0; k < hmm_numBGStates; k++) {
        bgFrequencyStates.add(k);
      }
    }
    

    //start by filling in the gamma values just in the denominator
    //array to hold the terms to be summed for the denominator    
    double[][] denominatorExps = new double[bgFrequencyStates.size()][regionSizeTotal];

    int regionOffset = 0;
    for (int i = 0; i < analysisRegions.size(); i++) {
      HMMAnalysisRegion currentRegion = analysisRegions.get(i);
      int width = currentRegion.getWidth();
      DoubleMatrix2D logGamma = currentRegion.getLogGamma();
      
      for (int k = 0; k < bgFrequencyStates.size(); k++) {
        int stateIndex = bgFrequencyStates.get(k);
        System.arraycopy(logGamma.viewColumn(stateIndex).toArray(), 0, denominatorExps[k], regionOffset, width);
      }
      
      regionOffset = regionOffset + width;
    }
  
    //copy the whole denominator array to the numerator exps array
    double[][] numeratorExps = new double[bgFrequencyStates.size()][regionSizeTotal];
    for (int k = 0; k < bgFrequencyStates.size(); k++) {
      System.arraycopy(denominatorExps[k], 0, numeratorExps[k], 0, regionSizeTotal);
    }
    
    //and add the read quantities to the numerator terms
    //also figure out max. exps.
    regionOffset = 0;
    double numeratorMaxExps[] = new double[bgFrequencyStates.size()];
    double denominatorMaxExps[] = new double[bgFrequencyStates.size()];
    Arrays.fill(numeratorMaxExps, Double.NEGATIVE_INFINITY);
    Arrays.fill(denominatorMaxExps, Double.NEGATIVE_INFINITY);
    for (int i = 0; i < analysisRegions.size(); i++) {
      HMMAnalysisRegion currentRegion = analysisRegions.get(i);
      for (int m = 0; m < currentRegion.getWidth(); m++) {
        double logNumReadsAtLoc = Math.log(currentRegion.numReadsAssignedTo(m));
        for (int k = 0; k < bgFrequencyStates.size(); k++) {
          //check the denominator max. exp.
          denominatorMaxExps[k] = Math.max(denominatorMaxExps[k], denominatorExps[k][regionOffset + m]);
          
          //add the read quantity
          numeratorExps[k][regionOffset + m] += logNumReadsAtLoc;
          
          //check the numerator max. exp.
          numeratorMaxExps[k] = Math.max(numeratorMaxExps[k], numeratorExps[k][regionOffset + m]);
        }
      }
    }
    
    double[] numeratorSums = new double[bgFrequencyStates.size()];
    double[] denominatorSums = new double[bgFrequencyStates.size()];
    for (int k = 0; k < bgFrequencyStates.size(); k++) {
      numeratorSums[k] = Numerical.log_sum_exp(numeratorExps[k], numeratorMaxExps[k]);
      denominatorSums[k] = Numerical.log_sum_exp(denominatorExps[k], denominatorMaxExps[k]);
      

      double newBindingFreq = Math.exp(numeratorSums[k] - denominatorSums[k]);
      
      //don't let the bindingFreq become 0...
      newBindingFreq = Math.max(newBindingFreq, minSignalFreq);
      //warn if the binding frequency is getting very low
//      if (newBindingFreq < 2) {
//        logger.warn("Binding Frequency updated to " + newBindingFreq + " for state " + bgFrequencyStates.get(k) + ".");
//      }
      hmm_binding.setQuick(bgFrequencyStates.get(k), newBindingFreq);
      
    }
  }

  
  /**
   * 
   */
  private void updateMotifBindingFrequencies(boolean coupled) {
    /**
     * Update the binding frequencies for the motif states
     */
    //start by filling in the gamma values just in the denominator
    //array to hold the terms to be summed for the denominator    
    double[][] denominatorExps;
//    if (coupled) {
//      denominatorExps = new double[bindingFrequencyStates.size()][2*regionSizeTotal];
//    }
//    else {
      denominatorExps = new double[2*hmm_numMotifs][regionSizeTotal];
//    }

    int regionOffset = 0;
    for (int i = 0; i < analysisRegions.size(); i++) {
      HMMAnalysisRegion currentRegion = analysisRegions.get(i);
      int width = currentRegion.getWidth();
      DoubleMatrix2D logGamma = currentRegion.getLogGamma();
      
      for (int n = 0; n < hmm_numMotifs; n++) {
        int fwdState = hmm_motifFwdStates[n];
        int revState = fwdState + 1;
        System.arraycopy(logGamma.viewColumn(fwdState).toArray(), 0, denominatorExps[n], regionOffset, width);
//        if (coupled) {
//          System.arraycopy(logGamma.viewColumn(revStateIndex).toArray(), 0, denominatorExps[k], regionSizeTotal + regionOffset, width);
//        }
//        else {
          System.arraycopy(logGamma.viewColumn(revState).toArray(), 0, denominatorExps[n + hmm_numMotifs], regionOffset, width);
//        }
      }
      
      regionOffset = regionOffset + width;
    }
  
    //copy the whole denominator array to the numerator exps array
    double[][] numeratorExps = new double[2*hmm_numMotifs][regionSizeTotal];
    for (int k = 0; k < numeratorExps.length; k++) {
      System.arraycopy(denominatorExps[k], 0, numeratorExps[k], 0, regionSizeTotal);
    }
    
    //and add the read quantities to the numerator terms
    //also figure out max. exps.
    regionOffset = 0;
    double numeratorMaxExps[] = new double[2*hmm_numMotifs];
    double denominatorMaxExps[] = new double[2*hmm_numMotifs];
    Arrays.fill(numeratorMaxExps, Double.NEGATIVE_INFINITY);
    Arrays.fill(denominatorMaxExps, Double.NEGATIVE_INFINITY);
    for (int i = 0; i < analysisRegions.size(); i++) {
      HMMAnalysisRegion currentRegion = analysisRegions.get(i);
      int[][] mkReadCounts = regionToMKReadCounts.get(currentRegion);
      for (int m = 0; m < currentRegion.getWidth(); m++) {        
        for (int k = 0; k < hmm_numMotifs; k++) {
          int readCount = mkReadCounts[m][hmm_motifFwdStates[k]];
          double logNumReadsAtLoc;
          if (readCount >= 0) {
            logNumReadsAtLoc = Math.log(readCount);
          }
          else {
            logNumReadsAtLoc = Double.NEGATIVE_INFINITY;
          }
          //check the denominator max. exp.
          denominatorMaxExps[k] = Math.max(denominatorMaxExps[k], denominatorExps[k][regionOffset + m]);
          denominatorMaxExps[k + hmm_numMotifs] = Math.max(denominatorMaxExps[k + hmm_numMotifs], 
              denominatorExps[k + hmm_numMotifs][regionOffset + m]);
          
          //add the read quantity
          numeratorExps[k][regionOffset + m] += logNumReadsAtLoc;
          numeratorExps[k + hmm_numMotifs][regionOffset + m] += logNumReadsAtLoc;
          
          //check the numerator max. exp.
          numeratorMaxExps[k] = Math.max(numeratorMaxExps[k], numeratorExps[k][regionOffset + m]);
          numeratorMaxExps[k + hmm_numMotifs] = Math.max(numeratorMaxExps[k + hmm_numMotifs], 
              numeratorExps[k + hmm_numMotifs][regionOffset + m]);
        }
      }
    }
    
    double[] numeratorSums = new double[2*hmm_numMotifs];
    double[] denominatorSums = new double[2*hmm_numMotifs];
    for (int k = 0; k < hmm_numMotifs; k++) {
      numeratorSums[k] = Numerical.log_sum_exp(numeratorExps[k], numeratorMaxExps[k]);
      denominatorSums[k] = Numerical.log_sum_exp(denominatorExps[k], denominatorMaxExps[k]);

      numeratorSums[k + hmm_numMotifs] = 
        Numerical.log_sum_exp(numeratorExps[k + hmm_numMotifs], 
            numeratorMaxExps[k + hmm_numMotifs]);
      
      denominatorSums[k + hmm_numMotifs] = 
        Numerical.log_sum_exp(denominatorExps[k + hmm_numMotifs], 
            denominatorMaxExps[k + hmm_numMotifs]);

      
      if (coupled) {
        double combinedNumerator = Numerical.log_add(numeratorSums[k], numeratorSums[k + hmm_numMotifs]);
        double combinedDenominator = Numerical.log_add(denominatorSums[k], denominatorSums[k + hmm_numMotifs]);
        double newBindingFreq = Math.exp(combinedNumerator - combinedDenominator);
        
        //don't let the bindingFreq become 0...
        newBindingFreq = Math.max(newBindingFreq, minSignalFreq);
        //warn if the binding frequency is getting very low
        if (newBindingFreq < 2) {
          logger.warn("Binding Frequency updated to " + newBindingFreq 
              + " for motif " + k + ".");
        }
        hmm_binding.setQuick(hmm_motifFwdStates[k], newBindingFreq);
        hmm_binding.setQuick(hmm_motifFwdStates[k] + 1, newBindingFreq);
        
      }
      else {
        double newBindingFreqFwd = Math.exp(numeratorSums[k] - denominatorSums[k]);
        double newBindingFreqRev = 
          Math.exp(numeratorSums[k + hmm_numMotifs] - denominatorSums[k + hmm_numMotifs]);
     
        //don't let the bindingFreq become 0...
        newBindingFreqFwd = Math.max(newBindingFreqFwd, minSignalFreq);
        newBindingFreqRev = Math.max(newBindingFreqRev, minSignalFreq);

        //warn if the binding frequency is getting very low
        if (newBindingFreqFwd < 2) {
          logger.warn("Binding Frequency updated to " + newBindingFreqFwd + " for state " + hmm_motifFwdStates[k] + ".");
        }
        if (newBindingFreqRev < 2) {
          logger.warn("Binding Frequency updated to " + newBindingFreqRev + " for state " + (hmm_motifFwdStates[k] + 1) + ".");
        }
      
        hmm_binding.setQuick(hmm_motifFwdStates[k], newBindingFreqFwd);
        hmm_binding.setQuick(hmm_motifFwdStates[k] + 1, newBindingFreqRev);
      }
    }
  }
  
  /**
   * FIXME needs to eb updated for GHMM state model
   */
  private void batchUpdateReadAssignments() {
    int model_min = model.getMin();
    int model_max = model.getMax();
    
    for (int i = 0; i < analysisRegions.size(); i++) {
      HMMAnalysisRegion currentRegion = analysisRegions.get(i);
      int currentRegionWidth = currentRegion.getWidth();
      int currentRegionStart = currentRegion.getGenomicRegion().getStart();
      DoubleMatrix2D logGamma = currentRegion.getLogGamma();
      
      List<ReadHit> hits = currentRegion.getSignalReads();
      int[] newAssignments = new int[hits.size()];
      for (int n = 0; n < hits.size(); n++) {
        ReadHit currentHit = hits.get(n);
        //calculate the range of locations to try
        int mStart, mEnd;
        if (currentHit.getStrand() == '+') {
          mStart = currentHit.getFivePrime() - currentRegionStart + model_min;
          mEnd = currentHit.getFivePrime() - currentRegionStart + model_max;
        }
        else {
          //subtract the min and max in order to keep mStart as the lesser value
          mStart = currentHit.getFivePrime() - currentRegionStart - model_max;
          mEnd = currentHit.getFivePrime() - currentRegionStart - model_min;
        }
        //ensure the values are valid
        mStart = Math.max(0, mStart);
        mEnd = Math.min(currentRegionWidth - 1, mEnd);
        
        int maxLoc = -1;
        double maxLocSum = Double.POSITIVE_INFINITY;
        int assignedLoc = currentRegion.getReadAssignment(currentHit);

        for (int m = mStart; m <= mEnd; m++) {
          double[] sumExps = new double[hmm_numStates];
          double maxExp = Double.NEGATIVE_INFINITY;
          double logReadProb = currentRegion.computeLogReadProb(currentHit, m);
          for (int k = 0; k < hmm_numStates; k++) {
            double logBindingProb;
            if (assignedLoc != m) {
              logBindingProb = this.computeLogBindingProb(currentRegion, m, k, currentRegion.numReadsAssignedTo(m) + 1);
            }
            else {
              logBindingProb = this.computeLogBindingProb(currentRegion, m, k, currentRegion.numReadsAssignedTo(m));
            }
            /**
             * logReadProb and logBindingProb will both be negative, so we can't
             * take the log of them directly as would be required to keep logGamma
             * in log-space. Instead, we take the log of the negative of them,
             * but the sign of the sum of the components of sumExps be positive
             * instead of negative. So, we need to find the minimum sum instead
             * of the maximum sum.             
             */
            sumExps[k] = logGamma.getQuick(m, k) + Math.log(-(logReadProb + logBindingProb));
            maxExp = Math.max(maxExp, sumExps[k]);
          }
          double sum = Numerical.log_sum_exp(sumExps, maxExp);
          if (sum < maxLocSum) {
            maxLoc = m;
            maxLocSum = sum;
          }
        }
        newAssignments[n] = maxLoc;
      }
      /**
       * update all the assignments to the new assignments
       * 
       * Note: If only a few read assignments change it should be quickest to
       * update them individually. If many read assignments are changing
       * it would be quicker to do some kind of batch update, but that will
       * probably only happen on the first handful of iterations
       */
      logger.debug(Arrays.toString(newAssignments));      
      
      for (int n = 0; n < hits.size(); n++) {
        ReadHit currentHit = hits.get(n);
        currentRegion.assignRead(currentHit, newAssignments[n]);
      }      
      currentRegion.printNumReadsByLoc();
    }
  }
    
   
  /**
   * 
   */
  private void nonBatchUpdateReadAssignments(boolean isStochastic) {
    int model_min = model.getMin();
    int model_max = model.getMax();
    
    for (int i = 0; i < analysisRegions.size(); i++) {
      HMMAnalysisRegion currentRegion = analysisRegions.get(i);
      int currentRegionWidth = currentRegion.getWidth();
      int currentRegionStart = currentRegion.getGenomicRegion().getStart();
      DoubleMatrix2D logGamma = currentRegion.getLogGamma();
      int[][] mkReadCounts = regionToMKReadCounts.get(currentRegion);
      
      
      ArrayList<ReadHit> hits = new ArrayList<ReadHit>(currentRegion.getSignalReads());

      Collections.shuffle(hits, randomSrc);
      int[] newAssignments = new int[hits.size()];
      for (int n = 0; n < hits.size(); n++) {
        ReadHit currentHit = hits.get(n);
        //calculate the range of locations to try
        int mStart, mEnd;
        if (currentHit.getStrand() == '+') {
          mStart = currentHit.getFivePrime() - currentRegionStart + model_min;
          mEnd = currentHit.getFivePrime() - currentRegionStart + model_max;
        }
        else {
          //subtract the min and max in order to keep mStart as the lesser value
          mStart = currentHit.getFivePrime() - currentRegionStart - model_max;
          mEnd = currentHit.getFivePrime() - currentRegionStart - model_min;
        }
        //ensure the values are valid
        mStart = Math.max(0, mStart);
        mEnd = Math.min(currentRegionWidth - 1, mEnd);
        
        double[] weights = new double[mEnd - mStart + 1];
        double[] weightRunningSum = new double[mEnd - mStart + 1];        
        int maxLoc = -1;
        double maxLocWeight = Double.NEGATIVE_INFINITY;        
          
        
        double[] sumExps = new double[hmm_numEffectiveStates];
        double maxExp = Double.NEGATIVE_INFINITY;
        
        double[] testSumExps = new double[hmm_numEffectiveStates];  
        double testMaxExp = Double.NEGATIVE_INFINITY;
        int testMaxLoc = -1;
        double testMaxLocWeight = Double.POSITIVE_INFINITY;

        
        int assignedLoc = currentRegion.getReadAssignment(currentHit);
        double logReadProb;
        double logBindingProb;
        
        for (int m = mStart; m <= mEnd; m++) {
          
          Arrays.fill(sumExps, Double.NEGATIVE_INFINITY);
          maxExp = Double.NEGATIVE_INFINITY;
          testMaxExp = Double.NEGATIVE_INFINITY; 
          logReadProb = currentRegion.computeLogReadProb(currentHit, m);
          int count = 0;          
          for (int k = 0; k < hmm_numStates; k++) {
          	/**
          	 * include all instances of this state (starting at prior positions)
          	 * that could account for reads at this position
          	 */
          	int pStart = Math.max(0, (m - stateDurations[k] + 1));
          	for (int p = pStart; p <= m; p++) {
          		if (logGamma.get(p, k) != Double.NEGATIVE_INFINITY) {
            		if ((assignedLoc < p) || (assignedLoc > (p + stateDurations[k] - 1))) {
            			logBindingProb = this.computeLogBindingProb(currentRegion, m, k, mkReadCounts[p][k] + 1);
            		}
            		else {
            			logBindingProb = this.computeLogBindingProb(currentRegion, m, k, mkReadCounts[p][k]);
            		}
            		
            		sumExps[count] = logGamma.getQuick(p, k) + logReadProb + logBindingProb;
            		maxExp = Math.max(maxExp, sumExps[count]);
            		
            		testSumExps[count] = logGamma.getQuick(p, k) + Math.log(-(logReadProb + logBindingProb));
            		testMaxExp = Math.max(testMaxExp, testSumExps[count]);
            		
            		//old code
//                /**
//                 * logReadProb and logBindingProb will both be negative, so we can't
//                 * take the log of them directly as would be required to keep logGamma
//                 * in log-space. Instead, we take the log of the negative of them,
//                 * but the sign of the sum of the components of sumExps becomes positive
//                 * instead of negative (since they're log probabilities). So, we 
//                 * need to find the minimum sum instead of the maximum sum.             
//                 */
//                sumExps[count] = logGamma.getQuick(p, k) + Math.log(-(logReadProb + logBindingProb));
//                double test = Math.exp(logGamma.getQuick(p,k)) * (logReadProb + logBindingProb);
//                maxExp = Math.max(maxExp, sumExps[count]);
          		}
          		else {
          			sumExps[count] = Double.NEGATIVE_INFINITY;
          			testSumExps[count] = Double.NEGATIVE_INFINITY;
          		}
              count++;          		
          	}
          }
          weights[m - mStart] = Math.exp(Numerical.log_sum_exp(sumExps, maxExp));         
          if ((m - mStart) > 0) {
          	weightRunningSum[m - mStart] = weightRunningSum[m - mStart - 1] + weights[m - mStart];
          }
          else {
          	weightRunningSum[m - mStart] = weights[m - mStart];
          }
          if (weights[m - mStart] > maxLocWeight) {
            maxLoc = m; //note: m _is_ a 0-based offset from the start of the region
            maxLocWeight = weights[m - mStart];
          }
          
          double testWeight = Numerical.log_sum_exp(testSumExps, testMaxExp);          
          if (testWeight < testMaxLocWeight) {
            testMaxLoc = m;
            testMaxLocWeight = testWeight;
          }
          
        }
//        if (testMaxLoc != maxLoc) {
//          logger.error("max loc and test max loc are not equal");
//        }
        if (isStochastic) {
        	double locRand = randomSrc.nextDouble() * weightRunningSum[mEnd - mStart];
        	int newLoc = 0;
        	while ((weightRunningSum[newLoc] < locRand) && (newLoc < weightRunningSum.length)) {
        	  newLoc++;
        	}
          newAssignments[n] = newLoc + mStart; 
          
        }
        else {
        	newAssignments[n] = maxLoc;
        }
        
        if (assignedLoc != newAssignments[n]) {
          this.updateMKReadCounts(mkReadCounts, assignedLoc, newAssignments[n]);
          currentRegion.assignRead(currentHit, newAssignments[n]);
        }
      }
 
      logger.debug(Arrays.toString(newAssignments));      
      
//      for (int n = 0; n < hits.size(); n++) {
//        ReadHit currentHit = hits.get(n);
//        currentRegion.assignRead(currentHit, newAssignments[n]);
//      }      
      currentRegion.printNumReadsByLoc();
    }
  }
  
  
  /**
   * 
   * @param stateStartPos the 0-based start position for the next state
   * @param regionSize the number of bases in the region
   * @return
   */
  private DoubleMatrix2D getAdjustedTransProbs(int stateStartPos, int regionSize) {
    int maxStateLength = regionSize - stateStartPos;
    if (maxStateLength >= maxStateDuration) {
      return hmm_LogA;
    }
    else {
      /**
       * This Map is set up so that the available duration should be used as the
       * key to find the appropriate adjusted probabilities.
       */
      return hmm_LogAdjustedTransProb.ceilingEntry(maxStateLength).getValue();
    }    
  }
  
  
  private void checkLikelihood(DoubleMatrix2D logAlpha, DoubleMatrix2D logBeta) {
  	double logLikelihood = Numerical.log_sum_exp(logAlpha.viewRow(0).copy().assign(logBeta.viewRow(0), Functions.plus));
  	TreeMap<Integer, Double> mismatches = new TreeMap<Integer, Double>();
  	for (int m = 1; m < logAlpha.rows(); m++) {
  		double currentLL = Numerical.log_sum_exp(logAlpha.viewRow(m).copy().assign(logBeta.viewRow(m), Functions.plus));
  		for (int j = 0; j < hmm_numStates; j++) {
  			if (isStateFixedDuration[j] && (stateDurations[j] > 1)) {
  				int start = Math.max(0, m - stateDurations[j] + 1);
  				for (int n = start; n < m; n++) {
  					currentLL = Numerical.log_add(currentLL, (logAlpha.getQuick(n, j) + logBeta.getQuick(n,j)));
  				}
  			}
  		}
  		if (Math.abs(logLikelihood - currentLL) > (double)1E-6) {
  			mismatches.put(m, currentLL);
  		}
  	}
  	
  	if (!mismatches.isEmpty()) {
  		StringBuffer errormsg = new StringBuffer();
  		errormsg.append("Log-Likelihood calculations from forward-backward algorithm are NaN or inconsistent\n" + "0: " + logLikelihood);
  		for (int pos : mismatches.keySet()) {
  			errormsg.append(pos + ": " + mismatches.get(pos) + "\n");
  		}
  		logger.error(errormsg.toString());
  		System.exit(1);
  	}
  }
  
  /**
   * Compute the log probability of a sequence of bases according to the
   * specified pwm
   * @param pwm The matrix of non-log base probabilities 
   * @param sequenceAsInts The complete sequence for a region represented as ints
   * @param firstBaseIndex The index into the sequence of the base at which to
   * start computing the probability under the pwm
   * @return
   */
  public double computeLogSequenceProb(DenseDoubleMatrix2D pwm, int[] sequenceAsInts, int firstBaseIndex) {
  	double logSequenceProb = 0;
  	for (int i = 0; i < pwm.rows(); i++) {
  		logSequenceProb += Math.log(pwm.getQuick(i, sequenceAsInts[firstBaseIndex + i]));  		
  	}
  	return logSequenceProb;
  }
  
  
  /**
   * Compute the poisson probability of a binding event with the specified 
   * number of hits
   * @param numhits
   * @param bindingFreq
   * @return
   */
  public double computeLogBindingProb(double bindingFreq, int numhits) {
    double logBindingProb = Numerical.poissonLogPDF(bindingFreq, numhits); 
    return logBindingProb;
  }
  

  /**
   * Compute the poisson probability of a binding event at the specified 
   * location with the specified number of hits, taking into account the 
   * specified state
   * @param state
   * @param loc
   * @param numhits
   * @return
   */
  public double computeLogBindingProb(HMMAnalysisRegion region, int loc, int state, int numhits) {
    if (hmm_useWCEforBGBinding && isBGState[state]) {
      double bgBindingFreq = region.getWCEBindingFreq(loc);
      double logBindingProb = Numerical.poissonLogPDF(bgBindingFreq, numhits);
      return logBindingProb;
    }
    else {
      return this.computeLogBindingProb(hmm_binding.get(state), numhits);
    }
  }
  
  /**
   * Update the cache of how many reads are assigned within the range of each 
   * state
   * @param mkReadCounts
   * @param oldLoc
   * @param newLoc
   */
  private void updateMKReadCounts(int[][] mkReadCounts, int oldLoc, int newLoc) {
    for (int k = 0; k < hmm_numStates; k++) {
      if (!isStateFixedDuration[k] || (stateDurations[k] == 1)) {
        if (mkReadCounts[oldLoc][k] > 0) {
          mkReadCounts[oldLoc][k]--;
        }
        else if (mkReadCounts[oldLoc][k] == 0) {
          logger.error("decrementing a 0 read count");
        }
        if (mkReadCounts[newLoc][k] >= 0) {
          mkReadCounts[newLoc][k]++;
        }
      }
      else {
        int mStart = Math.max(0, (oldLoc - stateDurations[k] + 1));
        for (int m = mStart; m <= oldLoc; m++) {
          if (mkReadCounts[m][k] > 0) {
            mkReadCounts[m][k]--;
          }
          else if (mkReadCounts[m][k] == 0) {
            logger.error("decrementing a 0 read count");
          }
        }
        mStart = Math.max(0, (newLoc - stateDurations[k] + 1));
        for (int m = mStart; m <= newLoc; m++) {
          if (mkReadCounts[m][k] >= 0) {
            mkReadCounts[m][k]++;         
          }
        }        
      }
    }
  }
  
  
  public List<Feature> getFeaturesFromViterbi() {
    ArrayList<Feature> features = new ArrayList<Feature>();
    
    for (HMMAnalysisRegion region : analysisRegions) {
      DoubleMatrix2D viterbiLogProb = region.getViterbiLogProb();
      int[][] viterbiTraceback = region.getViterbiTraceback();
      int[] viterbiPath = new int[region.getWidth()];
      
      /**
       * find the end state with the highest path probability, for states that
       * have duration > 1 we need to check a different offset in the 
       * viterbi matrix...
       */
      int maxEndIdx = -1;
      int maxOffset = 1;
      double maxEndProb = Double.NEGATIVE_INFINITY;      
      for (int k = 0; k < hmm_numStates; k++) {
        int offset = 1;
        if (isStateFixedDuration[k] && (stateDurations[k] > 1)) {
          offset = stateDurations[k];
        }
        double endProb = viterbiLogProb.getQuick(region.getWidth() - offset, k);
        if (endProb > maxEndProb) {
          maxEndProb = endProb;
          maxEndIdx = k;
          maxOffset = offset;
        }
      }
      viterbiPath[region.getWidth() - 1] = maxEndIdx;
      //if it's a state with duration > 1 then skip back that far
      for (int n = region.getWidth() - 2; n >= region.getWidth() - maxOffset; n--) {
        viterbiPath[n] = viterbiPath[region.getWidth()-1];
      }

      /**       
       * trace back and create a feature wherever the motif state was entered
       */
      int m = region.getWidth() - maxOffset;
      while (m > 0) {
        viterbiPath[m-1] = viterbiTraceback[m][viterbiPath[m]];
        int stateDuration = stateDurations[viterbiPath[m-1]];
        
        //if it's a motif state add a feature
        if (isMotifState[viterbiPath[m-1]]) {
          int motifIndex = stateToMotifMap.get(viterbiPath[m-1]);
          int motifRegionStart = region.getGenomicRegion().getStart() + m - stateDuration;
          int motifRegionEnd = motifRegionStart + hmm_motifLengths[motifIndex] - 1;
          
          EnrichedFeature motifPeakFeature = 
            new EnrichedFeature(new Region(region.getGenomicRegion().getGenome(), 
                region.getGenomicRegion().getChrom(),
                motifRegionStart,
                motifRegionEnd));
          
          int signalHits = 0;
          for (int i = 0; i < hmm_motifLengths[motifIndex] ; i++) {
            signalHits += region.numReadsAssignedTo(m - stateDuration + i);
          }
          motifPeakFeature.signalTotalHits = signalHits;
          if (isMotifFwdState[viterbiPath[m-1]]) {
            motifPeakFeature.strand = '+';
          }
          else {
            motifPeakFeature.strand = '-';
          }
          features.add(motifPeakFeature);
        }
        
        //if it's a state with duration > 1 then skip back that far
        for (int n = m - 2; n >= m - stateDuration; n--) {
          viterbiPath[n] = viterbiPath[m-1];
        }
        m = m - stateDuration;        
      }
      //logger.debug("Region: " + region.getGenomicRegion().toString() + ", Viterbi Path: " + Arrays.toString(viterbiPath));      
    }
    
    logger.debug(features.size() + " Features:\n" + features.toString());
    return features;
  }
     
  
  /**
   * 
   * @param base
   * @return
   */
  public static int baseToIndex(char base) {
    switch(base) {
      case 'a': return 0;
      case 'A': return 0;
      case 'c': return 1;
      case 'C': return 1;
      case 'g': return 2;
      case 'G': return 2;
      case 't': return 3;
      case 'T': return 3;
      default: throw new IllegalArgumentException("Invalid Base");
    }
  }
  
  
  /**
   * 
   * @param base
   * @return
   */
  public static int baseIndexToRevCompIndex(int baseIndex) {
    switch(baseIndex) {
    	case 0: return 3;
    	case 1: return 2;
    	case 2: return 1;
    	case 3: return 0;
      default: throw new IllegalArgumentException("Invalid Base Index");
    }
  }
  
  
  private void printPi(){
    logger.debug("logPi:\n" + hmm_LogPi.toString() + "\n\n");
  }
  
  private void printTransitionProbabilities(){
    logger.debug("log tranistion prob.\n" + hmm_LogA.toString() + "\n\n");    
  }
  
  private void printSequenceEmission(){
    for (int i = 0; i < hmm_phi.length; i++) {
      logger.debug("Sequence Emission for state: " + i + "\n" + hmm_phi[i].toString() + "\n\n");
    }
  }
  
  private void printBindingFrequencies(){
    logger.debug("Binding Frequencies:\n" + hmm_binding.toString() + "\n\n");
  }
  
  private void printReadAssignments(){
    StringBuffer readAssignSB = new StringBuffer();
    readAssignSB.append("Read Assignments:\n");
    for (HMMAnalysisRegion region : analysisRegions) {
      readAssignSB.append(region.getGenomicRegion().toString() + "\n");
      String chrom = region.getGenomicRegion().getChrom();
      int regionStart = region.getGenomicRegion().getStart();
      TreeMap<Integer, Integer> numReadsByLoc = region.getNumReadsByLoc();
      Set<Integer> keys = numReadsByLoc.keySet();
      for (Integer loc : keys) {
        int genLoc = loc + regionStart;
        readAssignSB.append(chrom + ":" + genLoc + ":" + numReadsByLoc.get(loc) + "\n");
      }
      readAssignSB.append("\n\n");
    }
  }
  
  private void printMLStateSequence(){
  //TODO
  }
  
  
  /*
   * (non-Javadoc)
   * 
   * @see edu.mit.csail.cgs.shaun.deepseq.discovery.FeatureFinder#printError()
   */
  @Override
  public void printError() {
    // TODO Auto-generated method stub

  }

}

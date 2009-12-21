package edu.mit.csail.cgs.utils.hmm;

import java.text.*;
import java.io.*;
import java.util.*;

/** 
 * Implements a hidden Markov model and its algorithms including
 * Forward-Backward, Viterbi, K-means, Baum-Welch, and
 * Kullback-Leibler distance measure.
 *
 * @author		Kanav Kahol
 * @author		Troy L. McDaniel
 */
public class HMM
{
  /** number of states */
  public int numStates;

  /** size of output vocabulary */
  public int sigmaSize;
  
  /** number of dimensions */
  public int dimensions;
  
  /** output vocabulary */
  public double V[][];

  /** initial state probabilities */
  public double pi[];

  /** transition probabilities */
  public double a[][];

  /** emission probabilities */
  public double b[][];

  /**
   * Class constructor. Initializes an HMM given
   * number of states and output vocabulary size.
   *
   * @param numStates 		number of states
   * @param sigmaSize 		size of output vocabulary
   */
  public HMM(int numStates, int sigmaSize)
  {
  	// make sure that the number of states and size of
  	// vocab is at least one
  	if(numStates < 1)
  		numStates = 1;
  	if(sigmaSize < 1)
  		sigmaSize = 1;
    this.numStates = numStates;
    this.sigmaSize = sigmaSize;
    // If an HMM is created using this constructor, i.e., without
    // specifying an output vocab, the number of dimenions for
    // the observations symbols is unknown and thus is set to -1.
    // Basically, this info is not used and ignore in this case.
    this.dimensions = -1;

    pi = new double[numStates];
    a = new double[numStates][numStates];
    b = new double[numStates][sigmaSize];
  }

  /**
   * Class constructor. Initializes an HMM a given
   * number of states, output vocabulary size, and
   * number of dimensions for its observation symbols.
   *
   * @param numStates 		number of states
   * @param sigmaSize 		size of output vocabulary
   * @param dimensions		number of dimensions
   */
  public HMM(int numStates, int sigmaSize, int dimensions)
  {
  	// make sure that the number of states, size of
  	// vocab, and number of dimensions are at least one
  	if(numStates < 1)
  		numStates = 1;
  	if(sigmaSize < 1)
  		sigmaSize = 1;
  	if(dimensions < 1)
  		dimensions = 1;
    this.numStates = numStates;
    this.sigmaSize = sigmaSize;
    this.dimensions = dimensions;

    pi = new double[numStates];
    a = new double[numStates][numStates];
    b = new double[numStates][sigmaSize];
    V = new double[sigmaSize][dimensions];
  }

  /**
   * Baum-Welch algorithm for hidden Markov models. Given an
   * observation sequence o, it will train this HMM using 
   * o to increase the probability of this HMM generating o.
   *
   * @param o 		the observation sequence
   * @param steps	the number of iterations performed
   */
  public void baumWelch(int[] o, int steps)
  { 
    int T = o.length;
    double[][] fwd;
    double[][] bwd;

    double pi1[] = new double[numStates];
    double a1[][] = new double[numStates][numStates];
    double b1[][] = new double[numStates][sigmaSize];

    for (int s = 0; s < steps; s++) {
      // calculation of Forward and Backward Variables from the current model
      fwd = forwardProc(o);
      bwd = backwardProc(o);

      // re-estimation of initial state probabilities
      for (int i = 0; i < numStates; i++)
			pi1[i] = gamma(i, 0, o, fwd, bwd);

      // re-estimation of transition probabilities
      for (int i = 0; i < numStates; i++) {
		for (int j = 0; j < numStates; j++) {
		  double num = 0;
		  double denom = 0;
		  for (int t = 0; t <= T - 1; t++) {
		    num += xi(t, i, j, o, fwd, bwd);
		    denom += gamma(i, t, o, fwd, bwd);
		  }
		  a1[i][j] = divide(num, denom);
		}
      }

      // re-estimation of emission probabilities
      for (int i = 0; i < numStates; i++) {
		for (int k = 0; k < sigmaSize; k++) {
		  double num = 0;
		  double denom = 0;
	
		  for (int t = 0; t <= T - 1; t++) {
		    double g = gamma(i, t, o, fwd, bwd);
		    num += g * (k == o[t] ? 1 : 0);
		    denom += g;
		  }
		  b1[i][k] = divide(num, denom);
		}
      }
      pi = pi1;
      a = a1;
      b = b1;
    }
  }


  /**
   * Calculation of Forward variables f(i,t) for state i at time
   * t for sequence o with the current HMM parameters.
   *
   * @param o		the observation sequence
   * @return 		a 2d array containing Forward variables
   *				f(i,t) over states and time
   */
  public double[][] forwardProc(int[] o)
  {
    int T = o.length;
    double[][] fwd = new double[numStates][T];

    // initialization (time 0)
    for (int i = 0; i < numStates; i++)
      fwd[i][0] = pi[i] * b[i][o[0]];

    // induction
    for (int t = 0; t <= T-2; t++) {
      for (int j = 0; j < numStates; j++) {
		fwd[j][t+1] = 0;
		for (int i = 0; i < numStates; i++)
		  fwd[j][t+1] += (fwd[i][t] * a[i][j]);
		fwd[j][t+1] *= b[j][o[t+1]];
      }
    }

    return fwd;
  }

  /**
   * Calculation of Backward variables b(i,t) for state i at time
   * t for sequence o with the current HMM parameters.
   *
   * @param o		the observation sequence
   * @return 		a 2d array containing Backward variables
   *				b(i,t) over states and time
   */
  public double[][] backwardProc(int[] o)
  {
    int T = o.length;
    double[][] bwd = new double[numStates][T];

    // initialization (time 0)
    for (int i = 0; i < numStates; i++)
      bwd[i][T-1] = 1;

    // induction
    for (int t = T - 2; t >= 0; t--) {
      for (int i = 0; i < numStates; i++) {
		bwd[i][t] = 0;
		for (int j = 0; j < numStates; j++)
		  bwd[i][t] += (bwd[j][t+1] * a[i][j] * b[j][o[t+1]]);
      }
    }

    return bwd;
  }

  /**
   * Calculation of xi_t(i, j), which is the probability
   * P(i_t = s_i, i_t+1 = s_j | o, hmm), that is, 
   * the probability of being in state i at time t and state j 
   * at time t+1 given observation sequence o and this HMM.
   *
   * @param t		the time
   * @param i		the number of state s_i
   * @param j		the number of state s_j
   * @param o		the observation sequence
   * @param fwd		the Forward variables for o
   * @param bwd		the Backward variables for o
   * @return		P(i_t = s_i, i_t+1 = s_j | o, hmm)
   */
  public double xi(int t, int i, int j, int[] o, double[][] fwd, double[][] bwd)
  {
    double num, denom = 0.0;
    
    // numerator
    if (t == o.length - 1)
      num = fwd[i][t] * a[i][j];
    else
      num = fwd[i][t] * a[i][j] * b[j][o[t+1]] * bwd[j][t+1];

	// denominator
    for (int k = 0; k < numStates; k++)
      denom += (fwd[k][t] * bwd[k][t]);

    return divide(num, denom);
  }

  /**
   * Calculation of gamma_t(i), which is the probability
   * P(i_t = s_i | o, hmm), that is, the probability
   * of being in state i at time given observation sequence 
   * o and this HMM.
   *
   * @param i		the number of state s_i
   * @param t		the time
   * @param o		the observation sequence
   * @param fwd		the Forward variables for o
   * @param bwd		the Backward variables for o
   * @return		P(i_t = s_i | o, hmm)
   */
  public double gamma(int i, int t, int[] o, double[][] fwd, double[][] bwd)
  {
    double num, denom = 0.0;
    
    // numerator
    num = fwd[i][t] * bwd[i][t];

	// denominator
    for (int j = 0; j < numStates; j++)
      denom += fwd[j][t] * bwd[j][t];

    return divide(num, denom);
  }
  
  
  /**
   * Viterbi algorithm for hidden Markov models. Given an
   * observation sequence o, Viterbi finds the best possible
   * state sequence for o on this HMM, along with the state
   * optimized probability.
   *
   * @param o 		the observation sequence
   * @return		a 2d array consisting of the minimum cost, U,
   *				at position (0,0) and the best possible path 
   *				beginning at (1,0). The state optimized 
   *				probability is Euler's number raised to the 
   *				power of -U.
   */
  public double[][] viterbi(int[] o)
  {
  	int T = o.length;
  	int min_state;
  	double min_weight, weight;
  	double stateOptProb = 0.0;
  	int[] Q = new int[T];
  	int[][] sTable = new int[numStates][T];
  	double[][] aTable = new double[numStates][T];
  	double[][] answer = new double[2][T];
  	
  	// calulate accumulations and best states for time 0
  	for(int i = 0; i < numStates; i++) {
  		aTable[i][0] = -1*Math.log(pi[i]) - Math.log(b[i][o[0]]);
  		sTable[i][0] = 0;
  	}
  	
  	// fill up the rest of the tables
  	for(int t = 1; t < T; t++) {
  		for(int j = 0; j < numStates; j++) {
  			min_weight = aTable[0][t-1] - Math.log(a[0][j]);
  			min_state = 0;
  			for(int i = 1; i < numStates; i++) {
  				weight = aTable[i][t-1] - Math.log(a[i][j]);
  				if(weight < min_weight) {
  					min_weight = weight;
  					min_state = i;
  				}
  			}
  			aTable[j][t] = min_weight - Math.log(b[j][o[t]]);
  			sTable[j][t] = min_state;
  		}
  	}
  	
  	// find minimum value for time T-1
  	min_weight = aTable[0][T-1];
  	min_state = 1;
  	for(int i = 1; i < numStates; i++) {
  		if(aTable[i][T-1] < min_weight) {
  			min_weight = aTable[i][T-1];
  			min_state = i;
  		}
  	}

  	// trace back to find the optimized state sequence
  	Q[T-1] = min_state;
  	for(int t = T-2; t >= 0; t--)
  		Q[t] = sTable[Q[t+1]][t+1];
  	
  	// store answers and return them
  	answer[0][0] = min_weight;
  	for(int i = 0; i < T; i++)
  		answer[1][i] = Q[i];
  	return answer;
  }
  
  /** 
   * K-means algorithm for hidden Markov models. Given training
   * data, this algorithm will classify the data into a specified
   * number of states, then using the final classification, calculate
   * initial, transition, and bias probabilities.
   *
   * @param data		  the filename of training data file
   * @param N			  the number of states
   * @param steps		  the number of iterations
   * @throws IOException  If training data file is invalid
   */
  public static HMM kmeansLearner(String data, int N, int steps) throws IOException {
  	HMM newHmm = null;
  	File trainingFile = new File(data);
  	FileReader in = new FileReader(trainingFile);
  	int D, w, T, numSymbols, symbolIndex = 0;
  	int[] count;
  	int[][] classification;
  	double[] current;
  	double[][] means, new_means, symbolSet;
  	double[][][] trainingData;
  	
  	/* read training data starting with header information */
  	// D
  	char c;
  	String entry = "";
  	c = (char)in.read();
	do {
		entry = entry.concat("" + c);
		c = (char)in.read();
	}while(!Character.isWhitespace(c));
	D = (new Integer(entry)).intValue();
	entry = "";
	
  	// T
  	c = (char)in.read();
	do {
		entry = entry.concat("" + c);
		c = (char)in.read();
	}while(!Character.isWhitespace(c));
	T = (new Integer(entry)).intValue();
	entry = "";
  	
  	// w
  	c = (char)in.read();
	do {
		entry = entry.concat("" + c);
		c = (char)in.read();
	}while(!Character.isWhitespace(c));
	w = (new Integer(entry)).intValue();
	entry = "";
	
	// read training data
	trainingData = new double[w][T][D];
  	for(int i = 0; i < w; i++) {
  		for(int j = 0; j < T; j++) {
  			for(int d = 0; d < D; d++) {
  				c = (char)in.read();
				do {
					entry = entry.concat("" + c);
					c = (char)in.read();
				}while(!Character.isWhitespace(c));
				trainingData[i][j][d] = (new Double(entry)).doubleValue();
				entry = "";
			}
		}
	}
	in.close();	// close file
	
	// set up matrices
	means = new double[N][D];
	current = new double[D];
	new_means = new double[N][D];
	count = new int[N];
	classification = new int[w][T];
	symbolSet = new double[w*T][D];
	
	// set initial means
	setInitialMeans(means, trainingData, D, w, T, N);
			
	// classify training data
	boolean done = false;
	int iterations = 0, identifiedClass;
	double dist, min_dist;
	while(!done) {
		for(int i = 0; i < w; i++) {
			for(int j = 0; j < T; j++) {
				for(int d = 0; d < D; d++)
					current[d] = trainingData[i][j][d];
				min_dist = distance(current, means, 0);
				identifiedClass = 0;
				for(int k = 1; k < N; k++) {
					dist = distance(current, means, k);
					if(dist < min_dist) {
						min_dist = dist;
						identifiedClass = k;
					}
				}
				classification[i][j] = identifiedClass;
			}
		}
				
		// calculate the new means
		for(int i = 0; i < N; i++) {
			for(int d = 0; d < D; d++)
				new_means[i][d] = 0;
			count[i] = 0;
		}
		for(int i = 0; i < w; i++) {
			for(int j = 0; j < T; j++) {
				for(int d = 0; d < D; d++)
					new_means[classification[i][j]][d] += trainingData[i][j][d];
				count[classification[i][j]]++;
			}
		}
		for(int k = 0; k < N; k++)
			for(int d = 0; d < D; d++)
				new_means[k][d] = divide(new_means[k][d], count[k]);

		// calculate the RMS error between the old and new means
		double total_error = 0.0;
		for(int i = 0; i < N; i++) {
			double RMS = 0.0;
			
			for(int d = 0; d < D; d++)
				RMS += (means[i][d] - new_means[i][d]) * (means[i][d] - new_means[i][d]);
			total_error += Math.sqrt(RMS);
		}
		
		// check error
		if(total_error < 0.01 || ++iterations >= steps)
			done = true;
		else {
			// zero out arrays, copy new means into means, and repeat
			for(int i = 0; i < N; i++)
				for(int d = 0; d < D; d++)
					means[i][d] = 0;
			for(int i = 0; i < N; i++)
				for(int d = 0; d < D; d++)
					means[i][d] = new_means[i][d];
			for(int i = 0; i < w; i++)
				for(int j = 0; j < T; j++)
					classification[i][j] = 0;
		}
	}

	// find symbols
	symbolIndex = getSymbols(symbolSet, trainingData, D, w, T);
	
	// create a new HMM
	newHmm = new HMM(N, symbolIndex, D);
	
	// copy symbol data into V
	for(int i = 0; i < symbolIndex; i++)
		for(int d = 0; d < D; d++)
			newHmm.V[i][d] = symbolSet[i][d];
	
	// calculate probabilities using final classification
	probabilitiesFromKmeans(newHmm, classification, trainingData, N, D, w, T);
	
	return newHmm;
  }
  
  /* Calculates and sets the initial means for the K-means 
   * algorithm.
   */
  private static void setInitialMeans(double[][] means, double[][][] trainingData, 
                                      int D, int w, int T, int N)
  {
	double[][] sortedarray = new double[w*T][3];
	int baseIndex;
	Random rand = new Random();
	
	/* Construct the initial sortedarray (currently unsorted) by storing each
	   training data point as a square root of a sum of its elements squared.
	   Note that its position within the file is stored. */
	baseIndex = 0;	// multiples of T
	for(int i = 0; i < w; i++) {
  		for(int j = 0; j < T; j++) {
  			double val = 0.0;
  			for(int d = 0; d < D; d++)
  				val += (trainingData[i][j][d]) * (trainingData[i][j][d]);
  			sortedarray[baseIndex + j][0] = Math.sqrt(val);
  			sortedarray[baseIndex + j][1] = i;
  			sortedarray[baseIndex + j][2] = j;
  		}
  		baseIndex += T;
  	}
    
    // Use insertion sort to sort the training data points. Coordinates
  	// are kept track of during the sorting process.
	for(int j = 1; j < w*T; j++) {
		int i;
		double coord1, coord2, key;
		key = sortedarray[j][0];
		coord1 = sortedarray[j][1];
		coord2 = sortedarray[j][2];
		i = j-1;
		while( i > -1 && sortedarray[i][0] > key) {
			sortedarray[i+1][0] = sortedarray[i][0];
			sortedarray[i+1][1] = sortedarray[i][1];
			sortedarray[i+1][2] = sortedarray[i][2];
			i = i - 1;
		}
		sortedarray[i+1][0] = key;
		sortedarray[i+1][1] = coord1;
		sortedarray[i+1][2] = coord2;
	}
	
	/* Using sortedarray, randomly select values within segments of length
	   w. When a value is selected, we extract the coordinates of this value
	   to find the original data point within our training data file. This
	   point is selected as an initial mean. */
	baseIndex = 0;	// multiples of w
	for(int i = 0; i < N; i++) {
		int randVal = rand.nextInt(w);
		int c1 = (int)sortedarray[baseIndex + randVal][1];
		int c2 = (int)sortedarray[baseIndex + randVal][2];
		for(int d = 0; d < D; d++)
			means[i][d] = trainingData[c1][c2][d];
		baseIndex += w;
		if(baseIndex > w*T-1)
			baseIndex = 0;
	}
  }
  	
  /* Calculates the initial, transition, and bias probabilities for
   * the K-means algorithm.
   */
  private static void probabilitiesFromKmeans(HMM newHmm, int[][] classification, 
	                                            double[][][] trainingData, int N, int D, 
	                                            int w, int T)
  {
  	// calculate initial probabilities
	for(int i = 0; i < N; i++) {
		int curr_count = 0;
		for(int j = 0; j < w; j++)
			if(classification[j][0] == i)
				curr_count++;
		newHmm.pi[i] = divide(curr_count, w);
	}
	
	// calculate tranisition probabilities
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			int curr_count = 0;
			int total_count = 0;
			for(int k = 0; k < w; k++) {
				for(int t = 0; t < T-1; t++) {
					if(classification[k][t] == i && classification[k][t+1] == j)
						curr_count++;
					if(classification[k][t] == i)
						total_count++;
				}
				if(classification[k][T-1] == i)
					total_count++;
			}
			newHmm.a[i][j] = divide(curr_count, total_count);
		}
	}
	
	// calculate emission probabilities
	int[] count = new int[newHmm.sigmaSize];
	double[] O = new double[D];
	for(int i = 0; i < N; i++) {
		int total_count = 0;
		for(int a = 0; a < newHmm.sigmaSize; a++)
			count[a] = 0;
		for(int j = 0; j < w; j++) {
			for(int t = 0; t < T; t++) {
				if(classification[j][t] == i) {
					int k = 0;
					boolean found = false;
					for(int d = 0; d < D; d++)
						O[d] = trainingData[j][t][d];
					for(int b = 0; b < newHmm.sigmaSize && !found; b++)
						if(euclidDistance(O, newHmm.V, D, b) < 0.0001) {
							k = b;
							found = true;
						}
					if(found)
						count[k]++;
					total_count++;
				}
			}
		}
		for(int k = 0; k < newHmm.sigmaSize; k++)
			newHmm.b[i][k] = divide(count[k], total_count);
	}
  }
  
  /* Determines all (distinct) symbols used by this HMM given a training set.
   * Symbols are stored in a symbol set, which is used to construct V, the
   * output vocabulary of this HMM. The set is initialized with negative ones,
   * so if a symbol in your model is represented by all negatives ones, you
   * should change the initialization value.
   */
  private static int getSymbols(double[][] symbolSet, double[][][] trainingData,
								  int D, int w, int T)
  {
	double[] symbol = new double[D];
	int symbolIndex = 0;
	
	// initialize the symbol set
	// change if a symbol can consist of all negative ones
	for(int i = 0; i < w*T; i++)
		for(int d = 0; d < D; d++)
			symbolSet[i][d] = -1.0;
	// search for new symbols and each time a new symbol is found, add
	// it to the symbol set and increment the number of symbols
	for(int i = 0; i < w; i++) {
		for(int j = 0; j < T; j++) {
			for(int d = 0; d < D; d++)
				symbol[d] = trainingData[i][j][d];
			if(!symbolCheck(symbol, symbolSet, D, w, T)) {
				for(int d = 0; d < D; d++)
					symbolSet[symbolIndex][d] = symbol[d];
				symbolIndex++;
			}
		}
	}
	return symbolIndex;	// return the number of symbols
  }
  
  /**
   * Calculates the distance between two hidden Markov models using
   * the Kullback-Leibler distance measure. The models should have 
   * the same size vocabulary.
   *
   * @param hmm2		the hidden Markov model that will be compared
   *					against this HMM
   * @return			the symmetric distance between this HMM and hmm2
   *					(if their vocabulary size differs, -1 is returned)
   */
  public double hmmDistance(HMM hmm2)
  {
  	int[] o1;
  	int[] o2;
  	int s = sigmaSize;
  	double[][] weight11, weight21, weight12, weight22;
  	double d12, d21;
  	Random rand = new Random();
  	
  	// vocabulary size must be equal since sequences are generated from
  	// one HMM and checked on the other
  	if(sigmaSize != hmm2.sigmaSize)
  		return -1;
  	
  	/* Generate observation sequences from both HMMs. Note that the
  	   length is 4000 since research has found that the distance
  	   stabilizes at a length of 4000.*/
  	o1 = observationGeneration(4000);
  	o2 = hmm2.observationGeneration(4000);
  	
  	// calculate minimum weights using Viterbi
  	weight11 = viterbi(o1);
  	weight21 = hmm2.viterbi(o1);
  	weight12 = viterbi(o2);
  	weight22 = hmm2.viterbi(o2);
  	
  	// compute and return the symmetric distance
  	d12 = (1.0 / 4000.0) * Math.abs((-1*weight11[0][0]) - (-1*weight21[0][0]));
  	d21 = (1.0 / 4000.0) * Math.abs((-1*weight22[0][0]) - (-1*weight12[0][0]));
  	
  	return 0.5 * (d12 + d21);
  }
  
  /**
   * Generates an observation sequence of a given length using
   * this HMM.
   *
   * @param length		the length of the observation sequence to be generated
   * @return			the generated observation sequence
   */
  public int[] observationGeneration(int length)
  {
  	Random rand = new Random();
  	double randVal, total;
  	int currentState = 0;
  	boolean done = false;
  	int[] o;
  	
  	// make sure length is valid
  	if(length < 1)
  		length = 1;
  	
  	o = new int[length];
  	
  	// determine the initial state based on intial probabilities
  	randVal = rand.nextDouble();
  	if(randVal <= pi[0]) {
  		currentState = 0;
  		done = true;
  	}
  	total = pi[0];
  	for(int i = 1; i < numStates && !done; i++) {
  		if(randVal > total && randVal <= (total + pi[i])) {
  			currentState = i;
  			done = true;
  		}
  		total += pi[i];
  	}

  	for(int i = 0; i < length; i++) {
  		// for currentState, determine the symbol seen based on 
  		// bias probabilities
  		randVal = rand.nextDouble();
  		total = 0.0;
  		done = false;
  		if(randVal <= b[currentState][0]) {
  			o[i] = 0;
  			done = true;
  		}
  		total = b[currentState][0];
  		for(int j = 1; j < sigmaSize && !done; j++) {
  			if(randVal > total && randVal <= (total + b[currentState][j])) {
  				o[i] = j;
  				done = true;
  			}
  			total += b[currentState][j];
  		}
  		
  		// determine which state to transition to based on 
  		// transition probabilities
  		randVal = rand.nextDouble();
  		total = 0.0;
  		done = false;
	  	if(randVal <= a[currentState][0]) {
	  		currentState = 0;
	  		done = true;
	  	}
	  	total = a[currentState][0];
	  	for(int j = 1; j < numStates && !done; j++) {
	  		if(randVal > total && randVal <= (total + a[currentState][j])) {
	  			currentState = j;
	  			done = true;
	  		}
	  		total += a[currentState][j];
	  	}
  	}
  	
  	return o;
  }
  
  /* Calculates the Euclidean distance between two means: the current
   * mean and mean k.
   */
  private static double distance(double[] current, double[][] means, int k)
  {
  	int D = current.length;
  	double sum = 0;
  	
  	// calculate the euclidean distance between current and mean k
  	for(int d = 0; d < D; d++)
  		sum += (current[d] - means[k][d]) * (current[d] - means[k][d]);
  	return Math.sqrt(sum);
  }
  
 /*
  * Calculates the Euclidean distance between two symbols.
  */
  private static double euclidDistance(double[] O, double[][] V, int dim, int k)
  {
	int i;
	double sum=0;

	for(i = 0; i < dim; i++)
		sum = sum + (O[i] - V[k][i])*(O[i] - V[k][i]);
	return Math.sqrt(sum);
  }
  
  /* Checks if a symbol is within the symbol set. This is used to extract
   * distinct symbols from the training set.
   */
  private static boolean symbolCheck(double[] symbol, double[][] symbolSet,
	                                 int D, int w, int T)
  {
	boolean inSet = false;
	boolean match;
	
	// determine if the symbol is already in symbol set
	for(int i = 0; i < w*T && !inSet; i++) {
		match = true;
		for(int d = 0; d < D; d++)
			if(symbol[d] != symbolSet[i][d])
				match = false;
		if(match)
			inSet = true;
	}
	
	return inSet;
  }
  
  /**
   * Converts a sequence of actually values into their corresponding
   * indices based on the vocabulary of this HMM. Any time a hidden
   * Markov model algorithm whose input involves an observation sequence
   * is used, the observation sequence should be converted to indices
   * if it's not already.
   *
   * @param sequence		the observation sequence (actual values)
   * @param D				the number of dimensions of sequence
   * @param T				the length of sequence
   * @return				the observation sequence as indices
   */
  public int[] convert(double[][] sequence, int D, int T)
  {
  	double[] o = new double[D];
  	int[] new_o = new int[T];
  	double min_dist, dist;
  	boolean found;
  	
  	for(int i = 0; i < T; i++) {
  		// extract a symbol consisting of actual values
  		for(int d = 0; d < D; d++)
  			o[d] = sequence[i][d];
  		// match the symbol with it's index
  		found = false;
  		min_dist = euclidDistance(o, V, D, 0);
		new_o[i] = 0;
		for(int k = 1; k < sigmaSize && !found; k++) {
			dist = euclidDistance(o, V, D, k);
			if(dist < min_dist) {
				min_dist = dist;
				new_o[i] = k;
			}
			if(min_dist == 0)
				found = true;
		}
  	}
  	return new_o; 
  }
  
  /* Divides two doubles.
   * 0 / 0 = 0!
   */
  private static double divide(double n, double d)
  {
    if (n == 0)
      return 0;
    else
      return n / d;
  }
  
  /**
   * Prints all the parameters of this HMM. Parameters include initial, 
   * transition, and output probabilities.
   */
  public void print()
  {
    DecimalFormat fmt = new DecimalFormat();
    fmt.setMinimumFractionDigits(3);
    fmt.setMaximumFractionDigits(3);

	// print the intial probabilities
    for (int i = 0; i < numStates; i++)
      System.out.println("pi(" + i + ") = " + fmt.format(pi[i]));
    System.out.println();

	// print the transition probabilities
    for (int i = 0; i < numStates; i++) {
      for (int j = 0; j < numStates; j++)
		System.out.print("a(" + i + "," + j + ") = " + fmt.format(a[i][j]) + "  ");
      System.out.println();
    }

	// print the emission probabilities
    System.out.println();
    for (int i = 0; i < numStates; i++) {
      for (int k = 0; k < sigmaSize; k++)
		System.out.print("b(" + i + "," + k + ") = " + fmt.format(b[i][k]) + "  ");
      System.out.println();
    }
  }
  
  /**
   * Loads an HMM from a file. A new HMM is created and the read values
   * are stored as the HMM parameters. Returned is the new HMM.
   *
   * @param filename		name of file to read the HMM from
   * @returns				HMM loaded from the given file
   * @throws IOException	if the file format is incorrect
   */
  public static HMM load(String filename) throws IOException
  {
  	HMM newHmm;
  	File inputFile = new File(filename);
  	FileReader in = new FileReader(inputFile);
  	char c;
	int i, j, numStatesTemp = 0, sigmaSizeTemp = 0, dimensionsTemp = -1;
  	String temp = "";
  	
  	// read the number of states, vocab size, and dimensions 
  	// (dimensions will be equal to -1 if it wasn't used)
  	for(i = 0; i < 3; i++) {
	  	c = (char)in.read();
		do {
			temp = temp.concat("" + c);
			c = (char)in.read();
		}while(!Character.isWhitespace(c));
		switch(i) {
			case 0:
				numStatesTemp = (new Integer(temp)).intValue();
				break;
			case 1:
				sigmaSizeTemp = (new Integer(temp)).intValue();
				break;
			case 2:
				dimensionsTemp = (new Integer(temp)).intValue();
				break;
			default:
				break;
		}
		temp = "";
	}
	
	// create a new HMM
	if(dimensionsTemp == -1)
		newHmm = new HMM(numStatesTemp, sigmaSizeTemp);
	else
		newHmm = new HMM(numStatesTemp, sigmaSizeTemp, dimensionsTemp);
	
	// read the initial probabilities
	for(i = 0; i < newHmm.numStates; i++) {
		c = (char)in.read();
		do {
			temp = temp.concat("" + c);
			c = (char)in.read();
		}while(!Character.isWhitespace(c));
		newHmm.pi[i] = (new Double(temp)).doubleValue();
		temp = "";
	}
	
	// read the transition probabilities
	for(i = 0; i < newHmm.numStates; i++) {
		for(j = 0; j < newHmm.numStates; j++) {
			c = (char)in.read();
			do {
				temp = temp.concat("" + c);
				c = (char)in.read();
			}while(!Character.isWhitespace(c));
			newHmm.a[i][j] = (new Double(temp)).doubleValue();
			temp = "";
		}
	}
	
	// read the bias probabilities
	for(i = 0; i < newHmm.numStates; i++) {
		for(j = 0; j < newHmm.sigmaSize; j++) {
			c = (char)in.read();
			do {
				temp = temp.concat("" + c);
				c = (char)in.read();
			}while(!Character.isWhitespace(c));
			newHmm.b[i][j] = (new Double(temp)).doubleValue();
			temp = "";
		}
	}
	
	// read the output vocabulary if it was used
	// (dimensions == -1 means that the vocab wasn't used)
	if(newHmm.dimensions != -1) {
		for(i = 0; i < newHmm.sigmaSize; i++) {
			for(j = 0; j < newHmm.dimensions; j++) {
				c = (char)in.read();
				do {
					temp = temp.concat("" + c);
					c = (char)in.read();
				}while(!Character.isWhitespace(c));
				newHmm.V[i][j] = (new Double(temp)).doubleValue();
				temp = "";
			}
		}
	}
	
	in.close();		// close stream
	return newHmm;
  }
  
  /**
   * Writes this HMM to a file. Information written includes the
   * number of states, output vocabulary size, number of dimensions
   * (-1 if not used), initials, transitions, biases, and, if it was
   * used, the output vocabulary.
   *
   * @param filename		name of file to write this HMM to
   * @throws IOException	if a problem occurred while writing to the file
   */
  public void write(String filename) throws IOException
  {
  	File outputFile = new File(filename);
  	FileWriter out = new FileWriter(outputFile);
  	BufferedWriter bw = new BufferedWriter(out);
  	int i, j;
  	
  	// write the number of states, vocab size, and dimensions 
  	// (dimensions will be equal to -1 if it wasn't used)
  	bw.write(numStates + " " + sigmaSize + " " + dimensions);
  	bw.newLine();
  	
  	// write initial probabilities
  	for(i = 0; i < numStates-1; i++)
  		bw.write(pi[i] + " ");
  	bw.write(pi[i] + "");
  	bw.newLine();
  	
  	// write transition probabilities
  	for(i = 0; i < numStates; i++) {
  		for(j = 0; j < numStates; j++) {
  			if(i == numStates-1 && j == numStates-1) {
  				bw.write(a[i][j] + "");
  				bw.newLine();
  			}
  			else
  				bw.write(a[i][j] + " ");
  		}
  	}
  	
  	// write bias probabilities
  	for(i = 0; i < numStates; i++) {
  		for(j = 0; j < sigmaSize; j++) {
  			if(i == numStates-1 && j == sigmaSize-1) {
  				bw.write(b[i][j] + "");
  				bw.newLine();
  			}
  			else
  				bw.write(b[i][j] + " ");
  		}
  	}
  	
  	// write actual vocab if it was used
  	// (dimensions == -1 means that the vocab wasn't used)
  	if(dimensions != -1) {
  		for(i = 0; i < sigmaSize; i++) {
	  		for(j = 0; j < dimensions; j++) {
	  			if(i == sigmaSize-1 && j == dimensions-1) {
	  				bw.write(V[i][j] + "");
	  				bw.newLine();
	  			}
	  			else
	  				bw.write(V[i][j] + " ");
  			}
  		}
  	}
  	
  	bw.close();	// close stream
  }
}

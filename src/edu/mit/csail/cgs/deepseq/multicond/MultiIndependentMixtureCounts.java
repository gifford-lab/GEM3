package edu.mit.csail.cgs.deepseq.multicond;

import java.io.File;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.TreeSet;

import edu.mit.csail.cgs.deepseq.BindingModel;
import edu.mit.csail.cgs.utils.Utils;
import edu.mit.csail.cgs.utils.stats.StatUtil;
import edu.mit.csail.cgs.utils.Pair;

/**
 * Data points are inserted as <tt>(position, count)</tt> pairs.   <br>
 * In that respect, data points are assumed unique in terms of their position.  <br>
 * In other words, we assume that reads corresponding to the same position will be the
 * result of the same binding event.											<br>
 * Also, here components are declared specifically instead of being all the components
 * corresponding to the region of length <tt>M</tt>.
 * @author gio_fou
 *
 */
public class MultiIndependentMixtureCounts {
	
	/**************
	 **  PARAMS  **
	 *************/

	/** Total number of data points. */
	private int Ntot;
	
	/** Number of data points for each condition. */
	private int[] N;
	
	/** Number of unique positions for each condition. */
	private int[] V;	
	
	/** Components' positions */
	private int[] compPos;
	
	/** non-zero components' positions  */
	private int[] nonzeroCompIdx;
	
	/** Number of data points for control condition (if any). */
	//private int Nc;
	
	/** Number of conditions. Should be equal to N.length . */
	private int C;
	
	/** Number of components. Should be the same for all conditions. */
	private int M;
	
	/** Number of possible values that each read can take. */
	private int O;
	
	/** Determines which version to run; 													   <br>
	 * 1: HARD SERIAL, global EM is run first and then ML EM on each condition, separately,
	 * taking into account only the non-zero components from the global solution.   <br>
	 * 2: SOFT SERIAL, global EM is run first and then MAP EM on each condition, separately,
	 * taking into account the global prior weights as pseudo-counts.           <br>
	 * PARALLEL, global and condition-specific EM are run at the same time.
	 */
	private int trainVersion = 1;
	
	/** Global log-likelihood for all the conditions. */
	private double glob_loglik;
	
	/** Global log-likelihoods for each iteration. */
	private List<Double> glob_logliks;
	
	private double[] loglik;
	
	private List<Double>[] logliks;
	
	/** Current number of iterations. */
	private int numIters;
	
	/** Shows that the algorithm has converged w.r.t. to the global independent mixture. */
	private boolean hasConverged;
	
	/** The component with the minimum summed responsibility and its position */
	private Pair<Double, TreeSet<Integer>> minResp_minIndex;
	
    /**
     * - <tt>curr_alpha = 0</tt>  up to iteration <tt>ML_maxIters</tt>. 		<br>
     *   That is: <tt>iter < ML_maxIters</tt>		<br>
     * - <tt>curr_alpha = alpha*(iter - ML_maxIters)/(anneal_maxIters - ML_maxIters)</tt>  up to iteration <tt>anneal_maxIters</tt>.  <br>
     *   That is: <tt>ML_maxIters <= iter < anneal_maxIters</tt>					    	<br>
     * - <tt>curr_alpha = alpha</tt>   for <tt>iter > anneal_maxIters</tt>
     */
    private double curr_alpha;
    
    /** The initialization of the global prior weighting. */
    private double[] glob_prior_weight_0;
    
    /** The initialization of the prior weighting for each condition. */
    private double[][] prior_weight_0;

	 ///////////////////////
	 //  SET BY THE USER  //
	 ///////////////////////
    
    /** The OxM (transposed) emission probability matrix: p(x_n = o | z_n = j). <br>
     * Each column sums up to 1.												<br>
     * Here: O = M, that is, the number (and essentially the range of) possible
     * values that observed x_n can take equals the number of mixture components, M. */
    private double[][] emit_mat_transp;   // O x M = M x M  emission matrix
    
    private BindingModel bm;

    /** Printing step step for outputting information during training. */
    private int print_step = 50;
	
	/**
	 * If <tt>check_increased == true</tt>, it also checks whether the log-likelihood 
	 * (compared to that of the previous iteration) has been decreased. */
	private boolean check_increased = true;

	/** Maximum number of iterations of the EM. */
	private int maxIters = 500;

	/** Up to this value (<tt>ML_maxIters</tt>), we run ML EM (without the sparse prior). */
	private int ML_maxIters = 10;

	/** We run MAP EM up to this step with the alpha given at each iteration by the formula:   <br>
	 *  curr_alpha = alpha*(iter - ML_maxIters)/(anneal_maxIters - ML_maxIters)  */
	private int anneal_maxIters = 80;

	/** Flag that specifies if all the sparse prior will act on all components or
	 *  only in the worst (in terms of lowest responsibility) one.  			 */
	private boolean isBatchElimOn = true;
	
	private boolean mask_glob_prior_weight = false;

	/**
	 * If the ratio of the log-likelihood (compared to that of the previous iteration)
	 * is lower than the convergence threshold <tt>thres</tt>, the algorithm is assumed to have converged. <br>
	 * Good values are positive near to zero. E.g. +1e-8.  						 */
	private double thres = 1e-8; 
		
	/**
	 * If the log-likelihood (compared to that of the previous iteration) has been
	 * decreased to a ratio lower than the decreasing threshold <tt>decr_thres</tt>,
	 * the algorithm has surely not converged.									<br>
	 * Good values are negative near to zero. E.g. -1e-1. 						 */
    private double decr_thres = -1e-6;
    
	/** Each prior weight is considered non-zero (significant) if it has value more than <tt>prior_weight_thres</tt> */
	private double prior_weight_thres = 1e-4;
    
    /** Parameter of the negative Dirichlet-type prior.	The larger it is, the more sparseness it imposes. */
    private double alpha = 10.0;    
    
    ////  PRINTING PREFERENCES  ////
	private String   prior_weight_path = null; //"/Users/gio_fou/Desktop/multi_toy_example/";
	private BufferedWriter   glob_prior_weight_bw  = null;
	private BufferedWriter[] cond_prior_weight_bw  = null;
	
	private String suffix                = null;       // Suffix for the names of the files that are printed
	
	private String log_lik_path          = null;       //"/Users/gio_fou/Desktop/multi_toy_example/logliks.txt";
	private BufferedWriter   log_lik_bw  = null;
	

	/**********************
	 **  LEARNED PARAMS  **
	 **********************/

    /** Map containing C entries (the number of conditions).					 <br> 
     *  Each entry holds as a value, a V x M matrix representing 
     *  the responsibilities of each position (for that condition)  
     *  under each component: p(z_k = j | x_k = l)  							 <br>
     *  (j = 0...M-1), (k = 0...V-1)											 <br>
     *  Obviously, rows sum up to 1.
     */
	private Map<Integer, double[][]> ga;
	
    /** 1 x M array representing the sums of the posterior probabilities        <br>
     *  Sum_{n} p(z_n = j | x_n)                (j = 0...M-1)                    */
	private double[] glob_sum_ga;
	
	/** 1 x M array (summing to 1) representing the weighting of the components. <br>
	 *  Usually, high weights correspond to potential binding events.           */
	private double[] glob_prior_weight;
	
    /** C x M matrix representing the sums of the posterior probabilities        <br>
     *  Sum_{n} p(z_n = j | x_n)  (j = 0...M-1) for each condition.             */
	private double[][] sum_ga;
	
	/** C x M matrix (summing to 1) representing the weighting of the components for each condition. <br>
	 *  Usually, high weights correspond to potential binding events.           */
	private double[][] prior_weight;

	public MultiIndependentMixtureCounts() {
		
	}
	
	public void exec(int[][] pos, int[][] count, double[][] emit_mat) {
		char[][] strand = new char[0][0];
		exec(pos, count, strand, emit_mat);
	}
	
	public void exec(int[][] pos, int[][] count, char[][] strand, double[][] emit_mat) {
		int M = emit_mat.length;
		int[] compPos = new int[M];
		for(int j = 0; j < M; j++) { compPos[j] = j; }
		exec(pos, count, strand, compPos, emit_mat);
	}

	
	public void exec(int[][] pos, int[][] count, int[] compPos, double[][] emit_mat) { char[][] strand = new char[0][0]; exec(pos, count, strand, compPos, emit_mat); }
	
	/**
	 * In this implementation of exec, we use an emission matrix instead of the binding model.
	 * @param pos C x V matrix representing the (unique) positions of the data points
	 * @param count C x V matrix representing the counts in the corresponding (unique) positions
	 * @param strand C x V matrix representing the strands at these positions
	 * @param compPos 1 x M array representing the positions of the components         <br>
	 * For the sake of simplicity, we assume that the components' positions are sorted.
	 * @param emit_mat M x O emission matrix
	 */
	public void exec(int[][] pos, int[][] count, char[][] strand, int[] compPos, double[][] emit_mat) {
		/* Exceptions body */
		if(pos.length == 0) 			 					   { throw new IllegalArgumentException("An empty data set has been entered."); }
		if(strand.length != 0 && pos.length != strand.length) { throw new IllegalArgumentException("data and strand should have the same number of rows (conditions)."); }
		for(int t = 0; t < pos.length; t++)
			if(strand.length != 0 && pos[t].length != strand[t].length)
				throw new IllegalArgumentException("Each condition must have the same number of data and strand points.");
		if(compPos.length != emit_mat.length)
			throw new IllegalArgumentException("The length of the compPos array should equal the number of rows of the emit_mat matrix.");
		
		/* Set values */
		bm = null;
		C  = pos.length;
		this.compPos = compPos;
		M  = emit_mat.length;
		O  = emit_mat[0].length;
		V  = new int[C];
		N  = new int[C];
		for(int t = 0; t < C; t++) { V[t] = pos[t].length; }
		for(int t = 0; t < C; t++) { for(int k = 0; k < V[t]; k++) { N[t] = count[t][k]; } }
		Ntot = 0;
		for(int t = 0; t < C; t++) { Ntot += N[t]; }
		emit_mat_transp = transpose(emit_mat);

		/* Start EM */
		em(pos, count, strand);
	}//end of exec method
	
	public void exec(int[][] pos, int[][] count, int M, BindingModel bm) {
		char[][] strand = new char[0][0];
		exec(pos, count, strand, M, bm);
	}
	
	public void exec(int[][] pos, int[][] count, char[][] strand, int M, BindingModel bm) {
		int[] compPos = new int[M];
		for(int j = 0; j < M; j++) { compPos[j] = j; }
		exec(pos, count, strand, compPos, bm);
	}
	
	public void exec(int[][] pos, int[][] count, int[] compPos, BindingModel bm) { char[][] strand = new char[0][0]; exec(pos, count, strand, compPos, bm); }
	
	/**
	 * In this implementation of exec, we use the binding model instead of an emission matrix.
	 * @param pos C x V matrix representing the (unique) positions of the data points
	 * @param count C x V matrix representing the counts in the corresponding (unique) positions
	 * @param strand C x V matrix representing the strands at these positions
	 * @param compPos 1 x M array representing the positions of the components. 	<br>
	 * For the sake of simplicity, we assume that the components' positions are sorted.
	 * @param bm The binding model
	 */
	public void exec(int[][] pos, int[][] count, char[][] strand, int[] compPos, BindingModel bm) {
		/* Exceptions body */
		if(pos.length == 0) 			 					   { throw new IllegalArgumentException("An empty data set has been entered."); }
		if(strand.length != 0 && pos.length != strand.length) { throw new IllegalArgumentException("data and strand should have the same number of rows (conditions)."); }
		for(int t = 0; t < pos.length; t++)
			if(strand.length != 0 && pos[t].length != strand[t].length)
				throw new IllegalArgumentException("Each condition must have the same number of data and strand points.");

		/* Set values */
		emit_mat_transp = null;
		C       = pos.length;
		this.compPos = compPos;
		M = this.compPos.length;
		this.bm = bm;
		V  = new int[C];
		N  = new int[C];
		for(int t = 0; t < C; t++) { V[t] = pos[t].length; }
		for(int t = 0; t < C; t++) { for(int k = 0; k < V[t]; k++) { N[t] += count[t][k]; } }
		Ntot = 0;
		for(int t = 0; t < C; t++) { Ntot += N[t]; }

		/* Start EM */
		em(pos, count, strand);
	}
	
	/**
	 * EM algorithm for multi-condition data.
	 * @param pos C x V matrix representing the (unique) positions of the data points
	 * @param count C x V matrix representing the counts in the corresponding (unique) positions
	 * @param strand C x V matrix representing the strands at these positions
	 */
	private void em(int[][] pos, int[][] count, char[][] strand) {
		long start_em = System.currentTimeMillis();
		glob_prior_weight_0 = new double[M];
		prior_weight_0      = new double[C][M];
		init(glob_prior_weight_0, prior_weight_0);
		train(pos, count, strand);
		long end_em = System.currentTimeMillis();
//<----------		System.out.printf("Duration of EM: %s.%n", get_duration_info(end_em - start_em));
	}//end of em method
	
	/**
	 * Initializes the global prior weighting and the prior weighting for each
	 * condition if there are more than one conditions.
	 * @param glob_prior_weight_0 initialization of the 1 x M global prior weighting
	 * @param prior_weight_0 initialization of the C x M  prior weighting for each condition
	 */
	private void init(double[] glob_prior_weight_0, double[][] prior_weight_0) {
		// Put 0.1 to make sure that no components are rejected from the beginning
		for(int j = 0; j < M; j++) { glob_prior_weight_0[j] = 0.1 + Math.random(); }		
		StatUtil.normalize(glob_prior_weight_0);
		
		// If there are more than one conditions, initialize the prior weighting as well
		if(C > 1) {
			for(int t = 0; t < C; t++)
				for(int j = 0; j < M; j++) 
					prior_weight_0[t][j] = 0.1 + Math.random();
			
			StatUtil.normalize(prior_weight_0, 1);
		}
	}//end of init method
	
	/**
	 * Method where the EM training takes place.
	 * @param pos C x V matrix representing the (unique) positions of the data points
	 * @param count C x V matrix representing the counts in the corresponding (unique) positions
	 * @param strand C x V matrix representing the strands at these positions
	 */
	private void train(int[][] pos, int[][] count, char[][] strand) {
		
		/***************************************
		 **			    CHECKS                **
		 ***************************************/
		if(thres < 0 || thres > 1 || decr_thres > 0 || decr_thres < -1)
			throw new IllegalArgumentException("thres must be a positive number between 0 and 1 and particularly close to 0.\n" +
											   "decr_thres must be a negative number between -1 and 0 and particularly close to 0.");
		
		if(prior_weight_path != null) {
			try {
				String glob_prior_weight_file = prior_weight_path + "glob_prior" + (suffix == null ? "" : "_" + suffix) + ".txt";
				glob_prior_weight_bw = new BufferedWriter(new FileWriter(glob_prior_weight_file));
				if(C > 1) {
					String[] cond_prior_weight_file = new String[C];
					for(int t = 0; t < C; t++) { cond_prior_weight_file[t] = prior_weight_path + "prior" + (t+1) + (suffix == null ? "" : "_" + suffix) + ".txt"; }
					cond_prior_weight_bw = new BufferedWriter[C];
					for(int t = 0; t < C; t++) { cond_prior_weight_bw[t] = new BufferedWriter(new FileWriter(cond_prior_weight_file[t])); }
				}
			}
			catch (IOException e) { e.printStackTrace(); }
		}
		
		if(log_lik_path != null) {
			try { log_lik_bw = new BufferedWriter(new FileWriter(log_lik_path + "logliks" + (suffix == null ? "" : "_" + suffix) + ".txt")); } 
			catch (IOException e) { e.printStackTrace(); }
		}
		
		
		/****************************************
		 **			INITIALIZATIONS            **
		 ***************************************/		
		numIters           = 0;
		hasConverged       = false;
	
		// global initialization
		glob_loglik       = 0;
		glob_logliks      = new ArrayList<Double>();
		glob_sum_ga       = new double[M];
		glob_prior_weight = glob_prior_weight_0.clone();
		
		//condition-wise initialization
		if(C > 1) {
			loglik       = new double[C];
			logliks		 = new ArrayList[C];
			for(int t = 0; t < C; t++) { logliks[t] = new ArrayList<Double>(); }
			sum_ga       = new double[C][M];
			prior_weight = prior_weight_0.clone();
		}
		
		long start_curr_print_step = System.currentTimeMillis();
		/***********************************
		 **			 TRAINING             **
		 **********************************/
		
		AggregatedData ad = new AggregatedData(pos, count, strand);
		
		/*****************************************
		 **		    SERIAL VERSION              **
		 **        First, global EM,  	    	** 
		 **  then each condition, separately.   **
		 ****************************************/
		if(trainVersion == 1)      // Hard Serial
			serial_hard_train(pos, count, strand);
		else if(trainVersion == 2) // Soft Serial
			serial_soft_train(pos, count, strand);

		/*****************************************
		 **		    PARALLEL VERSION            **
		 **  Global and condition-specific EMs  ** 
		 **      run at the same time.			**
		 ****************************************/
		else  // Parallel Version
			parallel_train(pos, count, strand);
				
		boolean areAllZeroComponents = true;
		for(int j = 0; j < M; j++) {
			if(glob_prior_weight[j] > 0.0) {
				areAllZeroComponents = false;
				break;
			}
		}
		
		if(areAllZeroComponents)
			for(int t = 0; t < C; t++)
				prior_weight[t] = new double[M];
		
		// Keep only the components that exceed a minimum weight threshold
		if(C>1 && !areAllZeroComponents) {		
			for(int t = 0; t < C; t++)
				for(int j = 0; j < M; j++)
					if(prior_weight[t][j] < prior_weight_thres)
						prior_weight[t][j] = 0.0;
			
			StatUtil.normalize(prior_weight, 1);
			
			for(int j = 0; j < M; j++) {
				double sum = 0.0;
				for(int t = 0; t < C; t++)
					sum += prior_weight[t][j];
				
				if(sum < prior_weight_thres) { glob_prior_weight[j] = 0.0; }
			}
			StatUtil.normalize(glob_prior_weight);			
		}
		else if(!areAllZeroComponents) {
			for(int j = 0; j < M; j++)
				if(glob_prior_weight[j] < prior_weight_thres)
					glob_prior_weight[j] = 0.0;
			
			StatUtil.normalize(glob_prior_weight);
		}
		
		// Evaluate the global likelihood (of all data)
		int[] glob_pos     = AggregatedData.get_glob_pos();
		int[] glob_count   = AggregatedData.get_glob_count();
		char[] glob_strand = AggregatedData.get_glob_strand();
		glob_loglik = eval_loglik(glob_pos, glob_count, glob_strand, glob_prior_weight);
		
		// Construct (final) responsibilities
		ga = new HashMap<Integer, double[][]>();
		for(int t = 0; t < C; t++) {
			double[][] cond_ga  = new double[V[t]][M];
			char[] curr_strand = strand.length != 0 ? strand[t] : new char[0];
			for(int k = 0; k < V[t]; k++) {
				char str = curr_strand.length == 0 ? '+' : curr_strand[k];
				if(C > 1) {
					for(int j = 0; j < M; j++) { cond_ga[k][j] = eval_emit_prob(pos[t][k], str, j)*prior_weight[t][j]; }
				}
				else {
					for(int j = 0; j < M; j++) { cond_ga[k][j] = eval_emit_prob(pos[t][k], str, j)*glob_prior_weight[j]; }
				}
				StatUtil.normalize(cond_ga, 1);
			}
			ga.put(t, cond_ga);
		}			

		
		/****************************************
		 **	     PRINT USEFUL INFO             **
		 ***************************************/
		if(log_lik_path != null) { 
			print_double_array(log_lik_bw, Utils.ref2prim(glob_logliks.toArray(new Double[0])));
			if(C > 1) {
				for(int t = 0; t < C; t++) { print_double_array(log_lik_bw, Utils.ref2prim(logliks[t].toArray(new Double[0]))); }
			}
		}

/* --------->		
   if(numIters%print_step != 0) { 
			System.out.printf("%nIteration %d, duration of %d iters: %s, glob_loglik = %.2f%n", numIters, numIters%print_step, get_duration_info(System.currentTimeMillis()-start_curr_print_step), glob_loglik);
			if(C > 1) {
				for(int t = 0; t < C; t++)
					System.out.printf("Condition %d, iteration %d, duration of %d iters: %s, loglik = %.2f%n", t, numIters, numIters%print_step, get_duration_info(System.currentTimeMillis()-start_curr_print_step), loglik[t]);
			}
		}
<---------		
*/		
		close_all_files();
	}//end of train method
	
	private void serial_hard_train(int[][] pos, int[][] count, char[][] strand) {
		int[] glob_pos     = AggregatedData.get_glob_pos();
		int[] glob_count   = AggregatedData.get_glob_count();
		char[] glob_strand = AggregatedData.get_glob_strand();
		double prev_glob_loglik = Double.NEGATIVE_INFINITY;

/***************/
//		System.out.println("Global Prior weights initially");
//		for(int j = 0; j < M; j++)
//				System.out.printf("%.3f\t", glob_prior_weight[j]);
//			System.out.print("\n");
//		
/***************/

		while(numIters < maxIters && !hasConverged) {
			// Determine value of alpha for current iteration
			if(isBatchElimOn) {
				if(numIters < ML_maxIters)           { curr_alpha = 0; }
				else if(numIters < anneal_maxIters)  { curr_alpha = alpha*(numIters - ML_maxIters)/(anneal_maxIters - ML_maxIters); }
				else								 { curr_alpha = alpha; }
			}
			else {
				if(numIters < ML_maxIters) { curr_alpha = 0; }
				else                       { curr_alpha = Math.max(minResp_minIndex.car(), alpha/2); }
			}
			
			
		    /**	    		GLOBAL				   **/
			
			/****************************************
			 **			GLOBAL E STEP              **
			 ***************************************/
			E_step(glob_pos, glob_count, glob_strand, glob_prior_weight, glob_sum_ga);
			
			/****************************************
			 **		    GLOBAL M STEP              **
			 ***************************************/
			M_step(glob_prior_weight, glob_sum_ga, curr_alpha);			
							
			/****************************************
			 **			EVAL STATISTICS            **
			 ***************************************/
			
			// Store log-likelihood of the current iteration
			glob_loglik = eval_loglik(glob_pos, glob_count, glob_strand, glob_prior_weight);
			glob_logliks.add(glob_loglik);
			
			// Check if the thresholding criterion has been satisfied (after the maximum number of annealing iterations has been exceeded)
			if(numIters >= anneal_maxIters || !isBatchElimOn)
				hasConverged = convergence_check(glob_loglik, prev_glob_loglik, thres, check_increased, decr_thres);
			
			prev_glob_loglik = glob_loglik;
			numIters++;
			
			/****************************************
			 **	     PRINT USEFUL INFO             **
			 ***************************************/
//			for(int j = 0; j < M; j++)
//				System.out.printf("%.3f\t", glob_prior_weight[j]);
//			System.out.printf("\n");
/* --------->
  			if(prior_weight_path != null) {
				print_double_array(glob_prior_weight_bw, glob_prior_weight);
			}
			
			if(numIters == 1 || numIters%print_step == 0) {
				long end_curr_print_step = System.currentTimeMillis();
				System.out.printf("%nIteration %d, duration of %d iters: %s, glob_loglik = %.2f%n", numIters, numIters < print_step ? numIters:print_step, get_duration_info(end_curr_print_step-start_curr_print_step), glob_loglik);
				start_curr_print_step    = end_curr_print_step;
			}
<---------
*/			
		}//end of while loop
		
//		System.out.println("numIters\t" + numIters + "\t MLIters\t" + ML_maxIters + "\t annealIters\t" + anneal_maxIters + "\tmaxIters\t" + maxIters);
		
		if(C > 1 && glob_prior_weight.length > 0) {
			int M_temp = M;
			int[] compPos_temp = compPos.clone();
			double[][] emit_mat_transp_temp = null;
			if(emit_mat_transp != null) { emit_mat_transp_temp = emit_mat_transp.clone(); }
			
			List<Integer> nonzeroCompsIdxList = new ArrayList<Integer>();
			for(int j = 0; j < M; j++)
				if(glob_prior_weight[j] >= prior_weight_thres)
					nonzeroCompsIdxList.add(j);
			nonzeroCompIdx = new int[nonzeroCompsIdxList.size()];
			for(int k = 0; k < nonzeroCompIdx.length; k++) { nonzeroCompIdx[k] = nonzeroCompsIdxList.get(k); }
						
			// Temporarily change the true values to speed up inference
			M = nonzeroCompIdx.length;
			compPos = new int[M];
			for(int j = 0; j < M; j++) { compPos[j] = compPos_temp[nonzeroCompIdx[j]]; }
			for(int t = 0; t < C; t++) { prior_weight[t] = new double[M]; sum_ga[t] = new double[M]; }
			
			for(int t = 0; t < C; t++) { for(int j = 0; j < M; j++) { prior_weight[t][j] = 1.0/M; } }
			StatUtil.normalize(prior_weight, 1);
			
/********/
//			System.out.println("\n\nCondition-specific EM has started. \n\nNon zero indices");
//			for(int k = 0; k < nonzeroCompIdx.length; k++)
//				System.out.print(nonzeroCompIdx[k] + "\t");
//			
//			System.out.println("Prior weights initially");
//			for(int t = 0; t < C; t++) {
//				for(int j = 0; j < M; j++)
//					System.out.printf("%.3f\t", prior_weight[t][j]);
//				System.out.print("\n");
//			}
/********/
			
			if(emit_mat_transp != null) {
				for(int l = 0; l < emit_mat_transp.length; l++) { 
					emit_mat_transp[l] = new double[M];
					for(int j = 0; j < M; j++)
						emit_mat_transp[l][j] = emit_mat_transp_temp[l][nonzeroCompIdx[j]];
				}				
			}
			
			double[] prev_loglik = new double[C];
			Arrays.fill(prev_loglik, Double.NEGATIVE_INFINITY);
			
			for(int t = 0; t < C; t++) {				
				numIters     = 0;
				hasConverged = false;
				long start_curr_print_step = System.currentTimeMillis();
///* --------->*/				System.out.println("\nCondition " + (t+1) + " prior weights...");
				while(numIters < maxIters && !hasConverged) {
				    /**        CONDITION-WISE           	**/
					
					/****************************************
					 **			CONDITION-WISE E STEP      **
					 ***************************************/
					char[] curr_strand = strand.length != 0 ? strand[t] : new char[0];
					E_step(pos[t], count[t], curr_strand, prior_weight[t], sum_ga[t]);
			
					/****************************************
					 **		    CONDITION-WISE M STEP      **
					 ***************************************/
					M_step(prior_weight[t], sum_ga[t], 0.0);
					
					/****************************************
					 **			EVAL STATISTICS            **
					 ***************************************/
					
					// Store log-likelihood of the current iteration						
					loglik[t] = eval_loglik(pos[t], count[t], curr_strand, prior_weight[t]);
					logliks[t].add(loglik[t]);
					
					// Check if the thresholding criterion has been satisfied
					hasConverged = convergence_check(loglik[t], prev_loglik[t], thres, check_increased, decr_thres);
		
					prev_loglik[t] = loglik[t];
					numIters++;
					
					/****************************************
					 **	     PRINT USEFUL INFO             **
					 ***************************************/
//					for(int j =0; j < M; j++)
//						System.out.printf("%.3f\t", prior_weight[t][j]);
//					System.out.print("\n");
		/* --------->
		  			if(prior_weight_path != null) {
		  				print_double_array(cond_prior_weight_bw[t], prior_weight[t]);
					}
					
					if(numIters == 1 || numIters%print_step == 0) {
						long end_curr_print_step = System.currentTimeMillis();
						System.out.printf("Condition %d, iteration %d, duration of %d iters: %s, loglik = %.2f%n", t, numIters, numIters < print_step ? numIters:print_step, get_duration_info(end_curr_print_step-start_curr_print_step), loglik[t]);
						start_curr_print_step = end_curr_print_step;
					}
		<---------
		*/			
				}//end of while loop
			}
			
			// Change back to true (initial values)
			M = M_temp;
			compPos = compPos_temp;
			emit_mat_transp = emit_mat_transp_temp;
			for(int t = 0; t < C; t++) {
				double[] prior_weight_temp = prior_weight[t].clone();
				double[] sum_ga_temp       = sum_ga[t].clone();
				prior_weight[t] = new double[M]; 
				sum_ga[t] = new double[M];
				for(int k = 0; k < nonzeroCompIdx.length; k++) {
					prior_weight[t][nonzeroCompIdx[k]] = prior_weight_temp[k];
					sum_ga[t][nonzeroCompIdx[k]]       = sum_ga_temp[k];
				}
			}
		}
	}//end of serial_hard_train method
	
	private void serial_soft_train(int[][] pos, int[][] count, char[][] strand) {
		int[] glob_pos     = AggregatedData.get_glob_pos();
		int[] glob_count   = AggregatedData.get_glob_count();
		char[] glob_strand = AggregatedData.get_glob_strand();
		double prev_glob_loglik = Double.NEGATIVE_INFINITY;
		
		while(numIters < maxIters && !hasConverged) {
			// Determine value of alpha for current iteration
			if(isBatchElimOn) {
				if(numIters < ML_maxIters)           { curr_alpha = 0; }
				else if(numIters < anneal_maxIters)  { curr_alpha = alpha*(numIters - ML_maxIters)/(anneal_maxIters - ML_maxIters); }
				else								 { curr_alpha = alpha; }
			}
			else {
				if(numIters < ML_maxIters) { curr_alpha = 0; }
				else                       { curr_alpha = Math.max(minResp_minIndex.car(), alpha/2); }
			}
			
		    /**	    		GLOBAL				   **/
			
			/****************************************
			 **			GLOBAL E STEP              **
			 ***************************************/
			E_step(glob_pos, glob_count, glob_strand, glob_prior_weight, glob_sum_ga);
			
			/****************************************
			 **		    GLOBAL M STEP              **
			 ***************************************/
			M_step(glob_prior_weight, glob_sum_ga, curr_alpha);			
							
			/****************************************
			 **			EVAL STATISTICS            **
			 ***************************************/
			
			// Store log-likelihood of the current iteration
			glob_loglik = eval_loglik(glob_pos, glob_count, glob_strand, glob_prior_weight);
			glob_logliks.add(glob_loglik);
			
			// Check if the thresholding criterion has been satisfied (after the maximum number of annealing iterations has been exceeded)
			if(numIters >= anneal_maxIters || !isBatchElimOn)
				hasConverged = convergence_check(glob_loglik, prev_glob_loglik, thres, check_increased, decr_thres);
			
			prev_glob_loglik = glob_loglik;
			numIters++;
			
			/****************************************
			 **	     PRINT USEFUL INFO             **
			 ***************************************/
/* --------->
  			if(prior_weight_path != null) {
				print_double_array(glob_prior_weight_bw, glob_prior_weight);
			}
			
			if(numIters == 1 || numIters%print_step == 0) {
				long end_curr_print_step = System.currentTimeMillis();
				System.out.printf("%nIteration %d, duration of %d iters: %s, glob_loglik = %.2f%n", numIters, numIters < print_step ? numIters:print_step, get_duration_info(end_curr_print_step-start_curr_print_step), glob_loglik);
				start_curr_print_step    = end_curr_print_step;
			}
<---------
*/			
		}//end of while loop
		
		if(C > 1) {
			double[] gamma_pseudocounts = new double[M];
			double[] prev_loglik = new double[C];
			Arrays.fill(prev_loglik, Double.NEGATIVE_INFINITY);
			
			for(int t = 0; t < C; t++) {
				
				for(int j = 0; j < M; j++)
					gamma_pseudocounts[j] = -(N[t]*glob_prior_weight[j] +1 -1); // add 1 as a pseudo-count
				
				numIters     = 0;
				hasConverged = false;
				long start_curr_print_step = System.currentTimeMillis();
				while(numIters < maxIters && !hasConverged) {
				    /**        CONDITION-WISE           	**/
					
					/****************************************
					 **			CONDITION-WISE E STEP      **
					 ***************************************/
					char[] curr_strand = strand.length != 0 ? strand[t] : new char[0];
					E_step(pos[t], count[t], curr_strand, prior_weight[t], sum_ga[t]);
					
					/****************************************
					 **		    CONDITION-WISE M STEP      **
					 ***************************************/
					M_step(prior_weight[t], sum_ga[t], gamma_pseudocounts);
					
					/****************************************
					 **			EVAL STATISTICS            **
					 ***************************************/
					
					// Store log-likelihood of the current iteration						
					loglik[t] = eval_loglik(pos[t], count[t], curr_strand, prior_weight[t]);
					logliks[t].add(loglik[t]);
					
					// Check if the thresholding criterion has been satisfied
					hasConverged = convergence_check(loglik[t], prev_loglik[t], thres, check_increased, decr_thres);
		
					prev_loglik[t] = loglik[t];
					numIters++;
					
					/****************************************
					 **	     PRINT USEFUL INFO             **
					 ***************************************/
		/* --------->
		  			if(prior_weight_path != null) {
		  				print_double_array(cond_prior_weight_bw[t], prior_weight[t]);
					}
					
					if(numIters == 1 || numIters%print_step == 0) {
						long end_curr_print_step = System.currentTimeMillis();
						System.out.printf("Condition %d, iteration %d, duration of %d iters: %s, loglik = %.2f%n", t, numIters, numIters < print_step ? numIters:print_step, get_duration_info(end_curr_print_step-start_curr_print_step), loglik[t]);
						start_curr_print_step = end_curr_print_step;
					}
		<---------
		*/			
				}//end of while loop
			}
		}
	}//end of serial_soft_train method
	
	private void parallel_train(int[][] pos, int[][] count, char[][] strand) {
		int[] glob_pos     = AggregatedData.get_glob_pos();
		int[] glob_count   = AggregatedData.get_glob_count();
		char[] glob_strand = AggregatedData.get_glob_strand();
		double prev_glob_loglik = Double.NEGATIVE_INFINITY;
		
		while(numIters < maxIters && !hasConverged) {
			
			// Determine value of alpha for current iteration
			if(isBatchElimOn) {
				if(numIters < ML_maxIters)           { curr_alpha = 0; }
				else if(numIters < anneal_maxIters)  { curr_alpha = alpha*(numIters - ML_maxIters)/(anneal_maxIters - ML_maxIters); }
				else								 { curr_alpha = alpha; }
			}
			else {
				if(numIters < ML_maxIters) { curr_alpha = 0; }
				else                       { curr_alpha = Math.max(minResp_minIndex.car(), alpha/2); }
			}
			
		    /**	    		GLOBAL				   **/
			
			/****************************************
			 **			GLOBAL E STEP              **
			 ***************************************/
			E_step(glob_pos, glob_count, glob_strand, glob_prior_weight, glob_sum_ga);
			
			/****************************************
			 **		    GLOBAL M STEP              **
			 ***************************************/
			M_step(glob_prior_weight, glob_sum_ga, curr_alpha);			
			
			
			/**	    	CONDITION-WISE			   **/
			if(C > 1) {
				double[] glob_prior_weight_factor;
				if(mask_glob_prior_weight) { glob_prior_weight_factor = binarize_data(glob_prior_weight); }
				else                       { glob_prior_weight_factor = glob_prior_weight; }
				
				for(int t = 0; t < C; t++) {
					/***********************************************************
					 **     CONDITION-WISE PRIOR WEIGHTED BY GLOBAL PRIOR     **
					 **********************************************************/
					for(int j = 0; j < M; j++) { prior_weight[t][j] *= glob_prior_weight_factor[j]; }
					StatUtil.normalize(prior_weight[t]);
					
					/***********************************************************
					 **              CONDITION-WISE E STEP                    **
					 **********************************************************/
					char[] curr_strand = strand.length != 0 ? strand[t] : new char[0];
					E_step(pos[t], count[t], curr_strand, prior_weight[t], sum_ga[t]);
			
					/***********************************************************
					 **              CONDITION-WISE M STEP                    **
					 **********************************************************/
					M_step(prior_weight[t], sum_ga[t], 0.0);
				}
			}
			
			/****************************************
			 **			EVAL STATISTICS            **
			 ***************************************/
			
			// Store log-likelihood of the current iteration
			glob_loglik = eval_loglik(glob_pos, glob_count, glob_strand, glob_prior_weight);
			glob_logliks.add(glob_loglik);
			
			if(C > 1) {
				for(int t = 0; t < C; t++) {
					char[] curr_strand = strand.length != 0 ? strand[t] : new char[0];
					loglik[t] = eval_loglik(pos[t], count[t], curr_strand, prior_weight[t]);
					logliks[t].add(loglik[t]);
				}	
			}
			
			// Check if the thresholding criterion has been satisfied (after the maximum number of annealing iterations has been exceeded)
			if(numIters >= anneal_maxIters || !isBatchElimOn)
				hasConverged = convergence_check(glob_loglik, prev_glob_loglik, thres, check_increased, decr_thres);
			
			prev_glob_loglik = glob_loglik;
			numIters++;
			
			
			/****************************************
			 **	     PRINT USEFUL INFO             **
			 ***************************************/
/* --------->
  			if(prior_weight_path != null) {
				print_double_array(glob_prior_weight_bw, glob_prior_weight);
				if( C > 1) {
					for(int t = 0; t < C; t++) { print_double_array(cond_prior_weight_bw[t], prior_weight[t]); }
				}
			}
			
			if(numIters == 1 || numIters%print_step == 0) {
				long end_curr_print_step = System.currentTimeMillis();
				System.out.printf("%nIteration %d, duration of %d iters: %s, glob_loglik = %.2f%n", numIters, numIters < print_step ? numIters:print_step, get_duration_info(end_curr_print_step-start_curr_print_step), glob_loglik);
				if(C > 1) {
					for(int t = 0; t < C; t++)
						System.out.printf("Condition %d, iteration %d, duration of %d iters: %s, loglik = %.2f%n", t, numIters, numIters < print_step ? numIters:print_step, get_duration_info(end_curr_print_step-start_curr_print_step), loglik[t]);
				}
				start_curr_print_step    = end_curr_print_step;
			}
<---------
*/			
		}//end of while loop
	}//end of parallel_train method
	
	/**
	 * E step of the EM algorithm												<br>
	 * The identification whether we have an emission matrix or a binding model is
	 * being done inside the <tt>eval_emit_prob</tt> method.
	 * @param pos 1 x V array representing the (unique) positions of the data points
	 * @param count 1 x V array representing the counts in the corresponding (unique) positions
	 * @param strand 1 x V array representing the strands at these positions
	 * @param prior_weight 1 x M array representing the prior weighting
	 * @param sum_ga 1 x M array representing the sums of the posterior probabilities
	 */
	private void E_step(int[] pos, int[] count, char[] strand, double[] prior_weight, double[] sum_ga) {
		for(int j = 0; j < M; j++) { sum_ga[j] = 0.0; }
		
		for(int k = 0; k < pos.length; k++) {
			char str    = strand.length == 0 ? '+' : strand[k];
			double[] ga = new double[M];
			for(int j = 0; j < M; j++) { ga[j] = eval_emit_prob(pos[k], str, j)*prior_weight[j]; }
			StatUtil.normalize(ga);
			for(int j = 0; j < M; j++) { sum_ga[j] += count[k]*ga[j]; }
		}
	}//end of E_step method
	
	/**
	 * M step of the EM algorithm
	 * @param prior_weight 1 x M array representing the prior weighting
	 * @param sum_ga 1 x M array representing the sums of the posterior probabilities
	 * @param alpha negative Dirichlet-type prior
	 */
	private void M_step(double[] prior_weight, double[] sum_ga, double alpha) {	
		double[] alpha_arr = new double[prior_weight.length];
		Arrays.fill(alpha_arr, alpha);
		M_step(prior_weight, sum_ga, alpha_arr);
	}//end of M_step method
	
	/**
	 * 
	 * @param prior_weight 1 x M array representing the prior weighting
	 * @param sum_ga 1 x M array representing the sums of the posterior probabilities
	 * @param pseudocounts pseudo-counts added/subtracted from the posterior probabilities. <br>
	 * If pseudo-counts are positive, then they are being subtracted by the posterior probabilities. <br>
	 * If they are negative, then they are being added by the posterior probabilities. <br>
	 * -- IT IS COUNTER-INTUITIVE, BUT I WON'T CHANGE IT FOR THE TIME BEING TO AVOID BUGS --
	 */
	private void M_step(double[] prior_weight, double[] sum_ga, double[] pseudocounts) {		
		if(isBatchElimOn) {
			for(int j = 0; j < M; j++) { prior_weight[j] = Math.max(0.0, sum_ga[j] - pseudocounts[j]); }	
		}
		else {
			double[] sum_ga_temp = sum_ga.clone();
			for(int j = 0; j < M; j++)
				if(sum_ga_temp[j] < 1e-100)
					sum_ga_temp[j] = Double.POSITIVE_INFINITY;
			
			minResp_minIndex = StatUtil.findMin(sum_ga_temp);
			// (minimum_responsibility-alpha) > zero => Perform sparse prior to all components
			if(minResp_minIndex.car() - pseudocounts[0] > 0) { 
				for(int j = 0; j < M; j++) { prior_weight[j] = Math.max(0.0, sum_ga[j] - pseudocounts[j]); } 
			}
			// There is at least one component with (responsibility-alpha) <= zero => Eliminate only the ones with the minimum value
			else {
				for(int j = 0; j < M; j++) { prior_weight[j] = sum_ga[j]; }
				for(int j_prime:minResp_minIndex.cdr()) { prior_weight[j_prime] = 0.0; } 
			}	
		}
		StatUtil.normalize(prior_weight);
	}//end of M_step method
	
	/**
	 * Evaluates the log-likelihood of the observed data.						<br>
	 * The identification whether we have an emission matrix or a binding model is
	 * being done inside the <tt>eval_emit_prob</tt> method. 
	 * @param pos 1 x V array representing the (unique) positions of the data points
	 * @param count 1 x V array representing the counts in the corresponding (unique) positions
	 * @param strand 1 x V array representing the strands at these positions
	 * @param prior_weight 1 x M array representing the prior weighting
	 * @return
	 */
	private double eval_loglik(int[] pos, int[] count, char[] strand, double[] prior_weight) {
		double loglik = 0.0;
		for(int k = 0; k < pos.length; k++) {
			char str   = strand.length == 0 ? '+' : strand[k];
			double sum = 0.0;
			for(int j = 0; j < M; j++) { sum    += eval_emit_prob(pos[k], str, j)*prior_weight[j]; }
			if(sum > 1e-20)            { loglik += count[k]*Math.log(sum); }
		}
		return loglik;
	}// end of eval_loglik method
	
	/**
	 * This method tests if the EM algorithm has converged. <br>
	 * It first examines if the likelihood has fallen below a decreasing threshold
	 * (<tt>decr_thres</tt>) and if no, it continues on by checking if the absolute
	 * likelihood difference is below than a convergence threshold (<tt>thres</tt>).  <br>
	 * We define as: <br>
	 * likelihood difference = loglik-prev_loglik								<br>
	 * average likelihood    = (abs(loglik) + abs(prev_loglik))/2				<br>
	 * 																			<br>
	 * If we are doing MAP estimation (using priors), the likelihood can decrease,
	 * even though the mode of the posterior is increasing.
	 * @param loglik the log-likelihood at the current iteration
	 * @param prev_loglik the log-likelihood at the previous iteration
	 * @param thres Acceptance threshold. If the likelihood absolute difference
	 * is less than <tt>thres</tt> times the average likelihood, where we accept.
	 * @param check_increased A boolean variable which allows to check if the
	 * log-likelihood has increased (compared to the previous step)
	 * @param decr_thres Decreasing threshold. If the likelihood difference has
	 * decreased more than <tt>decr_thres</tt> times the average likelihood, then it is not yet converged. <br>
	 * Otherwise, we let one further comparison where we check if the likelihood absolute difference
	 * is less than <tt>thres</tt> times the average likelihood, where we accept.
	 * @return <tt>true</tt> if it has converged, <tt>false</tt> otherwise.
	 */
	public boolean convergence_check(double loglik, double prev_loglik, double thres, boolean check_increased, double decr_thres) {
		
		if(thres < 0 || thres > 1 || decr_thres > 0 || decr_thres < -1)
			throw new IllegalArgumentException("thres must be a positive number between 0 and 1 and particularly close to 0.\n" +
					                           "decr_thres must be a negative number between -1 and 0 and particularly close to 0.");
		
		boolean hasConverged = false;
		
		double delta_loglik = loglik-prev_loglik;
		double avg_abs_loglik = (Math.abs(loglik) + Math.abs(prev_loglik) + Double.MIN_VALUE)/2;
		
		if(delta_loglik == Double.POSITIVE_INFINITY)
			return false;		
		
		if(loglik - prev_loglik >= 0.0 && loglik - prev_loglik < 1e-20)
			return true;
		
		if(check_increased)
			if(delta_loglik/avg_abs_loglik < decr_thres)
				return false;
		
		if( Math.abs(delta_loglik)/avg_abs_loglik < thres )
			hasConverged = true;
		
		return hasConverged;		
	}//end of convergence_check method

	
	/*********************
	 **   GET METHODS   **
	 ********************/

	public double[] get_glob_prior_weight()   { return glob_prior_weight; }
	
	public double[][] get_cond_prior_weight() { return prior_weight; }
	
	public Map<Integer, double[][]> get_resp() { return ga; }
	
	public double get_glob_loglik() { return glob_loglik; }
	
	public List<Double> get_glob_loglik_history() { return glob_logliks; }
	
	public List<Double>[] get_cond_loglik_history() { return logliks; }
	
	/*********************
	 **   SET METHODS   **
	 ********************/
	
	public void set_print_step(int new_print_step) { 
		if(new_print_step < 0) { throw new IllegalArgumentException("new_print_step must be a positive integer."); }
		print_step = new_print_step;
	}
	
	public void set_check_increased(boolean new_check_increased) { check_increased = new_check_increased; }
	
	public void set_maxIters(int new_maxIters) {
		if(new_maxIters < 0) { throw new IllegalArgumentException("new_maxIters must be a positive integer."); }
		if(new_maxIters < ML_maxIters || new_maxIters < anneal_maxIters) { throw new IllegalArgumentException("new_maxIters cannot be less than ML_maxIters or anneal_maxIters."); }
		maxIters = new_maxIters;
	}
	
	public void set_ML_maxIters(int new_ML_maxIters) {
		if(new_ML_maxIters < 0) { throw new IllegalArgumentException("new_ML_maxIters must be a positive integer."); }
		if(new_ML_maxIters > maxIters || new_ML_maxIters > anneal_maxIters) { throw new IllegalArgumentException("new_ML_maxIters cannot be greater than maxIters or anneal_maxIters."); }
		ML_maxIters = new_ML_maxIters;
	}

	public void set_anneal_maxIters(int new_anneal_maxIters) {
		if(new_anneal_maxIters < 0) { throw new IllegalArgumentException("new_anneal_maxIters must be a positive integer."); }
		if(new_anneal_maxIters < ML_maxIters || new_anneal_maxIters > maxIters) { throw new IllegalArgumentException("new_anneal_maxIters cannot be less than ML_maxIters or greater than maxIters."); }
		anneal_maxIters = new_anneal_maxIters;
	}
	
	public void set_isBatchElimOn(boolean new_isBatchElimOn) { isBatchElimOn = new_isBatchElimOn; }
	
	/** @param new_trainVersion 1 for Hard Serial, 2 for Soft Serial, 3 for Parallel */
	public void set_trainVersion(int new_trainVersion) { 
		if(new_trainVersion < 1 || new_trainVersion > 3)
			throw new IllegalArgumentException("The only valid values for trainVersion variable are 1:Hard Serial, 2:Soft Serial, 3:Parallel");
		trainVersion = new_trainVersion; 
	}

	/**  Sets a new value for the convergence threshold. Good values are positive values near 0.  */
	public void set_convergence_thres(double new_thres) {
		if(new_thres < 0 || new_thres > 1) { throw new IllegalArgumentException("new_thres must be within 0 and 1."); }
		thres = new_thres; 
	}
	
	/** Sets a new value for the decreasing threshold. Good values are negative values near 0. */
	public void set_decreasing_thres(double new_decr_thres) {
		if(new_decr_thres < -1 || new_decr_thres > 0) { throw new IllegalArgumentException("new_decr_thres must be within -1 and 0."); }
		decr_thres = new_decr_thres; 
	}
	
	public void set_mask_glob_prior_weight(boolean new_mask_glob_prior_weight) { mask_glob_prior_weight = new_mask_glob_prior_weight; }
	
	public void set_alpha(double new_alpha) {
		if(new_alpha < 0) { throw new IllegalArgumentException("new_alpha must be non-negative."); }
		alpha = new_alpha;
	}
	
	public void set_prior_weight_path(String new_prior_weight_path) { 
		if(new_prior_weight_path == null)
			throw new IllegalArgumentException("You cannot enter a null common path for the prior weights.");
		if( (!new_prior_weight_path.equals("")) && (new_prior_weight_path.lastIndexOf(File.separator) != new_prior_weight_path.length()-1) )
			throw new IllegalArgumentException("You must enter a valid path. One that ends with '/' in Unix systems or '\\\\' in Windows.");
		
		prior_weight_path = new_prior_weight_path; 
	}
	
	public void set_log_lik_file(String new_log_lik_path) {
		if(new_log_lik_path == null)
			throw new IllegalArgumentException("You cannot enter a null path for where the log-likelihood file will be located.");
		if( (!new_log_lik_path.equals("")) && (new_log_lik_path.lastIndexOf(File.separator) != new_log_lik_path.length()-1) )
			throw new IllegalArgumentException("You must enter a valid path. One that ends with '/' in Unix systems or '\\\\' in Windows.");
	
		log_lik_path = new_log_lik_path; 
	}	
	
	public void set_suffix(String new_suffix) { suffix = new_suffix; }
	
	/***********************************
	 **   AUXILIARY PRIVATE METHODS   **
	 ***********************************/
	private void close_all_files() {
		if(prior_weight_path != null) {
			try {
				glob_prior_weight_bw.close();
				if(C > 1) {
					for(int t = 0; t < C; t++) { cond_prior_weight_bw[t].close(); }	
				}
				
			} 
			catch (IOException e) { e.printStackTrace(); }			
		}
		
		if(log_lik_path != null) {
			try                   { log_lik_bw.close();  } 
			catch (IOException e) { e.printStackTrace(); }
		}
	}
	
	private void print_double_array(BufferedWriter bw, double[] f) {
		StringBuilder sb = new StringBuilder();
		
		if(f.length != 0) { sb.append(String.format("%.3f", f[0])); }
		if(f.length > 1) {
			for(int i = 1; i < f.length; i++)
				sb.append("\t" + String.format("%.3f", f[i]));
		}
		sb.append("\n");
		try {
			bw.write(sb.toString());
		} catch (IOException e) { e.printStackTrace(); }
	}
	
	/** binarize_data takes a double array and returns an array with ones wherever f is positive.*/
	private double[] binarize_data(double[] f) {
		double[] binarized_f = new double[f.length];
		for(int i = 0; i < binarized_f.length; i++) { 
			if(f[i] > 0.0) { binarized_f[i] = 1.0; }
		}
		return binarized_f;
	}

	private static void gather_data(List<Integer> glob_pos_list, List<Integer> glob_count_list, List<Character> glob_strand_list, int[][] pos, int[][] count, char[][] strand) {
		Map<Integer, Integer> aggr_posCount_plus  = new LinkedHashMap<Integer, Integer>();
		Map<Integer, Integer> aggr_posCount_minus = new LinkedHashMap<Integer, Integer>();
		
		int C   = pos.length;
		int[] V = new int[C];
		for(int t = 0; t < C; t++) { V[t] = pos[t].length; }
		
		for(int t = 0; t < C; t++) {
			for(int k = 0; k < V[t]; k++) {
				if(strand[t][k] == '+') {
					if(!aggr_posCount_plus.containsKey(pos[t][k])) { aggr_posCount_plus.put(pos[t][k], 0); }
					aggr_posCount_plus.put(pos[t][k], aggr_posCount_plus.get(pos[t][k]) + count[t][k]);
				}
				else {
					if(!aggr_posCount_minus.containsKey(pos[t][k])) { aggr_posCount_minus.put(pos[t][k], 0); }
					aggr_posCount_minus.put(pos[t][k], aggr_posCount_minus.get(pos[t][k]) + count[t][k]);					
				}	
			}
		}
				
		glob_pos_list.addAll(aggr_posCount_plus.keySet());
		glob_pos_list.addAll(aggr_posCount_minus.keySet());

		glob_count_list.addAll(aggr_posCount_plus.values());
		glob_count_list.addAll(aggr_posCount_minus.values());
		
		for(int i = 0; i < aggr_posCount_plus.size(); i++)  { glob_strand_list.add('+'); }
		for(int i = 0; i < aggr_posCount_minus.size(); i++) { glob_strand_list.add('-'); }
	}//end of gather_data method
	
	/**
	 * Returns the emission probability <tt>p(x=pos|z=compPos[event_idx])</tt> of the data point with
	 * position <tt>pos</tt> assigned to event <tt>compPos[event_idx]</tt>.
	 * @param pos position of the data point
	 * @param str the strand at that position
	 * @param event_idx index of the event (index of the compPos array)
	 * @return
	 */
	private double eval_emit_prob(int pos, char str, int event_idx) {
		int rel_pos = find_rel_pos(pos, str, event_idx);
		if(emit_mat_transp != null && bm == null)      { return emit_mat_transp[rel_pos][event_idx]; }
		else if(bm != null && emit_mat_transp == null) { return bm.probability(rel_pos); 	          }
		else                                           { throw new IllegalArgumentException("You cannot have an emission matrix and a binding model at the same time."); }
	}//end of eval_emit_prob method
	
	/**
	 * The <tt>find_rel_pos</tt> method finds the relative position of a data point
	 * based on two factors: its strand and whether we have been provided with an
	 * emission matrix or a binding model.										<br>
	 * 																			<br>
	 * Let's first consider the case where we have been provided with an emission matrix:   <br>
	 * - strand = '+': the relative position stays the same.					<br>
	 * - strand = '-': the relative position is on the symmetric position w.r.t.
	 * the y-axis located on the position of the event.
	 * 																			<br>
	 * When we are provided with a binding model instead, we have:				<br>
	 * - strand = '+': the relative position is the subtraction of the position
	 * of the event from the data point.										<br>
	 * - strand = '-': the relative position is the subtraction of the data point
	 * from the position of the event.											<br>
	 * That is, if a read (on the + strand) is upstream of an event, then its
	 * relative position to the event is negative. If it is downstream, it is positive.									 
	 * @param pos position of the data point
	 * @param strand the strand at that position
	 * @param event_idx index of the event (index of the compPos array)
	 * @return
	 */
	private int find_rel_pos(int pos, char str, int event_idx) {
		if(emit_mat_transp != null && bm == null) {
			int rel_data = str == '+' ? pos : compPos[event_idx] - (pos - compPos[event_idx]);
			if(rel_data < 0)   { rel_data = 0;   }
			if(rel_data > O-1) { rel_data = O-1; }
			return rel_data;
		}
		else if(bm != null && emit_mat_transp == null) {
			int rel_data = str == '+' ? pos - compPos[event_idx] : compPos[event_idx] - pos;
			return rel_data;
		}
		else { throw new IllegalArgumentException("You cannot have an emission matrix and a binding model at the same time."); }
	}//end of find_rel_pos method
	
     private double[][] transpose(double[][] a) { 
    	 double[][] b = new double[a[0].length][a.length];
		 for(int i = 0; i < b.length; i++) { for(int j = 0; j < b[i].length; j++)  b[i][j] = a[j][i];}
		 return b;
	 }//end of transpose method
     
     private static class AggregatedData {
    	 private static int[] glob_pos;
    	 private static int[] glob_count;
    	 private static char[] glob_strand;
    	 
    	 public AggregatedData(int[][] pos, int[][] count, char[][] strand) {
    		 List<Integer> glob_pos_list      = new ArrayList<Integer>();
    		 List<Integer> glob_count_list    = new ArrayList<Integer>();
    		 List<Character> glob_strand_list = new ArrayList<Character>();	
    		 gather_data(glob_pos_list, glob_count_list, glob_strand_list, pos, count, strand);
    		 glob_pos    = Utils.ref2prim(glob_pos_list.toArray(new Integer[0]));
    		 glob_count  = Utils.ref2prim(glob_count_list.toArray(new Integer[0]));
    		 glob_strand = Utils.ref2prim(glob_strand_list.toArray(new Character[0]));
    	 }
    	 
    	 public static int[]  get_glob_pos()     { return glob_pos;    }
    	 public static int[]  get_glob_count()   { return glob_count;  }
    	 public static char[] get_glob_strand()  { return glob_strand; }

     }//end of AggregatedData class

     private String get_duration_info(long msec_duration) {
       	if( msec_duration < 0) { throw new IllegalArgumentException("msec_duration must be a non-negative integer."); }
       	
       	StringBuilder sb = new StringBuilder();
       	long quot;
       	
       	long sec_duration_long   = msec_duration/1000;
       	double sec_duration = msec_duration/1000.0;
       	sec_duration -= sec_duration_long;
       	
       	// Get hour information
       	if( (quot = sec_duration_long/3600) != 0) { 
       		sb.append(quot + " h, "); 
       		sec_duration_long -= quot*3600;
       	}
       	
       	// Get minute information
       	if( (quot = sec_duration_long/60) != 0) {
       		sb.append(quot + " m, ");
       		sec_duration_long -= quot*60;
       	}
       	
       	sec_duration += sec_duration_long;
       	
       	// Get second information
       	sb.append(String.format("%.2f s", sec_duration));

       	return sb.toString();
       }//end of get_duration_info method
     
     
     public static void main(String[] args) {
    	 int[][] pos = {{2, 10, 15, 12, 15}, 
    			        {8, 15, 2},
    			        {16, 10, 15, 15}};
    	 
    	 int[][] count = {{3, 7, 1, 2, 2}, 
			        	  {2, 2, 5},
			        	  {4, 7, 6, 3}};
    	 
    	 char[][] str = {{'+', '+', '+', '+', '-'}, 
	        	  		 {'-', '-', '+'},
	        	  		 {'+', '+', '+', '-'}};
    	 
    	 MultiIndependentMixtureCounts mim = new MultiIndependentMixtureCounts();
    	 mim.C = 3;
    	 mim.V = new int[]{5,3,4};
   	  
    	 List<Integer> glob_pos_list = new ArrayList<Integer>();
    	 List<Integer> glob_count_list = new ArrayList<Integer>();
    	 List<Character> glob_strand_list = new ArrayList<Character>();
    	 mim.gather_data(glob_pos_list, glob_count_list, glob_strand_list, pos, count, str);
    	 
    	 int foo = 3;
    	 
     }
     
	
}//end of MultiIndependentMixtureCounts class

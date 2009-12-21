package edu.mit.csail.cgs.utils.randomgenerators;

import java.util.Arrays;
import java.util.Random;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.File;
import java.io.IOException;

import edu.mit.csail.cgs.utils.stats.StatUtil;

/**
 * This class will implement samplers of various distributions
 * @author geopapa
 *
 */
public class Sampler {
	
	private static Random r;
	
	/***********************
	 **  GAUSSIAN SAMPLES **
	 ***********************/
	
	public static double[][] gauss_rnd(int[] data_size, double[] mu) { return gauss_rnd(data_size, mu, new double[data_size.length]); }
	
	public static double[][] gauss_rnd(int[] data_size, double mu) { return gauss_rnd(data_size, mu); }
	
	public static double[][] gauss_rnd(int num_rows, int num_cols, double mu) { return gauss_rnd(num_rows, num_cols, mu, 1.0); }
	
	public static double[] gauss_rnd(int num_cols, double mu) { return gauss_rnd(num_cols, mu, 1.0); }

	/**
	 * It generates gaussian samples for <tt>data_size.length</tt> data sets.    <br>
	 * Each data set has data size <tt>data_size[t]</tt> with mean <tt>mu[t]</tt> and standard deviation <tt>std[t]</tt>.
	 * (where t in {0, ..., data_size.length-1}
	 * @param data_size the size of each data set
	 * @param mu the mean of the Gaussian
	 * @param std the standard deviation of the Gaussian
	 * @return
	 */
	public static double[][] gauss_rnd(int[] data_size, double[] mu, double[] std) {
		double[][] samples = new double[data_size.length][];
		for(int t = 0; t < samples.length; t++)
			samples[t] = gauss_rnd(data_size[t], mu[t], std[t]);
		
		return samples;
	}//end of gauss_rnd method
	
	/**
	 * It generates gaussian samples of mean <tt>mu</tt> and standard deviation <tt>std</tt>
	 * for <tt>data_size.length</tt> data sets, each having data_size[t] points. 
	 * (where t in {0, ..., data_size.length-1}
	 * @param data_size the size of each data set
	 * @param mu the mean of the Gaussian
	 * @param std the standard deviation of the Gaussian
	 * @return
	 */
	public static double[][] gauss_rnd(int[] data_size, double mu, double std) {
		double[][] samples = new double[data_size.length][];
		for(int t = 0; t < samples.length; t++)
			samples[t] = gauss_rnd(data_size[t], mu, std);
		
		return samples;
	}//end of gauss_rnd method
	
	public static double[][] gauss_rnd(int num_rows, int num_cols, double mu, double std) {
		double[][] samples = new double[num_rows][];
		for(int t = 0; t < samples.length; t++)
			samples[t] = gauss_rnd(num_cols, mu, std);
		
		return samples;
	}//end of gauss_rnd method
	
	public static double[] gauss_rnd(int num_cols, double mu, double std) {
		double[] samples = new double[num_cols];
		for(int n = 0; n < samples.length; n++)
			samples[n] = gauss_rnd(mu, std);
		
		return samples;
	}//end of gauss_rnd method
	
	public static double gauss_rnd(double mu) { return gauss_rnd(mu, 1.0); }
	
	/**
	 * It generates a gaussian sample of mean <tt>mu</tt> and standard deviation <tt>std</tt>.
	 * @param mu the mean of the gaussian
	 * @param std the standard deviation of the gaussian
	 * @return
	 */
	public static double gauss_rnd(double mu, double std) {
		if(std < 0) { throw new IllegalArgumentException("Standard deviation must be non-negative."); }
		if (r == null) initRNG();
		double sample = std*r.nextGaussian() + mu; 
		return sample;
	}//end of gauss_rnd method
	
	 private static synchronized void initRNG() {
		 if (r == null) 
			 r = new Random();
	 }//end of initRNG method
	
	

	/**********************************
	 **  CONTINUOUS UNIFORM SAMPLES **
	 **********************************/
	 
	 public static double[][] unif_rnd(int[] data_size, double[] lower_bound, double[] upper_bound) {
		 double[][] samples = new double[data_size.length][];
		 for(int t = 0; t < samples.length; t++)
			 samples[t] = unif_rnd(data_size[t], lower_bound[t], upper_bound[t]);

		 return samples;
	 }//end of unif_rnd method
	 
	 /**
	  * Gives <tt>num_rows</tt> data sets of continuous uniform samples, each having
	  * lower bound <tt>lower_bound[t]</tt> (inclusive) and upper bound <tt>upper_bound</tt> (exclusive). 
	  * @param num_rows
	  * @param num_cols
	  * @param lower_bound
	  * @param upper_bound
	  * @return
	  */
	 public static double[][] unif_rnd(int num_rows, int num_cols, double[] lower_bound, double[] upper_bound) {
		 double[][] samples = new double[num_rows][];
		 for(int t = 0; t < samples.length; t++)
			 samples[t] = unif_rnd(num_cols, lower_bound[t], upper_bound[t]);

		 return samples;
	 }//end of unif_rnd method

	 public static double[][] unif_rnd(int[] data_size, double lower_bound, double upper_bound) {
		 double[][] samples = new double[data_size.length][];
		 for(int t = 0; t < samples.length; t++)
			 samples[t] = unif_rnd(data_size[t], lower_bound, upper_bound);

		 return samples;
	 }//end of unif_rnd method
	 	 
	 public static double[][] unif_rnd(int num_rows, int num_cols, double lower_bound, double upper_bound) {
		 double[][] samples = new double[num_rows][];
		 for(int t = 0; t < samples.length; t++)
			 samples[t] = unif_rnd(num_cols, lower_bound, upper_bound);

		 return samples;
	 }//end of unif_rnd method
	 
	 public static double[][] unif_rnd(int[] data_size) {
		 double[][] samples = new double[data_size.length][];
		 for(int t = 0; t < samples.length; t++)
			 samples[t] = unif_rnd(data_size[t]);

		 return samples;
	 }//end of unif_rnd method
	 	 
	 public static double[] unif_rnd(int num_cols) {
		 return unif_rnd(num_cols, 0.0, 1.0);
	 }//end of unif_rnd method
 
	 public static double[] unif_rnd(int num_cols, double lower_bound, double upper_bound) {
		 double[] samples = new double[num_cols];
		 for(int n = 0; n < samples.length; n++)
			 samples[n] = unif_rnd(lower_bound, upper_bound);

		 return samples;
	 }//end of unif_rnd method

	 /**
	  * Gives a continuous uniform sample in the interval between <tt>lower_bound</tt> (inclusive)
	  * and <tt>upper_bound</tt> (exclusive).
	  * @param lower_bound
	  * @param upper_bound
	  * @return
	  */
	public static double unif_rnd(double lower_bound, double upper_bound) {
		if(upper_bound < lower_bound)
			throw new IllegalArgumentException("upper_bound has to be greater lower_bound.");
	
		double sample = Math.random();	
		double scale = upper_bound - lower_bound;
		sample       = scale*sample + lower_bound;
		return sample;
	}//end of unif_rnd method
	
	
	
	/*******************************
	 **  DISCRETE UNIFORM SAMPLES **
	 *******************************/
	 
	 public static int[][] unif_rnd(int[] data_size, int[] lower_bound, int[] upper_bound) {
		 int[][] samples = new int[data_size.length][];
		 for(int t = 0; t < samples.length; t++)
			 samples[t] = unif_rnd(data_size[t], lower_bound[t], upper_bound[t]);

		 return samples;
	 }//end of unif_rnd method
	 
	 /**
	  * Gives <tt>num_rows</tt> data sets of continuous uniform samples, each having
	  * lower bound <tt>lower_bound[t]</tt> and upper bound <tt>upper_bound</tt>. 
	  * @param num_rows
	  * @param num_cols
	  * @param lower_bound
	  * @param upper_bound
	  * @return
	  */
	 public static int[][] unif_rnd(int num_rows, int num_cols, int[] lower_bound, int[] upper_bound) {
		 int[][] samples = new int[num_rows][];
		 for(int t = 0; t < samples.length; t++)
			 samples[t] = unif_rnd(num_cols, lower_bound[t], upper_bound[t]);

		 return samples;
	 }//end of unif_rnd method

	 public static int[][] unif_rnd(int[] data_size, int lower_bound, int upper_bound) {
		 int[][] samples = new int[data_size.length][];
		 for(int t = 0; t < samples.length; t++)
			 samples[t] = unif_rnd(data_size[t], lower_bound, upper_bound);

		 return samples;
	 }//end of unif_rnd method
	 	 
	 public static int[][] unif_rnd(int num_rows, int num_cols, int lower_bound, int upper_bound) {
		 int[][] samples = new int[num_rows][];
		 for(int t = 0; t < samples.length; t++)
			 samples[t] = unif_rnd(num_cols, lower_bound, upper_bound);

		 return samples;
	 }//end of unif_rnd method
 
	 public static int[] unif_rnd(int num_cols, int lower_bound, int upper_bound) {
		 int[] samples = new int[num_cols];
		 for(int n = 0; n < samples.length; n++)
			 samples[n] = unif_rnd(lower_bound, upper_bound);

		 return samples;
	 }//end of unif_rnd method

	 /**
	  * Gives a discrete uniform sample in the interval between <tt>lower_bound</tt> (inclusive)
	  * and <tt>upper_bound</tt> (exclusive).
	  * @param lower_bound
	  * @param upper_bound
	  * @return
	  */
	public static int unif_rnd(int lower_bound, int upper_bound) {
		if(upper_bound < lower_bound)
			throw new IllegalArgumentException("upper_bound has to be greater lower_bound.");
	
		if (r == null) initRNG();
		int width = upper_bound - lower_bound;
		int sample = r.nextInt(width) + lower_bound;
		return sample;
	}//end of unif_rnd method
	
	
	
	/*******************************
	 **      BINOMIAL SAMPLES     **
	 *******************************/

	 public static int[][] binom_rnd(int[] data_size, int[] N, double[] p) {
		 int[][] samples = new int[data_size.length][];
		 for(int t = 0; t < samples.length; t++)
			 samples[t] = binom_rnd(data_size[t], N[t], p[t]);

		 return samples;
	 }//end of binom_rnd method
	
	 public static int[][] binom_rnd(int num_rows, int num_cols, int[] N, double[] p) {
		 int[][] samples = new int[num_rows][];
		 for(int t = 0; t < samples.length; t++)
			 samples[t] = binom_rnd(num_cols, N[t], p[t]);

		 return samples;
	 }//end of binom_rnd method
	
	 public static int[][] binom_rnd(int[] data_size, int N, double p) {
		 int[][] samples = new int[data_size.length][];
		 for(int t = 0; t < samples.length; t++)
			 samples[t] = binom_rnd(data_size[t], N, p);

		 return samples;
	 }//end of binom_rnd method
	
	 public static int[][] binom_rnd(int num_rows, int num_cols, int N, double p) {
		 int[][] samples = new int[num_rows][];
		 for(int t = 0; t < samples.length; t++)
			 samples[t] = binom_rnd(num_cols, N, p);

		 return samples;
	 }//end of binom_rnd method
	
	 public static int[] binom_rnd(int num_cols, int N, double p) {
		 int[] samples = new int[num_cols];
		 for(int n = 0; n < samples.length; n++)
			 samples[n] = binom_rnd(N, p);

		 return samples;
	 }//end of binom_rnd method
	 
	/**
	 * This method creates a sample from 0 to N-1 following a Binomial distribution, Bin(N, p), using the algorithm
	 * by Devroye taken from: 													<br>
	 * <a href="http://delivery.acm.org/10.1145/50000/42381/p216-kachitvichyanukul.pdf?key1=42381&key2=4588540521&coll=GUIDE&dl=GUIDE&CFID=48844090&CFTOKEN=30501264">
	 * Binomial Random Variate Generation</a>
	 * <br>
	 * Time Complexity: O(Np)
	 * @param N
	 * @param p
	 * @return
	 */
	public static int binom_rnd(int N, double p) {
		if(N <= 0 || p < 0 || p > 1)
			throw new IllegalArgumentException("n must be a non-negative integer, while p must be between 0.0 and 1.0");
		
		int x, y;
		double c;
		
		x = 0; y = 0; c = Math.log(1-p);
		
		if(c == 0 || N == 0) { return x; }
		
		boolean first_pass = true;
		while( y < N ) {
			if(!first_pass) { x++; }
			first_pass = false; 
			
			double unif_sample = Math.random();
			y += (int) Math.floor(Math.log(unif_sample)/c) + 1;
		}
		
		return x;
	}//end of binom_rnd method

	
	/**************************
	 **  MULTINOMIAL SAMPLES **
	 **************************/

	/**
	 * Creates <tt>data_size.length</tt> data sets, each of size <tt>data_size[t]</tt> 
	 * following a multinomial distribution <tt>prob[t]</tt>.
	 */
	public static int[][] multinomial_rnd(int[] data_size, double[][] prob) {
		int[][] samples = new int[data_size.length][];
		for(int t = 0; t < samples.length; t++)
			samples[t] = multinomial_rnd(data_size[t], prob[t]);
		
		return samples;
	}//end of multinomial_rnd method

	public static int[][] multinomial_rnd(int[] data_size, double[] prob) {
		int[][] samples = new int[data_size.length][];
		for(int t = 0; t < samples.length; t++)
			samples[t] = multinomial_rnd(data_size[t], prob);
		
		return samples;
	}//end of multinomial_rnd method
	
	
	public static int[][] multinomial_rnd(int num_rows, int num_cols, double[] prob) {
		int[][] samples = new int[num_rows][num_cols];
		for(int t = 0; t < num_rows; t++)
			samples[t] = multinomial_rnd(num_cols, prob);
		
		return samples;
	}//end of multinomial_rnd method
	
	public static int[] multinomial_rnd(int num_cols, double[] prob) {
		int[] samples = new int[num_cols];
		for(int n = 0; n < num_cols; n++)
			samples[n] = multinomial_rnd(prob);
		
		return samples;
	}//end of multinomial_rnd method
	
	/**
	 * 
	 * @param prob multinomial probability. Has to sum up to 1.
	 * @return
	 */
	public static int multinomial_rnd(double[] prob) {
  		StatUtil.normalize(prob);
		
		double[] cumsum = new double[prob.length-1];
		double sum = 0.0;
		for(int i = 0; i < prob.length-1; i++) { 
			  sum         += prob[i];
			  cumsum[i] = sum;
		}
		
		double unif_sample = Math.random();
		int i;
		for(i = 0; i < cumsum.length; i++)
			if(unif_sample < cumsum[i])
				break;
			
		return i;		
	}//end of multinomial_rnd method
	
	
	
	/***************************
	 **  MARKOV CHAIN SAMPLES **
	 **************************/
	
	
	/**
	 * Creates samples from a Markov Chain for <tt>data_size.length</tt> data sets, each having
	 * <tt>data_size[t]</tt> samples and having a <tt>init_prob[t]</tt> initial probability
	 * and a transition matrix <tt>trans_mat[t]</tt>. 
	 */
	public static int[][] mc_rnd(int[] data_size, double[][] init_prob, double[][][] trans_mat) {
		int[][] samples = new int[data_size.length][];
		for(int t = 0; t < samples.length; t++)
			samples[t] = mc_rnd(data_size[t], init_prob[t], trans_mat[t]);
		
		return samples;
	}//end of mc_rnd method
	
	public static int[][] mc_rnd(int[] data_size, double[][] init_prob, double[][][][] trans_mat) {
		int[][] samples = new int[data_size.length][];
		for(int t = 0; t < samples.length; t++)
			samples[t] = mc_rnd(data_size[t], init_prob[t], trans_mat[t]);
		
		return samples;
	}//end of mc_rnd method
	
	/**
	 * For inhomogeneous MCs
	 * @param num_rows
	 * @param num_cols
	 * @param init_prob
	 * @param trans_mat
	 * @return
	 */
	public static int[][] mc_rnd(int num_rows, int num_cols, double[] init_prob, double[][][] trans_mat) {
		int[][] samples = new int[num_rows][];
		for(int t = 0; t < samples.length; t++)
			samples[t] = mc_rnd(num_cols, init_prob, trans_mat);
		
		return samples;
	}//end of mc_rnd method

	
	/**
	 * Creates samples from a Markov Chain for <tt>num_rows</tt> data sets, each having
	 * <tt>num_cols</tt> samples and having a <tt>init_prob</tt> initial probability
	 * and a transition matrix <tt>trans_mat</tt>.
	 * @param num_rows
	 * @param num_cols
	 * @param init_prob
	 * @param trans_mat
	 * @return
	 */
	public static int[][] mc_rnd(int num_rows, int num_cols, double[] init_prob, double[][] trans_mat) {
		int[][] samples = new int[num_rows][];
		for(int t = 0; t < samples.length; t++)
			samples[t] = mc_rnd(num_cols, init_prob, trans_mat);
		
		return samples;
	}//end of mc_rnd method

	public static int[][] mc_rnd(int[] data_size, double[] init_prob, double[][] trans_mat_homog_part, double[][] trans_mat_inhomog_part) {
		int[][] samples = new int[data_size.length][];
		for(int t = 0; t < samples.length; t++)
			samples[t] = mc_rnd(data_size[t], init_prob, trans_mat_homog_part, trans_mat_inhomog_part);
		
		return samples;
	}//end of mc_rnd method
	
	public static int[][] mc_rnd(int num_rows, int num_cols, double[] init_prob, double[][] trans_mat_homog_part, double[][] trans_mat_inhomog_part) {
		int[][] samples = new int[num_rows][];
		for(int t = 0; t < samples.length; t++)
			samples[t] = mc_rnd(num_cols, init_prob, trans_mat_homog_part, trans_mat_inhomog_part);
		
		return samples;
	}//end of mc_rnd method
	
	/**
	 * Returns samples from an inhomogeneous Markov Chain.						<br>
	 * @param num_cols
	 * @param init_prob
	 * @param trans_mat_homog_part the constant (homogeneous) part of the transition matrix
	 * @param trans_mat_inhomog_part the non-constant (inhomogeneous) part of the transition matrix (last row)
	 * @return
	 */
	public static int[] mc_rnd(int num_cols, double[] init_prob, double[][] trans_mat_homog_part, double[][] trans_mat_inhomog_part) {
		int[] samples = new int[num_cols];
		samples[0] =  multinomial_rnd(init_prob);
		for(int n = 1; n < samples.length; n++) {
			//construct transition matrix
			double[][] trans_mat = new double[trans_mat_homog_part.length+1][trans_mat_homog_part[0].length];
			for(int k = 0; k < trans_mat.length-1; k++)
				trans_mat[k] = trans_mat_homog_part[k];
			
			trans_mat[trans_mat.length-1] = trans_mat_inhomog_part[n];
			
			samples[n] = multinomial_rnd(trans_mat[samples[n-1]]);
		}
		
		return samples;
	}//end of mc_rnd method 
	
	public static int[][] mc_rnd(int[] data_size, int[][] close_read_inds, double[] init_prob, double[][] trans_mat1_homog_part, double[][] trans_mat1_inhomog_part, double[][] trans_mat2) {
		int[][] samples = new int[data_size.length][];
		for(int t = 0; t < data_size.length; t++)
			samples[t] = mc_rnd(data_size[t], close_read_inds[t], init_prob, trans_mat1_homog_part, trans_mat1_inhomog_part, trans_mat2);
		
		return samples;
	}//end of mc_rnd method
	
	public static int[][] mc_rnd(int num_rows, int num_cols, int[][] close_read_inds, double[] init_prob, double[][] trans_mat1_homog_part, double[][] trans_mat1_inhomog_part, double[][] trans_mat2) {
		int[][] samples = new int[num_rows][];
		for(int t = 0; t < num_rows; t++)
			samples[t] = mc_rnd(num_cols, close_read_inds[t], init_prob, trans_mat1_homog_part, trans_mat1_inhomog_part, trans_mat2);
		
		return samples;
	}//end of mc_rnd method

	
	/**
	 * 
	 * @param num_cols N
	 * @param close_read_inds 1xK (K<N) indices for reads that are within distance <tt>d</tt> from their previous read
	 * @param init_prob 1xM
	 * @param trans_mat1_homog_part (M-1)xM Homogeneous transition matrix part when adjacent reads are very close together (<=d)
	 * @param trans_mat1_inhomog_part KxM Inomogeneous transition matrix part when adjacent reads are very close together (<=d)
	 * @param trans_mat2 Homogeneous transition matrix when adjacent reads are far away (>d)
	 * @return
	 */
	public static int[] mc_rnd(int num_cols, int[] close_read_inds, double[] init_prob, double[][] trans_mat1_homog_part, double[][] trans_mat1_inhomog_part, double[][] trans_mat2) {
		int[] samples = new int[num_cols];
		samples[0] =  multinomial_rnd(init_prob);
		
		Arrays.sort(close_read_inds);
		int n = 1;
		for(int i = 0; i < close_read_inds.length; i++) {
			
			// adjacent reads are far away
			while(n < close_read_inds[i] && n < num_cols) {
				samples[n] = multinomial_rnd(trans_mat2[samples[n-1]]);
				n++;
			}
			
			if( n < num_cols ) {
				// When n stops being less than close_read_inds[i], it becomes certainly
				// equal to it. So, we take trans_mat1
				//construct transition matrix 1 (for close adjacent reads)
				double[][] trans_mat1 = new double[trans_mat1_homog_part.length+1][trans_mat1_homog_part[0].length];
				for(int k = 0; k < trans_mat1.length-1; k++)
					trans_mat1[k] = trans_mat1_homog_part[k];
				
				trans_mat1[trans_mat1.length-1] = trans_mat1_inhomog_part[i];
				
				samples[n] = multinomial_rnd(trans_mat1[samples[n-1]]);
				n++;				
			}
			else {
				break;
			}
			
		}
		
		while(n < num_cols) {
			samples[n] = multinomial_rnd(trans_mat2[samples[n-1]]);
			n++;
		}
		
		return samples;
	}//end of mc_rnd method 

	/*
	 * Creates <tt>num_cols<tt> samples of an inhomogeneous Markov Chain, with an initial probability
	 * <tt>init_prob</tt> and a transition matrix <tt>trans_mat</tt>.
	 */
	public static int[] mc_rnd(int num_cols, double[] init_prob, double[][][] trans_mat) {
		int[] samples = new int[num_cols];
		
		samples[0] =  multinomial_rnd(init_prob);
		for(int n = 1; n < samples.length; n++) { samples[n] = multinomial_rnd(trans_mat[n][samples[n-1]]); }
		
		return samples;
	}//end of mc_rnd method

	
	/*
	 * Creates <tt>num_cols<tt> samples of a Markov Chain, with an initial probability
	 * <tt>init_prob</tt> and a transition matrix <tt>trans_mat</tt>.
	 */
	public static int[] mc_rnd(int num_cols, double[] init_prob, double[][] trans_mat) {
		int[] samples = new int[num_cols];
		
		samples[0] =  multinomial_rnd(init_prob);
		for(int n = 1; n < samples.length; n++) { samples[n] = multinomial_rnd(trans_mat[samples[n-1]]); }
		
		return samples;
	}//end of mc_rnd method

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		/*
		double[] prob = {.2, .5, .3};
		int[] samples = Sampler.multinomial_sample(prob, 1000);
		for(int i = 0; i < samples.length; i++)
			System.out.print(samples[i]);
		
		double count0 = 0;
		double count1 = 0;
		double count2 = 0;
		for(int i = 0; i < samples.length; i++) {
			if(samples[i] == 0) {count0++;}
			else if(samples[i] == 1) {count1++;}
			else {count2++;}
		}
		System.out.println();
		System.out.printf("0: %.2f, 1: %.2f, 2: %.2f\n", count0/samples.length, count1/samples.length, count2/samples.length);
		*/
		
		/*
		double[] init_prob =   {.8, .2};
		double[][] trans_mat = {{.6, .4},
				               {.3, .7}};
		
		int num_rows = 1000; int num_cols = 1000;
		int[][] samples = Sampler.mc_sample(num_rows, num_cols, init_prob, trans_mat);
		double count_init_0 = 0;
		double count00 = 0; double count01 = 0; double count10=0; double count11 = 0;
		for(int t = 0; t < num_rows; t++) {
			if(samples[t][0] == 0) { count_init_0++; }
			for(int n = 1; n < num_cols; n++) {
				if(samples[t][n-1] == 0) {
					if(samples[t][n] == 0)
						count00++;
					else
						count01++;
				}
				else {
					if(samples[t][n] == 0)
						count10++;
					else
						count11++;
				}
			}
		}
		
		double sum0 = count00+count01;
		double sum1 = count10+count11;
		
		System.out.printf("init_0: %.2f, init_1: %.2f\n\n", count_init_0/num_rows, (1-count_init_0/num_rows));
		
		System.out.printf("a0->0: %.2f, a0->1: %.2f\n" +
				          "a1->0: %.2f, a1->1: %.2f\n", count00/sum0, count01/sum0, count10/sum1, count11/sum1);
	 */
	
	/*
	 double mu = 10.0;
	 double std = 4.0;
	 double[] samples = Sampler.gauss_rnd(20, mu, std);
	 */
		
	 /*
	 double[] mu  = {5, 50, 1000};
	 double[] std = {1, 10, 200};
	 int[] num_samples = {5, 10, 50};
	 
	 //double[][] samples = Sampler.gauss_rnd(num_samples, mu, std);
	 double[][] prob = {
			            {.1, .9},
			            {.7, .3},
			            {.5, .5}
	 					};
	// int[][] sampler = Sampler.multinomial_rnd(num_samples, prob);
	 
	 
	 double[][] init_prob = {
			 			   {.2, .8},
			 			   {.5, .5},
			 			   {.35, .65}
	                      };
	 double[][][] trans_mat = {
			 					{{.8, .2},
			 					{.6, .4}}
			 					,
			 					{{.6, .4},
				 				{.5, .5}}
			 					,
								{{.1, .9},
								{.2, .8}}
	 						  };
	 
	 
	 int[][] sampler = Sampler.mc_rnd(num_samples, init_prob, trans_mat);
	 */
		
	int[] unif_samples = unif_rnd(20, -3, 1);
	int[] binom_samples = binom_rnd(1000000, 25, .99);
	
	File f = new File("/Users/gio_fou/Desktop/binom_samples2.txt");
	BufferedWriter bw = null;
	try {
		bw = new BufferedWriter(new FileWriter(f));
		
		for(int n = 0; n < binom_samples.length; n++)
			bw.write(binom_samples[n] + "\t");
	}
	catch (IOException e) { e.printStackTrace(); }
	finally {
		try {
			bw.close();
		} 
		catch (IOException e) { e.printStackTrace(); }
	}


		
	
	int num_cols = 20;
	
	int[] close_read_inds = {2, 4, 6};
	
	double[] init_prob = {.4, .5, .1};
	
	double[][] trans_mat1_homog_part = {
										{.4, .3, .3},
										{.2, .2, .6}
								  	   };
	
	double[][] trans_mat1_inhomog_part = {
									 	  {.3, .6, .1},
									 	  {.8, .1, .1},
									 	  {.2, .1, .8}
			  							};
	
	double[][] trans_mat2 = {
							 {.1, .7, .2},
							 {.1, .6, .3},
							 {.7, .2, .1}
							};
 
	
	int[] samples1 = mc_rnd(num_cols, close_read_inds, init_prob, trans_mat1_homog_part, trans_mat1_inhomog_part, trans_mat2);
	
	 System.out.println("Fin");

	}//end of main method

}

package edu.mit.csail.cgs.utils.probability;

public class NormalDistributionGio {
	
	private double mean, variance, stddev;
	private double coeff, logcoeff;
    private double denom;
    
    /**
     * constructor
     * @param m mean
     * @param v variance
     */
    public NormalDistributionGio(Double m, Double v)
    {
    	mean = m;
    	variance = v;
    	coeff = 1.0 / Math.sqrt(2.0 * Math.PI * variance);
        logcoeff = (Math.log(Math.PI) + Math.log(variance) + Math.log(2.0)) / -2.0;
        denom = 2.0 * variance;
    	stddev = Math.sqrt(variance);
    }//end of constructor method
    
    /**
     * <tt>value</tt> is following a normal distribution <tt>N(m,v)</tt>
     * @param value the value whose z-score will be calculated
     * @return the z-score of this <tt>value</tt>
     */
    public double calcZScore(Double value) { 
    	double diff = value - mean;
    	return Math.abs(diff / stddev);
    }//end of calcZScore method
    
    /**
     * <tt>value</tt> is following a normal distribution <tt>N(m,v)</tt>
     * @param value the value whose likelihood will be calculated
     * @return the likelihood of this <tt>value</tt>
     */
    public double calcProbability(Double value) {
    	return Math.exp(calcLogProbability(value));
    }//end of calcProbability method
    
    /**
     * <tt>value</tt> is following a normal distribution <tt>N(m,v)</tt>
     * @param value the value whose log-likelihood will be calculated
     * @return the log-likelihood of this <tt>value</tt>
     */
    public double calcLogProbability(Double value) { 
        double diff = value - mean;
        double exponent = - (diff*diff) / denom;
        return exponent + logcoeff;
    }//end of calcLogProbability method
    
    /**
     * static method when we don't initialize the constructor<br>
     * <tt>value</tt> is following a normal distribution <tt>N(mean,var)</tt> 
     * @param value the last value of a set of <tt>n</tt> data points whose (the set) total log-likelihood will be calculated
     * @param n # of data points
     * @param prevMean previous mean (of <tt>n-1</tt> data points)
     * @param prevVar previous variance (of <tt>n-1</tt> data points)
     * @return the likelihood of the set of <tt>n</tt> data points
    */
    public double dynCalcTotalProbability(Double value, int n, Double prevMean, Double prevVar) {
    	return Math.exp(dynCalcTotalLogProbability(value, n, prevMean, prevVar));
    }//end of dynCalcProbability method
    
    /**
     * <b>Total Log-likelihood</b><br>
     * <pre>
     *                  1    n
     * mean_n = mu_n = ---* Sum {y_i}
     *                  n   i=1
     * 
     *                        1    n
     * var_n = (sigma_n)^2 = ---* Sum {(y_i - mean_n)^2} , where: sigma_n: standard deviation
     *                        n   i=1
     *           
     *                      n                         1
     * total-log-lik_n = - --- * log(2*pi*var_n) - ------- * [(n-1)*( var_(n-1) + ( mean_n - mean_(n-1) )^2 ) + (y_n-mean_n)^2]
     *                      2                      2*var_n
     * </pre>
     * where:
     * <pre>
     *                    n                 1            (y_i - mean_n)^2
     * total-log-lik_n = Sum {  log(-----------------) - -----------------  }
     *                   i=1         sqrt(2*pi*var_n)         2*var_n
     * </pre>
     *  
     * static method when we don't initialize the constructor<br>
     * <tt>value</tt> is following a normal distribution <tt>N(mean,var)</tt> 
     * @param value the last value of a set of <tt>n</tt> data points whose (the set) total log-likelihood will be calculated
     * @param n # of data points
     * @param prevMean previous mean (of <tt>n-1</tt> data points)
     * @param prevVar previous variance (of <tt>n-1</tt> data points)
     * @return the log-likelihood of the set of <tt>n</tt> data points
     */
    public double dynCalcTotalLogProbability(Double value, int n, Double prevMean, Double prevVar) { 
        return -((double)n/2)*Math.log(2.0*Math.PI*variance) - (1/(2.0*variance))*((n-1)*(prevVar + (mean-prevMean)*(mean-prevMean)) + (value-mean)*(value-mean));
    }//end of dynCalcLogProbability method
    
      
    /**********************
     *** STATIC METHODS ***
     **********************/
    
    /**
     * static method when we don't initialize the constructor<br>
     * <tt>value</tt> is following a normal distribution <tt>N(mean,var)</tt> 
     * @param value value the value whose z-score will be calculated
     * @param mean mean
     * @param var variance
     * @return the z-score of this <tt>value</tt>
     */
    public static double calcZScore(Double value, Double mean, Double var)
    {
    	double diff = value - mean;
    	double stddev = Math.sqrt(var);
    	return Math.abs(diff / stddev);
    }//end of static calcZScore method
    
    /**
     * static method when we don't initialize the constructor<br>
     * <tt>value</tt> is following a normal distribution <tt>N(mean,var)</tt>
     * @param value the value whose likelihood will be calculated
     * @param mean mean
     * @param var variance
     * @return the likelihood of this <tt>value</tt>
     */
    public static double calcProbability(Double value, Double mean, Double var) {
    	return Math.exp(calcLogProbability(value, mean, var));
    }//end of static calcProbability method
    
    /**
     * static method when we don't initialize the constructor<br>
     * <tt>value</tt> is following a normal distribution <tt>N(mean,var)</tt> 
     * @param value the value whose log-likelihood will be calculated
     * @param mean mean
     * @param var variance
     * @return the log-likelihood of this <tt>value</tt>
     */
    public static double calcLogProbability(Double value, Double mean, Double var) {
    	double logcoeff = (Math.log(Math.PI) + Math.log(var) + Math.log(2.0)) / -2.0;
    	double denom = 2.0*var;
    	double diff = value - mean;
    	double exponent = - (diff*diff) / denom;
    	return exponent + logcoeff;
    }//end of static calcLogProbability method
    
   /**
     * static method when we don't initialize the constructor<br>
     * <tt>value</tt> is following a normal distribution <tt>N(mean,var)</tt> 
     * @param value the last value of a set of <tt>n</tt> data points whose (the set) total log-likelihood will be calculated
     * @param n # of data points
     * @param mean current mean (of <tt>n</tt> data points)
     * @param var current variance (of <tt>n</tt> data points)
     * @param prevMean previous mean (of <tt>n-1</tt> data points)
     * @param prevVar previous variance (of <tt>n-1</tt> data points)
     * @return the likelihood of the set of <tt>n</tt> data points
    */
    public static double dynCalcTotalProbability(Double value, int n, Double mean, Double var, Double prevMean, Double prevVar) {
    	return Math.exp(dynCalcTotalLogProbability(value, n, mean, var, prevMean, prevVar));
    }//end of static dynCalcProbability method
    
    /**
     * <b>Total Log-likelihood</b><br>
     * <pre>
     *                  1    n
     * mean_n = mu_n = ---* Sum {y_i}
     *                  n   i=1
     * 
     *                        1    n
     * var_n = (sigma_n)^2 = ---* Sum {(y_i - mean_n)^2} , where: sigma_n: standard deviation
     *                        n   i=1
     *           
     *                      n                         1
     * total-log-lik_n = - --- * log(2*pi*var_n) - ------- * [(n-1)*( var_(n-1) + ( mean_n - mean_(n-1) )^2 ) + (y_n-mean_n)^2]
     *                      2                      2*var_n
     * </pre>
     * where:
     * <pre>
     *                    n                 1            (y_i - mean_n)^2
     * total-log-lik_n = Sum {  log(-----------------) - -----------------  }
     *                   i=1         sqrt(2*pi*var_n)         2*var_n
     * </pre>
     * 
     * static method when we don't initialize the constructor<br>
     * <tt>value</tt> is following a normal distribution <tt>N(mean,var)</tt> 
     * @param value the last value of a set of <tt>n</tt> data points whose (the set) total log-likelihood will be calculated
     * @param n # of data points
     * @param mean current mean (of <tt>n</tt> data points)
     * @param var current variance (of <tt>n</tt> data points)
     * @param prevMean previous mean (of <tt>n-1</tt> data points)
     * @param prevVar previous variance (of <tt>n-1</tt> data points)
     * @return the log-likelihood of the set of <tt>n</tt> data points
     */
    public static double dynCalcTotalLogProbability(Double value, int n, Double mean, Double var, Double prevMean, Double prevVar) { 
        return -((double)n/2)*Math.log(2.0*Math.PI*var) - (1/(2.0*var))*((n-1)*(prevVar + (mean-prevMean)*(mean-prevMean)) + (value-mean)*(value-mean));
    }//end of static dynCalcLogProbability method
	
}//end of NormalDistributionGio class

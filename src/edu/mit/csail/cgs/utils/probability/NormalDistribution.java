/*
 * Created on Feb 1, 2008
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.utils.probability;

import cern.jet.random.Normal;
import cern.jet.random.engine.DRand;

public class NormalDistribution {
    
    private double mean, variance, stddev;
    private double coeff, logcoeff;
    private double denom;
    
    public static double cdf(double mean, double std, double x){
    	DRand engine = new DRand();
    	Normal normal = new Normal(mean, std, engine);
    	return normal.cdf(x);
    }
    
    public NormalDistribution(double m, double v) { 
        mean = m;
        variance = v;
        coeff = 1.0 / Math.sqrt(2.0 * Math.PI * variance);
        logcoeff = (Math.log(Math.PI) + Math.log(variance) + Math.log(2.0)) / -2.0;
        denom = 2.0 * variance;
        stddev = Math.sqrt(variance);
    }
    
    public double calcZScore(Double value) { 
    	double diff = value - mean;
    	return Math.abs(diff / stddev);
    }

    public double calcProbability(Double value) { 
        double diff = value - mean;
        double exponent = - (diff*diff) / denom;
        return coeff * Math.exp(exponent);
    }

    public double calcLogProbability(Double value) { 
        double diff = value - mean;
        double exponent = - (diff*diff) / denom;
        return exponent + logcoeff;
    }
    
    public String toString() { return String.format("[Gaussian, mean:%.3f, sd:%.3f]", mean, stddev); }
}

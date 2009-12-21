/*
 * Created on Feb 1, 2008
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.utils.probability;

import cern.jet.random.Exponential;
import cern.jet.random.engine.DRand;

public class ExponentialDistribution {
    
    private double lambda;
    private Exponential expt;
    
    public ExponentialDistribution(double lmbda) { 
        lambda = lmbda;
        expt = new Exponential(lambda, new DRand());
    }
    
    public double calcProbability(Double value) {
        return expt.pdf(value);
    }

    public double calcLogProbability(Double value) {
        return Math.log(expt.pdf(value));
    }
    
    public String toString() { 
        return String.format("[Exponential, lambda:%f]", lambda); 
    }
    
    public double calculateCumulProbability(Double value) { 
        return expt.cdf(value);
    }
}

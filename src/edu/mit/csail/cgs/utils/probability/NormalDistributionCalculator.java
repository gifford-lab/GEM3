/*
 * Created on Feb 7, 2008
 */
package edu.mit.csail.cgs.utils.probability;

import java.util.*;

public class NormalDistributionCalculator {

    private LinkedList<Double> values;
    private double sum;
    
    public NormalDistributionCalculator() { 
        values = new LinkedList<Double>();
        sum = 0.0;
    }
    
    public NormalDistribution calculateDistribution() {
        double divisor = (double)Math.max(1, values.size());
        double mean = sum / divisor;
        double varsum = 0.0;
        for(double v : values) { 
            double diff = v - mean;
            varsum += (diff*diff);
        }
        double var = values.size() >= 1 ? varsum / divisor : 1.0;
        return new NormalDistribution(mean, var);
    }
    
    public void reset() { 
        values.clear();
        sum = 0.0;
    }
    
    public double getMean() { 
        return sum / (double)Math.max(1, values.size());
    }
    
    public void addValues(Collection<Double> vs) { 
        for(double v : vs) { addValue(v); }
    }
    
    public void addValue(double v) { 
        if(Double.isNaN(v) || Double.isInfinite(v)) {
            String msg = Double.isNaN(v) ? "NaN" : "Infinite";
            throw new IllegalArgumentException(msg);
        }
        
        values.add(v);
        sum += v;
    }
}

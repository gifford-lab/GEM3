/*
 * Created on Mar 5, 2006
 */
package edu.mit.csail.cgs.utils.probability;

import java.util.*;
import java.text.*;

/**
 * @author tdanford
 */
public class FiniteDistribution {
    
    private static Random rand;
    private static NumberFormat nf;
    
    static { 
        rand = new Random();
        nf = DecimalFormat.getInstance();
        nf.setMaximumFractionDigits(4);
    }
    
    private double[] vals;

    public FiniteDistribution(int s) {
        if(s <= 0) { throw new IllegalArgumentException("non-zero length req: " + s); }
        vals = new double[s];
        clear(1.0);
        normalize();
    }
    
    public FiniteDistribution(Collection<Integer> w) { 
        vals = new double[w.size()];
        int i = 0; 
        for(int v : w) { 
            vals[i++] = (double)v;
        }
        normalize();
    }
    
    public FiniteDistribution(double[] probs) { 
    	vals = probs.clone();
    	normalize();
    }
    
    public FiniteDistribution(int s, int v) { 
    	this(s);
    	clear(0.0);
    	vals[v] = 1.0;
    }
    
    public int size() { return vals.length; }
    public double getProb(int i) { return vals[i]; }
    
    public boolean equals(Object o) { 
        if(!(o instanceof FiniteDistribution)) { 
            return false;
        }
        FiniteDistribution fd = (FiniteDistribution)o;
        if(vals.length != fd.vals.length) { return false; }
        for(int i = 0; i < vals.length; i++) { 
            if(vals[i] != fd.vals[i]) { return false; }
        }
        
        return true;
    }
    
    public int hashCode() { 
        int code = 17;
        
        for(int i = 0; i < vals.length; i++) { 
            long bits = Double.doubleToLongBits(vals[i]); 
            code += (int)(bits >> 32); code *= 37;
        }
        
        return code;
    }
    
    public String toString() { 
        StringBuilder sb = new StringBuilder();
        sb.append("[");
        for(int i = 0; i < vals.length; i++) { 
            sb.append(nf.format(vals[i]));
            if(i < vals.length-1) { sb.append(" "); }
        }
        sb.append("]");
        return sb.toString();
    }
    
    public int sampleValue() { 
        double p = rand.nextDouble();
        for(int i = 0; i < vals.length; i++) { 
            p -= vals[i];
            if(p <= 0.0) { return i; }
        }
        return vals.length-1;
    }
    
    public void clear() { 
    	clear(0.0);
    }
    
    public void clear(double clearValue) { 
    	for(int i = 0; i < vals.length; i++) { 
    		vals[i] = clearValue;
    	}
    }
    
    public void addValue(Integer i) { 
    	addValue(i, 1.0);
    }
    
    public void addValue(Integer i, Double w) { 
    	vals[i] += w;
    }

    public void normalize() { 
        double sum = 0.0;
        for(int i = 0; i < vals.length; i++) {
            sum += vals[i];
        }
        
        for(int i = 0; i < vals.length; i++) { 
            if(sum > 0.0) { 
                vals[i] /= sum;
            } else { 
                vals[i] = 1.0 / (double)vals.length;
            }
        }
    }
}

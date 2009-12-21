/*
 * Created on Jun 9, 2005
 */
package edu.mit.csail.cgs.utils.probability;

import java.util.Vector;
import java.io.*;

/**
 * @author tdanford
 */
public class Hypergeometric {
    
    public static void main(String[] args) { 
        int birthdays = Integer.parseInt(args[0]);
        double prob = 0.0;
        for(int i = 1; prob < 0.5 && i < birthdays; i++) { 
            prob = birthday_prob(birthdays, i);
            System.out.println(i + ": " + prob);
        }
    }
    
    public static void old_main(String[] args) { 
        int total = Integer.parseInt(args[0]);
        int theta = Integer.parseInt(args[1]);
        int sample = Integer.parseInt(args[2]);
        int x = Integer.parseInt(args[3]);
        Hypergeometric h = new Hypergeometric();
        double log_prob = h.log_hypgeomPValue(total, theta, sample, x);
        double prob = Math.exp(log_prob);
        System.out.println("Prob: " + prob);
    }
    
    private Vector<Double> logFactorial;
    
    public Hypergeometric() { 
        logFactorial = new Vector<Double>();
		logFactorial.setSize(21);
        logFactorial.set(0, -Double.MAX_VALUE);
        double sum = 0.0;
        for(int i = 1; i <= 20; i++) { 
            logFactorial.set(i, sum);
            sum += Math.log((double)(i + 1));
        }
    }
    
    public static double birthday_prob(int birthdays, int people) { 
        double log_first = 0.0;        
        for(int i = birthdays; i > birthdays-people; i--) { 
            log_first += Math.log((double)i);
        }
        
        double log_second = (double)people * Math.log((double)birthdays);
        double total = log_first - log_second;
        double prob = 1.0 - Math.exp(total);
        return prob;
    }
    
    private void extendLogFactorial(int max) {
        int currentMax = logFactorial.size() - 1;
        if(max <= currentMax) { return; }
		logFactorial.setSize(max+1);
        double sum = logFactorial.get(currentMax);
        sum += Math.log((double)(currentMax + 1));
        for(int i = currentMax + 1; i <= max; i++) { 
            logFactorial.set(i, sum);
            sum += Math.log((double)(i + 1));            
        }
    }
    
    public double log_add(double v1, double v2) { 
        //if(v1 > v2) { 
            return v1 + Math.log((double)1.0 + Math.exp(v2 - v1));
        //} else {
            //return log_add(v2, v1);
        //}
    }
    
    public double log_factorial(int n) { 
        extendLogFactorial(n);
        return logFactorial.get(n);
    }
    
    public double log_choose(int N, int x) { 
        if(x < 0 || N < 0 || x > N) { 
			throw new IllegalArgumentException("params: " + N + " " + x); 
		}
        if(x == N || x == 0 || N == 0) { return 0.0; }
        extendLogFactorial(N);
        double numer = logFactorial.get(N);
        double denom = logFactorial.get(N - x) + logFactorial.get(x);
        return numer - denom;
    }
    
    public double log_hypgeom(int N, int theta, int n, int x) {
        /**
         * factors: 
         * numer1: (theta choose x)
         * numer2: (N - theta choose n - x)
         * denom: (N choose n)
         */
        double n1 = log_choose(theta, x);
        double n2 = log_choose(N-theta, n-x);
        double d = log_choose(N, n);
        return (n1 + n2) - d;
    }
    
    public double log_hypgeomPValue(int N, int theta, int n, int x) {
        try { 
            int max = n; 
            if(theta < max) { max = theta; }
            double log_sum = 0.0;
            for(int i = x; i <= max; i++) { 
                if(i == x) { 
                    log_sum = log_hypgeom(N, theta, n, i);
                } else { 
                    log_sum = log_add(log_sum, log_hypgeom(N, theta, n, i));
                }
            }
            return log_sum;
        } catch(IllegalArgumentException iae) {
            String error = "Error occurred with parameters: " + N + ", " + theta + ", " + n + ", " + x;
            throw new IllegalArgumentException(error);
        }
    }
}



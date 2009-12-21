package edu.mit.csail.cgs.utils.probability;

import java.lang.*;
import java.io.*;
import java.util.*;
import java.text.*;
import java.math.*;

import edu.mit.csail.cgs.utils.numeric.Numerical;

public abstract class Binomial { 
    
    public static NumberFormat nf = new DecimalFormat("##0.##E0");

    public static void main(String[] args) { 
        int n = Integer.parseInt(args[0]);
        int r = Integer.parseInt(args[1]);
        double theta = Double.parseDouble(args[2]);
        double sig = log_binomial_significance(r, n, theta);
        System.out.println("N: " + n + ", R:" + r + ", Theta: " + theta);
        System.out.println("Sig: " + Math.exp(sig));
    }

	public static double log_sum(double lv1, double lv2) { 
        if (Double.isNaN(lv2) || Double.isNaN(lv1)) {
            throw new ArithmeticException("Can't have NaN args");
        }
		if(lv2 >= lv1) { 
            double d = lv2 - lv1;
            double ed = Math.exp(d);
            if (Double.isInfinite(ed)) {
                return Math.max(lv1,lv2);
            } else {
                return lv1 + Math.log(1.0 + ed);
            }
		} else {
			return log_sum(lv2, lv1);
		}
	}

	public static double log_binomial_significance(int r, int n, double theta) { 
		if(theta <= 0.0 || theta >= 1.0 ||
			n < 1 || r < 0 || r > n) { 
				System.err.println("n: " + n);
				System.err.println("r: " + r);
				System.err.println("theta: " + theta);
				throw new IllegalArgumentException(); 
		}
		if(r == 0) { return 1.0; }

		double lsum = log_binomial(r, n, theta);
		for(int x = r+1; x <= n; x++) { 
            double lb = log_binomial(x, n, theta);
			lsum = log_sum(lsum, lb);
		}
		return lsum;
	}

	public static double log_binomial(int r, int n, double theta) { 
		if(theta <= 0.0 || theta >= 1.0 ||
			n < 1 || r < 0) { 
			System.err.println("n: " + n);
			System.err.println("r: " + r);
			System.err.println("theta: " + theta);
			throw new IllegalArgumentException(); 
		}
		if(r > n) { return Double.NEGATIVE_INFINITY; }
		
		double coeff = Numerical.log_bico(n, r);
		double f1 = (double)r * Math.log(theta);
		double f2 = (double)(n-r) * Math.log(1.0 - theta);
		return coeff + f1 + f2;
	}
}

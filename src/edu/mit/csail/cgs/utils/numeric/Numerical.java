package edu.mit.csail.cgs.utils.numeric;

import java.lang.*;
import java.io.*;
import java.util.*;

import cern.colt.matrix.DoubleMatrix1D;
import cern.jet.math.Arithmetic;
import cern.jet.math.Functions;

import edu.mit.csail.cgs.utils.Pair;

/* 
 * Implementations of code from Numerical Recipes, for
 * general use.
 */

public abstract class Numerical {    
  
  public static class NumericalException extends RuntimeException {
    public NumericalException() { super(); }
    public NumericalException(String m) { super(m); }
  }

  public static double log_add(double v1, double v2) {
    if(Double.isNaN(v1)) { return v2; }
    if(Double.isNaN(v2)) { return v1; }
    if(v1 >= v2) { 
      return v1 + Math.log((double)1.0 + Math.exp(v2 - v1));
    } else {
      return log_add(v2, v1);
    }
  }
  
  public static double log_subtract(double v1, double v2) {
    return v1 + Math.log((double)1.0 - Math.exp(v2 - v1));
  }

  public static double log_sum_exp(DoubleMatrix1D exponents) {
    double maxExp = exponents.aggregate(Functions.max, Functions.identity);
    if (maxExp == Double.NEGATIVE_INFINITY) {
      return Double.NEGATIVE_INFINITY;
    }
    else {
      double expSum = exponents.aggregate(Functions.plus, Functions.chain(Functions.exp, Functions.minus(maxExp)));
      return maxExp + Math.log(expSum);
    }
  }

  public static double log_sum_exp(double[] exponents) {
    if (exponents.length >= 1) {
      double maxExponent = exponents[0];
      for (int i = 1; i < exponents.length; i++) {
        maxExponent = Math.max(maxExponent, exponents[i]);
      }
      return Numerical.log_sum_exp(exponents, maxExponent);
    }
    else {
      throw new IllegalArgumentException("array has zero-length");
    }
  }

  public static double log_sum_exp(double[] exponents, double maxExp) {      
    if (maxExp == Double.NEGATIVE_INFINITY) {
      return Double.NEGATIVE_INFINITY;
    }
    if (exponents.length >= 1) {
      double expSum = 0;
      for (int i = 0; i < exponents.length; i++) {
        expSum = expSum + Math.exp(exponents[i] - maxExp);        
      }
      return maxExp + Math.log(expSum);
    }
    else {
      throw new IllegalArgumentException("array has zero-length");
    }
  }
    
    
    private static double gammln_cof[] = {
        76.18009172947146, -86.50532032941677, 
        24.01409824083091, -1.231739572450155, 
        0.1208650973866179e-2, -0.5395239384953e-5 };
    
    public static final double LOG_ZERO = Math.log(1.0e-99);
    
    private static final double EPS1 = 0.001;
    private static final double EPS2 = 1.0e-8;
    private static final int ITMAX = 100;
    private static final int MAXIT = 100;
    private static final double EPS = 3.0e-7;
    private static final double FPMIN = 1.0e-30;
    
    private static double factln_a[];
    
    static { 
        factln_a = new double[101];
        for(int i = 0; i < factln_a.length; i++) { factln_a[i] = (float)0.0; }
    }
    
    // pg. 625 (Numerical Recipes in C)
    public static Pair<Double,Double> kstwo(double[] data1, double[] data2) { 
        int j1 = 1, j2 = 1;
        double fn1 = 0.0, fn2 = 0.0;
        
        Arrays.sort(data1);
        Arrays.sort(data2);
        
        double en1 = (double)data1.length;
        double en2 = (double)data2.length;
        double d = 0.0;
        double d1, d2, dt;
        
        while(j1 <= data1.length && j2 <= data2.length) {
            d1 = data1[j1-1]; 
			d2 = data2[j2-1];
            if(d1 <= d2) {  
                fn1 = (double)j1 / en1;
                j1++;
            }
            
            if(d2 <= d1) { 
                fn2 = (double)j2 / en2;
                j2++;
            }
            
            if((dt = Math.abs(fn2-fn1)) > d) { 
                d = dt;
            }
        }
        
        double en = Math.sqrt(en1 * en2 / (en1 + en2));
        double prob = probks((en + 0.12 + 0.11/en) * d);
        
        return new Pair<Double,Double>(d, prob);
    }
    
    // pg. 626 (Numerical Recipes in C)
    public static double probks(double alam) { 
        
        double a2 = -2.0 * alam * alam;
        double sum = 0.0, termbf = 0.0, fac = 2.0;
        
        for(int j = 1; j <= 100; j++) {
            double dj = (double)j;
            double term = fac * Math.exp(a2 * dj * dj);
            sum += term;
            
            if(Math.abs(term) <= EPS1 * termbf ||
                    Math.abs(term) <= EPS2 * sum) { 
                return sum;
            }
            
            fac = -fac;
            termbf = Math.abs(term);
        }
        
        return 1.0;
    }
    
    // pg. 218 (Numerical Recipes in C)
    public static Pair<Double,Double> gser(double a, double x) {
        double gamser, gln;
        int n;
        double sum, del, ap;
        
        gln = gammln(a);
        if(x <= 0.0) { 
            if(x < 0.0) { throw new IllegalArgumentException(); }
            gamser = 0.0;
            return new Pair<Double,Double>(gamser, gln);
        } else { 
            ap = a;
            del = sum = 1.0 / a;
            for(n = 1; n <= ITMAX; n++) { 
                ap += 1.0;
                del *= x / ap;
                sum += del;
                if(Math.abs(del) < Math.abs(sum) * EPS) { 
                    gamser = sum * Math.exp(-x + a + Math.log(x) - gln);
                    return new Pair<Double,Double>(gamser, gln);
                }
            }
            throw new NumericalException("a too large, ITMAX too small");
        }
        
    }
    
    // pg. 617 in Numerical Recipes in C
    public static Pair<Double,Double> avevar(Collection<Double> values) { 
        double s, ep;
        double ave, var;
        ave = 0.0;
        int n = values.size();
        for(Double v : values) { ave += v; }
        ave /= (double)n;
        var = ep = 0.0;
        for(Double v : values) { 
            s = v - ave;
            ep += s;
            var += (s * s);
        }
        var = (var - ep * ep / (double)n) / (double)(n - 1);
        return new Pair<Double,Double>(ave, var);
    }
    
    // pg. 617-618 in Numerical Recipes in C
    public static Pair<Double,Double> tutest(Collection<Double> data1, Collection<Double> data2) { 
        Pair<Double,Double> results = null;
        double var1, var2, df, ave1, ave2;
        int n1 = data1.size(), n2 = data2.size();
        
        Pair<Double,Double> p1 = avevar(data1), p2 = avevar(data2);
        ave1 = p1.getFirst(); var1 = p1.getLast();
        ave2 = p2.getFirst(); var2 = p2.getLast();

        double temp1 = var1 / (double)n1, temp2 = var2 / (double)n2;
        double temp = temp1 + temp2;
        double stemp = temp * temp, stemp1 = temp1 * temp1, stemp2 = temp2 * temp2;
        
        double t = (ave1 - ave2) / Math.sqrt(temp);
        df = stemp / (stemp1 / (double)(n1 - 1) + stemp2 / (double)(n2-1));
        double prob = betai(0.5 * df, 0.5, df/(df + (t * t)));
        
        results = new Pair<Double,Double>(t, prob);
        return results;
    }
    
    // pg. 219 in Numerical Recipes in C.
    public static Pair<Double,Double> gcf(double a, double x) { 
        Pair<Double,Double> ret = null;
        
        double gammcf, gln;
        int i; 
        double an, b, c, d, del, h;
        
        gln = gammln(a);
        b = x + 1.0 - a;
        c = 1.0 / FPMIN;
        d = 1.0 / b;
        h = d;
        for(i = 1; i <= ITMAX; i++) { 
            an = -((double)i) * (((double)i)-a);
            b += 2.0;
            d = an * d + b;
            if(Math.abs(d) < FPMIN) { d = FPMIN; }
            c = b + an / c;
            if(Math.abs(c) < FPMIN) { c = FPMIN; }
            d = 1.0 / d;
            del = d * c;
            h *= del;
            if(Math.abs(del-1.0) < EPS) {
                break;
            }
        }
        if(i > ITMAX) { throw new NumericalException("a too large, ITMAX too small"); }
        gammcf = Math.exp(-x + a * Math.log(x) - gln) * h;
        ret = new Pair<Double,Double>(gammcf, gln);
        return ret;
    }
    
    // pg. 220 (Numerical Recipes in C)
    public static double erff(double x) {
        return x < 0.0 ? -gammp(0.5, x*x) : gammp(0.5, x*x);
    }
    
    // pg. 220 (Numerical Recipes in C)
    public static double erffc(double x) { 
        return x < 0.0 ? 1.0 + gammp(0.5, x*x) : gammq(0.5, x*x);
    }
    
    // pg. 221 (Numerical Recipes in C) 
    public static double erfcc(double x) { 
        double t, z, ans;
        z = Math.abs(x);
        t = 1.0 / (1.0 + 0.5 * z);
        ans = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + 
                t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 +
                        t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
        return x >= 0.0 ? ans : 2.0-ans;
    }
    
    // pg. 218 (Numerical Recipes in C)
    public static double gammp(double a, double x) {
        if(x < 0.0 || x <= 0.0) { throw new IllegalArgumentException(); }
        
        if(x < (a + 1.0)) { 
            Pair<Double,Double> p = gser(a, x);
            return p.getFirst();
        } else { 
            Pair<Double,Double> p = gcf(a, x);
            return 1.0 - p.getFirst();
        }
    }
    
    // pg. 218 (Numerical Recipes in C);
    public static double gammq(double a, double x) { 
        if(x < 0.0 || x <= 0.0) { throw new IllegalArgumentException(); }
        
        if(x < (a + 1.0)) { 
            Pair<Double,Double> p = gser(a, x);
            return 1.0-p.getFirst();
        } else { 
            Pair<Double,Double> p = gcf(a, x);
            return p.getFirst();
        }
    }
    
    // pg. 214 (Numerical Recipes in C)
    public static double gammln(double xx) { 
        double x, y, tmp, ser;
        int j;
        y = x = xx;
        tmp = x + 5.5;
        tmp -= (x + 0.5) * Math.log(tmp);
        ser = 1.000000000190015;
        for(j = 0; j <= 5; j++) { ser += gammln_cof[j] / ++y; }
        return -tmp + Math.log(2.5066282746310005 * ser / x);
    }
    
    // pg 215
    public static double factln(int n) { 
        if(n<0) { 
            throw new IllegalArgumentException("negative factorial: " + n); 
        }
        
        if(n <= 1) { return 0.0; }
        if(n <= 100) { 
            if(factln_a[n] <= 0.0) {
                factln_a[n] = gammln((float)n + 1.0);
            }
            return factln_a[n];
        } else { 
            return gammln((float)n + 1.0);
        }
    }
    
    public static double binomial(int draws, int hits, double theta) {
        if(theta < 0.0 || theta > 1.0) { throw new IllegalArgumentException(); }
        if(draws < 0 || hits < 0 || hits > draws) { throw new IllegalArgumentException(); }
        double hit_prob = 1.0, nohit_prob = 1.0;
        double notheta = 1.0 - theta;
        for(int i = 0; i < hits; i++) { hit_prob *= theta; }
        for(int i = 0; i < draws-hits; i++) { nohit_prob *= notheta; } 
        return bico(draws, hits) * hit_prob * nohit_prob;
    }
    
    public static double log_binomial(int draws, int hits, double theta) { 
        if(theta <= 0.0 || theta >= 1.0) { throw new IllegalArgumentException("theta: " + theta); }
        if(draws < 0 || hits < 0 || hits > draws) { 
            throw new IllegalArgumentException("d/h: " + draws + "," + hits); 
        }
        double log_theta = Math.log(theta);
        double log_notheta = Math.log(1.0 - theta);
        double hit_prob = 0.0, nohit_prob = 0.0;
        for(int i = 0; i < hits; i++) { hit_prob += log_theta; }
        for(int i = 0; i < draws-hits; i++) { nohit_prob += log_notheta; } 
        return log_bico(draws, hits) + hit_prob + nohit_prob;        
    }
    
    public static double log_binomialPValue(int draws, int hits, double theta) { 
        double sum = log_binomial(draws, hits, theta);
        for(int i = hits+1; i < draws; i++) { 
            sum = log_add(sum, log_binomial(draws, i, theta));
        }
        return sum;
    }
    
    public static double log_bico(int n, int k) { 
        return factln(n) - factln(k) - factln(n-k);
    }
    
    public static double bico(int n, int k) { 
        return Math.floor(0.5 + Math.exp(factln(n) - factln(k) - factln(n-k)));
    }

    // pg. 641 (Numerical Recipes in C).
    public static SpearmanResult spear(double[] data1, double[] data2) { 
        int n = data1.length;
        if(data2.length != n) { 
            throw new IllegalArgumentException("arrays must be of same length"); 
        }
        
        double d, zd, probd, rs, probrs;
        int j;
        double vard, t, sg, sf, fac, en3n, en, df, aved; 
        Double[] wksp1, wksp2;
        SpearmanResult sr = null;

        wksp1 = new Double[n];
        wksp2 = new Double[n];
        for(j = 1; j <= n; j++) { 
            wksp1[j-1] = data1[j-1];
            wksp2[j-1] = data2[j-1];
        }
        Sorter<Double,Double> sorter = new Sorter<Double,Double>();
        sorter.sort2(wksp1, wksp2);
        sf = crank(wksp1);
        sorter.sort2(wksp2, wksp1);
        sg = crank(wksp2);
        d = 0.0;
        
        for(j = 1; j <= n; j++) { 
            double temp = wksp1[j-1] - wksp2[j-1];
            d += temp * temp;
        }
        
        en = (double)n;
        en3n = en * en * en - en;
        aved = en3n / 6.0 - (sf + sg) / 12.0;
        fac = (1.0 - sf / en3n) * (1.0-sg / en3n);
        
        double temp = en + 1.0;
        temp *= temp;
        vard = ((en-1.0)*en*en*temp / 36.0) * fac;
        zd = (d - aved) / Math.sqrt(vard);
        probd = erfcc(Math.abs(zd)/1.4142136);
        rs = (1.0 - (6.0/en3n)*(d + (sf + sg)/12.0))/Math.sqrt(fac);
        fac = (rs + 1.0) * (1.0-rs);
        if(fac > 0.0) { 
            t = rs * Math.sqrt((en-2.0)/fac);
            df = en - 2.0;
            probrs = betai(0.5 * df, 0.5, df/(df + t * t));
        } else { 
            probrs = 0.0;
        }
        
        sr = new SpearmanResult(d, zd, probd, rs, probrs);
        return sr;
    }
    
    // pg. 642 (Numerical Recipes in C)
    public static double crank(Double[] w) { 
        int n = w.length;
        int j = 1, ji, jt;
        double t, rank;
        double s = 0.0;
        
        while(j < n) { 
            if(w[j] != w[j-1]) {
                w[j-1] = (double)j;
                ++j;
            } else { 
                for(jt = j+1; jt <= n && w[jt-1] == w[j-1]; jt++)
                    ;
                rank = 0.5 * (double)(j + jt - 1);
                for(ji = j; ji <= (jt-1); ji++) { 
                    w[ji-1] = rank;
                }
                t = (double)(jt - j);
                s += (t * t * t - t);
                j = jt;
            }
        }
        if(j == n) { w[n-1] = (double)n; }
        return s;
    }
    
    // pg. 227 (Numerical Recipes in C)
    public static double betacf(double a, double b, double x) { 
        int m, m2;
        double aa, c, d, del, h, qab, qam, qap;
        
        qab = a + b;
        qap = a + 1.0;
        qam = a - 1.0;
        c = 1.0;
        d = 1.0 - qab * x / qap;
        if(Math.abs(d) < FPMIN) { d = FPMIN; }
        d = 1.0 / d;
        h = d;
        for(m = 1; m <= MAXIT; m++) { 
            m2 = 2 * m;
            double dm = (double)m, dm2 = (double)m2;
            aa = dm * (b-dm) * x /((qam + dm2)*(a + dm2));
            d = 1.0 + aa * d;
            if(Math.abs(d) < FPMIN) { d = FPMIN; }
            c = 1.0 + aa / c;
            if(Math.abs(c) < FPMIN) { c = FPMIN; }
            d = 1.0 / d;
            h *= d * c;
            aa = -(a+dm) * (qab + dm)*x / ((a+dm2)*(qap+dm2));
            d = 1.0 + aa * d;
            if(Math.abs(d) < FPMIN) { d = FPMIN; }
            c = 1.0 + aa / c;
            if(Math.abs(c) < FPMIN) { c = FPMIN; }
            d = 1.0/d;
            del = d * c;
            h *= del;
            if(Math.abs(del-1.0) < EPS) { 
                break;
            }
            
        }
        if(m > MAXIT) { throw new NumericalException("a or b too big or MAXIT too small."); }
        return h;
    }

    // pg. 227 (Numerical Recipes in C)
    public static double betai(double a, double b, double x) {
        double bt;
        if(x < 0.0 || x > 1.0) { throw new IllegalArgumentException(); }
        if(x == 0.0 || x == 1.0) { 
            bt = 0.0; 
        } else { 
            bt = Math.exp(gammln(a+b) - gammln(a) - gammln(b) + 
                    a * Math.log(x) + b * Math.log(1.0-x));
        }
        
        if(x < (a + 1.0) / (a + b + 2.0)) { 
            return bt * betacf(a, b, x) / a;
        } else { 
            return 1.0-bt * betacf(b, a, 1.0-x) / b;
        }
    }

    public static class SpearmanResult { 
        private double d, zd, probd, rs, probrs;
        
        public SpearmanResult(double _d, double _zd, double _probd,
                double _rs, double _probrs) { 
            d = _d;
            zd = _zd;
            probd = _probd;
            rs = _rs;
            probrs = _probrs;
        }
        
        public double getD() { return d; }
        public double getZD() { return zd; }
        public double getProbD() { return probd; }
        public double getRS() { return rs; }
        public double getProbRS() { return probrs; }
    }

    public static class Sorter<X extends Comparable, Y> { 
        
        public Sorter() {
        }
        
        public void sort2(X[] arr, Y[] brr) {
            if(arr.length != brr.length) { throw new IllegalArgumentException(); }
            Sortable<X, Y>[] array = new Sortable[arr.length];
            for(int i = 0; i < arr.length; i++) { 
                array[i] = new Sortable<X, Y>(arr[i], brr[i]);
            }
            Arrays.sort(array);
            for(int i = 0; i < array.length; i++) { 
                arr[i] = array[i].getValue();
                brr[i] = array[i].getData();
            }
        }
    }
    
    public static class Sortable<X extends Comparable,Y> implements Comparable<Sortable<X,Y>> {
        
        private X value;
        private Y data;
        
        public Sortable(X v, Y d) { 
            value = v;
            data = d;
        }
        
        public X getValue() { return value; }
        public Y getData() { return data; }
        
        public int hashCode() { 
            int code = 17;
            code += data.hashCode(); code *= 37;
            code += data.hashCode(); code *= 37;
            return code;
        }
        
        public boolean equals(Object o) { 
            if(!(o instanceof Sortable)) { return false; }
            Sortable s = (Sortable)o;
            if(!value.equals(s.value)) { return false; }
            if(!data.equals(s.data)) { return false; }
            return true;
        }

        /* (non-Javadoc)
         * @see java.lang.Comparable#compareTo(java.lang.Object)
         */
        public int compareTo(Sortable<X,Y> s) {
            return value.compareTo(s.value);
        } 
    }
    
    /**
     * Computes the log-space PDF of the Poisson distribution with the
     * specified mean. Code lifted from the COLT Poisson class so that it can
     * be implemented as a static method instead of having to instantiate a
     * COLT Poisson object
     * @param mean
     * @param k
     * @return
     */
    public static double poissonLogPDF(double mean, int k) {
      return k * Math.log(mean) - Arithmetic.logFactorial(k) - mean;
    }
    
	/**
	 * LambertW function
	 * default value for <tt>maxIters = 100</tt>
	 * @see edu.mit.csail.cgs.utils.numeric.Numerical#lambertW(int, double, double, int)
	 */
	public static double lambertW(int branch_index, double z, double thres) {
		return lambertW(branch_index, z, thres, 100);
	}//end of lambertW method
	
	/**
	 * LambertW function
	 * default value for <tt>thres = 1e-12</tt>
	 * @see edu.mit.csail.cgs.utils.numeric.Numerical#lambertW(int, double, double, int)
	 */
	public static double lambertW(int branch_index, double z, int maxIters) {
		return lambertW(branch_index, z, 1e-12, maxIters);
	}//end of lambertW method

	/**
	 * LambertW function
	 * default values for <tt>thres = 1e-12, maxIters = 100</tt>
	 * @see edu.mit.csail.cgs.utils.numeric.Numerical#lambertW(int, double, double, int)
	 */
	public static double lambertW(int branch_index, double z) {
		return lambertW(branch_index, z, 1e-12, 100);
	}//end of lambertW method
	
	/**
	 * This method finds the solution of the equation: W * exp(W) = z 			<br>
	 * for only the W_{0} and W_{-1} branches.
	 * This function is known as the Lambert W function.							<br>
	 * For more information see:												<br>
	 * <a href="http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.47.767">
	 * Structure Learning in Conditional Probability Models via an Entropic Prior and Parameter Extinction
	 * </a>
	 * <br>
	 * <a href="http://www.cs.uwaterloo.ca/research/tr/1993/03/W.pdf">
	 * On the Lambert W Function
	 * </a> 
	 * <br>
	 * <a href="http://en.wikipedia.org/wiki/Lambert_function">
	 * Lambert W function
	 * </a> 
	 * @param branch_index Valid values 0 or -1.								<br>
	 * 0 for principal (W_{0}), -1 for W_{-1} branch 
	 * @param z The value whose we wish to determine the W s.t. W * exp(W) = z	<br>
	 * Valid values for z: z > -exp(-1)
	 * @param thres precision threshold
	 * @param maxIters maximum number of iterations for searching for W
	 * @return
	 */
	public static double lambertW(int branch_index, double z, double thres, int maxIters) {
		
		if(branch_index != -1 && branch_index != 0) {
			throw new IllegalArgumentException("The only valid values for branch_index is -1 (W_{-1}) and 0 (W_{0}).");
		}
		
		double W;
		int numIters;
		
		if( z < -Math.exp(-1) ) {
		    System.out.printf("The real branches of Lambert W function are not defined for z < -exp(-1) = %.4f.\n", -Math.exp(-1));
		    W = Double.NaN;
		    return W;
		}

		numIters = 0;
		W = initW(branch_index, z);

		if(Double.isNaN(W)) { return W; }

		double delta;
		while(numIters < maxIters) {
		    delta = W*Math.exp(W) - z;
		    
		    if(delta == 0) { break; }
		    if(Math.abs(-delta/(Math.exp(W)*(W+1))) < thres) { break; } 
		   
		    W -= delta/(Math.exp(W)*(W+1) - (delta*(W+2)/(2*(W+1))));
		    numIters++;
		}
		
		return W;

	}//end of lambertW method

	
	private static double initW(int branch_index, double z) {
		
		double W_init; 
		
		// Positive z
		if(z >= 0) {
			
			// branch W_{0}
			if(branch_index == 0) {
				if( z <= 500.0 ) { W_init = 0.665*(1 + 0.0195*Math.log(z + 1.0))*Math.log(z + 1.0) + 0.04; }
				else 			 {W_init  = Math.log(z - 4.0) - (1.0 - 1.0/Math.log(z))*Math.log(Math.log(z)); }
			}
			
			// branch W_{-1} - It is NOT defined for z >= 0
			else { W_init = Double.NaN; }
			
			
		}
		
		// Negative z
		else {
			
			// z not too near -Math.exp(-1)
			if(Math.abs(z + Math.exp(-1)) >  0.01) {
				
				// branch W_{0}
				if(branch_index == 0)     { W_init = 0; }
				
				// branch W_{-1}
				else                      { W_init = Math.log(-z) - Math.log(-Math.log(-z)); }
			}
			
			// z too near -Math.exp(-1)
			else {
				if( z == -Math.exp(-1))   { W_init = -1; }
				
				else {
					// branch W_{0}
					if(branch_index == 0) { W_init = -1 + Math.sqrt(2*(Math.exp(1)*z +1)); }
					
					// branch W_{-1}
					else 				  { W_init = -1 - Math.sqrt(2*(Math.exp(1)*z +1)); }
				}
			
			}
		}
		
		return W_init;
	}//end initW method

    
}//end of Numerical class

/*
 * Created on May 24, 2006
 */
package edu.mit.csail.cgs.utils;

import java.text.*;
import java.io.*;
import java.util.*;

/**
 * @author tdanford
 */
public class Enrichment implements Comparable<Enrichment> {
    
    private static NumberFormat nf;
    private static DecimalFormat sci;
    
    static { 
        nf = DecimalFormat.getInstance();
        nf.setMaximumFractionDigits(6);
        nf.setMinimumFractionDigits(6);
        
        sci = new DecimalFormat("0.0000E0");
    }
    
    private String category;
    private int N, theta, n, x;
    private double log_pvalue;
    
    public Enrichment(String c, int _N, int th, int _n, int _x, double pv) { 
        category = c;
        N = _N;
        theta = th;
        n = _n; 
        x = _x; 
        log_pvalue = pv;
    }
    
    public boolean passedThreshold(double thresh) { return log_pvalue <= thresh; }
    public String getCategory() { return category; }
    public int getN() { return N; }
    public int getTheta() { return theta; }
    public int getn() { return n; }
    public int getx() { return x; }
    public double getLogPValue() { return log_pvalue; }
    
    public double getBGFrac() { return N > 0 ? (double)theta / (double)N : 0.0; }
    public double getFGFrac() { return n > 0 ? (double)x / (double)n : 0.0; }
    
    public String toString() { 
        return String.format("%.5f\t%s (%d %d)", Math.exp(log_pvalue), category, x, n);
    }
    public String reportLine() {
        return String.format("%s, pval= %.2e, %d / %d = %.3f (freq in GO %d / %d = %.3f)",
                             getCategory(),
                             Math.exp(getLogPValue()),
                             getx(),
                             getn(),
                             (double)((double)getx())/getn(),
                             getTheta(),
                             getN(),
                             (double)((double)getTheta())/getN());
    }
    
    public int compareTo(Enrichment e) { 
        if(log_pvalue < e.log_pvalue) { return -1; }
        if(log_pvalue > e.log_pvalue) { return 1; }
        if(x > e.x) { return -1; }
        if(x < e.x) { return 1; }
        return category.compareTo(e.category);
    }
    
    public int hashCode() { 
        int code = 17;
        code += category.hashCode(); code *= 37;
        return code;
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof Enrichment)) { return false; }
        Enrichment e = (Enrichment)o;
        if(!category.equals(e.category)) { return false; }
        //if(N != e.N || n != e.n || x != e.x || theta != e.theta) { return false; }
        return true;
    }
    
}

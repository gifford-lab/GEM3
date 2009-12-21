/*
 * Created on Oct 9, 2005
 */
package edu.mit.csail.cgs.utils;

/**
 * @author tdanford
 */
public class BinCounter {
    
    private int[] bins;
    private int totalCount;

    public BinCounter(int numBins) {
        bins = new int[numBins];
        for(int i = 0; i < bins.length; i++) { 
            bins[i] = 0;
        }
        totalCount = 0;
    }

    public int getNumBins() { return bins.length; }
    public int getTotalCount() { return totalCount; }
    
    public int getMaxBinCount(){
    	int max=0;
    	for(int x=0; x<bins.length; x++){
    		if(bins[x]>max){max=bins[x];}
    	}return(max);
    }
    public int getMinBinCount(){
    	int min=Integer.MAX_VALUE;
    	for(int x=0; x<bins.length; x++){
    		if(bins[x]<min){min=bins[x];}
    	}return(min);
    }
    
    public void addToBin(int b, int c) { 
        bins[b] += c;
        totalCount += c;
    }
    
    public int getBin(int b) { return bins[b]; }
    
    public double getBinFraction(int b) {
        if(totalCount == 0) { return 0.0; }
        return (double)bins[b] / (double)totalCount;
    }
    
    public double[] cumulativeFromLeft() { 
        double[] cdf = new double[bins.length];
        int seen = 0;
        for(int i = 0; i < cdf.length; i++) { 
            seen += bins[i];
            cdf[i] = 0.0;
            if(totalCount > 0) { 
                cdf[i] = (double)seen / (double)totalCount;
            }
        }
        return cdf;
    }

    public double[] cumulativeFromRight() { 
        double[] cdf = new double[bins.length];
        int seen = 0;
        for(int i = cdf.length-1; i >= 0; i--) { 
            seen += bins[i];
            cdf[cdf.length-1-i] = 0.0;
            if(totalCount > 0) { 
                cdf[cdf.length-1-i] = (double)seen / (double)totalCount;
            }
        }
        return cdf;
    }
}

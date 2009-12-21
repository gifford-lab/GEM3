package edu.mit.csail.cgs.utils;

import java.util.Arrays;
import java.util.List;

import edu.mit.csail.cgs.utils.Utils;

/**
 * Identical to MATLAB's <tt>hist</tt> and <tt>histc</tt> functions.			<br>
 * It counts the number of values in vector x that fall between the elements 
 * in the edges vector (which must contain monotonically nondecreasing values).	<br>
 * 																				<br>
 * The edges form essentially the boundaries of the bins. So, if we provide
 *  <tt>binEdges.length</tt> edges, then we will have <tt>binEdges.length-1</tt> bins. <br>
 *  Bins are inclusive from the left and exclusive from the right: [ , )			<br>
 *  Lastly, if values fall outside the range that the bins represent, they are
 *  assigned to the most extreme bins.											<br>
 *  In more detail:																<br>
 *  boundaries of bins[0]:    [binEdges[0], binEdges[1])							<br>
 *  boundaries of bins[1]:    [binEdges[1], binEdges[2])							<br>
 *  ...																			<br>
 *  boundaries of bins[bin.length-1]:    [binEdges[binEdges.length-2], binEdges[binEdges.length-1])			<br>
 *  																			<br>
 *  We can also provide directly bin centers. In that case, the number of bins
 *  will equal the number of bin centers.										<br>
 *  The boundaries of each bin will be:											<br>
 *  boundaries of bins[0]:    [binCenters[0] - (binCenters[1]-binCenters[0])/2.0, (binCenters[0]+binCenters[1])/2.0)							<br>
 *  boundaries of bins[1]:    [(binCenters[0]+binCenters[1])/2.0, (binCenters[1]+binCenters[2])/2.0)							<br>
 *  ...																			<br>
 *  boundaries of bins[bins.length-1]:    [(binCenters[binCenters.length-2]+binCenters[binCenters.length-1])/2.0, binCenters[binCenters.length-1] + (binCenters[binCenters.length-1]-binCenters[binCenters.length-2])/2.0)							<br>
 *  
 *  																			<br> 
 *  In addition, we can provide the number of bins (<tt>numBins</tt>), where in this case the bin boundaries
 *  will be automatically determined from the data.								<br>
 *  In this case, we create <tt>numBins</tt> equally sized bins each having width
 *  <tt>binwidth = (max-min)/range</tt>.										<br>
 *  In more detail:																<br>
 *  boundaries of bins[0]:    [min, min + binwidth)								<br>
 *  boundaries of bins[1]:    [min + binwidth, min + 2*binwidth)  				<br>
 *  ...																			<br>
 *   boundaries of bins[bins.length-1]:    [min + (bins.length-1)*binwidth, max)	<br>
 *   																			<br>
 *  Lastly, we can provide the width of the bins (<tt>binwidth</tt>), where in this case 
 *  the number of bins will be automatically determined from the data.			<br>
 *  In this case, we create <tt>numBins = Math.ceil((upper_bound-lower_bound)/binwidth)</tt> 
 *  equally sized bins each having width <tt>binwidth</tt>.						<br> 
 *  In more detail:																<br>
 *  boundaries of bins[0]:    [lower_bound, lower_bound + binwidth)								<br>
 *  boundaries of bins[1]:    [lower_bound + binwidth, lower_bound + 2*binwidth)  				<br>
 *  ...																			<br>
 *   boundaries of bins[bins.length-1]:    [lower_bound + (bins.length-1)*binwidth, upper_bound)

 * @author yguo
 */
public class HistogramCount{
    
    private double[] binEdges;
    private double[] binCenters;
    private int[] bins;
    private int totalCount;
    
    public HistogramCount(int numBins, List<Double> values) {
    	this(numBins, Utils.ref2prim(values.toArray(new Double[0])));
    }//end of HistogramCount constructor
    
    /**
     * Creates <tt>numBins</tt> equally sized bins starting from the minimum value
     * of <tt>values</tt> to the maximum value of <tt>values</tt>. 
     * @param numBins number of bins to be created
     * @param values values to be binned
     */
    public HistogramCount(int numBins, double[] values) {
    	if(numBins < 0) { throw new IllegalArgumentException("numBins has to be a positive integer."); }
    	
    	double[] values_temp = values.clone();
    	Arrays.sort(values_temp);
    	double min   = values_temp[0];
    	double max   = values_temp[values_temp.length-1];
    	double range = max-min;
    	double binwidth = range/numBins;
    	
    	binEdges = new double[numBins+1];
    	int i = 0;
    	double bound = min;
    	while(bound <= max) {
    		binEdges[i++] = bound;
    		bound += binwidth;
    	}
		binCenters = new double[binEdges.length-1];
		for(i = 0; i < binCenters.length; i++) { binCenters[i] = (binEdges[i]+binEdges[i+1])/2.0; }

        bins = new int[binEdges.length-1];
        totalCount = 0;

        addAllValues(values_temp);
    }//end of HistogramCount constructor
    
    
    /**
     * Creates <tt>numBins</tt> equally sized bins starting from the minimum value
     * of <tt>values</tt> to the maximum value of <tt>values</tt>. 
     * @param numBins number of bins to be created
     * @param values values to be binned
     */
    public HistogramCount(int numBins, int[] values) {
    	if(numBins < 0) { throw new IllegalArgumentException("numBins has to be a positive integer."); }
    	
    	int[] values_temp = values.clone();
    	Arrays.sort(values_temp);
    	int min   = values_temp[0];
    	int max   = values_temp[values_temp.length-1];
    	int range = max-min;
    	double binwidth = (double)range/numBins;
    	
    	binEdges = new double[numBins+1];
    	int i = 0;
    	double bound = min;
    	while(bound <= max) {
    		binEdges[i++] = bound;
    		bound += binwidth;
    	}
		binCenters = new double[binEdges.length-1];
		for(i = 0; i < binCenters.length; i++) { binCenters[i] = (binEdges[i]+binEdges[i+1])/2.0; }

        bins = new int[binEdges.length-1];
        totalCount = 0;

        addAllValues(values_temp);
    }//end of HistogramCount constructor
    
    public HistogramCount(double lower_bound, double upper_bound, int[] values) {
    	this(lower_bound, upper_bound);
    	addAllValues(values);
    }//end of HistogramCount constructor
    
    public HistogramCount(double lower_bound, double upper_bound, double[] values) {
    	this(lower_bound, upper_bound);
    	addAllValues(values);
    }//end of HistogramCount constructor
    
    public HistogramCount(double lower_bound, double upper_bound, List<Double> values) {
    	this(lower_bound, upper_bound);
    	addAllValues(values);
    }//end of HistogramCount constructor
   
    public HistogramCount(double lower_bound, double upper_bound) {
    	this(lower_bound, upper_bound, 1.0);
    }//end of HistogramCount constructor
    
    public HistogramCount(double lower_bound, double upper_bound, double binwidth, List<Double> values) {
    	this(lower_bound, upper_bound, binwidth);
    	addAllValues(values);
    }//end of HistogramCount constructor    
    
    public HistogramCount(double lower_bound, double upper_bound, double binwidth, double[] values) {
    	this(lower_bound, upper_bound, binwidth);
    	addAllValues(values);
    }//end of HistogramCount constructor

    /**
     * Default binwidth = 1
     */
    public HistogramCount(double lower_bound, double upper_bound, double binwidth, int[] values) {
    	this(lower_bound, upper_bound, binwidth);
    	addAllValues(values);
    }//end of HistogramCount constructor
    
    /**
     * This method forms a histogram, with each bin having width: <tt>binwidth</tt>
     * @param lower_bound the lower bound of the bins
     * @param upper_bound the upper bound of the bins
     * @param binwidth the bin width
     */
    public HistogramCount(double lower_bound, double upper_bound, double binwidth) {
    	if((upper_bound - lower_bound) < 0) { throw new IllegalArgumentException("upper_bound has to be greater than lower_bound."); }
    	
    	double range = upper_bound - lower_bound;
    	int numBins  = (int) Math.ceil(range/binwidth);
    	
    	binEdges = new double[numBins+1];
    	int i = 0;
    	double bound = lower_bound;
    	while(bound < upper_bound) {
    		binEdges[i++] = bound;
    		bound += binwidth; 		
    	}
    	binEdges[i] = upper_bound;
    	
		binCenters = new double[binEdges.length-1];
		for(i = 0; i < binCenters.length; i++) { binCenters[i] = (binEdges[i]+binEdges[i+1])/2.0; }
 	
        bins = new int[binEdges.length-1];
        totalCount = 0;    	
    }//end of HistogramCount constructor
        
    /** Except for just constructing the histogram, it adds values in it. */
    public HistogramCount(double[] binSpecs, char flag, int[] values) {
    	this(binSpecs, flag);
    	addAllValues(values);
    }//end of HistogramCount constructor
    
    /** Except for just constructing the histogram, it adds values in it. */
    public HistogramCount(double[] binSpecs, char flag, List<Double> values) {
    	this(binSpecs, flag);
    	addAllValues(values);
    }//end of HistogramCount constructor
    
    /** Except for just constructing the histogram, it adds values in it. */
    public HistogramCount(double[] binSpecs, char flag, double[] values) {
    	this(binSpecs, flag);
    	addAllValues(values);
    }//end of HistogramCount constructor
    
    /**
     * If the user wants to catch all values, use Double.MIN_VALUE and Double.MAX_VALUE as edge
     * @param binSpecs Either bin centers or bin edges (which are monotonically nondecreasing values)
     * @param flag 'c' means that <tt>binSpecs</tt> are bin centers, while 'e' that 
     * <tt>binSpecs</tt> are bin edges.											<br> <br>
     * Note: Either bin centers or bin edges are assumed to be sorted.
     */
    public HistogramCount(double[] binSpecs, char flag) {
    	if(flag == 'c') {
    		binCenters = binSpecs;
        	binEdges = new double[binCenters.length+1];
        	for(int i = 1; i < binCenters.length; i++) { binEdges[i] = (binCenters[i-1] + binCenters[i])/2.0; }
        	binEdges[0] = binCenters[0] - (binCenters[1] - binCenters[0])/2.0;
        	binEdges[binEdges.length-1] = binCenters[binCenters.length-1] + (binCenters[binCenters.length-1] - binCenters[binCenters.length-2])/2.0;
    	}
    	else if(flag == 'e') {
    		binEdges = binSpecs;
    		binCenters = new double[binEdges.length-1];
    		for(int i = 0; i < binCenters.length; i++) { binCenters[i] = (binEdges[i]+binEdges[i+1])/2.0; }
    	}
    	else {throw new IllegalArgumentException("flag can take either the values: 'c' representing bin centers or 'e' representing bin edges.\n" +
    											 "Do not include the single quotes in input."); }
    	
        bins = new int[binEdges.length-1];
        totalCount = 0;
    }//end of HistogramCount constructor

   
	/***********************************
	 **   GET HISTOGRAM STATISTICS    **
	 ***********************************/ 
    
    public int      getNumBins()     { return bins.length; }
    public int      getTotalCount()  { return totalCount;  }
    public double[] getBinEdges()    { return binEdges;    }
    public double[] getBinCenters()  { return binCenters;  }
    public int[]    getBins()        { return bins;        }
    /**
     * Returns the counts of the bin with this index
     * @param b bin index
     * @return
     */
    public int getBin(int b)      { return bins[b]; }


	public int getMaxBinCount(){
    	int max=0;
    	for(int x=0; x<bins.length; x++){ if(bins[x]>max){max=bins[x];} }
    	return max;
    }
	
    public int getMinBinCount(){
    	int min=Integer.MAX_VALUE;
    	for(int x=0; x<bins.length; x++){ if(bins[x]<min){min=bins[x];} }
    	return min;
    }
        
    public double getBinFraction(int b) {
        if(totalCount == 0) { return 0.0; }
        return (double)bins[b]/totalCount;
    }
    
    public double[] cumulativeFromLeft() { 
        double[] cdf = new double[bins.length];
        double seen = 0;
        for(int i = 0; i < cdf.length; i++) { 
            seen += bins[i];
            cdf[i] = 0.0;
            if(totalCount > 0) { 
                cdf[i] = seen / totalCount;
            }
        }
        return cdf;
    }

    public double[] cumulativeFromRight() { 
        double[] cdf = new double[bins.length];
        double seen = 0;
        for(int i = cdf.length-1; i >= 0; i--) { 
            seen += bins[i];
            cdf[cdf.length-1-i] = 0.0;
            if(totalCount > 0) { 
                cdf[cdf.length-1-i] = seen / totalCount;
            }
        }
        return cdf;
    }
    
    
	/***********************************
	 **    PERFORM BIN OPERATIONS     **
	 ***********************************/ 
   
    /**
     * If a value falls on a bin edge, it goes to the right bin
     * i.e. left inclusive
     * @param v value to be binned
     * @return
     */
    public int findBin(double v){
    	if(v < binEdges[0]) { return 0; }
    	else if (v >= binEdges[binEdges.length-1]) { return bins.length-1; }
    	else { 
    		int index = Arrays.binarySearch(binEdges, v); 
    		if (index < 0) { return -index-2; }
    		else { return index; }
    	}
    }
    
    /**
     * @param b bin index
     * @param c counts to be added in the bin with this index
     */
    public void addToBin(int b, int c) { 
        bins[b] += c;
        totalCount += c;
    } 

    /**
     * Adds the specified value to the corresponding bin one time
     * @param v value to be added
     */
    public void addValue(double v) { 
    	addValue(v, 1);
    }
    
    /**
     * Adds the specified value to the corresponding bin <tt>count</tt> times
     * @param v
     * @param count
     */
    public void addValue(double v, int count) { 
    	addToBin(findBin(v), count);
    } 
    
    public void addAllValues(double[] values) { 
    	for (double v:values){ addValue(v, 1); }
    }
    
    public void addAllValues(List<Double> values) { 
    	for (double v:values){ addValue(v, 1); }
    }
    
    public void addAllValues(int[] values) { 
    	for (int v:values){ addValue(v, 1); }
    }
    
    public void printContents(){
    	for(int i = 0; i < bins.length; i++) { 
            System.out.println(String.format("%d\t%f", i, bins[i]));
        }
    }
}//end of HistogramCount class

package edu.mit.csail.cgs.utils.stats;

import java.lang.reflect.Array;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;

import cern.jet.random.Beta;
import cern.jet.random.Binomial;
import cern.jet.random.Gamma;
import cern.jet.random.Uniform;
import cern.jet.random.engine.DRand;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KMAC;
import edu.mit.csail.cgs.deepseq.discovery.kmer.Kmer;
import edu.mit.csail.cgs.deepseq.discovery.kmer.mtree.MTree;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.probability.NormalDistribution;

public class StatUtil {
//	static cern.jet.random.engine.RandomEngine engine = new cern.jet.random.engine.MersenneTwister();
	static private double[] logFactorials;
	static private ConcurrentHashMap<String, Double> HGP_table = new ConcurrentHashMap<String, Double>();
	static private ConcurrentHashMap<String, Double> log10_HGP_table = new ConcurrentHashMap<String, Double>();
	static private ConcurrentHashMap<String, Double> HGPDF_table = new ConcurrentHashMap<String, Double>();
	static private ConcurrentHashMap<String, Double> log10_HGPDF_table = new ConcurrentHashMap<String, Double>();
	static public String getCacheSize(){
		return String.format("HGP_table=%s, log10_HGP_table=%s, HGPDF_table=%s, log10_HGPDF_table=%s", 
				HGP_table.size(), log10_HGP_table.size(), HGPDF_table.size(), log10_HGPDF_table.size());
	}
	static public int cacheAccessCount = 0;
	
	public static double uniform_rnd() {
		return Uniform.staticNextDouble();
	}	
	
	public static double beta_rnd(double a, double b) {
		return Beta.staticNextDouble(a,b);
	}
		
	public static void dirichlet_rnd(double[] x, double[] alpha, int vSize) {
		double sum = 0.0;
		double xt = 0.0;
		int i = 0;
		
		for(i=0;i<vSize;i++) {
			if (alpha[i] <= 0.0) {
				xt = 0.0;
			} else {
				xt = Gamma.staticNextDouble(alpha[i],1.0);
			}
			if (Double.isNaN(xt)) {
				xt = 0.0;
			}
			x[i] = xt;
			sum = sum + x[i];
		}
		
		for(i=0;i<vSize;i++) {
			x[i] = x[i]/sum;
		}
	}
	
	
	public static double logDirichletPDF(double[] x, double[] alpha, int vSize) {
		double alpha_sum = 0.0;
		int i = 0;
		double lik = 0.0;
		
		for (i=0;i<vSize;i++) {
			if (alpha[i] > 0.0) {
				if (x[i] > 0.0) {
					lik = lik + (alpha[i] - 1.0)*Math.log(x[i]);
				}
				lik = lik - cern.jet.stat.Gamma.logGamma(alpha[i]);
				alpha_sum = alpha_sum + alpha[i];
			}
		}
		
		lik = lik+cern.jet.stat.Gamma.logGamma(alpha_sum);
		
		return lik;
	}
	
	
	public static int multinomial_rnd(double[] p, int vSize) {
		double sum = 0.0;
		double mass = 0.0;
		int cc = 0;
		int i = 0;
		
		for (i=0;i<vSize;i++) {
			sum += p[i];
		}
		
		mass = Uniform.staticNextDouble() * sum;
		
		i = 0;
		while (true) {
			mass -= p[i];
			if (mass <= 0.0) break;
			i++;
			cc++;
		}
		return cc;
	}
	/**
	 * Round to the nearest integer.<br>
	 * This is different from Math.round() for negative number.
	 * i.e., 2.5->3, 2.499->2, -1.499->-1, -1.5-->-2
	 */
	public static int round(double n){
		return (int)((n>0?0.5:-0.5)+n);
	}
	
	/**
	 * Returns the median of this array
	 * @param x array of type <tt>double</tt>
	 * @return
	 */
	public static double median(double[] x) {
		if (x.length == 1) {
			return x[0];
		}
		
		double[] y = new double[x.length];
		System.arraycopy(x,0,y,0,x.length);
		Arrays.sort(y);
		int n = x.length;
		double m = ((double) n)/2.0;
		if (Math.floor(m) == m) {
			n = n / 2;
			return (y[n] + y[n-1])/2.0;
		} else {
			n = (n-1)/2;
			return y[n];
		}
	}
	
	/**
	 * Returns the median of this array
	 * @param x array of type <tt>double</tt>
	 * @return
	 */
	public static float median(float[] x) {
		if (x.length == 1) {
			return x[0];
		}
		
		float[] y = new float[x.length];
		System.arraycopy(x,0,y,0,x.length);
		Arrays.sort(y);
		int n = x.length;
		float m = ((float) n)/2.0f;
		if (Math.floor(m) == m) {
			n = n / 2;
			return (y[n] + y[n-1])/2.0f;
		} else {
			n = (n-1)/2;
			return y[n];
		}
	}
	
	/** 
	 * Count occurrences of the elements
	 */
	public static TreeMap<Integer, Integer> countOccurences (ArrayList<Integer> nums){
		TreeMap<Integer, Integer> map = new TreeMap<Integer, Integer>();
		for (int i:nums){
			if (!map.containsKey(i)){
				map.put(i, 1);
			}
			else{
				map.put(i, map.get(i)+1);
			}
		}
		return map;
	}	
	/** 
	 * Sort the elements by their occurrences, ascending order
	 * @return Pair of [elements, counts]
	 */
	public static Pair<int[], int[]> sortByOccurences (ArrayList<Integer> integerList){
		TreeMap<Integer, Integer> map = countOccurences(integerList);
		int[] elements = new int[map.keySet().size()];
		int[] counts = new int[elements.length];
		int i=0;
		for (int e:map.keySet()){
			elements[i]=e;
			counts[i]=map.get(e);
			i++;
		}
		int[] idx = findSort(counts);
		int[] sortedElements = new int[elements.length];
		for (int j=0;j<idx.length;j++)
			sortedElements[j] = elements[idx[j]];
		return new Pair<int[],int[]>(sortedElements, counts);
	}

	/**
	 * Sorts the array and returns the positions of the original array corresponding
	 * to the ordered elements <br>
	 * Accepts only primitive integers<br>
	 * The array is sorted in ascending order after the method is called
	 * @param a integer array to be sorted
	 * @return positions of the original array corresponding to the ordered elements
	 */
	public static int[] findSort(int[] a) {
		int[] sortedInds = new int[a.length];
		Map<Integer, ArrayList<Integer>> val2Index = new HashMap<Integer, ArrayList<Integer>>();
		
		for(int i = 0; i < a.length; i++)
		{
			if( !val2Index.containsKey(a[i])) 
				val2Index.put(a[i], new ArrayList<Integer>());
			
			val2Index.get(a[i]).add(i);
		}
		
		Arrays.sort(a);
	
		Set<Integer> uniqueEls = new LinkedHashSet<Integer>();
		for(int i = 0; i < a.length; i++) { uniqueEls.add(a[i]); }
		int count = 0;
		for(Object key:uniqueEls.toArray()){
			List<Integer> currVal_idxs = val2Index.get(key);
			for(Integer curr_idx:currVal_idxs)
				sortedInds[count++] = curr_idx;
		}
		return sortedInds;	
	}//end of findSort method 
	
	/**
	 * Sorts the double array and returns the positions of the original array corresponding
	 * to the ordered elements
	 * @param a double array to be sorted
	 * @return positions of the original array corresponding to the ordered elements
	 */
	public static int[] findSort(double[] a) {
		int[] sortedInds = new int[a.length];
		Map<Double, ArrayList<Integer>> val2Index = new HashMap<Double, ArrayList<Integer>>();
		
		for(int i = 0; i < a.length; i++)
		{
			if( !val2Index.containsKey(a[i])) 
				val2Index.put(a[i], new ArrayList<Integer>());
			val2Index.get(a[i]).add(i);
		}
		Arrays.sort(a);
	
		Set<Double> uniqueEls = new LinkedHashSet<Double>();
		for(int i = 0; i < a.length; i++) { uniqueEls.add(a[i]); }
		int count = 0;
		for(Object key:uniqueEls.toArray()){
			List<Integer> currVal_idxs = val2Index.get(key);
			for(Integer curr_idx:currVal_idxs)
				sortedInds[count++] = curr_idx;
		}
		return sortedInds;	
	}//end of findSort method - accepts only primitive double
		
	/**
	 * @see  edu.mit.csail.cgs.utils.stats.StatUtil#findSort(Object[], Comparator)
	 */
	public static <T> int[] findSort(T[] a) {
		return findSort(a, null);
	}
	
	
	/**
	 * Orders the elements of the array <tt>a</tt> in ascending order (based on
	 * the <tt>Comparator c</tt>) and returns the indexes of the sorted elements
	 * in the original array. <br>
	 * For reverse ordering set: <tt> c = java.util.Collections.reverseOrder() </tt>
	 * @param <T> type of array to be sorted
	 * @param a Array of type <tt>T</tt> to be sorted
	 * @param c Comparator
	 * @return the indexes of the sorted elements in the original array
	 */
	public static <T> int[] findSort(T[] a, Comparator<? super T> c) {
		int[] sortedInds = new int[a.length];
		
		Map<T, ArrayList<Integer>> val2Index = new HashMap<T, ArrayList<Integer>>();
		
		for(int i = 0; i < a.length; i++)
		{
			if( !val2Index.containsKey(a[i])) 
				val2Index.put(a[i], new ArrayList<Integer>());
			
			val2Index.get(a[i]).add(i);
		}
		
		Arrays.sort(a, c);
		
		Set<T> uniqueEls = new LinkedHashSet<T>(Arrays.asList(a));
		int count = 0;
		for(Object key:uniqueEls.toArray()){
			List<Integer> currVal_inds = val2Index.get(key);
			for(Integer currInd:currVal_inds)
				sortedInds[count++] = currInd;
		}
		
		return sortedInds;
	}//end of findSort method
	
	
	/**
	 * Returns the maximum of this array as well as all the positions in the array
	 * corresponding to the maximum
	 * @param x The array whose maximum is going to be determined
	 * </br>
	 * <tt>x</tt> should be a non-empty array
	 * @return The maximum of this array as well as the positions where these maximums are found (in ordered form)
	 */
	public static <T extends Comparable<T> >    Pair<T, TreeSet<Integer>> findMax(T[] x)
	{ 
		if( x.length == 0 )
			throw new IndexOutOfBoundsException("Hey dude, you cannot enter an empty array...");
		
		T maximum = x[0];
		TreeSet<Integer> max_index = new TreeSet<Integer>();   max_index.add(0);
		
		int count = 0;
		for(T el: x)
		{
			if(el.compareTo(maximum) > 0)
			{
				maximum = el;
				max_index.clear();
			    max_index.add(count);
			}
			else if(el.compareTo(maximum) == 0)
			{
				max_index.add(count);
			}
			count++;
		}
		
		return new Pair<T, TreeSet<Integer>>(maximum, max_index);
	}//end of findMax method 
	
	/**
	 * Returns the maximum of this double array as well as all the positions in the array
	 * corresponding to the maximum
	 * @param x The double array whose maximum is going to be determined
	 * </br>
	 * <tt>x</tt> should be a non-empty double array
	 * @return The maximum of this array as well as the positions where these maximums are found (in ordered form)
	 */
	public static Pair<Double, TreeSet<Integer>> findMax(double[] x)
	{ 
		if( x.length == 0 )
			throw new IndexOutOfBoundsException("Hey dude, you cannot enter an empty array...");
		
		Double maximum = x[0];
		TreeSet<Integer> max_index = new TreeSet<Integer>();   
		max_index.add(0);
		
		int count = 0;
		for(Double el: x)
		{
			if(el.compareTo(maximum) > 0)
			{
				maximum = el;
				max_index.clear();
			    max_index.add(count);
			}
			else if(el.compareTo(maximum) == 0)
			{
				max_index.add(count);
			}
			count++;
		}
		
		return new Pair<Double, TreeSet<Integer>>(maximum, max_index);
	}//end of findMax method 
	/**
	 * Returns the maximum of this double array as well as all the positions in the array
	 * corresponding to the maximum
	 * @param x The double array whose maximum is going to be determined
	 * </br>
	 * <tt>x</tt> should be a non-empty double array
	 * @return The maximum of this array as well as the positions where these maximums are found (in ordered form)
	 */
	public static Pair<Integer, TreeSet<Integer>> findMax(int[] x)
	{ 
		if( x.length == 0 )
			throw new IndexOutOfBoundsException("Hey dude, you cannot enter an empty array...");
		
		Integer maximum = x[0];
		TreeSet<Integer> max_index = new TreeSet<Integer>();   
		max_index.add(0);
		
		int count = 0;
		for(Integer el: x)
		{
			if(el.compareTo(maximum) > 0)
			{
				maximum = el;
				max_index.clear();
			    max_index.add(count);
			}
			else if(el.compareTo(maximum) == 0)
			{
				max_index.add(count);
			}
			count++;
		}
		
		return new Pair<Integer, TreeSet<Integer>>(maximum, max_index);
	}//end of findMax method 
	public static double getMax(double[] x){
		double max = Double.MIN_VALUE;
		for (int i=0;i<x.length;i++){
			if (max<x[i])
				max = x[i];
		}
		return max;
	}
	public static float getMax(float[] x){
		float max = Float.MIN_VALUE;
		for (int i=0;i<x.length;i++){
			if (max<x[i])
				max = x[i];
		}
		return max;
	}	
	public static int getMax(int[] x){
		int max = Integer.MIN_VALUE;
		for (int i=0;i<x.length;i++){
			if (max<x[i])
				max = x[i];
		}
		return max;
	}
	public static double getMin(double[] x){
		double min = Double.MAX_VALUE;
		for (int i=0;i<x.length;i++){
			if (min>x[i])
				min = x[i];
		}
		return min;
	}
	public static int getMin(int[] x){
		int min = Integer.MAX_VALUE;
		for (int i=0;i<x.length;i++){
			if (min>x[i])
				min = x[i];
		}
		return min;
	}
	/**
	 * Returns the minimum of this array as well as all the positions in the array
	 * corresponding to the minimum
	 * @param x The array whose minimum is going to be determined 
	 * </br>
	 * <tt>x</tt> should be a non-empty array
	 * @return The minimum of this array as well as the positions where these minimums are found (in ordered form)
	 */
	public static <T extends Comparable<T> >    Pair<T, TreeSet<Integer>> findMin(T[] x)
	{ 
		if( x.length == 0 )
			throw new IndexOutOfBoundsException("Hey dude, you cannot enter an empty array...");
		
		T minimum   = x[0];
		TreeSet<Integer> min_index = new TreeSet<Integer>();   min_index.add(0);
		
		int count = 0;
		for(T el: x)
		{
			if(el.compareTo(minimum) < 0)
			{
				minimum = el;
				min_index.clear();
				min_index.add(count);
			}
			else if(el.compareTo(minimum) == 0)
			{
				min_index.add(count);
			}
			count++;
		}
		
		return new Pair<T, TreeSet<Integer>>(minimum, min_index);
	}//end of findMin method 
	

	/**
	 * Returns the minimum of this double array as well as all the positions in the array
	 * corresponding to the minimum
	 * @param x The double array whose minimum is going to be determined 
	 * </br>
	 * <tt>x</tt> should be a non-empty double array
	 * @return The minimum of this array as well as the positions where these minimums are found (in ordered form)
	 */
	public static  Pair<Double, TreeSet<Integer>> findMin(double[] x)
	{ 
		if( x.length == 0 )
			throw new IndexOutOfBoundsException("Hey dude, you cannot enter an empty array...");
		
		double minimum   = x[0];
		TreeSet<Integer> min_index = new TreeSet<Integer>();   min_index.add(0);
		
		int count = 0;
		for(Double el: x)
		{
			if(el.compareTo(minimum) < 0)
			{
				minimum = el;
				min_index.clear();
				min_index.add(count);
			}
			else if(el.compareTo(minimum) == 0)
			{
				min_index.add(count);
			}
			count++;
		}
		
		return new Pair<Double, TreeSet<Integer>>(minimum, min_index);
	}//end of findMin method 
	
	/**
	 * Returns the position of the integer number <tt>f</tt> in the list. <br>
	 * If it is not there, <tt>-1</tt> is returned.
	 * @param x List of integers
	 * @param f number to be sought
	 * @return
	 */
	public static int find(ArrayList<Integer> x, int f) {
		int found = -1;
		int i = 0;
		
		while((i < x.size()) & (found == -1)) {
			if (x.get(i) == f) {
				found = i;
			}
			i++;
		}
		
		return found;
	}


	/**
	 * Returns all the positions (as a list) of the array that hold the satisfying
	 * relation with the element <tt>el</tt>.
	 * @param <T> the type of the elements of <tt>x</tt> array
	 * @param x The array whose positions satisfying the operation are going to be found
	 * @param operator the operator. Valid values: <tt>{ >, <, >=, <=, = }</tt>
	 * @param el the element which will be compared to the elements of <tt>x</tt> array
	 * according to the operation dictated by the operator
	 * @return The positions which satisfy the operation
	 */
	public static <T extends Comparable<T> >    ArrayList<Integer> find(T[] x, String operator, T el)
	{ 
		if( x.length == 0 )
			throw new IndexOutOfBoundsException("Hey dude, you cannot enter an empty array...");
		
		char operChar = 'n';
		if(operator.equals(">")) operChar = '>';
		else if(operator.equals("<")) operChar = '<';
		else if(operator.equals("<=")) operChar = 'L';
		else if(operator.equals(">=")) operChar = 'G';
		else if(operator.equals("=")) operChar = '=';
		
		ArrayList<Integer> pos = new ArrayList<Integer>();
		switch(operChar)
		{
		    case '>': {
				int count = 0;
				for(T curr_x: x)
				{
					if(curr_x.compareTo(el) > 0)
					{
						pos.add(count);
					}
					count++;
				}
				break;
	    	}
		
		    case '<': {
				int count = 0;
				for(T curr_x: x)
				{
					if(curr_x.compareTo(el) < 0)
					{
						pos.add(count);
					}
					count++;
				}
				break;				
	    	}
		    
		    case 'G': {
		    	int count = 0;
				for(T curr_x: x)
				{
					if(curr_x.compareTo(el) >= 0)
					{
						pos.add(count);
					}
					count++;
				}
				break;
	    	}
		    
		    case 'L': {
		    	int count = 0;
				for(T curr_x: x)
				{
					if(curr_x.compareTo(el) <= 0)
					{
						pos.add(count);
					}
					count++;
				}
				break;
	    	}
		    
		    case '=': {
		    	int count = 0;
				for(T curr_x: x)
				{
					if(curr_x.compareTo(el) == 0)
					{
						pos.add(count);
					}
					count++;
				}
				break;
	    	}
		    
		    default: {
		    	throw new IllegalArgumentException("You have entered an invalid value for the operator argument.\nValid values are: { >, <, >=, <=, = }");
				
	    	}
		}
		
		return pos;
	}//end of find method 
	
	
	/**
	 * Permutes an integer array  <br>
	 * <u>Note</u>: Assumes that <tt>permInds</tt> are in a valid form. That is,
	 * they are valued from 0 to a.length-1 (obviously, not necessarily in order). 
	 * @param an integer array to be permuted
	 * @param permInds permutation indices
	 * @return the permuted array
	 */
	public static int[] permute(int[] a, int[] permInds) {
		
		if(permInds.length != a.length)
			throw new IllegalArgumentException("a and permInds must have the same length.");
		
		int[] temp = new int[a.length];
		for(int i = 0; i < a.length; i++)
			temp[i] = a[permInds[i]];
		
		return temp;
	}//end of permute method
	
	/**
	 * Permutes a float array  <br>
	 * <u>Note</u>: Assumes that <tt>permInds</tt> are in a valid form. That is,
	 * they are valued from 0 to a.length-1 (obviously, not necessarily in order). 
	 * @param a float array to be permuted
	 * @param permInds permutation indices
	 * @return the permuted array
	 */
	public static float[] permute(float[] a, int[] permInds) {
		
		if(permInds.length != a.length)
			throw new IllegalArgumentException("a and permInds must have the same length.");
		
		float[] temp = new float[a.length];
		for(int i = 0; i < a.length; i++)
			temp[i] = a[permInds[i]];
		
		return temp;
	}//end of permute method
	
	
	/**
	 * Permutes a char array  <br>
	 * <u>Note</u>: Assumes that <tt>permInds</tt> are in a valid form. That is,
	 * they are valued from 0 to a.length-1 (obviously, not necessarily in order). 
	 * @param a char array to be permuted
	 * @param permInds permutation indices
	 * @return the permuted array
	 */
	public static char[] permute(char[] a, int[] permInds) {
		
		if(permInds.length != a.length)
			throw new IllegalArgumentException("a and permInds must have the same length.");
		
		char[] temp = new char[a.length];
		for(int i = 0; i < a.length; i++)
			temp[i] = a[permInds[i]];
		
		return temp;
	}//end of permute method

	
	/**
	 * Permutes an array
	 * @param <T> type of array to be sorted
	 * @param a Array of type <tt>T</tt> to be permuted
	 * @param permInds permutation indices
	 * @return the permuted array
	 */
	public static <T> T[] permute(T[] a, int[] permInds) {
		
		if(permInds.length != a.length)
			throw new IllegalArgumentException("a and permInds must have the same length.");
		
		int[] tempInds = permInds.clone(); Arrays.sort(tempInds);
		for(int i = 0; i < tempInds.length; i++) {
			if(tempInds[i] != i)
				throw new IllegalArgumentException("permInds must have all numbers from 0 to permInds.length-1.");
		}
		
		T[] temp = (T[]) Array.newInstance(a.getClass().getComponentType(), a.length);
		for(int i = 0; i < a.length; i++)
			temp[i] = a[permInds[i]];
		
		return temp;
	}//end of permute method
	
	
	/**
	 * Effective for an ORDERED integer array
	 * @see edu.mit.csail.cgs.utils.stats.StatUtil#searchFrom(Comparable[], String, Comparable, int)
	 */
	public static int searchFrom(int[] a, String operator, int value, int posIndex) {
		if( a.length == 0 )
			throw new IndexOutOfBoundsException("Hey dude, you cannot enter an empty array...");
		
		if(posIndex < 0 || posIndex > a.length)
			throw new IndexOutOfBoundsException("posIndex must be within 0 and a.length.");
		
		if(posIndex == a.length) { posIndex = a.length-1; }
		
		int ind;
		
		char operChar = 'n';
		if(operator.equals(">")) operChar = '>';
		else if(operator.equals("<")) operChar = '<';
		else if(operator.equals("<=")) operChar = 'L';
		else if(operator.equals(">=")) operChar = 'G';
		
		switch(operChar) {
		
			case 'G': {         // a >= value
				ind = a.length;
				int count = 0; int temp = -1;
				while(true) {
					if( (posIndex-count >= 0) && (a[posIndex-count] >= value ) )
						temp = posIndex-count;
					else if( (posIndex-count < 0) || (a[posIndex-count] < value ) )
						break;

					count++;
				}
				ind = temp!=-1 ? temp : ind;

				if( ind == a.length ) {
					count = 1; temp = -1;
					while(true) {
						if( (posIndex+count < a.length) && (a[posIndex+count] >= value ) )
							{ temp = posIndex+count; break; }
						else if( posIndex+count >= a.length )
							break;

						count++;				
					}
					ind = temp!=-1 ? temp : ind;	
				}

				break;
			}//end of case 'G'

			case '>': {          // a > value
				ind = a.length;
				int count = 0; int temp = -1;
				while(true) {
					if( (posIndex-count >= 0) && (a[posIndex-count] > value ) )
						temp = posIndex-count;
					else if( (posIndex-count < 0) || (a[posIndex-count] <= value ) )
						break;

					count++;
				}
				ind = temp!=-1 ? temp : ind;

				if( ind == a.length ) {
					count = 1; temp = -1;
					while(true) {
						if( (posIndex+count < a.length) && (a[posIndex+count] > value ) )
							{ temp = posIndex+count; break; }
						else if( posIndex+count >= a.length )
							break;

						count++;				
					}
					ind = temp!=-1 ? temp : ind;	
				}

				break;
			}//end of case '>'
			
			case 'L': {       // a <= value
				ind = -1;
				int count = 0; int temp = -1;
				while(true) {
					if( (posIndex+count < a.length) && (a[posIndex+count] <= value ) )
						temp = posIndex+count; 
					else if( (posIndex+count >= a.length) || (a[posIndex+count] > value ) )
						break;

					count++;
				}
				ind = temp!=-1 ? temp : ind;

				if( ind == -1 ) {
					count = 1; temp = -1;
					while(true) {
						if( (posIndex-count >= 0) && (a[posIndex-count] <= value ) )
							{ temp = posIndex-count; break; }
						else if( posIndex-count < 0 )
							break;

						count++;				
					}
					ind = temp!=-1 ? temp : ind;	
				}

				break;
			}//end of case 'L'

	
			case '<': {       // a < value
				ind = -1;
				int count = 0; int temp = -1;
				while(true) {
					if( (posIndex+count < a.length) && (a[posIndex+count] < value ) )
						temp = posIndex+count; 
					else if( (posIndex+count >= a.length) || (a[posIndex+count] >= value ) )
						break;

					count++;
				}
				ind = temp!=-1 ? temp : ind;

				if( ind == -1 ) {
					count = 1; temp = -1;
					while(true) {
						if( (posIndex-count >= 0) && (a[posIndex-count] < value ) )
							{ temp = posIndex-count; break; }
						else if( posIndex-count < 0 )
							break;

						count++;				
					}
					ind = temp!=-1 ? temp : ind;	
				}

				break;
			}//end of case '<'
			
		    default: {
		    	throw new IllegalArgumentException("You have entered an invalid value for the operator argument.\nValid values are: { >, <, >=, <= }");
				
	    	}

		}//end of switch statement		
		
		return ind; 
	}
	
	
	/**
	 * Assumes that <tt>a</tt> is an <b>ordered</b> array. <br>
	 * If <b>flag: >= </b> or <b>flag: > </b>, it will return the index <tt>ind</tt> of the smallest element  
	 * that is greater equal or strictly greater than <tt>value</tt>.   <br>
	 * That is, <tt>a[ind]</tt> to <tt>a[a.length-1]</tt> will be greater equal or
	 * strictly greater than <tt>value</tt>  <br>
	 * If nothing is found, it will return the length of <tt>a</tt> (<tt>a.length</tt>) <p>
     *
	 * If <b>flag: <= </b> or <b>flag: < </b>, it will return the index <tt>ind</tt> of the largest element  
	 * that is less equal or strictly less than <tt>value</tt>   <br>
	 * That is, <tt>a[0]</tt> to <tt>a[ind]</tt> will be less equal or
	 * strictly less than <tt>value</tt>  <br>
	 * If nothing is found, it will return <tt>-1</tt>  <p>
	 * @param <T> type of array to be sorted
	 * @param a array to be sorted
	 * @param operator operator of comparison <br>
	 * Acceptable values: ">=", "<=", ">", "<"
	 * @param value value where we want: <tt> a :operator: value </tt> <br>
	 * E.g. if <tt>operator: > </tt>, we want: <tt> a > value </tt>
	 * @param posIndex index where search will begin
	 * @return
	 */
	public static <T extends Comparable<T> > int searchFrom(T[] a, String operator, T value, int posIndex) {
		
		if( a.length == 0 )
			throw new IndexOutOfBoundsException("Hey dude, you cannot enter an empty array...");
		
		if(posIndex < 0 || posIndex > a.length)
			throw new IndexOutOfBoundsException("posIndex must be within 0 and a.length.");
		
		if(posIndex == a.length) { posIndex = a.length-1; }
		
		int ind;
		
		char operChar = 'n';
		if(operator.equals(">")) operChar = '>';
		else if(operator.equals("<")) operChar = '<';
		else if(operator.equals("<=")) operChar = 'L';
		else if(operator.equals(">=")) operChar = 'G';
		
		switch(operChar) {
		
			case 'G': {         // a >= value
				ind = a.length;
				int count = 0; int temp = -1;
				while(true) {
					if( (posIndex-count >= 0) && (a[posIndex-count].compareTo(value) >= 0 ) )
						temp = posIndex-count;
					else if( (posIndex-count < 0) || (a[posIndex-count].compareTo(value) < 0 ) )
						break;

					count++;
				}
				ind = temp!=-1 ? temp : ind;

				if( ind == a.length ) {
					count = 1; temp = -1;
					while(true) {
						if( (posIndex+count < a.length) && (a[posIndex+count].compareTo(value) >= 0 ) )
							{ temp = posIndex+count; break; }
						else if( posIndex+count >= a.length )
							break;

						count++;				
					}
					ind = temp!=-1 ? temp : ind;	
				}

				break;
			}//end of case 'G'

			case '>': {          // a > value
				ind = a.length;
				int count = 0; int temp = -1;
				while(true) {
					if( (posIndex-count >= 0) && (a[posIndex-count].compareTo(value) > 0 ) )
						temp = posIndex-count;
					else if( (posIndex-count < 0) || (a[posIndex-count].compareTo(value) <= 0 ) )
						break;

					count++;
				}
				ind = temp!=-1 ? temp : ind;

				if( ind == a.length ) {
					count = 1; temp = -1;
					while(true) {
						if( (posIndex+count < a.length) && (a[posIndex+count].compareTo(value) > 0 ) )
							{ temp = posIndex+count; break; }
						else if( posIndex+count >= a.length )
							break;

						count++;				
					}
					ind = temp!=-1 ? temp : ind;	
				}

				break;
			}//end of case '>'
			
			case 'L': {       // a <= value
				ind = -1;
				int count = 0; int temp = -1;
				while(true) {
					if( (posIndex+count < a.length) && (a[posIndex+count].compareTo(value) <= 0 ) )
						temp = posIndex+count; 
					else if( (posIndex+count >= a.length) || (a[posIndex+count].compareTo(value) > 0 ) )
						break;

					count++;
				}
				ind = temp!=-1 ? temp : ind;

				if( ind == -1 ) {
					count = 1; temp = -1;
					while(true) {
						if( (posIndex-count >= 0) && (a[posIndex-count].compareTo(value) <= 0 ) )
							{ temp = posIndex-count; break; }
						else if( posIndex-count < 0 )
							break;

						count++;				
					}
					ind = temp!=-1 ? temp : ind;	
				}

				break;
			}//end of case 'L'

	
			case '<': {       // a < value
				ind = -1;
				int count = 0; int temp = -1;
				while(true) {
					if( (posIndex+count < a.length) && (a[posIndex+count].compareTo(value) < 0 ) )
						temp = posIndex+count; 
					else if( (posIndex+count >= a.length) || (a[posIndex+count].compareTo(value) >= 0 ) )
						break;

					count++;
				}
				ind = temp!=-1 ? temp : ind;

				if( ind == -1 ) {
					count = 1; temp = -1;
					while(true) {
						if( (posIndex-count >= 0) && (a[posIndex-count].compareTo(value) < 0 ) )
							{ temp = posIndex-count; break; }
						else if( posIndex-count < 0 )
							break;

						count++;				
					}
					ind = temp!=-1 ? temp : ind;	
				}

				break;
			}//end of case '<'
			
		    default: {
		    	throw new IllegalArgumentException("You have entered an invalid value for the operator argument.\nValid values are: { >, <, >=, <= }");
				
	    	}

		}//end of switch statement		
		
		return ind; 
	}//end of searchFrom method
	
	public static Double mean(Double[] a) {
		if(a.length == 0) { throw new IllegalArgumentException("The entered array is empty."); }
		Double sum = 0.0;
		for(int n = 0; n < a.length; n++) { sum += a[n]; }
		return sum/a.length;
	}//end of mean method
	
	public static double mean(double[] a) {
		if(a.length == 0) { throw new IllegalArgumentException("The entered array is empty."); }
		double sum = 0.0;
		for(int n = 0; n < a.length; n++) { sum += a[n]; }
		return sum/a.length;
	}//end of mean method

	
	public static Double mean(Integer[] a) {
		if(a.length == 0) { throw new IllegalArgumentException("The entered array is empty."); }
		Double[] a_d = new Double[a.length];
		for(int n = 0; n < a_d.length; n++) { a_d[n] = Double.valueOf(a[n]); }
		return mean(a_d);
	}//end of mean method
	
	public static double mean(int[] a) {
		if(a.length == 0) { throw new IllegalArgumentException("The entered array is empty."); }
		double[] a_d = new double[a.length];
		for(int n = 0; n < a_d.length; n++) { a_d[n] = a[n]; }
		return mean(a_d);
	}//end of mean method
	
	
	public static double std(double[] a, double mean, boolean isUnbiasedEstimator) {
		if(a.length == 0) { throw new IllegalArgumentException("The entered array is empty."); }
		int N = a.length;
		int divisor = isUnbiasedEstimator ? N-1 : N;
		
		double sum = 0.0;
		for(int n = 0; n < N; n++) { sum += (a[n]-mean)*(a[n]-mean); }
		return Math.sqrt(sum/divisor);
	}//end of std method
	
	/** Default is the unbiased estimator std. The one where the divisor is N-1 **/
	public static double std(double[] a, double mean) { return std(a, mean, true); }
	
	public static double std(double[] a, boolean isUnbiasedEstimator) { return std(a, mean(a), isUnbiasedEstimator); }
	
	public static double std(double[] a) { return std(a, mean(a), true); }
	
	public static double std(int[] a, double mean, boolean isUnbiasedEstimator) {
		if(a.length == 0) { throw new IllegalArgumentException("The entered array is empty."); }
		int N = a.length;
		int divisor = isUnbiasedEstimator ? N-1 : N;
		
		double sum = 0.0;
		for(int n = 0; n < N; n++) { sum += (a[n]-mean)*(a[n]-mean); }
		return Math.sqrt(sum/divisor);
	}//end of std method
	
	/** Default is the unbiased estimator std. The one where the divisor is N-1 **/
	public static double std(int[] a, double mean) { return std(a, mean, true); }
	
	public static double std(int[] a, boolean isUnbiasedEstimator) { return std(a, mean(a), isUnbiasedEstimator); }
	
	public static double std(int[] a) { return std(a, mean(a), true); }
	
	public static Double std(Double[] a, Double mean, boolean isUnbiasedEstimator) {
		if(a.length == 0) { throw new IllegalArgumentException("The entered array is empty."); }
		int N = a.length;
		int divisor = isUnbiasedEstimator ? N-1 : N;
		
		Double sum = 0.0;
		for(int n = 0; n < N; n++) { sum += (a[n]-mean)*(a[n]-mean); }
		return Math.sqrt(sum/divisor);
	}//end of std method
	
	/** Default is the unbiased estimator std. The one where the divisor is N-1 **/
	public static Double std(Double[] a, Double mean) { return std(a, mean, true); }
	
	public static Double std(Double[] a, boolean isUnbiasedEstimator) { return std(a, mean(a), isUnbiasedEstimator); }
	
	public static Double std(Double[] a) { return std(a, mean(a), true); }

	public static Double std(Integer[] a, Double mean, boolean isUnbiasedEstimator) {
		if(a.length == 0) { throw new IllegalArgumentException("The entered array is empty."); }
		int N = a.length;
		int divisor = isUnbiasedEstimator ? N-1 : N;
		
		Double sum = 0.0;
		for(int n = 0; n < N; n++) { sum += (a[n]-mean)*(a[n]-mean); }
		return Math.sqrt(sum/divisor);
	}//end of std method
	
	/** Default is the unbiased estimator std. The one where the divisor is N-1 **/
	public static Double std(Integer[] a, Double mean) { return std(a, mean, true); }
	
	public static Double std(Integer[] a, boolean isUnbiasedEstimator) { return std(a, mean(a), isUnbiasedEstimator); }
	
	public static Double std(Integer[] a) { return std(a, mean(a), true); }
	
	public static double correlation(double[] x,double[] y) {
		double v = 0.0;
		double vx = 0.0;
		double vy = 0.0;
		double muX = mean(x);
		double muY = mean(y);
		int i = 0;
		for (i=0;i<x.length;i++) {
			v += (x[i] - muX)*(y[i] - muY);
			vx += Math.pow(x[i] - muX,2.0);
			vy += Math.pow(y[i] - muY,2.0);
		}
		v = v/(Math.sqrt(vx)*Math.sqrt(vy));
		return v;
	}
	
	/**
	 * Returns the odds ratio (effect size) of a group of samples with <tt>p</tt> positives
	 * and <tt>n</tt> negatives, the size of the total positive set <tt>P</tt>
	 * and the total negative size <tt>N</tt>.<br>
	 * https://en.wikipedia.org/wiki/Effect_size#Odds_ratio
	 * https://en.wikipedia.org/wiki/Odds_ratio
	 * @param P total positive set size
	 * @param N total negative set size
	 * @param p # of samples in positive
	 * @param n # of samples in negative
	 * @return
	 */
	public static double odds_ratio(int P, int N, double p, double n, double posPseudocount, double negPseudocount) {
		double or = (p+posPseudocount)/(P-p+posPseudocount) / ((n+negPseudocount)/(N-n+negPseudocount));
		return or;	
	}

	
	/**
	 * Returns the log-likelihood of an array of numbers (<tt>x</tt>) distributed
	 * normally.
	 * @param mu mean value for each element of array <tt>x</tt>
	 * @param s common variance
	 * @param x array of type <tt>double</tt>
	 * @return
	 */
	public static double logNormPDF(double[] mu, double s, double[] x) {
		int i = 0;
		double l = 0.0;
		double t = 0.0;
		for (i=0;i<x.length;i++) {
			t = x[i];
			if (!Double.isNaN(t)) {
				l = l + -0.5*Math.log(2.0*Math.PI)-Math.log(s)-Math.pow(t-mu[i],2.0)/(2.0*Math.pow(s,2.0));
			}
		}
		return l;		
	}
	
	
	/**
	 * Returns the hypergeometric cumulative probability of number <tt>x</tt> 
	 * when the sample size is <tt>n</tt>, the size of the positive set <tt>s</tt>
	 * and the population size <tt>N</tt>.
	 * The precision of the method is up to 1E-8, determined empirically
	 * This method use cern.jet.stat.Gamma.logGamma to compute factorial values. 
	 * It is more appropriate for one time on-the-fly calculation.
	 * @param x # observed successes in the sample
	 * @param N population size
	 * @param s # of successes in population (e.g., size of positive set in the population)
	 * @param n sample size
	 * @return
	 */
	public static double hyperGeometricCDF(int x, int N, int s, int n) {
		double v = 0.0;
		int k = 0;
		for (k=0;k<=x;k++) {
			v += hyperGeometricPDF(k,N,s,n);
		}
		
		if (v > 1.0)
			v = 1.0;
		return v;	
	}
	
	
	/**
	 * Returns the hypergeometric density probability of number <tt>x</tt> 
	 * when the sample size is <tt>n</tt>, the size of the positive set <tt>s</tt>
	 * and the population size <tt>N</tt>.  <br>
	 * It works with the gamma function and in log space in an attempt to avoid
	 * numerical problems.
	 * @param x # observed successes in the sample
	 * @param N population size
	 * @param s # of successes in population (e.g., size of positive set in the population)
	 * @param n sample size
	 * @return
	 */
	public static double hyperGeometricPDF(int x, int N, int s, int n) {
		double xd = (double) x;
		double Nd = (double) N;
		double sd = (double) s;
		double nd = (double) n;
		double kx = cern.jet.stat.Gamma.logGamma(sd+1)-cern.jet.stat.Gamma.logGamma(x+1)-cern.jet.stat.Gamma.logGamma(sd-x+1);
		double mn = cern.jet.stat.Gamma.logGamma(Nd+1)-cern.jet.stat.Gamma.logGamma(nd+1)-cern.jet.stat.Gamma.logGamma(Nd-nd+1);
		double mknx = cern.jet.stat.Gamma.logGamma(Nd-sd+1)-cern.jet.stat.Gamma.logGamma(nd-x+1)-cern.jet.stat.Gamma.logGamma(Nd-sd-(nd-x)+1);
		return Math.exp(kx + mknx - mn);
	}
	
	/**
	 * Returns the hypergeometric cumulative probability of number <tt>x</tt> 
	 * when the sample size is <tt>n</tt>, the size of the positive set <tt>s</tt>
	 * and the population size <tt>N</tt>.<br>
	 * If the CDF is close to 1, the precision of the method is up to 1E-8, determined empirically
	 * This method use a cache to store the factorial values. 
	 * So it is more suitable to compute many hgp values with a fix population.
	 * @param x # observed successes in the sample
	 * @param N population size
	 * @param s # of successes in population (e.g., size of positive set in the population)
	 * @param n sample size
	 * @return
	 */
	public static double hyperGeometricCDF_cache(int x, int N, int s, int n) {
		String key = String.format("%d_%d_%d_%d", x, N, s, n);
		if (HGP_table.containsKey(key)){
			cacheAccessCount++;
			return HGP_table.get(key);
		}
		else{
			double v = 0.0;
			int k = 0;
			for (k=0;k<=x;k++) {
				double v1= hyperGeometricPDF_cache(k,N,s,n);
				v += v1;
			}
			
			if (v > 1.0)
				v = 1.0;
			HGP_table.put(key, v);
			return v;	
		}
	}
	/**
	 * Returns the hypergeometric cumulative probability
	 * High precision implementation with BigDecimal, but very slow	
	 * @param x # observed successes in the sample
	 * @param N population size
	 * @param s # of successes in population (e.g., size of positive set in the population)
	 * @param n sample size
	 */
	public static double log10_hyperGeometricCDF_cache(int x, int N, int s, int n) {
		BigDecimal v = BigDecimal.ZERO;
		for (int k=0;k<=x;k++) {
			v = v.add(hyperGeometricPDF_cache_BIG(k,N,s,n));
//			v = v.setScale(v.scale()-v.precision()+20, BigDecimal.ROUND_DOWN);	
		}
		int newScale  = v.scale()-v.precision()+10;
		BigDecimal v2 = v.setScale(newScale, BigDecimal.ROUND_DOWN);		
		double d = v2.unscaledValue().doubleValue();
		return Math.log10(d)-v2.scale();	
	}
	
	/** Returns the hypergeometric cumulative probability
	 *  Pretty high precision implementation of log10_hyperGeometricCDF, much faster than BigDecimal version 
	 *  Summing the PDF on log10 space, avoid losing precision 
	 * @param x # observed successes in the sample
	 * @param N population size
	 * @param s # of successes in population (e.g., size of positive set in the population)
	 * @param n sample size
	 **/
	public static double log10_hyperGeometricCDF_cache_appr(int x, int N, int s, int n) {
		String key = String.format("%d_%d_%d_%d", x, N, s, n);
		if (log10_HGP_table.containsKey(key)){
			cacheAccessCount++;
			return log10_HGP_table.get(key);
		}
		else{
			double Lx = log10_hyperGeometricPDF_cache(x,N,s,n);
			double sum = 1;
			for (int k=x-1;k>=0;k--) {
				sum += Math.pow(10, log10_hyperGeometricPDF_cache(k,N,s,n)-Lx);
			}
			Lx += Math.log10(sum);
			HGP_table.put(key, Lx);
			return Lx;
		}
	}
	/**
	 * Returns the hypergeometric density probability of number <tt>x</tt> 
	 * when the sample size is <tt>n</tt>, the size of the positive set <tt>s</tt>
	 * and the population size <tt>N</tt>.  <br>
	 * This method use a cache to store the factorial values. 
	 * So it is more suitable to compute many hgp values with a fix population.
	 * @param x # observed successes in the sample
	 * @param N population size
	 * @param s # of successes in population (e.g., size of positive set in the population)
	 * @param n sample size
	 * @return
	 */
	public static double hyperGeometricPDF_cache(int x, int N, int s, int n) {
		if (x+N-s-n<0)
			return 0;
		String key = String.format("%d_%d_%d_%d", x, N, s, n);
		if (HGPDF_table.containsKey(key)){
			cacheAccessCount++;
			return HGPDF_table.get(key);
		}
		else{
	
			// extend the cache if N is large than existing cache
			int len = (logFactorials==null)?0:logFactorials.length;
			if (logFactorials==null || logFactorials.length < N+1){
				double[] old = logFactorials;
				logFactorials = new double[N+1];
				if (len!=0)
					System.arraycopy(old, 0, logFactorials, 0, len);
				else{
					logFactorials[0]=0;
					len++;
				}
				for (int i=len;i<=N;i++)
					logFactorials[i] = logFactorials[i-1]+Math.log(i);
			}
			// compute
			double kx = logFactorials[s]-logFactorials[x]-logFactorials[s-x];
			double mknx = logFactorials[N-s]-logFactorials[n-x]-logFactorials[N-s-(n-x)];
			double mn = logFactorials[N]-logFactorials[n]-logFactorials[N-n];
			double result = Math.exp(kx + mknx - mn);
			HGPDF_table.put(key, result);
			return result;
		}
	}
	
	public static double log10_hyperGeometricPDF_cache (int x, int N, int s, int n) {
		if (x+N-s-n<0)
			return Double.NEGATIVE_INFINITY;
		String key = String.format("%d_%d_%d_%d", x, N, s, n);
		if (log10_HGPDF_table.containsKey(key)){
			cacheAccessCount++;
			return log10_HGPDF_table.get(key);
		}
		else{
		// extend the cache if N is large than existing cache
		int len = (logFactorials==null)?0:logFactorials.length;
		if (logFactorials==null || logFactorials.length < N+1){
			double[] old = logFactorials;
			logFactorials = new double[N+1];
			if (len!=0)
				System.arraycopy(old, 0, logFactorials, 0, len);
			else{
				logFactorials[0]=0;
				len++;
			}
			for (int i=len;i<=N;i++)
				logFactorials[i] = logFactorials[i-1]+Math.log(i);
		}
		// compute
		double kx = logFactorials[s]-logFactorials[x]-logFactorials[s-x];
		double mknx = logFactorials[N-s]-logFactorials[n-x]-logFactorials[N-s-(n-x)];
		double mn = logFactorials[N]-logFactorials[n]-logFactorials[N-n];
		double result = (kx + mknx - mn)*Math.log10(Math.exp(1));
		log10_HGPDF_table.put(key, result);

		return result;
		}
	}
	
	/**Use BigDecimal to compute hyperGeometric, high precision, but very slow, use for verification
	 * @param x # observed successes in the sample
	 * @param N population size
	 * @param s # of successes in population (e.g., size of positive set in the population)
	 * @param n sample size
	 */
	public static BigDecimal hyperGeometricPDF_cache_BIG (int x, int N, int s, int n) {
		if (x+N-s-n<0)
			return BigDecimal.ZERO;
		// extend the cache if N is large than existing cache
		int len = (logFactorials==null)?0:logFactorials.length;
		if (logFactorials==null || logFactorials.length < N+1){
			double[] old = logFactorials;
			logFactorials = new double[N+1];
			if (len!=0)
				System.arraycopy(old, 0, logFactorials, 0, len);
			else{
				logFactorials[0]=0;
				len++;
			}
			for (int i=len;i<=N;i++)
				logFactorials[i] = logFactorials[i-1]+Math.log(i);
		}
		// compute
		double kx = logFactorials[s]-logFactorials[x]-logFactorials[s-x];
		double mknx = logFactorials[N-s]-logFactorials[n-x]-logFactorials[N-s-(n-x)];
		double mn = logFactorials[N]-logFactorials[n]-logFactorials[N-n];
		BigDecimal v = exp_BIG(kx + mknx - mn);
		return v;
	}
	
	public static BigDecimal exp_BIG(double v){
		int p = (int) Math.floor(Math.log10(v<0?-v:v));
		double unscaledValue = v/Math.pow(10, p);
		BigDecimal exp = BigDecimal.valueOf(Math.exp(unscaledValue)).pow((int)Math.round(Math.pow(10, p)));
		return exp;
	}

	/**
	 * 
	 * @param n
	 * @return
	 */
	public static int[] randPerm(int n) {
		int[] order = new int[n];
		ArrayList<Integer> nm = new ArrayList<Integer>();
		int i = 0;
		for (i=0;i<n;i++) {
			nm.add(i);
		}
		
		int c = 0;
		for (i=n-1;i>=0;i--) {
			if (i > 0) {
				c = Uniform.staticNextIntFromTo(0,i);
			} else {
				c = 0;
			}
			order[i] = nm.get(c);
			nm.remove(c);
		}
		return order;
		//Uniform.staticNextIntFromTo(0,numClusters-1);
	}
	
	
	/* 
	 * smooth the curve in y[]
	 */
	public static double[] cubicSpline(double[] y, int stepSize, int averageWidth){
		// averageWidth should be smaller than stepSize
		if (stepSize<averageWidth)
			averageWidth = stepSize;
		
		// sample some data points to create the spline
		double[]x= new double[y.length/stepSize];
		double[]yx = new double[y.length/stepSize];
		// first sample point, average from first data points
		x[0]=0;
		yx[0]=0;
		for (int j=0;j<averageWidth;j++){
			yx[0]+=y[j]/averageWidth;
		}
		// middle points
		for (int i=1; i<x.length-1; i++){
			x[i]=i*stepSize+stepSize/2;
			yx[i]=0;
			for (int j=0;j<averageWidth;j++){
				yx[i]+=y[(int)x[i]+(j-averageWidth/2)]/averageWidth;
			}
		}
		// last sample point, average from last data points
		x[x.length-1] = y.length-1;
		yx[x.length-1] = 0;
		for (int j=0;j<averageWidth;j++){
			yx[x.length-1]+=y[y.length-1-j]/averageWidth;
		}
		
		CubicSpline cs =new CubicSpline(x, yx);
		
		//read smoothed data from the spline
		double[] yy = new double[y.length];
		for (int i=0; i<y.length; i++){
			yy[i]=cs.interpolate(i);
		}
		return yy;
	}
	
	// take a prob. density dist., smooth it using Gaussian kernel density
	// assume Y is evenly spaced.
	//TODO: boundary effect
	public static double[] gaussianSmoother( double[]Y, double width){
		int length = Y.length;
		double[] yy = new double[length];
		double[] yy_weight = new double[length];
		// width --> stdev, width*width --> variance
		NormalDistribution gaussian = new NormalDistribution(0, width*width);
		double total=0;
		for (int i=0;i<length;i++){
			yy[i]=2.0E-300;		// init with very small number
			yy_weight[i]=2.0E-300;	
			for (int j=0;j<length;j++){
				double gaussianProb = gaussian.calcProbability((double)j-i);
				yy[i] += Y[j]*gaussianProb;
				yy_weight[i] += gaussianProb;
			}
			yy[i] /= yy_weight[i];
			total += yy[i];
		}
		// normalize
		for (int i=0;i<length;i++){
			yy[i] /= total;
		}
		return yy;
	}
	/** take a prob. density dist., smooth it using specified kernel density
	 * 
	 * @param Y should be evenly spaced.
	 * @param kernel kernel density should be symmetrical (e.g. Gaussian).
	 * @return normalized smooth prob. density dist.
	 */
	public static double[] symmetricKernelSmoother( double[]Y, double[]kernel){
		int length = Y.length;
		int kernel_length = kernel.length;
		double[] yy = new double[length];
		double total=0;
		for (int i=0;i<length;i++){
			double v=kernel[0]*Y[i] + 2.0E-300;		// init with very small number
            double weight=kernel[0];
            for (int j = 1; j < kernel_length && i+j < length; j++) {
                v+=Y[i+j]*kernel[j];
                weight += kernel[j];                
            }
            for (int j = 1; j < kernel_length && i-j >= 0; j++) {
                v+=Y[i-j]*kernel[j];
                weight += kernel[j];                
            }
			v = v / weight;
            yy[i] = v;
			total+=v;
		}
		for (int i=0;i<length;i++){
			yy[i]=yy[i]/total;
		}
		return yy;
	}	
	
	// K-L divergence for discrete probability distributions P and Q 
	//	http://en.wikipedia.org/wiki/Kullback-Leibler_divergence
	public static double KL_Divergence( double[]P, double[]Q){
		double d=0;
		double[] p=P.clone();	// clone so that the input array is not changed
		double[] q=Q.clone();
		// Make sure that p and q are all proper probability values
		mutate_normalize(p);
		mutate_normalize(q);
		for (int i=0;i<p.length;i++){
			d+=p[i]*(Math.log(p[i])-Math.log(q[i]));
		}
		return d;
	}
	// log of K-L divergence 
	public static double log_KL_Divergence( double[]P, double[]Q){
		return Math.log(KL_Divergence(P,Q));
	}
	public static double log10_KL_Divergence( double[]reference, double[]Q){
		return Math.log10(KL_Divergence(reference,Q));
	}
	
	// Uses COLT binomial test
	// Binomial CDF assuming scaled control. k=scaled control, n=scaled control+signal
	public static double binomialPValue(double k, double n){
        return binomialPValue(k,n,.5);
    }
	/**
	 * Binomial test (COLT lib)
	 * @param k 
	 * @param n	Total
	 * @param p	Prob
	 * @return
	 */
	
    public static double binomialPValue(double k, double n, double p){
		if (n==0)
			return Double.NaN;
		double pval=1;
		Binomial b = new Binomial((int)Math.ceil(n), p, new DRand());
		pval = b.cdf((int) Math.ceil(k));
		return(pval);		
	}
	
	/** this method will mutate the input array */
	public static void mutate_normalize(double[] dist){
		double total=0;
		for (int i=0;i<dist.length;i++){
			if (dist[i]<=0)
				dist[i]=1e-20;
			total += dist[i];
		}
		for(int i=0;i<dist.length;i++){
			dist[i] /=total;
		}
	}
	
	public static double entropy(double[] p) {
		double sum = 0.0;
		for(int j = 0; j < p.length; j++) { sum += p[j]*Math.log(p[j]); }
		sum *= -1;
		return sum;
	}//end of entropy method
	

	/**
	 * This method normalizes the given array (vector) and returns the normalizing constant.  <br>
	 * @param a matrix to be normalized
	 * @return
	 */
	public static double normalize(double[] a) {
		double sum = 0.0;
		for(int i = 0; i < a.length; i++){ sum += a[i]; }
		if(sum == 0) { sum = 1;}
		for(int i = 0; i < a.length; i++){ a[i] /= sum; }
		return sum;
	}// end of normalize method

	
	/**
	 * This method normalizes the given matrix and returns the normalizing constant.  <br>
	 * However, it treats the whole matrix as an array (vector).
	 * @param a matrix to be normalized
	 * @return
	 */
	public static double normalize(double[][] a) {
		double sum = 0.0;
		for(int i = 0; i < a.length; i++){ for(int j = 0; j < a[i].length; j++) { sum += a[i][j]; } }
		if(sum == 0) { sum = 1;}
		for(int i = 0; i < a.length; i++){ for(int j = 0; j < a[i].length; j++) { a[i][j] /= sum; } }
		return sum;
	}// end of normalize method
	
	
	/**
	 * This method will normalize the matrix on the specified dimension and 
	 * returns the normalizing constants.
	 * @param a matrix to be normalized
	 * @param dim dimension on which the normalization will take place. <br>
	 * dim = 1, normalize w.r.t. to the first dimension, aka rows.  <br>
	 * dim = 2, normalize w.r.t. to the second dimension, aka columns.
	 * @return
	 */
	public static double[] normalize(double[][] a, int dim) {
		double[] sum ;
		if(dim == 1)
			sum = new double[a.length];
		else if(dim == 2)
			sum = new double[a[0].length];
		else
			throw new IllegalArgumentException("The only valid values for dim are: 1 (rows), 2 (columns).");
		
		if(dim == 1) {
			for(int i = 0; i < a.length; i++) {
				for(int j = 0; j < a[i].length; j++)
					sum[i] += a[i][j];
				
				if(sum[i] == 0) { sum[i] = 1; };
			}
			
			for(int i = 0; i < a.length; i++) {
				for(int j = 0; j < a[i].length; j++)
					 a[i][j] /= sum[i];
			}
		}
		else {
			for(int j = 0; j < a[0].length; j++) {
				for(int i = 0; i < a.length; i++)
					sum[j] += a[i][j];
				
				if(sum[j] == 0) { sum[j] = 1; };
			}
			
			for(int j = 0; j < a[0].length; j++) {
				for(int i = 0; i < a.length; i++)
					 a[i][j] /= sum[j];
			}			
		}
		
		return sum;		
	}// end of normalize method
	
	
	/**
	 * This method will normalize the matrix on the specified dimension and 
	 * returns the normalizing constants.
	 * When normalizing w.r.t. rows or columns, it does so for each of element of
	 * the third dimension separately.
	 * @param a matrix to be normalized
	 * @param dim dimension on which the normalization will take place. <br>
	 * dim = 1, normalize w.r.t. to the first dimension. <br>
	 * dim = 2, normalize w.r.t. to the second dimension, aka rows.  <br>
	 * dim = 3, normalize w.r.t. to the third dimension, aka columns.
	 * @return
	 */
	public static double[][] normalize(double[][][] a, int dim) {
		double[][] sum ;
		if(dim == 1)
			sum = new double[a[0].length][a[0][0].length];
		else if(dim == 2)
			sum = new double[a.length][a[0].length];
		else if(dim == 3)
			sum = new double[a.length][a[0][0].length];
		else
			throw new IllegalArgumentException("The only valid values for dim are: 1 (rows), 2 (columns).");
		
		if(dim == 2) {
			for(int k = 0; k < a.length; k++) {
				for(int i = 0; i < a[k].length; i++) {
					for(int j = 0; j < a[k][i].length; j++)
						sum[k][i] += a[k][i][j];
					
					if(sum[k][i] == 0) { sum[k][i] = 1; }
				}		
			}
			
			for(int k = 0; k < a.length; k++) {
				for(int i = 0; i < a[k].length; i++) {
					for(int j = 0; j < a[k][i].length; j++)
						 a[k][i][j] /= sum[k][i];
				}		
			}			
		}
		else if(dim == 3){
			for(int k = 0; k < a.length; k++) {
				for(int j = 0; j < a[k][0].length; j++) {
					for(int i = 0; i < a[k].length; i++)
						sum[k][j] += a[k][i][j];
					
					if(sum[k][j] == 0) { sum[k][j] = 1; }
				}		
			}
			
			for(int k = 0; k < a.length; k++) {
				for(int j = 0; j < a[k][0].length; j++) {
					for(int i = 0; i < a[k].length; i++)
						 a[k][i][j] /= sum[k][j];
				}		
			}
		}
		else {
			for(int i = 0; i < a[0].length; i++) {
				for(int j = 0; j < a[0][i].length; j++) {
					for(int k = 0; k < a.length; k++)
						sum[i][j] += a[k][i][j];
				
				if(sum[i][j] == 0) { sum[i][j] = 1; }					
				}
			}		
			
			for(int i = 0; i < a[0].length; i++) {
				for(int j = 0; j < a[0][i].length; j++)
					for(int k = 0; k < a.length; k++)
						 a[k][i][j] /= sum[i][j];
				}		
		}
		
		return sum;		
	}// end of normalize method

	/** data structure to store the points for density clustering
	 */
	public class DensityClusteringPoint implements Comparable<DensityClusteringPoint>{
		public int id;		// the original id of the data point
		public float density=0;		// density as total count by the neighborhood
		public float densitySxN=0;		// sqrt (self density * neighborhood density)
		public float delta=0;
		public float gamma=0;
		public int delta_id;
		public ArrayList<DensityClusteringPoint> members = new ArrayList<DensityClusteringPoint>();
		public TreeSet<Integer> memberIds = new TreeSet<Integer>();
		
		public String toString(){
			return String.format("id=%d\tdensity=%.1f\tdelta=%.1f\tgamma=%.1f\tdelta_id=%d\tmembers=%d", id, density, delta, gamma, delta_id, members.size());
		}
		/**
		 * Sort DensityClusteringPoint by descending density value (default)
		 */		
		public int compareTo(DensityClusteringPoint p) {
			if(density>p.density){return(-1);}
			else if(density<p.density){return(1);}
			else return(0);
		}
		/**
		 * Sort DensityClusteringPoint by descending gamma value
		 */
		public int compareToByGamma(DensityClusteringPoint p) {
			if(gamma>p.gamma){return(-1);}
			else if(gamma<p.gamma){return(1);}
			else return(0);
		}
		/**
		 * Sort DensityClusteringPoint by ascending id
		 */
		public int compareToById(DensityClusteringPoint p) {
			if(id<p.id){return(-1);}
			else if(id>p.id){return(1);}
			else return(0);
		}
		/** Sort DensityClusteringPoint by descending S*N density value
		 */
		public int compareToByDensitySxN(DensityClusteringPoint p) {
			if(densitySxN>p.densitySxN){return(-1);}
			else if(densitySxN<p.densitySxN){return(1);}
			else return(0);
		}
		
	}
	


	private static void printHGP(int pos, int neg){
		System.out.println(String.format("%d\t%d\t%.2f\t%.2f", pos, neg, 		 
				KMAC.computeHGP(5000, 5000, pos, neg),
				odds_ratio(5000, 5000, pos, neg, 3, 3)));
	}
	
	private static void printHGP(int total, int motif1, int motif2, int overlap){
		System.out.println(String.format("%d\t%d\t%d\t%d\t%.2f\t%.2f", total, motif1, motif2, overlap,		 
				KMAC.computeHGP(motif1, total-motif1, overlap, motif2-overlap),
				odds_ratio(motif1, total-motif1, overlap, motif2-overlap, 10, 10)));
	}
	 public static void main(String[] args){
//		 System.out.println( Math.log10(binomialPValue(0.0, 11.0+0.0)));
//		 System.out.println( Math.log10(binomialPValue(3.3, 24.0+3.3)));
//		 Poisson poisson = new Poisson(0, new DRand());
//		 poisson.setMean(2.7);System.out.println(1-poisson.cdf(56));
//		 System.out.println(poisson.pdf(1));
//		System.out.println(log10_hyperGeometricCDF_cache_appr(4,10000,5000,4+1));
//		System.out.println(hyperGeometricCDF_cache(2405,41690+40System.out.println(hyperGeometricCDF(3298,10000,5000,3298+2));506,41690,2405+2));
//		 printHGP(1480,465);
//		 printHGP(1325,328);
//		 printHGP(1360,327);
//		 printHGP(1517,467);
		 printHGP(3456, 741,919,300);
//		 for (int i=0;i<=5000;i+=100){
//			 for (int j=0;j<=5000;j+=100){
//				 printHGP(i, j);
//			 }
//		 }

//		System.out.println(hyperGeometricCDF(2405,41690+40506,41690,2405+2));
//		System.out.println(hyperGeometricCDF_cache(2,41690+40506,40506,3298+2));
//		System.out.println(hyperGeometricCDF_cache(2,41690+40506,40506,2405+2));
//		System.out.println(hyperGeometricCDF_cache(3,8+7,8,3+2));
//		System.out.println(hyperGeometricCDF_cache(2,8+7,7,3+2));
//		 System.out.println(log10_hyperGeometricCDF_cache_appr(1,41690+40506,1,5000));
//		for (int i=10700;i<=10800;i=i+10){
////			System.out.println(Math.log10(hyperGeometricPDF_cache_BIG(99,41690+40506,40506,i+100).doubleValue()));
////			System.out.println(log10_hyperGeometricPDF_cache(i,41690+40506,40506,i+100));
//			System.out.println(i+"\t"+log10_hyperGeometricCDF_cache_appr(i,40876+40873,40873,i+31761));
////			System.out.println(log10_hyperGeometricCDF_cache(99,41690+40506,40506,i+5000));
//			System.out.println("-----------");
//		}
//		System.out.println(Long.MAX_VALUE);
//		 System.out.println(correlation(new double[]{1.0, 2,3,4},new double[]{2.0, -5,6,-8}));
	}
}//end of StatUtil class 41690 / 40506

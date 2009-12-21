/*
 * Author: tdanford
 * Date: Sep 8, 2008
 */
package edu.mit.csail.cgs.utils;

import java.lang.reflect.Array;
import java.util.*;

public class ArrayUtils {
	
	public static Integer[] range(int start, int end) {
		if(start > end) { 
			throw new IllegalArgumentException(String.format("[%d, %d) illegal coordinates.", 
					start, end));
		}
		Integer[] array = new Integer[end-start];
		for(int i = start, j = 0; i < end; i++, j++) { 
			array[j] = i;
		}
		return array;
	}
	
	/** 
	 * Returns an array composed of only the subset of items from the original 
	 * array which are accepted by the given predicate. 
	 * 
	 * @param <T>
	 * @param a
	 * @param pred
	 * @return
	 */
	public static <T> T[] mask(T[] a, Predicate<T> pred) { 
		TreeSet<Integer> inds = new TreeSet<Integer>();
		for(int i = 0; i < a.length; i++) { 
			if(pred.accepts(a[i])) { 
				inds.add(i);
			}
		}
		int len = inds.size();
		Class cls = a.getClass().getComponentType();
		T[] sarray = (T[])Array.newInstance(cls, len);
		int i = 0;
		for(Integer ind : inds) { 
			sarray[i++] = a[ind];
		}
		return sarray;
	}
	
	/**
	 * Reverses the given array -- returns a copy of the array in reversed order 
	 * (in other words, the ordering of the original array is left unmodified).  
	 * 
	 * @param <T>
	 * @param a
	 * @return
	 */
	public static <T> T[] reverse(T[] a) { 
		T[] array = a.clone();
		for(int i = 0; i < array.length; i++) { 
			array[i] = a[a.length-i-1];
		}
		return array;
	}
	
	public static <T> T[] append(T[] a, T last) { 
		if(a == null) { throw new IllegalArgumentException("Null array to ArrayUtils.append()"); }
		Class cls = a.getClass().getComponentType();
		T[] sarray = (T[])Array.newInstance(cls, a.length+1);
		sarray[sarray.length-1] = last;
		for(int i = 0; i < a.length; i++){ 
			sarray[i] = a[i];
		}
		return sarray;		
	}
	
	public static <T> T[] prepend(T first, T[] a) { 
		Class cls = first.getClass().getComponentType();
		T[] sarray = (T[])Array.newInstance(cls, a.length+1);
		sarray[0] = first;
		for(int i = 0; i < a.length; i++){ 
			sarray[i+1] = a[i];
		}
		return sarray;
	}
	
	/**
	 * Returns a new array whose elements are the elements of the arrays a1 and 
	 * a2, concatenated in order.  
	 * 
	 * @param <T>
	 * @param a1
	 * @param a2
	 * @return
	 */
	public static <T> T[] concat(T[] a1, T[] a2) {
		int len = a1.length + a2.length;
		Class cls = a1.getClass().getComponentType();
		int length = Math.max(0, len);
		T[] sarray = (T[])Array.newInstance(cls, length);
		for(int i = 0; i < a1.length; i++) { 
			sarray[i] = a1[i];
		}
		for(int i = 0; i < a2.length; i++) { 
			sarray[a1.length+i] = a2[i];
		}
		return sarray;
	}
	
	public static <T> T[] cat(T[]... as) { 
		T[] array = as.length > 0 ? as[0] : null;
		for(int i = 1; i < as.length; i++) { 
			array = concat(array, as[i]);
		}
		return array;
	}

	public static <T> Iterator<T> asIterator(T[] array) {
		return new ArrayIterator<T>(array);
	}
	
	public static <T> T[] subArray(T[] array, int start, int end) { 
		Class cls = array.getClass().getComponentType();
		int length = Math.max(0, end-start);
		T[] sarray = (T[])Array.newInstance(cls, length);
		for(int i = start; i < end; i++) { 
			sarray[i-start] = array[i];
		}
		return sarray;
	}

	public static <T> T[] tail(T[] c) {
		return subArray(c, 1, c.length);
	}

	public static <T> Collection<T> asCollection(T[] array) {
		ArrayList<T> lst = new ArrayList<T>();
		for(int i = 0; i < array.length; i++) { 
			lst.add(array[i]);
		}
		return lst;
	}

	public static <T> T[] asArray(T[] arr, Iterator<T> itr) {
		ArrayList<T> lst = new ArrayList<T>();
		while(itr.hasNext()) { 
			lst.add(itr.next());
		}
		Class tclass = arr.getClass().getComponentType();
		return lst.toArray(arr);
	}
}

class ArrayIterator<T> implements Iterator<T> {
	
	private T[] array;
	private int idx;
	
	public ArrayIterator(T[] a) { 
		array = a;
		idx = 0;
	}

	public boolean hasNext() {
		return idx < array.length;
	}

	public T next() {
		return array[idx++];
	}

	public void remove() {
		throw new UnsupportedOperationException();
	} 
}

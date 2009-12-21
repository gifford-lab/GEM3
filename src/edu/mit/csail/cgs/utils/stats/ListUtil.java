package edu.mit.csail.cgs.utils.stats;

import java.util.*;

public class ListUtil {

	public static double euclidean(List<Double> x,List<Double> y) {
		double v = 0;
		int i = 0;
		for (i=0;i<x.size();i++) {
			v += Math.pow(x.get(i)-y.get(i),2.0);
		}
		return Math.sqrt(v);
	}
	
	public static String toString(List l) {
		StringBuffer buf = new StringBuffer();
		for (int i=0; i<l.size()-1; i++) {
			buf.append(l.get(i).toString());
			buf.append("\t");
		}
		buf.append(l.get(l.size()-1).toString());
		return buf.toString();
	}
	
	public static double sum(List<Double> l) {
		double sum = 0.0d;
		for (Double d : l) {
			sum += d;
		}
		return sum;
	}
	
	public static double average(List<Double> arr) {
		double sum = 0.0d;
		int count = 0;
		for (int i=0; i<arr.size(); i++) {
			if (!Double.isNaN(arr.get(i))) {
				sum += arr.get(i);
				count++;
			}
		}
		return sum / ((double)count);
	}
	
}

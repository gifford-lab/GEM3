package edu.mit.csail.cgs.utils.stats;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Vector;

public class VectorUtil {
	public static void mult(Vector<Double> x, double m) {
		int i = 0;
		for (i=0;i<x.size();i++) {
			x.set(i,x.get(i)*m);
		}
	}
	
	public static void add(Vector<Double> x,double a) {
		int i = 0;
		for(i=0;i<x.size();i++) {
			x.set(i,x.get(i) + a);
		}
	}
	
	public static double mean(Vector<Double> x) {
		double m = 0;
		double n = 0.0;
		int i = 0;
		for (i=0;i<x.size();i++) {
			if (!Double.isNaN(x.get(i))) {
				m = m + x.get(i);
				n = n + 1.0;
			}
		}
		return m/n;
	}
	// standard deviation
	public static double stdev(Vector<Double> x) {
		double mean = mean(x);
		double m = 0;
		double n = 0.0;
		int i = 0;
		for (i=0;i<x.size();i++) {
			if (!Double.isNaN(x.get(i))) {
				m = m + (x.get(i)-mean)*(x.get(i)-mean);
				n = n + 1.0;
			}
		}
		return Math.sqrt(m/n);
	}
	
	public static double sum(Vector<Double> x) {
		double m = 0;
		int i = 0;
		for (i=0;i<x.size();i++) {
			if (!Double.isNaN(x.get(i))) {
				m = m + x.get(i);
			}
		}
		return m;
	}
	
	public static double max(Vector<Double> x) {
		double m = x.get(0);
		int i = 0;
		for (i=0;i<x.size();i++) {
			if (!Double.isNaN(x.get(i))) {
				if (x.get(i) > m) {
					m = x.get(i);
				}
			}
		}
		return m;
	}
	
	public static void zeros(Vector<Double> x,int n) {
		x.clear();
		int i = 0;
		for (i=0;i<n;i++) {
			x.add(new Double(0.0));
		}
	}
	
	public static double correlation(Vector<Double> x,Vector<Double> y) {
		double v = 0.0;
		double vx = 0.0;
		double vy = 0.0;
		double muX = mean(x);
		double muY = mean(y);
		int i = 0;
		for (i=0;i<x.size();i++) {
			v += (x.get(i) - muX)*(y.get(i) - muY);
			vx += Math.pow(x.get(i) - muX,2.0);
			vy += Math.pow(y.get(i) - muY,2.0);
		}
		
		v = v/(Math.sqrt(vx)*Math.sqrt(vy));
		return v;
	}
	
	public static double euclidean(Vector<Double> x,Vector<Double> y) {
		double v = 0;
		int i = 0;
		for (i=0;i<x.size();i++) {
			v += Math.pow(x.get(i)-y.get(i),2.0);
		}
		return Math.sqrt(v);
	}
	
	public static double percentile(Vector<Double> x,double p) {
		int rp = (int) Math.round(p*((double) x.size()));
		
		Vector<Double> y = new Vector<Double>(x);
		Collections.sort(y);
		return y.get(rp);
	}
	
	public static int[] sortOrder(Vector<Double> x) {
		int[] order = new int[x.size()];
		int i = 0;
		Vector<DoubleSort> v = new Vector<DoubleSort>();
		DoubleSort ds = null;
		for (i=0;i<x.size();i++) {
			ds = new DoubleSort();
			ds.pos = i;
			ds.value = x.get(i);
			v.add(ds);
		}
		
		Collections.sort(v);
		
		for (i=0;i<v.size();i++) {
			ds = v.get(i);
			order[i] = ds.pos;
		}
		
		return order;
	}
}

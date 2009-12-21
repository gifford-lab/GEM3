package edu.mit.csail.cgs.utils.stats;

// maximum likelihood fit to Dirichlet distribution
// uses the simple fixed-point iteration described in
// "Estimating a Dirichlet distribution" by T. Minka.
public class DirichletFit {
	public static double[] fit(double[][] p) {
		double[] a = momentMatch(p);
		
		
		
		return a;
	}
	
	public static double[] momentMatch(double[][] p) {
		int i = 0;
		int j = 0;
		int n = p.length;
		int m = p[0].length;
		double[] a = new double[m];
		double[] m2 = new double[m];
		
		for (i=0;i<n;i++) {
			for (j=0;j<m;j++) {
				a[j] = a[j] + p[i][j];
				m2[j] = m2[j] + Math.pow(p[i][j],2.0);
			}
		}
		
		for (j=0;j<m;j++) {
			a[j] = a[j]/((double) n);
			m2[j] = m2[j]/((double) n);
		}
		
		double[] s = new double[m];
		for (j=0;j<m;j++) {
			s[j] = (a[j] - m2[j])/(m2[j]-Math.pow(a[j],2.0));
		}
		
		// each dimension of p gives an independent estimate of s, so take the median.
		double sm = StatUtil.median(s);
		for (j=0;j<m;j++) {
			a[j] = a[j]*sm;
		}
		return a;
	}
}

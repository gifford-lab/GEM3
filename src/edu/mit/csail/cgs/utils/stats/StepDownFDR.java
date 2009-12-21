package edu.mit.csail.cgs.utils.stats;

// implements the step-down false discovery rate method
// for multiple hypothesis testing from Benjamini 95
public class StepDownFDR {
	// pvals = list of p-values (must be sorted from smallest to largest)
	// FDR = false discovery rate
	// returns an array of booleans, indicating whether the p-value
	// is significant
	public static boolean[] correctPvals(double[] pvals,double FDR) {
		int n = pvals.length;
		boolean[] acceptPval = new boolean[n];
		int i = 0;
		for (i=0;i<n;i++) {
			acceptPval[i] = false;
		}
		
		// now correct the p-values
		boolean contt = true;
		double ap = 0.0;
		i = n-1;
		while(i >= 0 && contt) {
			ap = ((double) (i + 1))*FDR/((double) n);
			if (pvals[i] <= ap) {
				contt = false;
			} else {
				i--;
			}
		}
		
		int j = 0;
		if (i>=0) {
			for(j=0;j<=i;j++) {
				acceptPval[j] = true;
			}
		}
		
		return acceptPval;
	}
}

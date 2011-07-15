/*
 * Author: tdanford
 * Date: May 22, 2008
 */
package edu.mit.csail.cgs.ewok.verbs.motifs;

import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;

public class WeightMatrixScoreProfile {

	private WeightMatrix matrix;
	private double[] forward, reverse;
	
	public WeightMatrixScoreProfile(WeightMatrix m, double[] f, double[] r) {
		matrix = m;
		if(f.length != r.length) { throw new IllegalArgumentException(); }
		forward = (double[])f.clone();
		reverse = (double[])r.clone();
	}
	
	public int length() { return forward.length; }
	public WeightMatrix getMatrix() { return matrix; }
	public double[] getForwardScores() { return forward; }
	public double[] getReverseScores() { return reverse; }
	
	public char getMaxStrand() { 
		int max = -1;
		char str = '+';
		double maxScore = 0.0;
		for(int i = 0; i < forward.length; i++) { 
			double ms = getMaxScore(i);
			if(max == -1 || ms > maxScore) { 
				maxScore= ms; 
				max = i;
				str = getMaxStrand(i);
			}
		}
		return str;
	}
	
	public char getMaxStrand(int i) { 
		if(forward[i] >= reverse[i]) { 
			return '+';
		} else { 
			return '-';
		}
	}
	/** max score (from 2 strands) at this position */
	public double getMaxScore(int i) { 
		return Math.max(forward[i], reverse[i]); 
	}
	/** max score over the whole sequence*/
	public double getMaxScore(){
		return(getMaxScore(getMaxIndex()));
	}
	public int getMaxIndex() { 
		int max = -1; 
		double maxScore = matrix.getMinScore();
		for(int i = 0; i < forward.length; i++) { 
			double ms = getMaxScore(i);
			if(max == -1 || ms > maxScore) { 
				maxScore= ms; 
				max = i;
			}
		}
		return max;
	}
}

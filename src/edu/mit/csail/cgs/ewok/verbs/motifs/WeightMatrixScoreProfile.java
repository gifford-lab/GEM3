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
			double ms = getHigherScore(i);
			if(max == -1 || ms > maxScore) { 
				maxScore= ms; 
				max = i;
				str = getHigherScoreStrand(i);
			}
		}
		return str;
	}
	
	public char getHigherScoreStrand(int i) { 
		char s = getMaxStrand_both(i);
		return s=='='?'+':s;
	}
	/** 
	 * Get the strand that gives larger score on position i<br>
	 * This is a more accurate version of getMaxStrand(i), considering palindromic motif may match on both strands
	 * @param i position in the sequence
	 * @return strand, '=' represents matches on both strand have the same scores
	 */
	public char getMaxStrand_both(int i) { 
		if(forward[i] > reverse[i]) { 
			return '+';
		} 
		else if(forward[i] < reverse[i]) { 
			return '-';
		}
		else{ 
			return '=';
		}
	}
	
	/** Higher score (from 2 strands) at this position */
	public double getHigherScore(int i) { 
		return Math.max(forward[i], reverse[i]); 
	}
	/** get score (from forward strand) at this position */
	public double getForwardScore(int i) { 
		return forward[i]; 
	}
	/** max score over the whole sequence*/
	public double getMaxScore(){
		return(getHigherScore(getMaxIndex()));
	}
	/** max score over the whole sequence*/
	public double getMaxForwordScore(){
		return(getForwardScore(getMaxIndex()));
	}
	public int getMaxIndex() { 
		int max = -1; 
		double maxScore = matrix.getMinScore();
		for(int i = 0; i < forward.length; i++) { 
			double ms = getHigherScore(i);
			if(max == -1 || ms > maxScore) { 
				maxScore= ms; 
				max = i;
			}
		}
		return max;
	}
	public int getForwardMaxIndex() { 
		int max = -1; 
		double maxScore = matrix.getMinScore();
		for(int i = 0; i < forward.length; i++) { 
			double ms = getForwardScore(i);
			if(max == -1 || ms > maxScore) { 
				maxScore= ms; 
				max = i;
			}
		}
		return max;
	}
	/**
	 * Find the best score and position/strand in the sequence
	 * @param strand_type
	 * @return
	 */
	public PwmMatch getBestMatch(int strand_type){
	      double maxSeqScore = Double.NEGATIVE_INFINITY;
	      int maxScoringShift = 0;
	      char maxScoringStrand = '+';
	      for (int j=0;j<forward.length;j++){
	    	  if (strand_type!=1){
	    		  double score = getHigherScore(j);
		    	  if (maxSeqScore<score || (maxSeqScore==score && maxScoringStrand=='-')){	// equal score, prefer on '+' strand
		    		  maxSeqScore = score;
		    		  maxScoringShift = j;
		    		  maxScoringStrand = getHigherScoreStrand(j);
		    	  }
	    	  }
	    	  else{
	    		  double score = getForwardScore(j);
	    		  if (maxSeqScore<score){
	    			  maxSeqScore = score;
		    		  maxScoringShift = j;
	    		  }
	    	  }
	      }
	      return new PwmMatch(maxSeqScore, maxScoringShift, maxScoringStrand);
	}
	
	public class PwmMatch{
		public PwmMatch(double score, int shift, char strand){
			this.score = score;
			this.shift = shift;
			this.strand = strand;
		}
		public double score;
		public int shift;
		public char strand;
	}

}



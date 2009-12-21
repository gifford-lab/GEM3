package edu.mit.csail.cgs.utils.stats;

public class InfoTheoryUtil {
	
	// the symbols in the sequence are coded as 0...numStates-1
	// the entropy if returned in nats
	// p is a vector to calculate probabilities
	public static double discreteEntropy(int[] sequence,int[] p) {
		double v = 0.0;
		double vt = 0.0;
		int i = 0;
		int numStates = p.length;
		
		for (i=0;i<numStates;i++) {
			p[i] = 0;
		}
		
		for (i=0;i<sequence.length;i++) {
			p[sequence[i]]++;
		}
		for (i=0;i<numStates;i++) {
			if (p[i] > 0) {
				vt = ((double) p[i])/((double) sequence.length);
				v = v - vt*Math.log(vt);
			}
		}
		
		if (v <= 0.0) {
			v = 0.0;
		}
		
		return v;
	}
	
	public static double discreteEntropy(int[] sequence,int numStates) {
		int[] p = new int[numStates];
		return discreteEntropy(sequence,p);
	}
	
	// the sequences must be of the same length
	// returns the mutual information in nats
	public static double discreteMutualInformation(int[] sequenceX,int[] sequenceY,int[] pX, int[] pY, int[][] pXY) {
		double v = 0.0;
		double vt = 0.0;
		int i = 0;
		int j = 0;
		int numStatesX = pX.length;
		int numStatesY = pY.length;
		
		v = v + discreteEntropy(sequenceX,pX);
		v = v + discreteEntropy(sequenceY,pY);
		
		for (i=0;i<numStatesX;i++) {
			for (j=0;j<numStatesY;j++) {
				pXY[i][j] = 0;
			}
		}
		
		for (i=0;i<sequenceX.length;i++) {
			pXY[sequenceX[i]][sequenceY[i]]++;
		}
		
		for (i=0;i<numStatesX;i++) {
			for (j=0;j<numStatesY;j++) {
				if (pXY[i][j] > 0) {
					vt = ((double) pXY[i][j])/((double) sequenceX.length);
					v = v + vt*Math.log(vt);
				}
			}
		}
		
		if (v <= 0.0) {
			v = 0.0;
		}
		
		return v;
	}
	
	public static double discreteMutualInformation(int[] sequenceX, int[] sequenceY,int numStates) {
		int[] pX = new int[numStates];
		int[] pY = new int[numStates];
		int[][] pXY = new int[numStates][numStates];
		return discreteMutualInformation(sequenceX,sequenceY,pX,pY,pXY);
	}
	
	// zeros are treated specially (essentially as missing values) - only entries for which
	// one sequence has non-zero entries are considered
	public static double discreteMutualInformationOverlap(int[] sequenceX,int[] sequenceY,int[] pX, int[] pY, int[][] pXY) {
		double v = 0.0;
		double vt = 0.0;
		int i = 0;
		int j = 0;
		int numStatesX = pX.length;
		int numStatesY = pY.length;
		int seqL = 0;
		
		for (i=0;i<sequenceX.length;i++) {
			if (sequenceX[i] > 0 | sequenceY[i] > 0) {
				seqL++;
				pX[sequenceX[i]]++;
				pY[sequenceY[i]]++;
				pXY[sequenceX[i]][sequenceY[i]]++;
			}
		}
		
		// calculate entropy of X
		for (i=0;i<numStatesX;i++) {
			if (pX[i] > 0) {
				vt = ((double) pX[i])/((double) seqL);
				v = v - vt*Math.log(vt);
			}
		}
		
		//  calculate entropy of Y
		for (i=0;i<numStatesY;i++) {
			if (pY[i] > 0) {
				vt = ((double) pY[i])/((double) seqL);
				v = v - vt*Math.log(vt);
			}
		}
		
		// calculate joint entropy
		for (i=0;i<numStatesX;i++) {
			for (j=0;j<numStatesY;j++) {
				if (pXY[i][j] > 0) {
					vt = ((double) pXY[i][j])/((double) seqL);
					v = v + vt*Math.log(vt);
				}
			}
		}
		
		if (v <= 0.0) {
			v = 0.0;
		}
		
		return v;
	}
	
	public static double discreteMutualInformationOverlap(int[] sequenceX, int[] sequenceY,int numStates) {
		int[] pX = new int[numStates];
		int[] pY = new int[numStates];
		int[][] pXY = new int[numStates][numStates];
		return discreteMutualInformationOverlap(sequenceX,sequenceY,pX,pY,pXY);
	}
}

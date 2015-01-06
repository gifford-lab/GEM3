package edu.mit.csail.cgs.utils.stats;

import java.util.HashMap;
import java.util.TreeSet;

import edu.mit.csail.cgs.utils.Pair;

public class ROC {
	double[] sortedScores;
	HashMap<Double, Integer> posScoreMap, negScoreMap;
	int posTotalCount, negTotalCount;
	double[] tprs, fprs; // true positive rates, false positive rates
	
	public ROC(double[] posScores, double[] negScores){
		posTotalCount = posScores.length;
		negTotalCount = negScores.length;
		
		posScoreMap = new HashMap<Double, Integer>();
		for (double d: posScores){
			if (posScoreMap.containsKey(d))
				posScoreMap.put(d, posScoreMap.get(d)+1);
			else
				posScoreMap.put(d, 1);
		}
		negScoreMap = new HashMap<Double, Integer>();
		for (double d: negScores){
			if (negScoreMap.containsKey(d))
				negScoreMap.put(d, negScoreMap.get(d)+1);
			else
				negScoreMap.put(d, 1);
		}
		
		TreeSet<Double> sScores = new TreeSet<Double>();		
		sScores.addAll(posScoreMap.keySet());
		sScores.addAll(negScoreMap.keySet());
		sortedScores = new double[sScores.size()];
		int idx=0;
		for (double s: sScores.descendingSet())	{	// descending sorting
			sortedScores[idx] = s;
			idx++;
		}
		
		for (double s: sortedScores){
			if (!posScoreMap.containsKey(s))
				posScoreMap.put(s,0);
			if (!negScoreMap.containsKey(s))
				negScoreMap.put(s,0);
		}
		
		tprs = new double[sortedScores.length+1];
		fprs = new double[sortedScores.length+1];
		tprs[0]=0;
		fprs[0]=0;
		for (int i=1;i<tprs.length;i++){
			double score = sortedScores[i-1];
			tprs[i] = tprs[i-1] + posScoreMap.get(score)/(double)posTotalCount;
			fprs[i] = fprs[i-1] + negScoreMap.get(score)/(double)negTotalCount;
		}
	}
	
	/**
	 * Compute partial ROC given a false positive rate
	 * @param falsePositiveRate #FP/#negatives, i.e. 1-specificity
	 * @return
	 */
	public double partialROC (double falsePositiveRate){
		if (falsePositiveRate > 1)
			falsePositiveRate = 1;
		
		double truePositiveRate = 0;			// interpolated TPR corresponding to the FPR
		int fprIdx = 0;			// the index for the smallest point that is larger than falsePositiveRate
		for (int i=1;i<fprs.length;i++){
			if (fprs[i]==falsePositiveRate){
				truePositiveRate=tprs[i];
				fprIdx = i;
				break;
			}
			else if(fprs[i]>falsePositiveRate){
				fprIdx = i;
				truePositiveRate = (tprs[i-1]-tprs[i])/(fprs[i-1]-fprs[i])*(falsePositiveRate-fprs[i])+tprs[i];
				break;
			}
			else
				continue;				
		}
		
		double area = 0;
		for (int i=1;i<fprIdx;i++)
			area += (tprs[i-1]+tprs[i])*(fprs[i]-fprs[i-1])/2;
		area += (tprs[fprIdx-1]+truePositiveRate)*(falsePositiveRate-fprs[fprIdx-1])/2;
		
		return area;		
	}
	
	/**
	 * Compute partial optimal cutoff point by maximizing Youden Index (i.e. sensitivity + specificity)
     * 
	 * @param falsePositiveRate #FP/#negatives, i.e. 1-specificity
	 * @return the optimal score and index for retrieving optimal sensitivity and false positive rate <br><br>
	 * 
     * Ref: http://www.medicalbiostatistics.com/roccurve.pdf
     * Youden index (J) that maximizes the vertical distance from line of equality to the point [x, y] as shown in Figure 3. The x represents (1- specificity) and y represents sensitivity. In other words, the Youden index J is the point on the ROC curve which is farthest from line of equality (diagonal line). The main aim of Youden index is to maximize the difference between TPR (sn) and FPR (1 – sp) and little algebra yields J = max[sn+sp]. The value of J for continuous test can be located by doing a search of plausible values where sum of sensitivity and specificity can be maximum. Youden index is more commonly used criterion because this index reflects the intension to maximize the correct classification rate and is easy to calculate. Many authors advocate this criterion.
	 */
	public Pair<Double,Integer> partialOptimalPoint (double falsePositiveRate){
		if (falsePositiveRate > 1)
			falsePositiveRate = 1;
		
    	double max = 1;
    	int maxIdx = 0;
    	for (int i=0;i<fprs.length;i++){
    		if (fprs[i]>falsePositiveRate){
    			if (max==1){
    				max = tprs[i] - fprs[i]+1;
    				maxIdx = i;
    			}
    			break;
    		}
    		double J = tprs[i] - fprs[i]+1;
    		if (J>max){
    			max = J;
    			maxIdx = i;
    		}
    	}
    	return new Pair<Double,Integer>(sortedScores[maxIdx-1], maxIdx);
	}
	
	public Pair<Double,Double> getPoint(int index){
		return new Pair<Double, Double>(tprs[index],fprs[index]);
	}
	public Pair<Integer,Integer> getHitCounts(double score){
		int posHit=0;
		int negHit=0;
		for (double s:sortedScores){
			if (s<score)
				break;
			posHit+=posScoreMap.get(s);
			negHit+=negScoreMap.get(s);
		}
		return new Pair<Integer, Integer>(posHit, negHit);
	}
	public double getScore (int index){
		return sortedScores[index-1];
	}
}

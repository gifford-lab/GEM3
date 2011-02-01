package edu.mit.csail.cgs.ewok.verbs.motifs;

import java.util.ArrayList;
import java.util.HashSet;

import edu.mit.csail.cgs.ewok.verbs.motifs.KmerEngine.KmerMatch;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;

public class Kmer implements Comparable<Kmer>{
	String kmerString;
	public String getKmerString() {	return kmerString;}
	int k;
	public int getK(){return k;}
	int count;	//all hit count
	int seqHitCount; //one hit at most for one sequence, to avoid simple repeat
	double weight;
	double hg;
	int negCount;
	
	Kmer reference;		//
	int shift;			// the shift to achieve best score
	public int getShift(){return shift;}
	int score;			// best possible number of matches (with or w/o shift) wrt reference Kmer
	public int getScore(){return score;}
	
	
	public Kmer(String kmerStr, int hitCount){
		this.kmerString = kmerStr;
		this.k = kmerString.length();
		this.seqHitCount = hitCount;
	}
	
	/** 
	 * Use reverse compliment to represent the kmer
	 */
	public void RC(){
		kmerString = getKmerRC();
	}
	public int getNegCount() {
		return negCount;
	}
	public void setNegCount(int negCount) {
		this.negCount = negCount;
	}

	// sort kmer by weight
	public int compareByWeight(Kmer o) {
		double diff = o.weight-weight;
		return diff==0?kmerString.compareTo(o.kmerString):(diff<0)?-1:1;  // descending
	}
	
	// default, sort kmer by seqHitCount
	public int compareTo(Kmer o) {
		double diff = o.seqHitCount-seqHitCount;
		return diff==0?kmerString.compareTo(o.kmerString):(diff<0)?-1:1; // descending
	}
	public String toString(){
		return kmerString+"\t"+seqHitCount+"\t"+negCount+"\t"+String.format("%.1f", Math.log10(hg));
	}
	public static String toHeader(){
		return "EnrichedKmer\tPosCt\tNegCt\tHGP_10";
	}

	public int getSeqHitCount() {
		return seqHitCount;
	}
	public void setSeqHitCount(int count) {
		seqHitCount=count;
	}
	
	/** 
	 * Set the reference Kmer
	 * Find the best score, strand(RC), shift for this Kmer to align with reference Kmer 
	 * @param ref
	 */
	public void setReference(Kmer ref){
		reference = ref;
		byte[] thisBytes = kmerString.getBytes();
		byte[] refBytes = ref.kmerString.getBytes();
		score = 0;
		shift = -99;
		for (int s=-2;s<=+2;s++){
			int count = 0;
			for (int i=-2;i<refBytes.length+2;i++){
				if (i<0 || i>refBytes.length-1 ||i+s<0 || i+s>refBytes.length-1 )
					continue;
				if (refBytes[i]==thisBytes[i+s]){
					count ++;
				}
			}
			if (count>score){
				score = count;
				shift = s;
			}
		}
		// try RC
		byte[] rcBytes = getKmerRC().getBytes();
		boolean useRC=false;
		for (int s=-2;s<=+2;s++){
			int count = 0;
			for (int i=-2;i<refBytes.length+2;i++){
				if (i<0 || i>refBytes.length-1 ||i+s<0 || i+s>refBytes.length-1 )
					continue;
				if (refBytes[i]==rcBytes[i+s]){
					count ++;
				}
			}
			if (count>score){
				score = count;
				shift = s;
				useRC = true;
			}
		}
		if (useRC)
			RC();
	}
	
	/**
	 * calculate the best shift for input kmer to align with this kmer
	 * allow for 2 mismatches, or 1 shift + 1 mismatch, or 2 shift
	 * @param kmer
	 * @return best shift for input kmer
	 */
//	public int shift(String kmer){
//		byte[] thisBytes = this.kmerString.getBytes();
//		byte[] kmerBytes = kmer.getBytes();
//		int maxScore = 0;
//		int maxScoreShift = -99;
//		for (int s=-2;s<=+2;s++){
//			int score = 0;
//			for (int i=-2;i<thisBytes.length+2;i++){
//				if (i<0 || i>thisBytes.length-1 ||i+s<0 || i+s>thisBytes.length-1 )
//					continue;
//				if (thisBytes[i]==kmerBytes[i+s]){
//					score ++;
//				}
//			}
//			if (score>k*0.8){
//				if (score>maxScore){
//					maxScore = score;
//					maxScoreShift = s;
//				}
//			}
//		}
//		return maxScoreShift;
//	}
	
	public String getKmerRC(){
		return SequenceUtils.reverseComplement(kmerString);
	}
	
//	public static void main(String[] args){
//		Kmer k1 = new Kmer("CCAGAAGAGGGC", 32);
//		Kmer k2 = new Kmer("CCCTCTTCTGGC", 3);
//		k2.setReference(k1);
//		System.out.println(k2.getShift());
//	}
	
	
	
	
	/** 
	 ***************************** Some old code *****************************
	 */
	ArrayList<KmerMatch> hits;
	HashSet<Integer> seqHits;
	public Kmer(String kmerStr, ArrayList<KmerMatch> hits){
		this.kmerString = kmerStr;
		this.hits = hits;
		count = hits.size();
		k = kmerStr.length();
		seqHits = new HashSet<Integer>();
		for (KmerMatch hit:hits){
			seqHits.add(hit.seqId);
		}
		seqHitCount = seqHits.size();
	}	
	public double bias(int pos){
		double sum=0;
		for (KmerMatch hit:hits)
			sum+=Math.abs(hit.pos-pos);
		return sum/hits.size();
	}
	public double bias2(int pos){
		double sum=0;
		for (KmerMatch hit:hits)
			sum+=(hit.pos-pos)*(hit.pos-pos);
		return Math.sqrt(sum/hits.size());
	}
	public float[] getPositionCounts(int seqLength){
		float counts[] = new float[seqLength];
		for (KmerMatch hit:hits)
			counts[hit.pos]++;
		return counts;
	}
	public void calcWeight(double[] seqProbs, double[] positionProbs){
		weight = 0;
		for (KmerMatch hit:hits){
			int id = hit.seqId;
			if (seqHits.contains(id)||seqHits.contains(-id))
				weight += seqProbs[id]*positionProbs[hit.pos];
		}
	}
}

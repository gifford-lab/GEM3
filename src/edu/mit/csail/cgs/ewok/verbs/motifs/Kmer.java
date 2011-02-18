package edu.mit.csail.cgs.ewok.verbs.motifs;

import edu.mit.csail.cgs.utils.sequence.SequenceUtils;

public class Kmer implements Comparable<Kmer>{
	String kmerString;
	public String getKmerString() {	return kmerString;}
	int k;
	public int getK(){return k;}
	int seqHitCount; //one hit at most for one sequence, to avoid simple repeat
	double strength;	// the total read counts from all events support this kmer
	public double getStrength(){return strength;}
	public void setStrength(double strength){this.strength = strength;}
	double hg;
	int negCount;
	
	int cluster=-1;			// non-negative integer, type of clustered motif 
	public void setCluster(int c){cluster=c;}
	public int getCluster(){return cluster;}
	Kmer reference;		//
	public Kmer getRef(){return reference;}
	int shift;			// the shift to achieve best score
	public int getShift(){return shift;}
	int score;			// best possible number of matches (with or w/o shift) wrt reference Kmer
	public int getScore(){return score;}	
	int globalShift;			// the shift to the first seed
	public int getGlobalShift(){return globalShift;}
	public void setGlobalShift(int s){globalShift=s;}
	
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
		double diff = o.strength-strength;
		return diff==0?kmerString.compareTo(o.kmerString):(diff<0)?-1:1;  // descending
	}
	
	// default, sort kmer by seqHitCount
	public int compareTo(Kmer o) {
		double diff = o.seqHitCount-seqHitCount;
		return diff==0?kmerString.compareTo(o.kmerString):(diff<0)?-1:1; // descending
	}
	public String toString(){
		return kmerString+"\t"+seqHitCount+"\t"+negCount+"\t"+
			   String.format("%.1f", Math.log10(hg))+"\t"+String.format("%.1f", strength);
	}
	public static String toHeader(){
		return "EnrichedKmer\tPosCt\tNegCt\tHGP_10\tWeight";
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
	 * Extend this kmer from reference kmer, if they are only offset base off
	 * @param ref
	 */
	public boolean extendKmer(Kmer ref, int offset){
		reference = ref;
		byte[] thisBytes = kmerString.getBytes();
		byte[] refBytes = ref.kmerString.getBytes();
		score = 0;
		shift = -99;
		for (int s=-offset;s<=offset;s++){
			for (int i=-offset;i<refBytes.length+offset;i++){
				if (i<0 || i>refBytes.length-1 ||i+s<0 || i+s>refBytes.length-1 )
					continue;
				if (refBytes[i]!=thisBytes[i+s])	// if mismatch
					break;
			}
			score = refBytes.length-Math.abs(s);
			shift = s;
		}
		// try RC
		byte[] rcBytes = getKmerRC().getBytes();
		boolean useRC=false;
		for (int s=-offset;s<=offset;s++){
			for (int i=-offset;i<refBytes.length+offset;i++){
				if (i<0 || i>refBytes.length-1 ||i+s<0 || i+s>refBytes.length-1 )
					continue;
				if (refBytes[i]!=rcBytes[i+s])	// if mismatch
					break;
			}
			int thisScore = refBytes.length-Math.abs(s);
			if (thisScore>score){
				score = thisScore;
				useRC = true;
			}
			shift = s;
		}
		
		if (useRC)
			RC();
		
		return score>0;
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
	
}

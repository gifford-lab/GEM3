package edu.mit.csail.cgs.deepseq.discovery.kmer;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;

import edu.mit.csail.cgs.utils.Pair;

/**
 * This KmerGroup class is used for recording the overlapping kmer instances mapped to the same binding position in a sequence<br>
 * This version is for KMAC1
 * @author yuchun
 */
public class KmerGroup implements Comparable<KmerGroup>{
	ArrayList<Kmer> kmers;
	int bs = 999;
	int clusterId = -1;
	int posHitGroupCount;
	int negHitGroupCount;
	String coveredSequence=null;
	boolean isKmersSorted=false;
	
	/** KmerGroup significance score [ OR or -log10(hgp) ] based on config.use_odds_ratio  */
	double kg_score;
	public double getScore() {return kg_score;	}
	public void setScore(double score) {this.kg_score = score;	}
	
//	public KmerGroup(int[] posCoveredWidth, int[] negCoveredWidth, ArrayList<Kmer> kmers, int bs){
//		if (kmers.isEmpty())
//			return;
//		this.bs = bs;
//		this.kmers = kmers;
//		
//		BitSet b_pos = new BitSet(kmers.get(0).posBits.length());
//		BitSet b_neg = new BitSet(kmers.get(0).negBits.length());
// 		for (Kmer km:kmers){
// 			b_pos.or(km.posBits);
// 			b_neg.or(km.negBits);
//		}
// 		
// 		// adjust hit count by requiring the KG hit should be equal or better than the training sequence to count it
// 		if (posCoveredWidth!=null){
//	 		int width = getCoveredWidth();
//	 		for (int i = b_pos.nextSetBit(0); i >= 0; i = b_pos.nextSetBit(i+1))
//	 			if(posCoveredWidth[i]>width)	// don't count those seqs that expect a better (wider k-mer coverage) hit
//	 				b_pos.clear(i);
//	 		for (int i = b_neg.nextSetBit(0); i >= 0; i = b_neg.nextSetBit(i+1))
//	 			if(negCoveredWidth[i]>width)	// don't count those seqs that expect a better (wider k-mer coverage) hit
//	 				b_neg.clear(i);
// 		}
// 		
// 		posHitGroupCount = b_pos.cardinality();
//		negHitGroupCount = b_neg.cardinality();
//	}	
	
	public KmerGroup(int[] posCoveredWidth, int[] negCoveredWidth, ArrayList<Kmer> kmers, int bs, double[]weights){
		if (kmers.isEmpty())
			return;
		this.bs = bs;
		this.kmers = kmers;			
		BitSet b_pos = new BitSet(kmers.get(0).posBits.length());
		BitSet b_neg = new BitSet(kmers.get(0).negBits.length());
 		for (Kmer km:kmers){
 			b_pos.or(km.posBits);
 			b_neg.or(km.negBits);
		}
 		
 		// adjust hit count by requiring that the KG hit should be equal or better than the training sequence to count it
 		if (posCoveredWidth!=null){
	 		int width = getCoveredWidth();
	 		for (int i = b_pos.nextSetBit(0); i >= 0; i = b_pos.nextSetBit(i+1))
	 			if(posCoveredWidth[i]>width)	// don't count those seqs that expect a better (wider k-mer coverage) hit
	 				b_pos.clear(i);
	 		for (int i = b_neg.nextSetBit(0); i >= 0; i = b_neg.nextSetBit(i+1))
	 			if(negCoveredWidth[i]>width)	// don't count those seqs that expect a better (wider k-mer coverage) hit
	 				b_neg.clear(i);
 		}

		if (weights==null){
    		posHitGroupCount = b_pos.cardinality();
		}
		else{
    		double weight=0;
    		for (int i = b_pos.nextSetBit(0); i >= 0; i = b_pos.nextSetBit(i+1))
    			weight+=weights[i];
    		posHitGroupCount = (int)(weight);
		}
		negHitGroupCount = b_neg.cardinality();
	}


	public KmerGroup(String[] posMatchStr, String[] negMatchStr, ArrayList<Kmer> kmers, int bs, String matchedSequence, double[]weights){
		if (kmers.isEmpty())
			return;
		this.bs = bs;
		this.kmers = kmers;	
		this.coveredSequence = matchedSequence;
		BitSet b_pos = new BitSet(kmers.get(0).posBits.length());
		BitSet b_neg = new BitSet(kmers.get(0).negBits.length());
 		for (Kmer km:kmers){
 			b_pos.or(km.posBits);
 			b_neg.or(km.negBits);
		}
 		
 		// adjust hit count by requiring that the KG matched string should contains the training hit sequence
 		if (posMatchStr!=null){
	 		int width = getCoveredWidth();
	 		for (int i = b_pos.nextSetBit(0); i >= 0; i = b_pos.nextSetBit(i+1))
	 			if(!matchedSequence.contains(posMatchStr[i]))	// don't count those seqs that are not sub-string of match
	 				b_pos.clear(i);
	 		for (int i = b_neg.nextSetBit(0); i >= 0; i = b_neg.nextSetBit(i+1))
	 			if(!matchedSequence.contains(negMatchStr[i]))	// don't count those seqs that are not sub-string of match
	 				b_neg.clear(i);
 		}

		if (weights==null){
    		posHitGroupCount = b_pos.cardinality();
		}
		else{
    		double weight=0;
    		for (int i = b_pos.nextSetBit(0); i >= 0; i = b_pos.nextSetBit(i+1))
    			weight+=weights[i];
    		posHitGroupCount = (int)(weight);
		}
		negHitGroupCount = b_neg.cardinality();
	}
	
	/** get the KG match end indices, relative to the sequence<br>
	 * 	return the start of the leftmost k-mer, and the end of the rightmost k-mer
	 * */
	public Pair<Integer,Integer> getMatchEndIndices(){
		int leftIdx = 999;
		int rightIdx = -999;
		for (Kmer km:kmers){
			if (km.getKmerStartOffset()<leftIdx)
				leftIdx = km.getKmerStartOffset();
			if (km.getKmerStartOffset()+km.kmerString.length()>rightIdx)
				rightIdx = km.getKmerStartOffset()+km.kmerString.length();
		}
		return new Pair<Integer,Integer>(leftIdx+bs, rightIdx+bs);
	}
	
	public String getCoveredSequence(){
		if (coveredSequence==null)
			coveredSequence = makeCoveredSequence();
		return coveredSequence;
	}
	
	public int getCoveredWidth(){
		int i=0;
		for (char c : getCoveredSequence().toCharArray())
			if (c!='N')
				i++;
		return i;
	}

	private String makeCoveredSequence(){
		int leftShift = 999;			// left most shift position
		int rightShift = -999;			// right most shift position
		int longest = 0;
		for (Kmer km:kmers){
			if (km.getKmerStartOffset() < leftShift)
				leftShift = km.getKmerStartOffset();
			if (km.getKmerStartOffset() > rightShift)
				rightShift = km.getKmerStartOffset();
			if (km.k > longest)
				longest = km.k;
		}
		char[] letters = new char[rightShift-leftShift+longest];
		int rightIdx = 0;
		for (Kmer km:kmers){
			char[] kChars = km.isSeedOrientation?km.kmerString.toCharArray():km.kmerRC.toCharArray();
			for (int i=0;i<kChars.length;i++){
				int idx = i+km.getKmerStartOffset()-leftShift;
				if (kChars[i]!='N'){
					letters[idx] = kChars[i];
					if (rightIdx<idx)
						rightIdx = idx;							
				}
			}
		}
		StringBuilder sb = new StringBuilder();
		for (int i=0;i<=rightIdx;i++){
			if (letters[i]!=0)
				sb.append(letters[i]);
			else
				sb.append('N');
		}
		return sb.toString();
	}
	
	public int getClusterId(){return clusterId;}
	public void setClusterId(int id){clusterId = id;}
	
	public ArrayList<Kmer> getKmers(){
		return kmers;
	}
	public Kmer getBestKmer(){
		if (kmers.isEmpty())
			return null;
		if (!isKmersSorted){
			Collections.sort(this.kmers);
			isKmersSorted = true;
		}
		return kmers.get(0);
	}
	public String getAllKmerString(){
		if (kmers.isEmpty())
			return null;
		if (!isKmersSorted){
			Collections.sort(this.kmers);
			isKmersSorted = true;
		}
		StringBuilder sb = new StringBuilder();
		for (Kmer km:kmers){
			sb.append(km.kmerString).append(",");
		}
		return sb.toString();
	}
//	public int getTotalKmerCount(){
//		int kmerCountSum = 0;
//		for (Kmer kmer:kmers){
//    		kmerCountSum+=kmer.getPosHitCount();	
//		}
//		return kmerCountSum;
//	}
	/** Get the number of sequences hit by any kmer in the group */
	public int getGroupHitCount(){
		return posHitGroupCount;
	}
	public int getGroupNegHitCount(){
		return negHitGroupCount;
	}
	public double getTotalKmerStrength(){
		double total = 0;
		for (Kmer kmer:kmers){
			// on first kpp round, kmers do not have strength value, use count here
			total+=kmer.getStrength()>1?kmer.getStrength():kmer.getPosHitCount();	
		}
		return total;
	}	
	/** Get the weighted kmer strength<cr>
	 *  The weight is 1 for top kmer, 1/k for other kmer 
	 *  Note: this is only approximate */
	public double getWeightedKmerStrength(){
		double total = kmers.get(0).getWeightedHitCount();
		double k = kmers.get(0).getK();
		for (int i=1;i<kmers.size();i++){
			total+=kmers.get(i).getWeightedHitCount()/k;	
		}
		return total;
	}			
	public int getPosBS(){
		return bs;
	}
	public int compareToByPosHitCount(KmerGroup kg) {		// descending pos hit count
		if(posHitGroupCount>kg.getGroupHitCount()){return(-1);}
		else if(posHitGroupCount<kg.getGroupHitCount()){return(1);}
		else return(0);
	}
	public int compareTo(KmerGroup kg) {					// descending score, ascending hgp
		if(kg_score>kg.getScore()){return(-1);}
		else if(kg_score<kg.getScore()){return(1);}
		else return(0);
	}
	public String toString(){
		Kmer best = getBestKmer();
		return String.format("%s|%b %d: %d+/%d-, kg_score=%.2f", best!=null?best.getKmerStrRC():"NULL", best!=null?getBestKmer().isSeedOrientation():false, bs, posHitGroupCount, negHitGroupCount, kg_score);
	}
}

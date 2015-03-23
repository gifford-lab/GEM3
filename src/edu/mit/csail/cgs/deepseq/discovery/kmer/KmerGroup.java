package edu.mit.csail.cgs.deepseq.discovery.kmer;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;

/**
 * This KmerGroup class is used for recording the overlapping kmer instances mapped to the same binding position in a sequence
 * @author yuchun
 */
public class KmerGroup implements Comparable<KmerGroup>{
	ArrayList<Kmer> kmers;
	int clusterId = -1;
	int bs = 999;
	int posHitGroupCount;
	int negHitGroupCount;
	/** hgp (log10) using the positive/negative sequences */
	double hgp;

	public KmerGroup(ArrayList<Kmer> kmers, int bs, int posSeqCount, int negSeqCount){
		this.bs = bs;
		this.kmers = kmers;
		Collections.sort(this.kmers);
		BitSet b_pos = new BitSet(posSeqCount);
		BitSet b_neg = new BitSet(negSeqCount);
 		for (Kmer km:kmers){
 			b_pos.or(km.posBits);
 			b_neg.or(km.negBits);
		}
 		posHitGroupCount = b_pos.cardinality();
		negHitGroupCount = b_neg.cardinality();
	}		
	public KmerGroup(ArrayList<Kmer> kmers, int bs, double[]weights, int posSeqCount, int negSeqCount){
		this.bs = bs;
		this.kmers = kmers;	
		Collections.sort(this.kmers);
		BitSet b_pos = new BitSet(posSeqCount);
		BitSet b_neg = new BitSet(negSeqCount);
 		for (Kmer km:kmers){
 			b_pos.or(km.posBits);
 			b_neg.or(km.negBits);
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
	public int getClusterId(){return clusterId;}
	public void setClusterId(int id){clusterId = id;}
	
	/** hgp (log10) using the positive/negative sequences */
	public double getHgp() {return hgp;	}
	/** hgp (log10) using the positive/negative sequences */
	public void setHgp(double hgp) {this.hgp = hgp;	}
	
	public ArrayList<Kmer> getKmers(){
		return kmers;
	}
	public Kmer getBestKmer(){
		return kmers.get(0);
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
	/** 
	 * get the expected binding position		
	 */
	public int getPosBS(){
		return bs;
	}
	public int compareToByPosHitCount(KmerGroup kg) {		// descending pos hit count
		if(posHitGroupCount>kg.getGroupHitCount()){return(-1);}
		else if(posHitGroupCount<kg.getGroupHitCount()){return(1);}
		else return(0);
	}
	public int compareTo(KmerGroup kg) {					// ascending hgp
		if(hgp<kg.getHgp()){return(-1);}
		else if(hgp>kg.getHgp()){return(1);}
		else return(0);
	}
	public String toString(){
		return String.format("%s %d: %d+/%d-, hpg=%.2f", getBestKmer().getKmerString(), bs, posHitGroupCount, negHitGroupCount, hgp);
	}
}
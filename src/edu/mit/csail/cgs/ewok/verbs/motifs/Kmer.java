package edu.mit.csail.cgs.ewok.verbs.motifs;

import java.util.ArrayList;
import java.util.HashSet;

public class Kmer implements Comparable<Kmer>{
	String kmerString;
	int k;
	int count;	//all hit count
	int seqHitCount; //one hit at most for one sequence, to avoid simple repeat
	ArrayList<KmerHit> hits;
	HashSet<Integer> seqHits;
	double weight;
	double hg;
	int negCount;
	
	public Kmer(String kmerStr, ArrayList<KmerHit> hits){
		this.kmerString = kmerStr;
		this.hits = hits;
		count = hits.size();
		k = kmerStr.length();
		seqHits = new HashSet<Integer>();
		for (KmerHit hit:hits){
			seqHits.add(hit.seqId);
		}
		seqHitCount = seqHits.size();
	}
	// sort kmer by weight
	public int compareByWeight(Kmer o) {
		double diff = o.weight-weight;
		return diff==0?kmerString.compareTo(o.kmerString):(diff<0)?-1:1;  // descending
	}
	
	// sort kmer by HG p-value
	public int compareTo(Kmer o) {
		double diff = o.hg-hg;
		return diff==0?kmerString.compareTo(o.kmerString):(diff<0)?1:-1; // ascending
	}
	public String toString(){
		return kmerString+"\t"+seqHitCount+"\t"+negCount+"\t"+String.format("%.1f", Math.log(hg));
	}
	public static String toHeader(){
		return "kmer\tSeqHitCount\tKmerCt\tWeight";
	}
	public double bias(int pos){
		double sum=0;
		for (KmerHit hit:hits)
			sum+=Math.abs(hit.pos-pos);
		return sum/hits.size();
	}
	public double bias2(int pos){
		double sum=0;
		for (KmerHit hit:hits)
			sum+=(hit.pos-pos)*(hit.pos-pos);
		return Math.sqrt(sum/hits.size());
	}
	public float[] getPositionCounts(int seqLength){
		float counts[] = new float[seqLength];
		for (KmerHit hit:hits)
			counts[hit.pos]++;
		return counts;
	}
	public void calcWeight(double[] seqProbs, double[] positionProbs){
		weight = 0;
		for (KmerHit hit:hits){
			int id = hit.seqId;
			if (seqHits.contains(id)||seqHits.contains(-id))
				weight += seqProbs[id]*positionProbs[hit.pos];
		}
	}
	public int getSeqHitCount() {
		return seqHitCount;
	}
}

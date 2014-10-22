package edu.mit.csail.cgs.deepseq.discovery.kmer;

import java.util.*;

import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;

public class WildcardKmer extends Kmer{
	private ArrayList<Kmer> kmers = new ArrayList<Kmer>();
	private int pos;
	private char[] s;

	public WildcardKmer(Kmer kmer, int pos){
		kmers.add(kmer);
		this.k = kmer.getK();
		this.pos = pos;
		this.s = kmer.getKmerString().toCharArray();
	}
	private WildcardKmer(ArrayList<Kmer> kmers, int pos, char[] s){
		this.kmers = kmers;
		this.pos = pos;
		this.s = s;		
	}
	/** 
	 * add the kmer to the sub-kmers for the wildcard<br>
	 * it is ok to RC() the kmer after adding it because the kmerString is not used anymore<br>
	 * we only care about the pos/neg hits
	 * @param kmer
	 */
	public void addKmer (Kmer kmer){
		kmers.add(kmer);
	}
	public ArrayList<Kmer> getKmers (){
		return kmers;
	}
	public boolean make(){
		if (kmers.size()==1)
			return false;
		s[pos]='N';
		kmerString = String.valueOf(s);
		kmerRC = SequenceUtils.reverseComplement(kmerString);
		for (Kmer km:kmers)
			posHits.addAll(km.getPosHits());
		return true;
	}
	public void setNegHits(){
		for (Kmer km:kmers)
			negHits.addAll(km.getNegHits());
	}
	public int getPosHitCount() {
			return posHits.size();
	}
	
	public WildcardKmer clone(){
		ArrayList<Kmer> newKmers = new ArrayList<Kmer>();
		for (Kmer km:kmers)
			newKmers.add(km.clone());
		WildcardKmer n = new WildcardKmer(newKmers, pos, s.clone());
		n.make();
		n.setNegHits();
		return n;
	}

	public static void printWildcardKmers(ArrayList<WildcardKmer> kmers, int posSeqCount, int negSeqCount, double score, 
			String filePrefix, boolean printShortFormat, boolean print_kmer_hits, boolean printKmersAtK){
		if (kmers==null || kmers.isEmpty())
			return;
		
		Collections.sort(kmers);
		
		StringBuilder sb = new StringBuilder();
		sb.append(String.format("#%d/%d\n", posSeqCount, negSeqCount));
		sb.append(String.format("#%.2f\n", score));
		if (printShortFormat)
			sb.append(Kmer.toShortHeader(kmers.get(0).getK()));
		else
			sb.append(Kmer.toHeader(kmers.get(0).getK()));
		sb.append("\n");
		for (Kmer kmer:kmers){
			if (printShortFormat)
				sb.append(kmer.toShortString()).append("\n");
			else{
				sb.append(kmer.toString());
				if (print_kmer_hits)
					sb.append("\t").append(hits2string(kmer.getPosHits())).append("\t").append(hits2string(kmer.getNegHits()));
				sb.append("\n");
			}
		}
		if (printKmersAtK)
			CommonUtils.writeFile(String.format("%s_kmers_k%d.txt", filePrefix, kmers.get(0).getK()), sb.toString());
		else
			CommonUtils.writeFile(String.format("%s_KSM.txt", filePrefix), sb.toString());
	}

}

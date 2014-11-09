package edu.mit.csail.cgs.deepseq.discovery.kmer;

import java.util.*;

import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;

public class GappedKmer extends Kmer{
	private HashMap<Kmer,Boolean> subKmers = new HashMap<Kmer,Boolean>();
	
	public GappedKmer(String wkString){
		k = wkString.length();
		kmerString = wkString;
		kmerRC = SequenceUtils.reverseComplement(kmerString);
	}

	/** 
	 * add the kmer to the sub-kmers for the gapped kmer<br>
	 * it is ok to RC() the sub-kmers after adding it because their kmerStrings are not used anymore<br>
	 * we only care about the pos/neg hits, and has stored the orientation
	 * @param kmer
	 */
	public void addSubKmer (Kmer kmer, boolean isSameOrientation){
		subKmers.put(kmer, isSameOrientation);
	}
	void linkSubKmers(){
		for (Kmer km:subKmers.keySet())
			km.addGappedKmer(this);
	}
	public Set<Kmer> getSubKmers (){
		return subKmers.keySet();
	}
	boolean getSubKmerOrientation(Kmer subkmer){
		return subKmers.get(subkmer);
	}
	public void mergePosHits(){
		posHits.clear();
		posBits.clear();
		for (Kmer km:subKmers.keySet()){
			posHits.addAll(km.getPosHits());
			posBits.or(km.posBits);
		}
		if (use_weighted_hit_count)
			setWeightedPosHitCount();
	}
	
	public void mergeNegHits(){
		negHits.clear();
		negBits.clear();
		for (Kmer km:subKmers.keySet()){
			negHits.addAll(km.getNegHits());
			negBits.or(km.negBits);
		}
	}
	
	public GappedKmer clone(){
		GappedKmer n = new GappedKmer(kmerString);
		for (Kmer km:subKmers.keySet())
			n.addSubKmer(km.clone(), subKmers.get(km));
		n.mergePosHits();
		n.mergeNegHits();
		return n;
	}
	
	public void addBasicKmersToSet(HashSet<Kmer> reg){
		for (Kmer km:subKmers.keySet())
			reg.add(km);
	}
	
	public static void printGappedKmers(ArrayList<GappedKmer> kmers, int k, int posSeqCount, int negSeqCount, double score, 
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
		for (GappedKmer kmer:kmers){
			if (printShortFormat)
				sb.append(kmer.toShortString()).append("\t").append(kmer.getSubKmers().size()).append("\n");
			else{
				sb.append(kmer.toString());
				if (print_kmer_hits)
					sb.append("\t").append(hits2string(kmer.getPosHits())).append("\t").append(hits2string(kmer.getNegHits()));
				sb.append("\n");
			}
		}
		if (printKmersAtK)
			CommonUtils.writeFile(String.format("%s_kmers_k%d.txt", filePrefix, k), sb.toString());
		else
			CommonUtils.writeFile(String.format("%s_KSM.txt", filePrefix), sb.toString());
	}

	public void removeSubkmers(Kmer km) {
		subKmers.remove(km);		
	}
}

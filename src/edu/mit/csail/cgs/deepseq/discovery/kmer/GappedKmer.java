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
	
	public void update(){
		mergePosHits();
		mergeNegHits();
		linkSubKmers();
	}
	
	public GappedKmer clone(){
		GappedKmer n = new GappedKmer(kmerString);
		for (Kmer km:subKmers.keySet())
			n.addSubKmer(km.clone(), subKmers.get(km));
		n.update();
		return n;
	}
	
	public void addBasicKmersToSet(HashSet<Kmer> reg){
		for (Kmer km:subKmers.keySet())
			reg.add(km);
	}
	
	/**
	 * Print a list of k-mers to a KSM file<br>
	 * It can print both gapped and ungapped kmers. For gapped kmers, the sub-kmers will also be printed.
	 * @param kmers
	 * @param k
	 * @param posSeqCount
	 * @param negSeqCount
	 * @param score
	 * @param filePrefix
	 * @param printShortFormat
	 * @param print_kmer_hits
	 * @param printKmersAtK
	 */
	public static void printGappedKmers(ArrayList<Kmer> kmers, int k, int posSeqCount, int negSeqCount, double score, 
			String filePrefix, boolean printShortFormat, boolean print_kmer_hits, boolean printKmersAtK){
		if (kmers==null || kmers.isEmpty())
			return;
		
		Collections.sort(kmers);
		
		int subKmerId = 0;
		HashMap<Kmer,Integer> allSubKmers = new HashMap<Kmer,Integer>();
		for (Kmer km: kmers){
			if (km instanceof GappedKmer){
				for(Kmer sk: ((GappedKmer) km).getSubKmers()){
					if (!allSubKmers.containsKey(sk)){
						allSubKmers.put(sk,subKmerId);
						subKmerId++;
					}
				}
			}
	    }
		Kmer[] subkmerList = new Kmer[allSubKmers.size()];
		for (Kmer km: allSubKmers.keySet())
			subkmerList[allSubKmers.get(km)] = km;
		
		StringBuilder sb = new StringBuilder();
		sb.append(String.format("#%d/%d\n", posSeqCount, negSeqCount));
		sb.append(String.format("#%.2f\n", score));
		if (printShortFormat){
			sb.append(GappedKmer.toShortHeader(k)).append("\n");
			for (Kmer kmer:kmers){
				if (kmer instanceof GappedKmer)
					sb.append(kmer.toShortString()).append("\t").append(((GappedKmer) kmer).getSubKmers().size()).append("\n");
				else
					sb.append(kmer.toShortString()).append("\n");
			}
		}
		else{
			sb.append(GappedKmer.toHeader(k)).append("\n");
			for (Kmer kmer:kmers){
				sb.append(kmer.toString2());				
				if (kmer instanceof GappedKmer){
					sb.append("\t");
					for (Kmer sk: ((GappedKmer) kmer).getSubKmers())
						sb.append(allSubKmers.get(sk)).append(",");
					sb.deleteCharAt(sb.length()-1);
					if (print_kmer_hits)
						sb.append("\tNA\tNA");	// the real hits are stored in sub-kmers
				}
				else{
					sb.append("\t").append(-1);
					if (print_kmer_hits)
						sb.append("\t").append(hitBits2string(kmer.posBits)).append("\t").append(hitBits2string(kmer.negBits));
				}
				sb.append("\n");
			}
			sb.append("\n");	// an empty line to signal that the following are the sub-kmers
			for (Kmer kmer:subkmerList){
				sb.append(kmer.toString2()).append("\t").append(allSubKmers.get(kmer));
				if (print_kmer_hits)
					sb.append("\t").append(hitBits2string(kmer.posBits)).append("\t").append(hitBits2string(kmer.negBits));
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
	public static String toHeader(int k){
		int length=2*k+1;
		String firstField = "# k-mer/r.c.";
		if (firstField.length()<length)
			firstField += CommonUtils.padding(length-firstField.length(), ' ');
		if (use_weighted_hit_count)
			return firstField+"\tOffset\tPosCt\twPosCt\tNegCt\tHGP_10\tCIDs\tPosHits\tNegHits";
		else
			return firstField+"\tOffset\tPosCt\tNegCt\tHGP_10\tCIDs\tPosHits\tNegHits";
	}

	public static String toShortHeader(int k){
		int length=2*k+1;
		String firstField = "# k-mer/r.c.";
		if (firstField.length()<length)
			firstField += CommonUtils.padding(length-firstField.length(), ' ');
		if (use_weighted_hit_count)
			return firstField+"\tPosCt\twPosCt\tNegCt\tHgp_lg10\tNumChildren";
		else
			return firstField+"\tPosCt\tNegCt\tHgp_lg10\tNumChildren";
	}
	
//	public static Kmer fromString(String str){
//		String[] f = str.trim().split("\t");
//		String[] f0f = f[0].split("/");
//		HashSet<Integer> posHits = new HashSet<Integer>();
//		HashSet<Integer> negHits = new HashSet<Integer>();
//		if (f.length==9){
//			String pos_id_string = f[7].trim();
//			if (!pos_id_string.equals("-1")){
//				String[] pos_ids = pos_id_string.split(" ");
//				for (String hit:pos_ids)
//					posHits.add(Integer.valueOf(hit));
//			}
//
//			String neg_id_string = f[8].trim();
//			if (!neg_id_string.equals("-1")){
//				String[] neg_ids = neg_id_string.split(" ");
//				for (String hit:neg_ids)
//					negHits.add(Integer.valueOf(hit));
//			}
//		}
//		else if (f.length==8){ // no wPos field
//			String pos_id_string = f[6].trim();
//			if (!pos_id_string.equals("-1")){
//				String[] pos_ids = pos_id_string.split(" ");
//				for (String hit:pos_ids)
//					posHits.add(Integer.valueOf(hit));
//			}
//
//			String neg_id_string = f[7].trim();
//			if (!neg_id_string.equals("-1")){
//				String[] neg_ids = neg_id_string.split(" ");
//				for (String hit:neg_ids)
//					negHits.add(Integer.valueOf(hit));
//			}
//		}
//		
//		Kmer kmer = new Kmer(f0f[0], posHits);	
//		kmer.clusterId = Integer.parseInt(f[1]);
//		kmer.kmerStartOffset = Integer.parseInt(f[2]);
//		kmer.setNegHits(negHits);
//		if (f.length==9){
//			kmer.hgp_lg10 = Double.parseDouble(f[6]);
//		}
//		else if (f.length==8){		// no wPos field
//			kmer.hgp_lg10 = Double.parseDouble(f[5]);
//		}
//
//		return kmer;
//	}
}

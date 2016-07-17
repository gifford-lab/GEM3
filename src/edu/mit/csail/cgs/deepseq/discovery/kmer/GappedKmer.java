package edu.mit.csail.cgs.deepseq.discovery.kmer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;

import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;

public class GappedKmer extends Kmer{
	private HashMap<Kmer,Boolean> baseKmers = new HashMap<Kmer,Boolean>();
	
	public GappedKmer(String wkString){
		k = wkString.length();
		kmerString = wkString;
		kmerRC = SequenceUtils.reverseComplement(kmerString);
	}

	/** 
	 * Add the base kmer to the base-kmers for the gapped kmer<br>
	 * SHOULD NOT RC() the base-kmers after adding it because base-kmers can map to multiple gapped k-mers, 
	 * their kmerStrings will be used in KSM output and later KSM scanning<br>
	 * in KMAC, we only care about the pos/neg hits, and has stored the orientation
	 * @param kmer
	 */
	public void addBaseKmer (Kmer kmer, boolean isSameOrientation){
		baseKmers.put(kmer, isSameOrientation);
	}
	void clearBaseKmers(){
		baseKmers.clear();
	}
	public Set<Kmer> getBaseKmers (){
		return baseKmers.keySet();
	}
	boolean getBaseKmerOrientation(Kmer subkmer){
		return baseKmers.get(subkmer);
	}
	void mergePosHits(double[] seq_weights ){
		posBits.clear();
		for (Kmer km:baseKmers.keySet()){
			posBits.or(km.posBits);
		}
		if (use_weighted_hit_count)
			setWeightedPosHitCount(seq_weights );
	}
	
	public void mergeNegHits(){
//		negHits.clear();
		negBits.clear();
		for (Kmer km:baseKmers.keySet()){
//			negHits.addAll(km.getNegHits());
			negBits.or(km.negBits);
		}
	}
	
	public void update(double[] seq_weights ){
		mergePosHits(seq_weights);
		mergeNegHits();
	}
	
	public void setMatrix(){
		matrix = new float[kmerString.length()][KMAC1.LETTERS.length];
		for (Kmer bk: getBaseKmers()){
			String ks = getBaseKmerOrientation(bk)?bk.kmerString:bk.kmerRC;
			for (int ii=0;ii<bk.getK();ii++){
				for (int j=0;j<KMAC1.LETTERS.length;j++){
					if (ks.charAt(ii) == KMAC1.LETTERS[j])
						matrix[ii][j] += bk.getPosHitCount();
				}
			}
		}
		// normalize at each position
		for (int ii=0;ii<matrix.length;ii++){
			float sum=0;
			for (int j=0;j<KMAC1.LETTERS.length;j++)
				sum += matrix[ii][j];
			for (int j=0;j<KMAC1.LETTERS.length;j++)
				matrix[ii][j] /= sum;
		}
		
		setMatrixRC();
	}
	/**
	 * Clone for output only. It does not set up correct GK-SK linkages.
	 */
	public GappedKmer clone(double[] seq_weights){
		GappedKmer n = new GappedKmer(getKmerStr());
		n.clusterId = clusterId;
		n.shift = shift;
		n.setNegBits((BitSet)negBits.clone());
		n.setPosBits((BitSet)posBits.clone());
		n.hgp_lg10 = hgp_lg10;
		n.kmerStartOffset = kmerStartOffset;
		n.isSeedOrientation = isSeedOrientation;
		for (Kmer km:baseKmers.keySet())
			n.addBaseKmer(km, baseKmers.get(km));
		n.update(seq_weights);
		return n;
	}
	
	public void addBaseKmersToSet(HashSet<Kmer> reg){
		for (Kmer km:baseKmers.keySet())
			reg.add(km);
	}
	
	public void removeBasekmers(Kmer km) {
		baseKmers.remove(km);		
	}

	/**
	 * Print a list of k-mers to a KSM file (KMAC1)<br>
	 * It can print both gapped and ungapped kmers. For gapped kmers, the base-kmers will also be printed.
	 * @param kmers
	 * @param kOriginal
	 * @param posSeqCount
	 * @param negSeqCount
	 * @param score
	 * @param filePrefix
	 * @param printShortFormat
	 * @param print_kmer_hits
	 * @param printKmersAtK
	 */
	public static void printKSM(ArrayList<Kmer> kmers, double[] seq_weights, int kOriginal, int gap, int posSeqCount, int negSeqCount, double score, 
			String filePrefix, boolean printShortFormat, boolean print_kmer_hits, boolean printKmersAtK){
		if (kmers==null || kmers.isEmpty())
			return;
		
		int k = kOriginal + gap;
		Collections.sort(kmers);
		
		int baseKmerId = 0;
		HashMap<Kmer,Integer> allBaseKmer2ID = new HashMap<Kmer,Integer>();
		for (Kmer km: kmers){
			if (km instanceof GappedKmer){
				for(Kmer sk: ((GappedKmer) km).getBaseKmers()){
					if (!allBaseKmer2ID.containsKey(sk)){
						allBaseKmer2ID.put(sk,baseKmerId);
						baseKmerId++;
					}
				}
			}
	    }
		Kmer[] basekmerList = new Kmer[allBaseKmer2ID.size()];
		for (Kmer km: allBaseKmer2ID.keySet())
			basekmerList[allBaseKmer2ID.get(km)] = km;
		
		StringBuilder sb = new StringBuilder();
		sb.append(String.format("#%d/%d\n", posSeqCount, negSeqCount));
		sb.append(String.format("#%.2f\n", score));
		if (printShortFormat){
			sb.append(GappedKmer.toShortHeader(k)).append("\n");
			for (Kmer kmer:kmers){
				if (kmer instanceof GappedKmer)
					sb.append(kmer.toShortString()).append("\t").append(((GappedKmer) kmer).getBaseKmers().size()).append("\n");
				else
					sb.append(kmer.toShortString()).append("\n");
			}
		}
		else{		// print KSM
			sb.append(GappedKmer.toHeader(k)).append("\n");
			for (Kmer kmer:kmers){
				sb.append(kmer.toString2());				
				if (kmer instanceof GappedKmer){
					sb.append("\t");
					for (Kmer bk: ((GappedKmer)kmer).getBaseKmers()){
						int id = allBaseKmer2ID.get(bk);
						// label the base kmer orientation w.r.t the seed orientation of the gapped k-mer, be consistent in the output KSM
						sb.append(((GappedKmer)kmer).getBaseKmerOrientation(bk)==kmer.isSeedOrientation ? "" : '-').append(id).append(",");
					}
					sb.deleteCharAt(sb.length()-1);
					if (print_kmer_hits)
						sb.append("\tN.A.\tN.A.");	// the real hits are stored in sub-kmers
				}
				else{
					sb.append("\t").append('*');
					if (print_kmer_hits)
						sb.append("\t").append(hitBits2string(kmer.posBits)).append("\t").append(hitBits2string(kmer.negBits));
				}
				sb.append("\n");
			}
			
			sb.append("$$$\n");	// $$$ to signal that the following are the sub-kmers
			
			for (Kmer kmer:basekmerList){
				sb.append(kmer.toString2()).append("\t").append(allBaseKmer2ID.get(kmer));
				if (print_kmer_hits)
					sb.append("\t").append(hitBits2string(kmer.posBits)).append("\t").append(hitBits2string(kmer.negBits));
				sb.append("\n");
			}
			
			if (seq_weights!=null){
				sb.append("%%%\n");	// %%% to signal that the following are the sequence weights
				for (double w: seq_weights)
					sb.append(String.format("%.4f\n", w));
			}
		}
		
		if (printKmersAtK)
			CommonUtils.writeFile(String.format("%s_kmers_%d+%d.txt", filePrefix, kOriginal, gap), sb.toString());
		else
			CommonUtils.writeFile(String.format("%s.KSM.txt", filePrefix), sb.toString());
	}
	

	public static String toHeader(int k){
		int length=2*k+1;
		String firstField = "# k-mer/r.c.";
		if (firstField.length()<length)
			firstField += CommonUtils.padding(length-firstField.length(), ' ');
		return firstField+"\tOffset\tPosCt\twPosCt\tNegCt\tHgpLg10\tCIDs\tPosHits (base85 encoding)\tNegHits (base85 encoding)";
	}

	public static String toShortHeader(int k){
		int length=2*k+1;
		String firstField = "# k-mer/r.c.";
		if (firstField.length()<length)
			firstField += CommonUtils.padding(length-firstField.length(), ' ');
		return firstField+"\tPosCt\twPosCt\tNegCt\tHgpLg10\tNumChildren";
	}
	/**
	 * Load KSM file and parse into gapped kmers and weights
	 * @param file
	 * @return
	 */
	public static KsmMotif loadKSM(File file){
		KsmMotif ksm = new KsmMotif();
		ArrayList<Kmer> kmers = new ArrayList<Kmer>();
		HashMap<Integer, Kmer> id2baseKmerMap = new HashMap<Integer, Kmer>();
		try {	
			BufferedReader bin = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
			String line = bin.readLine();
	        line = line.substring(1,line.length());			//remove # sign
	        String[] f = line.split("/");
	        ksm.posSeqCount = Integer.parseInt(f[0]);
	        ksm.negSeqCount = Integer.parseInt(f[1]);
		       
	        // threshold
	        line = bin.readLine();
	        line = line.substring(1,line.length());			//remove # sign
	        ksm.cutoff = Double.parseDouble(line);
	        
	        //load gapped k-mers
	        while((line = bin.readLine()) != null) { 
	        	if (line.startsWith("#"))
	        		continue;
	            if (line.equals("") || line.equals("$$$"))	// break at the empty line before the sub-kmers
	            	break;
	            line = line.trim();
	            Kmer kmer = GappedKmer.fromString(line);
	            kmers.add(kmer);
	        }	
	        
	        //load base k-mers
	        while((line = bin.readLine()) != null) { 
	            line = line.trim();
	            if (line.equals("") || line.equals("%%%"))	// break at the empty line between the base-kmers and sequence weights
	            	break;
	            Kmer kmer = GappedKmer.fromString(line);
	            id2baseKmerMap.put(Integer.parseInt(kmer.CIDs), kmer);		// for base-kmer, CIDs field is only one id
	        }	

	        // load sequence weights
	        ArrayList<Double> weights = new ArrayList<Double>();
	        while((line = bin.readLine()) != null) { 
	            line = line.trim();
	            if (line.equals(""))	// The end of file
	            	break;
	            weights.add(Double.parseDouble(line));
	        }
	        if (!weights.isEmpty()){
		        ksm.seq_weights = new double[weights.size()];
		        for (int i=0; i<weights.size(); i++){
		        	ksm.seq_weights[i] = weights.get(i);
		        }
	        }
	        if (bin != null) {
	            bin.close();
	        }
        } catch (IOException e) {
        	System.err.println("I/O Error when processing "+file.getName());
            e.printStackTrace(System.err);
        }
		for (Kmer km:kmers){
			if (km instanceof GappedKmer){
				GappedKmer gk = (GappedKmer) km;
				String[] f = gk.CIDs.split(",");
				for (String id: f){
					int idx = Integer.parseInt(id);
					boolean isSameStrand = true;
					if (idx<0){
						idx = -idx;
						isSameStrand = false;
					}
					else if (idx==0){
						isSameStrand = id.startsWith("0");		// if "-0", not same strand 
					}
					Kmer bk = id2baseKmerMap.get(idx);
					gk.addBaseKmer(bk, isSameStrand);
				}
				gk.update(ksm.seq_weights);
			}
			
		}
        kmers.trimToSize();
        ksm.kmers = kmers;
        
		return ksm;
	}
	

	public static Kmer fromString(String str){
		String[] f = str.trim().split("\t");
		String[] f0f = f[0].split("/");
		String kmerStr = f0f[0];
		Kmer kmer = null;
		if (kmerStr.contains("N")){// Gapped kmer
			kmer = new GappedKmer(kmerStr);
		}
		else{
			BitSet b_pos = BitSet.valueOf(CommonUtils.decodeAscii85StringToBytes(f[7]));
			kmer = new Kmer(kmerStr, b_pos);
			if (f.length>8)
				kmer.negBits = BitSet.valueOf(CommonUtils.decodeAscii85StringToBytes(f[8]));
		}
		kmer.isSeedOrientation = true;
		kmer.CIDs = f[6];
		kmer.kmerStartOffset = Integer.parseInt(f[1]);
		kmer.shift = kmer.kmerStartOffset;
		kmer.hgp_lg10 = Double.parseDouble(f[5]);

		return kmer;
	}
}

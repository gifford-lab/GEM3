package edu.mit.csail.cgs.ewok.verbs.motifs;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSPeak;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;

public class Kmer implements Comparable<Kmer>{
	static cern.jet.random.engine.RandomEngine randomEngine = new cern.jet.random.engine.MersenneTwister();
	
	String kmerString;
	public String getKmerString() {	return kmerString;}
	int k;
	public int getK(){return k;}
	int seqHitCount; //one hit at most for one sequence, to avoid simple repeat
	double strength;	// the total read counts from all events explained by this kmer
	public double getStrength(){return strength;}
	public void setStrength(double strength){this.strength = strength;}
	public void incrStrength(double strength){this.strength += strength;}
	double hgp = 1;
	/**  get hyper-geometric p-value of the kmer */
	public double getHgp() {
		return hgp;
	}
	/**  set hyper-geometric p-value of the kmer */
	public void setHgp(double hgp) {
		this.hgp = hgp;
	}
	int negCount;
	
	Kmer reference;		//
	public Kmer getRef(){return reference;}
	int shift;			// the position relative to seedKmer, after aligning this kmer to seedKmer
	/** get the position relative to seedKmer, after aligning this kmer to seedKmer */
	public int getShift(){return shift;}
	public void setShift(int s){shift=s;}
	int score;			// best possible number of matches (with or w/o shift) wrt reference Kmer
	public int getScore(){return score;}	
	String alignString;
	public void setAlignString(String str){alignString=str;}
	public String getAlignString(){return alignString;}
	
	/**
	 *  The offset of kmer start from the binding position of motif(PWM) (Pos_kmer-Pos_wm)
	 */
	int kmerStartOffset;			
	/**
	 *  Get the offset of kmer start from the binding position of motif(PWM) (Pos_kmer-Pos_wm)
	 */
	public int getKmerStartOffset(){return kmerStartOffset;}
	/**
	 *  Set the offset of kmer start from the binding position of motif(PWM) (Pos_kmer-Pos_wm)
	 */
	public void setKmerStartOffset(int s){kmerStartOffset=s;}
	int group=-1;			// the group of motif
	public int getGroup(){return group;}
	public void setGroup(int g){group=g;}
	
	public Kmer(String kmerStr, int hitCount){
		this.kmerString = kmerStr;
		this.k = kmerString.length();
		this.seqHitCount = hitCount;
//		if (randomEngine.nextDouble()>0.5)
//			this.kmerStartOffset = -(this.k-1-(this.k-1)/2);
//		else
//			this.kmerStartOffset = -(this.k-1)/2;
	}
	
	public Kmer clone(){
		Kmer n = new Kmer(getKmerString(), getSeqHitCount());
		n.strength = strength;
		n.shift = shift;
		n.negCount = negCount;
		n.hgp = hgp;
		n.alignString = alignString;
		n.kmerStartOffset = kmerStartOffset;
		return n;
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

	// sort kmer by strength
	public int compareByStrength(Kmer o) {
		double diff = o.strength-strength;
		return diff==0?kmerString.compareTo(o.kmerString):(diff<0)?-1:1;  // descending
	}
	// sort kmer by hgp
	public int compareByHGP(Kmer o) {
		double diff = o.hgp-hgp;
		return diff==0?this.compareTo(o):(diff<0)?1:-1;  // ascending HGP, descending seqHitCount
	}	
	// default, sort kmer by seqHitCount
	public int compareTo(Kmer o) {
		double diff = o.seqHitCount-seqHitCount;
		return diff==0?kmerString.compareTo(o.kmerString):(diff<0)?-1:1; // descending
	}
	public boolean hasString(String kmerString){
		return this.kmerString.equals(kmerString);
	}
	public String toString(){
		double hg_lg = Math.log10(hgp);
		if (hg_lg==Double.NEGATIVE_INFINITY)
			hg_lg=-400;
		return kmerString+"\t"+seqHitCount+"\t"+negCount+"\t"+String.format("%.1f", hg_lg)+"\t"+String.format("%.1f", strength)+"\t"+kmerStartOffset;
	}
	public String toNonOverlapString(){
		return kmerString+"\t"+seqHitCount+"\t"+String.format("%.1f", strength)+"\t"+kmerStartOffset;
	}
	public static String toHeader(){
		return "EnrichedKmer\tEventCt\tStrengt\tOffset";
	}
	public String toOverlapString(){
		double hg_lg = Math.log10(hgp);
		if (hg_lg==Double.NEGATIVE_INFINITY)
			hg_lg=-400;
		return kmerString+"\t"+seqHitCount+"\t"+negCount+"\t"+String.format("%.1f", hg_lg);
	}
	public static String toOverlapHeader(){
		return "OverlappedKmer\tPosCt\tNegCt\tHGP_10";
	}
	public static Kmer fromString(String str){
		String[] f = str.split("\t");
		Kmer kmer = new Kmer(f[0], Integer.parseInt(f[1]));
//		kmer.negCount = Integer.parseInt(f[2]);
//		kmer.hgp = Math.pow(10, Double.parseDouble(f[3]));
//		kmer.strength = Double.parseDouble(f[4]);
//		kmer.kmerStartOffset = Integer.parseInt(f[5]);
		return kmer;
	}

	public int getSeqHitCount() {
		return seqHitCount;
	}
	public void setSeqHitCount(int count) {
		seqHitCount=count;
	}
	public void incrSeqHitCount() {
		seqHitCount++;
	}
	public void mergeKmer(Kmer newKmer){
		if (kmerString.equals(newKmer.kmerString)){
			seqHitCount += newKmer.seqHitCount;
			strength += newKmer.strength;
		}
	}
//	
//	/** 
//	 * Set the reference Kmer
//	 * Find the best score, strand(RC), shift for this Kmer to align with reference Kmer 
//	 * @param ref
//	 */
//	public void setReference(Kmer ref){
//		reference = ref;
//		byte[] thisBytes = kmerString.getBytes();
//		byte[] refBytes = ref.kmerString.getBytes();
//		score = 0;
//		shift = -99;
//		for (int s=-2;s<=+2;s++){
//			int count = 0;
//			for (int i=-2;i<refBytes.length+2;i++){
//				if (i<0 || i>refBytes.length-1 ||i+s<0 || i+s>refBytes.length-1 )
//					continue;
//				if (refBytes[i]==thisBytes[i+s]){
//					count ++;
//				}
//			}
//			if (count>score){
//				score = count;
//				shift = s;
//			}
//		}
//		// try RC
//		byte[] rcBytes = getKmerRC().getBytes();
//		boolean useRC=false;
//		for (int s=-2;s<=+2;s++){
//			int count = 0;
//			for (int i=-2;i<refBytes.length+2;i++){
//				if (i<0 || i>refBytes.length-1 ||i+s<0 || i+s>refBytes.length-1 )
//					continue;
//				if (refBytes[i]==rcBytes[i+s]){
//					count ++;
//				}
//			}
//			if (count>score){
//				score = count;
//				shift = s;
//				useRC = true;
//			}
//		}
//		if (useRC)
//			RC();
//	}
//	
//	/**
//	 * Extend this kmer from reference kmer, if they are only offset base off
//	 * @param ref
//	 */
//	public boolean extendKmer(Kmer ref, int offset){
//		reference = ref;
//		byte[] thisBytes = kmerString.getBytes();
//		byte[] refBytes = ref.kmerString.getBytes();
//		score = 0;
//		shift = -99;
//		for (int s=-offset;s<=offset;s++){
//			for (int i=-offset;i<refBytes.length+offset;i++){
//				if (i<0 || i>refBytes.length-1 ||i+s<0 || i+s>refBytes.length-1 )
//					continue;
//				if (refBytes[i]!=thisBytes[i+s])	// if mismatch
//					break;
//			}
//			score = refBytes.length-Math.abs(s);
//			shift = s;
//		}
//		// try RC
//		byte[] rcBytes = getKmerRC().getBytes();
//		boolean useRC=false;
//		for (int s=-offset;s<=offset;s++){
//			for (int i=-offset;i<refBytes.length+offset;i++){
//				if (i<0 || i>refBytes.length-1 ||i+s<0 || i+s>refBytes.length-1 )
//					continue;
//				if (refBytes[i]!=rcBytes[i+s])	// if mismatch
//					break;
//			}
//			int thisScore = refBytes.length-Math.abs(s);
//			if (thisScore>score){
//				score = thisScore;
//				useRC = true;
//			}
//			shift = s;
//		}
//		
//		if (useRC)
//			RC();
//		
//		return score>0;
//	}
//	
	/**
	 * calculate the best shift for input kmer to align with this kmer
	 * allow for 2 mismatches, or 1 shift + 1 mismatch, or 2 shift
	 * @param kmerMatches
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
	
	public static void printKmers(ArrayList<Kmer> kmers, String filePrefix, boolean isOverlappedKmer){
		if (kmers==null || kmers.isEmpty())
			return;
		
		Collections.sort(kmers);
		
		StringBuilder sb = new StringBuilder();
		if (isOverlappedKmer)
			sb.append(Kmer.toOverlapHeader());
		else
			sb.append(Kmer.toHeader());
		sb.append("\n");
		for (Kmer kmer:kmers){
			if (isOverlappedKmer)
				sb.append(kmer.toOverlapString()).append("\n");
			else
				sb.append(kmer.toNonOverlapString()).append("\n");
		}
		CommonUtils.writeFile(String.format("%s_kmer_%d.txt",filePrefix, kmers.get(0).getK()), sb.toString());
	}
	
	public static ArrayList<Kmer> loadKmers(List<File> files){
		ArrayList<Kmer> kmers = new ArrayList<Kmer>();
		for (File file: files){
			try {	
				BufferedReader bin = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
				
		        String line;
		        bin.readLine();	// skip header
		        while((line = bin.readLine()) != null) { 
		            line = line.trim();
		            Kmer kmer = Kmer.fromString(line);
		            kmers.add(kmer);
		        }			
		        if (bin != null) {
		            bin.close();
		        }
	        } catch (IOException e) {
	        	System.err.println("Error when processing "+file.getName());
	            e.printStackTrace(System.err);
	        }
		}
		return kmers;
	}
	
}

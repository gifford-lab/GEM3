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
	int negCount;
	
	double strength;	// the total read counts from all events explained by this kmer
	public double getStrength(){return strength;}
	public void setStrength(double strength){this.strength = strength;}
	public void incrStrength(double strength){this.strength += strength;}
	
	double hgp_lg10 = 0;
	/**  get hyper-geometric p-value (log10) of the kmer */
	public double getHgp() {
		return hgp_lg10;
	}
	/**  set hyper-geometric p-value (log10) of the kmer */
	public void setHgp(double hgp) {
		this.hgp_lg10 = hgp;
	}
	
	int shift;			// the position relative to seedKmer, after aligning this kmer to seedKmer
	/** get the position relative to seedKmer, after aligning this kmer to seedKmer */
	public int getShift(){return shift;}
	public void setShift(int s){shift=s;}

	String alignString;
	public void setAlignString(String str){alignString=str;}
	public String getAlignString(){return alignString;}
	
	/** The offset of kmer start from the binding position of motif(PWM) (Pos_kmer-Pos_wm)*/
	int kmerStartOffset;			
	/** Get the offset of kmer start from the binding position of motif(PWM) (Pos_kmer-Pos_wm)*/
	public int getKmerStartOffset(){return kmerStartOffset;}
	/** Set the offset of kmer start from the binding position of motif(PWM) (Pos_kmer-Pos_wm)*/
	public void setKmerStartOffset(int s){kmerStartOffset=s;}
	
	public Kmer(String kmerStr, int hitCount){
		this.kmerString = kmerStr;
		this.k = kmerString.length();
		this.seqHitCount = hitCount;
	}
	
	public Kmer clone(){
		Kmer n = new Kmer(getKmerString(), getSeqHitCount());
		n.strength = strength;
		n.shift = shift;
		n.negCount = negCount;
		n.hgp_lg10 = hgp_lg10;
		n.alignString = alignString;
		n.kmerStartOffset = kmerStartOffset;
		return n;
	}
	
	/** Use reverse compliment to represent the kmer */
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
		double diff = o.hgp_lg10-hgp_lg10;
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
		return String.format("%s\t%d\t%d\t%d\t%.1f\t%.1f",kmerString, kmerStartOffset, seqHitCount, negCount, hgp_lg10, strength);
	}

	public static String toHeader(){
		return "EnrichedKmer\tOffset\tPosCt\tNegCt\tHGP_10\tStrength";
	}
	public String toShortString(){
		return kmerString+"\t"+seqHitCount+"\t"+negCount+"\t"+String.format("%.1f", hgp_lg10);
	}
	public static String toShortHeader(){
		return "OverlappedKmer\tPosCt\tNegCt\tHGP_10";
	}
	public static Kmer fromString(String str){
		String[] f = str.split("\t");
		Kmer kmer = new Kmer(f[0], Integer.parseInt(f[2]));
		kmer.kmerStartOffset = Integer.parseInt(f[1]);
		kmer.negCount = Integer.parseInt(f[3]);
		kmer.hgp_lg10 = Double.parseDouble(f[4]);
		kmer.strength = Double.parseDouble(f[5]);

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
	
	public String getKmerRC(){
		return SequenceUtils.reverseComplement(kmerString);
	}
	
	public static void printKmers(ArrayList<Kmer> kmers, String filePrefix, boolean printShortFormat){
		if (kmers==null || kmers.isEmpty())
			return;
		
		Collections.sort(kmers);
		
		StringBuilder sb = new StringBuilder();
		if (printShortFormat)
			sb.append(Kmer.toShortHeader());
		else
			sb.append(Kmer.toHeader());
		sb.append("\n");
		for (Kmer kmer:kmers){
			if (printShortFormat)
				sb.append(kmer.toShortString()).append("\n");
			else
				sb.append(kmer.toString()).append("\n");
		}
		CommonUtils.writeFile(String.format("%s_kmer_k%d.txt",filePrefix, kmers.get(0).getK()), sb.toString());
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

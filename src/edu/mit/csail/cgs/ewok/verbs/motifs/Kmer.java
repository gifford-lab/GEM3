package edu.mit.csail.cgs.ewok.verbs.motifs;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.TreeSet;

import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;

public class Kmer implements Comparable<Kmer>{
	static cern.jet.random.engine.RandomEngine randomEngine = new cern.jet.random.engine.MersenneTwister();
	
	String kmerString;
	public String getKmerString() {	return kmerString;}
	int k;
	public int getK(){return k;}
	
	HashSet<Integer> posHits = new HashSet<Integer>();
//	int posHitCount; //one hit at most for one sequence, to avoid simple repeat
	public int getPosHitCount() {return posHits.size();}
	public void setPosHits(HashSet<Integer> posHits) {
		this.posHits = posHits;
//		posHitCount = posHits.size();
	}
	public HashSet<Integer> getPosHits(){return posHits;}
	
//	int negHitCount;
	HashSet<Integer> negHits = new HashSet<Integer>();
	public int getNegHitCount() {return negHits.size();}
	public void setNegHits(HashSet<Integer> negHits) {
		this.negHits = negHits;
//		negHitCount = negHits.size();
	}
	public HashSet<Integer> getNegHits(){return negHits;}
	
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
	
	public Kmer(String kmerStr, HashSet<Integer> posHits ){
		this.kmerString = kmerStr;
		this.k = kmerString.length();
		setPosHits(posHits);
	}
	
	public Kmer clone(){
		Kmer n = new Kmer(getKmerString(), (HashSet<Integer>)posHits.clone());
		n.strength = strength;
		n.shift = shift;
		n.negHits = (HashSet<Integer>)negHits.clone();
		n.hgp_lg10 = hgp_lg10;
		n.alignString = alignString;
		n.kmerStartOffset = kmerStartOffset;
		return n;
	}
	
	/** Use reverse compliment to represent the kmer */
	public void RC(){
		kmerString = getKmerRC();
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
		double diff = o.getPosHitCount()-getPosHitCount();
		return diff==0?kmerString.compareTo(o.kmerString):(diff<0)?-1:1; // descending
	}
	public boolean hasString(String kmerString){
		return this.kmerString.equals(kmerString);
	}
	public String toString(){
		return String.format("%s/%s\t%d\t%d\t%d\t%.1f\t%.1f",kmerString,getKmerRC(),kmerStartOffset, getPosHitCount(), getNegHitCount(), hgp_lg10, strength);
	}

	public static String toHeader(){
		return "Enriched k-mer/r.c.\tOffset\tPosCt\tNegCt\tHGP_10\tStrength";
	}
	public String toShortString(){
		return kmerString+"/"+getKmerRC()+"\t"+getPosHitCount()+"\t"+getNegHitCount()+"\t"+String.format("%.1f", hgp_lg10);
	}
	public static String toShortHeader(){
		return "Enriched k-mer/r.c.\tPosCt\tNegCt\tHGP_10";
	}
	public static Kmer fromString(String str){
		String[] f = str.split("\t");
		String[] f0f = f[0].split("/");
		HashSet<Integer> posHits = new HashSet<Integer>();
		HashSet<Integer> negHits = new HashSet<Integer>();
		if (f.length==8){
			String f6 = f[6].trim();
			if (!f6.equals("")){
				String[] f6f = f6.split(" ");
				for (String hit:f6f)
					posHits.add(Integer.valueOf(hit));
			}
			String f7 = f[7].trim();
			if (!f7.equals("")){
				String[] f7f = f7.split(" ");
				for (String hit:f7f)
					negHits.add(Integer.valueOf(hit));
			}
		}
		Kmer kmer = new Kmer(f0f[0], posHits);	
		kmer.kmerStartOffset = Integer.parseInt(f[1]);
		kmer.hgp_lg10 = Double.parseDouble(f[4]);
		kmer.strength = Double.parseDouble(f[5]);
		kmer.setNegHits(negHits);

		return kmer;
	}

//	public void incrSeqHitCount() {
//		posHitCount++;
//	}
	public void mergeKmer(Kmer newKmer){
		if (kmerString.equals(newKmer.kmerString)){
			posHits.addAll(newKmer.posHits);
			negHits.addAll(newKmer.negHits);
//			posHitCount = posHits.size();
//			negHitCount = negHits.size();
			strength += newKmer.strength;
		}
	}
	
	public String getKmerRC(){
		return SequenceUtils.reverseComplement(kmerString);
	}
	
	public static void printKmers(ArrayList<Kmer> kmers, String filePrefix, boolean printShortFormat, boolean print_kmer_hits){
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
			else{
				sb.append(kmer.toString());
				if (print_kmer_hits)
					sb.append("\t").append(hits2string(kmer.getPosHits())).append("\t").append(hits2string(kmer.getNegHits()));
				sb.append("\n");
			}
		}
		CommonUtils.writeFile(String.format("%s_kmer_k%d.txt",filePrefix, kmers.get(0).getK()), sb.toString());
	}
	private static String hits2string(HashSet<Integer> ids){
		StringBuilder sb = new StringBuilder();
		TreeSet<Integer> sorted = new TreeSet<Integer>(ids);
		for (int id:sorted)
			sb.append(id).append(" ");
		return sb.toString();
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
	public static ArrayList<Kmer> loadKmers(File file){
		ArrayList<Kmer> kmers = new ArrayList<Kmer>();
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
		return kmers;
	}
	public static void main(String[] args){
		ArrayList<Kmer> kmers = Kmer.loadKmers(new File(args[0]));
	}
}

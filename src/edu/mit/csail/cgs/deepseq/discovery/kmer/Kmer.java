package edu.mit.csail.cgs.deepseq.discovery.kmer;

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
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public class Kmer implements Comparable<Kmer>{
	static cern.jet.random.engine.RandomEngine randomEngine = new cern.jet.random.engine.MersenneTwister();
	static boolean use_weighted_hit_count = false;
	static void set_use_weighted_hit_count(boolean weighted){
		use_weighted_hit_count = weighted;
	}
	static double[] seq_weights;
	static void set_seq_weights(double[] weights){
		seq_weights = weights;
	}	
	private String kmerString;
	public String getKmerString() {	return kmerString;}
	private String kmerRC;
	public String getKmerRC(){
		if (kmerRC==null)
			kmerRC=SequenceUtils.reverseComplement(kmerString);
		return kmerRC;
	}
	private int k;
	public int getK(){return k;}
	
	private HashSet<Integer> posHits = new HashSet<Integer>();
	/**	get posHitCount; one hit at most for one sequence, to avoid simple repeat<br>
	 * 	get a weighted version of hit count if use_weighted_hit_count is true
	 */
	public int getPosHitCount() {
		if (use_weighted_hit_count)
			return weightedPosHitCount;
		else
			return posHits.size();
	}
	public void setPosHits(HashSet<Integer> posHits) {
		this.posHits = posHits;
		if (use_weighted_hit_count)
			setWeightedPosHitCount();
		if (posHits.isEmpty())
			return;
		double[] ids = new double[posHits.size()];
		int count=0;
		for (int id:posHits){
			ids[count++]=id;
		}
		double mean = StatUtil.mean(ids);
		double median = StatUtil.median(ids);
		setTop(median-mean);
	}
	public HashSet<Integer> getPosHits(){return posHits;}
	
	private int weightedPosHitCount;
	public void setWeightedPosHitCount(){
		double weight=0;
		for (int i: posHits)
			weight+=seq_weights[i];
		weightedPosHitCount = (int)weight;
	}
//	public int getWeightedHitCount(){
//		return weightedPosHitCount;
//	}
	public double familyHgp;
	
	private double top;
	public void setTop(double top) {this.top = top;}
	public double getTop() {return top;	}
	
//	int negHitCount;
	private HashSet<Integer> negHits = new HashSet<Integer>();
	public int getNegHitCount() {return negHits.size();}
	public void setNegHits(HashSet<Integer> negHits) {
		this.negHits = negHits;
//		negHitCount = negHits.size();
	}
	public HashSet<Integer> getNegHits(){return negHits;}
	
	private double strength;	// the total read counts from all events explained by this kmer
	public double getStrength(){return strength;}
	public void setStrength(double strength){this.strength = strength;}
	public void incrStrength(double strength){this.strength += strength;}
	
	private double hgp_lg10 = 0;
	/**  get hyper-geometric p-value (log10) of the kmer */
	public double getHgp() {
		return hgp_lg10;
	}
	/**  set hyper-geometric p-value (log10) of the kmer */
	public void setHgp(double hgp) {
		this.hgp_lg10 = hgp;
	}
	
	private int shift;
	/** get the position relative to seedKmer, after aligning this kmer to seedKmer */
	public int getShift(){return shift;}
	public void setShift(int s){shift=s;}
	
	private int clusterId=-1;
	public int getClusterId(){return clusterId;}
	public void setClusterId(int id){clusterId=id;}

	private String alignString;
	public void setAlignString(String str){alignString=str;}
	public String getAlignString(){return alignString;}
	
	/** The offset of kmer start from the binding position of motif(PWM) (Pos_kmer-Pos_wm)*/
	private int kmerStartOffset;			
	/** Get the offset of kmer start from the binding position of motif(PWM) (Pos_kmer-Pos_wm)*/
	public int getKmerStartOffset(){return kmerStartOffset;}
	/** Set the offset of kmer start from the binding position of motif(PWM) (Pos_kmer-Pos_wm)*/
	public void setKmerStartOffset(int s){kmerStartOffset=s;}
	
	public Kmer(String kmerStr, Integer posHit ){
		this.kmerString = kmerStr;
		this.k = kmerString.length();
		posHits.add(posHit);
	}
	
	public Kmer(String kmerStr, HashSet<Integer> posHits ){
		this.kmerString = kmerStr;
		this.k = kmerString.length();
		setPosHits(posHits);
	}
	
	public Kmer clone(){
		HashSet<Integer> hits = new HashSet<Integer>();
		hits.addAll(this.posHits);
		Kmer n = new Kmer(getKmerString(), hits);
		n.clusterId = clusterId;
		n.strength = strength;
		n.shift = shift;
		n.negHits = new HashSet<Integer>();
		n.negHits.addAll(negHits);
		n.hgp_lg10 = hgp_lg10;
		n.alignString = alignString;
		n.kmerStartOffset = kmerStartOffset;
		return n;
	}
	
	/** Use reverse compliment to represent the kmer */
	public void RC(){
		String tmp = kmerString;
		kmerString = getKmerRC();
		kmerRC = tmp;
	}
	// sort kmer by strength
	public int compareByStrength(Kmer o) {
		double diff = o.strength-strength;
		return diff==0?kmerString.compareTo(o.kmerString):(diff<0)?-1:1;  // descending
	}
	// sort kmer by family kmer-group hgp
	public int compareByFamilyHGP(Kmer o) {
		double diff = o.familyHgp-familyHgp;
		return diff==0?this.compareTo(o):(diff<0)?1:-1;  // ascending HGP, descending seqHitCount
	}	
	// sort kmer by hgp
	public int compareByHGP(Kmer o) {
		double diff = o.hgp_lg10-hgp_lg10;
		return diff==0?this.compareTo(o):(diff<0)?1:-1;  // ascending HGP, descending seqHitCount
	}	
	// default, sort kmer by posHitCount
	public int compareTo(Kmer o) {
		double diff = o.getPosHitCount()-getPosHitCount();
		return diff==0?kmerString.compareTo(o.kmerString):(diff<0)?-1:1; // descending
	}
	// sort kmer by cluster ID then by HGP
	public int compareByClusterAndHGP(Kmer o) {
		double diff = o.clusterId-clusterId;
		return diff==0?this.compareByHGP(o):(diff<0)?1:-1;  // ascending clusterID, ascending HGP
	}	
	public boolean hasString(String kmerString){
		return this.kmerString.equals(kmerString);
	}
	public String toString(){
		if (use_weighted_hit_count)
			return String.format("%s/%s\t%d\t%d\t%d\t%d\t%d\t%.1f", 
				kmerString, getKmerRC(),clusterId, kmerStartOffset, posHits.size(), weightedPosHitCount, getNegHitCount(), hgp_lg10);
		else
			return String.format("%s/%s\t%d\t%d\t%d\t%d\t%.1f", 
				kmerString, getKmerRC(),clusterId, kmerStartOffset, posHits.size(), getNegHitCount(), hgp_lg10);
	}
	public static String toHeader(int k){
		int length=2*k+1;
		String firstField = "# k-mer/r.c.";
		if (firstField.length()<length)
			firstField += CommonUtils.padding(length-firstField.length(), ' ');
		if (use_weighted_hit_count)
			return firstField+"\tCluster\tOffset\tPosCt\twPosCt\tNegCt\tHGP_10";
		else
			return firstField+"\tCluster\tOffset\tPosCt\tNegCt\tHGP_10";
	}
	public String toShortString(){
		if (use_weighted_hit_count)
			return String.format("%s/%s\t%d\t%d\t%d\t%.1f", 
				kmerString, getKmerRC(), posHits.size(), weightedPosHitCount, getNegHitCount(), hgp_lg10);
		else
			return String.format("%s/%s\t%d\t%d\t%.1f", 
				kmerString, getKmerRC(), posHits.size(), getNegHitCount(), hgp_lg10);
	}
	public static String toShortHeader(int k){
		int length=2*k+1;
		String firstField = "# k-mer/r.c.";
		if (firstField.length()<length)
			firstField += CommonUtils.padding(length-firstField.length(), ' ');
		if (use_weighted_hit_count)
			return firstField+"\tPosCt\twPosCt\tNegCt\tHGP_10";
		else
			return firstField+"\tPosCt\tNegCt\tHGP_10";
	}
	public static Kmer fromString(String str){
		String[] f = str.trim().split("\t");
		String[] f0f = f[0].split("/");
		HashSet<Integer> posHits = new HashSet<Integer>();
		HashSet<Integer> negHits = new HashSet<Integer>();
		if (f.length==9){
			String pos_id_string = f[7].trim();
			if (!pos_id_string.equals("-1")){
				String[] pos_ids = pos_id_string.split(" ");
				for (String hit:pos_ids)
					posHits.add(Integer.valueOf(hit));
			}

			String neg_id_string = f[8].trim();
			if (!neg_id_string.equals("-1")){
				String[] neg_ids = neg_id_string.split(" ");
				for (String hit:neg_ids)
					negHits.add(Integer.valueOf(hit));
			}
		}
		else if (f.length==8){ // no wPos field
			String pos_id_string = f[6].trim();
			if (!pos_id_string.equals("-1")){
				String[] pos_ids = pos_id_string.split(" ");
				for (String hit:pos_ids)
					posHits.add(Integer.valueOf(hit));
			}

			String neg_id_string = f[7].trim();
			if (!neg_id_string.equals("-1")){
				String[] neg_ids = neg_id_string.split(" ");
				for (String hit:neg_ids)
					negHits.add(Integer.valueOf(hit));
			}
		}
		
		Kmer kmer = new Kmer(f0f[0], posHits);	
		kmer.clusterId = Integer.parseInt(f[1]);
		kmer.kmerStartOffset = Integer.parseInt(f[2]);
		kmer.setNegHits(negHits);
		if (f.length==9){
			kmer.hgp_lg10 = Double.parseDouble(f[6]);
		}
		else if (f.length==8){		// no wPos field
			kmer.hgp_lg10 = Double.parseDouble(f[5]);
		}

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

	
	public static void printKmers(ArrayList<Kmer> kmers, int posSeqCount, int negSeqCount, double score, 
			String filePrefix, boolean printShortFormat, boolean print_kmer_hits){
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
		CommonUtils.writeFile(String.format("%s_kmer_k%d.txt",filePrefix, kmers.get(0).getK()), sb.toString());
	}
	
	/**
	 * If the set of ids is empty (no negative hits), return "-1"
	 */
	private static String hits2string(HashSet<Integer> ids){
		StringBuilder sb = new StringBuilder();
		if (ids.isEmpty())
			return "-1";
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
			bin.readLine();			// skip header line, to be compatible with old file format
	        String line;
	        while((line = bin.readLine()) != null) { 
	        	if (line.startsWith("#"))
	        		continue;
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
        kmers.trimToSize();
		return kmers;
	}
	public static Pair<Integer, Integer> getTotalCounts(File file){
		int posSeqCount=-1;
		int negSeqCount=-1;
		try {	
			BufferedReader bin = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	        String line = bin.readLine();
	        line = line.substring(1,line.length());			//remove # sign
	        String[] f = line.split("/");
	        posSeqCount = Integer.parseInt(f[0]);
	        negSeqCount = Integer.parseInt(f[1]);
	        if (bin != null) {
	            bin.close();
	        }
        } catch (IOException e) {
        	System.err.println("Error when processing "+file.getName());
            e.printStackTrace(System.err);
        }
		return new Pair<Integer, Integer>(posSeqCount,negSeqCount);
	}
	
	public static void main(String[] args){
		ArrayList<Kmer> kmers = Kmer.loadKmers(new File(args[0]));
	}
}

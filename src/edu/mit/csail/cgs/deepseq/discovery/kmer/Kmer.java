package edu.mit.csail.cgs.deepseq.discovery.kmer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
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
	String CIDs = null;

	protected String kmerString;
	public String getKmerString() {	return kmerString;}
	public String getKmerStrRC() {	return kmerString+"/"+kmerRC;}
	public String getKmerRCStr() {	return kmerRC+"/"+kmerString;}
	protected String kmerRC;
	public String getKmerRC(){
		if (kmerRC==null)
			kmerRC=SequenceUtils.reverseComplement(kmerString);
		return kmerRC;
	}
	void setKmerString(String kmString){
		kmerString = kmString;
		kmerRC=SequenceUtils.reverseComplement(kmerString);
	}
	protected int k;
	public int getK(){return k;}
	private HashSet<GappedKmer> gappedKmers = null;
	void addGappedKmer(GappedKmer wk){
		if (gappedKmers == null)
			gappedKmers = new HashSet<GappedKmer>();
		gappedKmers.add(wk);	
	}
	void clearGappedKmers(){
		gappedKmers.clear();	
	}
	HashSet<GappedKmer> getGappedKmers(){
		return gappedKmers;
	}
	protected HashSet<Integer> posHits = new HashSet<Integer>();
	BitSet posBits = new BitSet();
	/**	get posHitCount; one hit at most for one sequence, to avoid simple repeat<br>
	 * 	get a weighted version of hit count if use_weighted_hit_count is true
	 */
	public int getPosHitCount() {
		if (use_weighted_hit_count)
			return weightedPosHitCount;
		else
			return posBits.cardinality();
	}
	public void setPosHits(HashSet<Integer> posHits) {
		this.posHits = posHits;
		if (posHits.isEmpty())
			return;
		
		posBits.clear();
		for (int id:posHits){
			posBits.set(id);
		}
		if (use_weighted_hit_count)
			setWeightedPosHitCount();
	}
	public HashSet<Integer> getPosHits(){return posHits;}
	
	private int weightedPosHitCount;
	protected void setWeightedPosHitCount(){
		double weight=0;
		if (seq_weights!=null){
			for (int i = posBits.nextSetBit(0); i >= 0; i = posBits.nextSetBit(i+1)) {
				weight+=seq_weights[i];
	 		}
			weightedPosHitCount = (int)weight;
		}
		else
			weightedPosHitCount = posBits.cardinality();
	}
	public int getWeightedHitCount(){
		return weightedPosHitCount;
	}
	public double familyHgp;
	
//	int negHitCount;
	protected HashSet<Integer> negHits = new HashSet<Integer>();
	BitSet negBits = new BitSet();
	public int getNegHitCount() {return negBits.cardinality();}
	public void setNegHits(HashSet<Integer> negHits) {
		this.negHits = negHits;
//		negHitCount = negHits.size();
		for (int id:negHits){
			negBits.set(id);
		}
	}
	public HashSet<Integer> getNegHits(){return negHits;}
	
	public int getNetHitCount(double posNegSeqRatio) {
		return getPosHitCount()-(int)(getNegHitCount()*posNegSeqRatio);
	}
	private double strength;	// the total read counts from all events explained by this kmer
	public double getStrength(){return strength;}
	public void setStrength(double strength){this.strength = strength;}
	public void incrStrength(double strength){this.strength += strength;}
	
	protected double hgp_lg10 = 0;
	/**  get hyper-geometric p-value (log10) of the kmer, this is usually a negative value */
	public double getHgp() {
		return hgp_lg10;
	}
	/**  set hyper-geometric p-value (log10) of the kmer */
	public void setHgp(double hgp) {
		this.hgp_lg10 = hgp;
	}
	
	protected int shift;
	/** get the position relative to seedKmer, after aligning this kmer to seedKmer<br>
	 * if isSeedOrientation=false, the position is kmerRC relative to seedKmer */
	public int getShift(){return shift;}
	public void setShift(int s){shift=s;}
	protected boolean isSeedOrientation=false;
	/** get whether this kmer is in the same orientation as the seedKmer, after aligning this kmer to seedKmer */
	public boolean isSeedOrientation(){return isSeedOrientation;}
	public void setSeedOrientation(boolean so){isSeedOrientation=so;}
	
	protected int clusterId=-1;
	public int getClusterId(){return clusterId;}
	public void setClusterId(int id){clusterId=id;}

	private String alignString;
	public void setAlignString(String str){alignString=str;}
	public String getAlignString(){return alignString;}
	
	/** The offset of kmer start from the binding position of motif(PWM) (Pos_kmer-Pos_wm)*/
	protected int kmerStartOffset;			
	/** Get the offset of kmer start from the binding position of motif(PWM) (Pos_kmer-Pos_wm)*/
	public int getKmerStartOffset(){return kmerStartOffset;}
	/** Set the offset of kmer start from the binding position of motif(PWM) (Pos_kmer-Pos_wm)*/
	public void setKmerStartOffset(int s){kmerStartOffset=s;}
	public Kmer(){}
	public Kmer(String kmerStr){
		this.kmerString = kmerStr;
		this.k = kmerString.length();
	}
	
	public Kmer(String kmerStr, HashSet<Integer> posHits ){
		this.kmerString = kmerStr;
		this.kmerRC = SequenceUtils.reverseComplement(kmerStr);
		this.k = kmerString.length();
		setPosHits(posHits);
	}
	public Kmer(String kmerStr, BitSet posBits ){
		this.kmerString = kmerStr;
		this.kmerRC = SequenceUtils.reverseComplement(kmerStr);
		this.k = kmerString.length();
		this.posBits = posBits;
	}
	
	public Kmer clone(){
		Kmer n = new Kmer(getKmerString(), (BitSet)(posBits.clone()));
		n.setWeightedPosHitCount();
		n.clusterId = clusterId;
		n.strength = strength;
		n.shift = shift;
		n.setNegBits((BitSet)(negBits.clone()));
		n.hgp_lg10 = hgp_lg10;
		n.alignString = alignString;
		n.kmerStartOffset = kmerStartOffset;
		if (gappedKmers!=null){
			n.gappedKmers = new HashSet<GappedKmer>();
			n.gappedKmers.addAll(gappedKmers);
		}
		n.isSeedOrientation = isSeedOrientation;
		return n;
	}
	
	protected void setNegBits(BitSet bitSet) {
		negBits = bitSet;
	}
	protected void setPosBits(BitSet bitSet) {
		posBits = bitSet;
	}
	/** Add this kmer into the register set*/
	public void addBasicKmersToSet(HashSet<Kmer> reg){
		reg.add(this);
	}
//	/** Use reverse compliment to represent the kmer */
//	public Kmer RC(){
//		String tmp = kmerString;
//		kmerString = getKmerRC();
//		kmerRC = tmp;
//		return this;
//	}
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
	/** sort kmer by ascending HGP, descending seqHitCount */
	public int compareByHGP(Kmer o) {
		double diff = o.hgp_lg10-hgp_lg10;
		return diff==0?this.compareTo(o):(diff<0)?1:-1;  // ascending HGP, descending seqHitCount
	}	
	// default, sort kmer by posHitCount, then by clusterID, then by sequence
	public int compareTo(Kmer o) {
		double diff = o.getPosHitCount()-getPosHitCount();
		double diff_id = o.clusterId-clusterId;
		return diff==0 ? 
				(diff_id==0? kmerString.compareTo(o.kmerString) : (diff_id<0)?1:-1) : 	// ascending clusterID
				((diff<0)?-1:1); // descending posHitCount
//		return diff==0 ? (diff_id==0? 0 : (diff<0)?1:-1): ((diff<0)?-1:1); // descending posHitCount
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
			return String.format("%s\t%d\t%d\t%d\t%d\t%d\t%.1f", 
				getKmerStrRC(),clusterId, shift, posBits.cardinality(), weightedPosHitCount, getNegHitCount(), hgp_lg10);
				// use posBits.cardinality() instead of getPosHitCount() to get raw pos hit count
		else
			return String.format("%s\t%d\t%d\t%d\t%d\t%.1f", 
				getKmerStrRC(),clusterId, shift, posBits.cardinality(), getNegHitCount(), hgp_lg10);
	}
	/**
	 * This method is used in the gapped k-mer printing. 
	 * <br>In KMAC1, a k-mer may be used in multiple KSM, therefore, the k-mer object store isSeedOrientation flag to make it in consistent orientation with the seed k-mer.
	 * @return
	 */
	public String toString2(){
		if (use_weighted_hit_count)
			return String.format("%s\t%d\t%d\t%d\t%d\t%.1f", 
				isSeedOrientation?getKmerStrRC():getKmerRCStr(), shift, posBits.cardinality(), weightedPosHitCount, getNegHitCount(), hgp_lg10);
				// use posBits.cardinality() instead of getPosHitCount() to get raw pos hit count
		else
			return String.format("%s\t%d\t%d\t%d\t%.1f", 
				isSeedOrientation?getKmerStrRC():getKmerRCStr(), shift, posBits.cardinality(), getNegHitCount(), hgp_lg10);
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
			return String.format("%s\t%d\t%d\t%d\t%.1f", 
					getKmerStrRC(), posBits.cardinality(), weightedPosHitCount, getNegHitCount(), hgp_lg10);
		else
			return String.format("%s\t%d\t%d\t%.1f", 
					getKmerStrRC(), posBits.cardinality(), getNegHitCount(), hgp_lg10);
	}
	public static String toShortHeader(int k){
		int length=2*k+1;
		String firstField = "# k-mer/r.c.";
		if (firstField.length()<length)
			firstField += CommonUtils.padding(length-firstField.length(), ' ');
		if (use_weighted_hit_count)
			return firstField+"\tPosCt\twPosCt\tNegCt\tHgp_lg10";
		else
			return firstField+"\tPosCt\tNegCt\tHgp_lg10";
	}
	public static Kmer fromString(String str){
		String[] f = str.trim().split("\t");
		String[] f0f = f[0].split("/");
		HashSet<Integer> posHits = new HashSet<Integer>();
		HashSet<Integer> negHits = new HashSet<Integer>();
		
		int numField = f.length;	// 9 with WPos field, 8 without
		int posIDIndex = numField-2;
		int negIDIndex = numField-1;
		int HgpIndex = numField-3;
		
		String pos_id_string = f[posIDIndex].trim();
		if (!pos_id_string.equals("-1")){
			String[] pos_ids = null;
			if (pos_id_string.charAt(0)=='{'){	// new BitSet format
				pos_id_string = pos_id_string.substring(1,pos_id_string.length()-1);
				pos_ids = pos_id_string.split(", ");
			}
			else{
				pos_ids = pos_id_string.split(" ");
			}
			for (String hit:pos_ids)
				posHits.add(Integer.valueOf(hit));
		}
		String neg_id_string = f[negIDIndex].trim();
		if (!neg_id_string.equals("-1")){
			String[] neg_ids = null;
			if (neg_id_string.charAt(0)=='{'){	// new BitSet format
				neg_id_string = neg_id_string.substring(1,neg_id_string.length()-1);
				neg_ids = neg_id_string.split(", ");
			}
			else{
				neg_ids = neg_id_string.split(" ");
			}
			for (String hit:neg_ids)
				negHits.add(Integer.valueOf(hit));
		}
		
		Kmer kmer = new Kmer(f0f[0], posHits);	
		kmer.setNegHits(negHits);
		kmer.setWeightedPosHitCount();
		kmer.clusterId = Integer.parseInt(f[1]);
		kmer.kmerStartOffset = Integer.parseInt(f[2]);
		kmer.hgp_lg10 = Double.parseDouble(f[HgpIndex]);
		
		return kmer;
	}

	public static void printKmers(ArrayList<Kmer> kmers, int posSeqCount, int negSeqCount, double score, 
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
					sb.append("\t").append(kmer.posBits.toString()).append("\t").append(kmer.negBits.toString());
				sb.append("\n");
			}
		}
		if (printKmersAtK)
			CommonUtils.writeFile(String.format("%s_kmers_k%d.txt", filePrefix, kmers.get(0).getK()), sb.toString());
		else
			CommonUtils.writeFile(String.format("%s_KSM.txt", filePrefix), sb.toString());
	}
	
	/**
	 * If the set of ids is empty (no negative hits), return "-1"
	 */
	protected static String hits2string(HashSet<Integer> ids){
		StringBuilder sb = new StringBuilder();
		if (ids.isEmpty())
			return "-1";
		TreeSet<Integer> sorted = new TreeSet<Integer>(ids);
		for (int id:sorted)
			sb.append(id).append(" ");
		return sb.toString();
	}
	/**
	 * If the set of ids is empty (no negative hits), return "-1"
	 */
	protected static String hitBits2string(BitSet bits){
		return CommonUtils.encodeBytesToAscii85(bits.toByteArray()).replace("\n", "");
	}
	/**
	 * Clone the list and clone the kmer elements, also remember the seedKmer and mapped it to the top<br>
	 * Deep clone, makes correct linking of gapped kmers and subkmers
	 * @param kmers_in
	 * @param seedKmer a kmer in kmers_in to be mapped to the cloned kmers
	 * @return a list of cloned kmers, with the top kmer being the seedKmer
	 */
	public static ArrayList<Kmer> deepCloneKmerList(ArrayList<Kmer> kmers_in, Kmer seedKmer){
		ArrayList<Kmer> kmers = new ArrayList<Kmer>();
		HashMap<Kmer,Kmer> subkmer2subkmer = new HashMap<Kmer,Kmer>(); 
		Kmer seedKmerClone = null;
		
		// clone the distinct subkmers
		for (Kmer km:kmers_in){
			if (km instanceof GappedKmer){
				GappedKmer gk = (GappedKmer) km;
				for (Kmer sk: gk.getBaseKmers()){
					if (!subkmer2subkmer.containsKey(sk)){
						Kmer sk2 = sk.clone();
						sk2.clearGappedKmers();
						subkmer2subkmer.put(sk, sk2);
					}
				}
			}
		}
		// linking gk with sk
		for (Kmer km:kmers_in){
			Kmer km2 = null;
			if (km instanceof GappedKmer){
				GappedKmer gk = (GappedKmer) km;
				GappedKmer gk2 = gk.clone();
				gk2.clearBaseKmers();
				for (Kmer sk: gk.getBaseKmers()){
					gk2.addBaseKmer(subkmer2subkmer.get(sk),gk.getSubKmerOrientation(sk));
				}
				gk2.update();
				kmers.add(gk2);
				km2 = gk2;
			}
			else{
				km2 = km.clone();
				kmers.add(km2);
			}
			if (km == seedKmer)
				seedKmerClone = km2;
		}// GappedKmers need special treatment because of the many to many connections between gappedkmers and their subkmers
		
		Collections.sort(kmers);
		
		// move seedKmerClone to the top
		ArrayList<Kmer> results = new ArrayList<Kmer>();
		results.add(seedKmerClone);
		for (Kmer km: kmers)
			if (km != seedKmerClone)
				results.add(km);
		results.trimToSize();
		
		return results;
	}
	
	/**
	 * Clone the list and clone the kmer elements<br>
	 * Shallow clone: do not link gappedKmers and their subkmers, copy for output only		
	 * @param kmers_in
	 * @return
	 */
	public static ArrayList<Kmer> shallowCloneKmerList(ArrayList<Kmer> kmers_in){
		ArrayList<Kmer> kmers = new ArrayList<Kmer>();
		for (Kmer km:kmers_in){
				kmers.add(km.clone());
		}
		kmers.trimToSize();
		return kmers;
	}
	
	/**
	 * Copy the kmer references to a new list<br>
	 * Do not clone kmer elements
	 * @param kmers_in
	 * @return
	 */
	public static ArrayList<Kmer> copyKmerList(ArrayList<Kmer> kmers_in){
		ArrayList<Kmer> kmers = new ArrayList<Kmer>();
		for (Kmer km:kmers_in)
			kmers.add(km);
		kmers.trimToSize();
		return kmers;
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
	public static void printKmerHashcode(ArrayList<Kmer> kmers){
		for (Kmer km: kmers)
			System.err.println(km.toShortString()+"\t"+km.hashCode());
		System.err.println();
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
	
	public static void main9(String[] args){
//		ArrayList<Kmer> kmers = Kmer.loadKmers(new File(args[0]));
	     BitSet bits1 = new BitSet(16);
	     BitSet bits2 = new BitSet(16);
	      
	     // set some bits
	     for(int i=0; i<16; i++) {
	        if((i%2) == 0) bits1.set(i);
	        if((i%5) != 0) bits2.set(i);
	     }
	     
	}
}

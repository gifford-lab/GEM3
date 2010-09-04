package edu.mit.csail.cgs.deepseq.utilities;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.List;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.StrandedBase;
import edu.mit.csail.cgs.utils.stats.StatUtil;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.datasets.general.Point;
/**
 * Modify from AlignmentFileReader.java
 * 
 * Here we use it as a memory cache to store the data.
 * The fivePrimes field for each chrom/strand will be distinct. 
 * Multiple reads mapped to the same bp position will be stored as counts.
 * Thus, the count field is different from AlignmentFileReader.java.
 * 
 * This class basically stores the hits coming from the correspoding files. <br>
 * We have made use of an unusual convention for reducing running time purposes. <br>
 * The hits are basically being represented by 3 main fields: <tt>fivePrimes, hitCounts</tt>
 * and <tt>hitIDs</tt>. <br>
 * Each of these fields are 3D arrays where:  <br>
 * - the first dimension corresponds to the chromosome that a hit belongs to (based on
 * the mapping from a chromosome as a <tt>String</tt> to an integer via the 
 * <tt>chrom2ID</tt> map).    <br>
 * - the second dimension corresponds to the strand. 0 for '+' (Watson), 1 for '-' (Crick). <br>
 * - the third dimension contains information for a hit (e.g. its fivePrimes, counts, ID).
 * 
 * @author Yuchun
 *
 */
public class ReadCache {

	private Genome gen;
	private int numChroms;
	private String name;
	
	private double totalHits;
	private int totalBases;
	private int[] binCounts;			// histogram bins, index: count at bases; value: number of bases.
	private int[] bin500Counts;			// histogram bins, index: count at bases; value: number of bases.
	static final int BINSIZE = 501;
	//Data structures for pre-loading
	
	/**
	 * Chromosomes are stored as IDs for efficiency (saving memory) 
	 */
	//protected int[] chrs=null;
	
	/**
	 * Five prime ends of the read hits. <br>
	 * First dimension represents the corresponding chromosome ID. <br>
	 * Second dimension represents the strand. 0 for '+', 1 for '-' <br>
	 * Third dimension contains the coordinates of the hits
	 */
	private int[][][] fivePrimes=null;
	
	private ArrayList<Integer>[][] fivePrimesList = null;
	
	/**
	 * Number of read hits that corresponds to the 5' position
	 * First dimension represents the corresponding chromosome ID. <br>
	 * Second dimension represents the strand. 0 for '+', 1 for '-' <br>
	 * Third dimension contains the number of hits at corresponding start position 
	 */
	private float[][][] hitCounts=null;
	
	private ArrayList<Float>[][] hitCountsList = null;
	
	/**
	 * Strands of the read hits
	 */
	//protected char[][] strands=null;
	
	private HashMap<String, Integer> chrom2ID=new HashMap<String,Integer>();
	
	private HashMap<Integer,String> id2Chrom=new HashMap<Integer,String>();
	
	public ReadCache(Genome g, String name){
		totalHits=0;
		totalBases=0;
		gen=g;
		List<String> chromList = g.getChromList();
		numChroms = chromList.size();
		this.name = name;
		
		//Initialize the chromosome name lookup tables
		int i=0; 
		for(String c:chromList){
			chrom2ID.put(c, i);
			id2Chrom.put(i, c);
			i++;
		}
	
		//Initialize the data structures
		fivePrimes    = new int[numChroms][2][];
		hitCounts = new float[numChroms][2][];
		
		fivePrimesList    = new ArrayList[numChroms][2];
		for(i = 0; i < fivePrimesList.length; i++) { for(int j = 0; j < fivePrimesList[i].length; j++) { fivePrimesList[i][j] = new ArrayList<Integer>(); } }
		
		hitCountsList = new ArrayList[numChroms][2];
		for(i = 0; i < hitCountsList.length; i++) { for(int j = 0; j < hitCountsList[i].length; j++) { hitCountsList[i][j] = new ArrayList<Float>(); } }
	}//end of ReadCache constructor
	
	public List<StrandedBase> getUnstrandedBases(Region r) {
		List<StrandedBase> bases = new ArrayList<StrandedBase>();
		bases.addAll(getStrandedBases(r,'+'));
		bases.addAll(getStrandedBases(r,'-'));
		return bases;
	}
	/*
	 * Return both stranded base counts
	 * Subtract IP count by Ctrl count (stranded) (scaled by ratio)
	 * We do not want to modify the readCaches, thus do it here on the fly
	 */
	public List<StrandedBase> getSubtractedBases(Region r, ReadCache ctrl, double ratio) {
		List<StrandedBase> bases = new ArrayList<StrandedBase>();
		
		for (int s=0;s<2;s++){
			List<StrandedBase> ip_bases = getStrandedBases(r,s==0?'+':'-');
			List<StrandedBase> ctrl_bases = ctrl.getStrandedBases(r,s==0?'+':'-');
			int ctrl_idx = 0;
			for (StrandedBase b: ip_bases){
				while(ctrl_idx<ip_bases.size() && 
						ctrl_bases.get(ctrl_idx).getCoordinate()<b.getCoordinate()){
					ctrl_idx++;
				}
				// if there is control reads at the same position, subtract it
				if (ctrl_bases.get(ctrl_idx).getCoordinate()==b.getCoordinate()){
					float count = (float) (b.getCount() - ctrl_bases.get(ctrl_idx).getCoordinate()*ratio);
					if (count>0){
						b.setCount(count);
						bases.add(b);
					}
				}
				else{
					bases.add(b);
				}
			}
		}
		return bases;
	}
	
	/**
	 * Loads hits in the region
	 * @param r
	 * @return
	 */
	public List<StrandedBase> getStrandedBases(Region r, char strand) {
		List<StrandedBase> bases = new ArrayList<StrandedBase>();
		String chr = r.getChrom();
		int chrID = chrom2ID.get(chr);
		int j = (strand=='+') ? 0 : 1;
		int[] tempStarts = fivePrimes[chrID][j];		
		if(tempStarts.length != 0) {
			int start_ind = Arrays.binarySearch(tempStarts, r.getStart());
			int end_ind   = Arrays.binarySearch(tempStarts, r.getEnd());
			
			if( start_ind < 0 ) { start_ind = -start_ind - 1; }
			if( end_ind < 0 )   { end_ind   = -end_ind - 1; }
			
			start_ind = StatUtil.searchFrom(tempStarts, ">=", r.getStart(), start_ind);
			end_ind   = StatUtil.searchFrom(tempStarts, "<=",   r.getEnd(), end_ind);
			
			for(int k = start_ind; k <= end_ind; k++) {
				bases.add(new StrandedBase(strand, tempStarts[k], hitCounts[chrID][j][k]));
			}	
		}	
		return bases;
	}//end of getStrandedBases method
	
	public float countHits(Region r) {
		return StrandedBase.countBaseHits(getUnstrandedBases(r));
	}

	/**
	 * Gets the stranded count of all hits (of all chromosomes) for the specified strand
	 * @param strand 
	 * @return
	 */
	protected double getStrandedTotalCount(char strand) {
		int strandInd = strand == '+' ? 0 : 1;
		double count = 0;
		for(int i = 0; i < hitCounts.length; i++) {
			float[] hitCountsTemp = hitCounts[i][strandInd];
			for(float el:hitCountsTemp)
				count += (double)el;
		}
		return count;
	}//end of getStrandedTotalCount method
	
	/**
	 * 	Add hits to data structure
	 * 	It is called for ReadDB loader, store data loaded from ReadDB
	 * 	It is called multiple times to retrieve all the data, then populateArrays() is called 
	 */
	public void addHits(String chrom, char strand, Collection<Integer> coords, Collection<Float> counts){
		int chrID   = chrom2ID.get(chrom);
		int strandInd = strand == '+' ? 0 : 1;
		fivePrimesList[chrID][strandInd].addAll(coords);
		hitCountsList[chrID][strandInd].addAll(counts);
		for (float c: counts)
			totalHits += c;
		totalBases += coords.size();
	}//end of addHits method	
	
	/**
	 * Add all hit starts from all chrom and strand
	 * It is called once for file reader that has loaded data into memory
	 * assuming the data structure of starts is same as file reader 
	 */
	public void addAllFivePrimes(ArrayList<int[][][]> allStarts){
		for(int i = 0; i < fivePrimesList.length; i++){			// chrom
			for(int j = 0; j < fivePrimesList[i].length; j++){	// strand
				int[][][] tmp = allStarts.get(0);
				int[] allPositions = tmp[i][j];
				tmp[i][j]=null;
				for (int k=1;k<allStarts.size();k++){		// for each of files/replicates
					tmp = allStarts.get(k);
					allPositions = mergeOrderedList(allPositions, tmp[i][j]);
					tmp[i][j]=null;
				}
				System.gc();
				if (allPositions.length==0)
					continue;
				// consolidate counts of same bp position
				int count = 1;
				int previous = allPositions[0];
				for (int m=1;m<allPositions.length;m++){
					if (allPositions[m]==previous){
						count++;
					}
					else{
						fivePrimesList[i][j].add(previous);				// now file reader stores 5' end
						hitCountsList[i][j].add((float)count);
						count=1;
						previous = allPositions[m];
					}
				}
				// add the last element
				fivePrimesList[i][j].add(previous);				// now file reader stores 5' end
				hitCountsList[i][j].add((float)count);

				// update stats
				totalBases += fivePrimesList[i][j].size();
				for (float c: hitCountsList[i][j])
					totalHits += c;
			}
		}
	}//end of addAllFivePrimes method
	
	private int[] mergeOrderedList(int[] a, int[] b){
		int[] result = new int[a.length+b.length];
		int ai=0; int bi=0;
		for (int i=0;i<result.length;i++){
			if (bi!=b.length && (ai==a.length || a[ai]>b[bi])){
				result[i]=b[bi];
				bi++;
			}
			else{
				result[i]=a[ai];
				ai++;
			}
		}
		return result;
	}//end of mergeOrderedList method
	
	/**
	 * Converts lists of Integers to integer arrays, deletes the lists for saving memory
	 * all array elements are ordered in terms of the array <tt>starts</tt>.
	 */
	public void populateArrays() {
		for(int i = 0; i < fivePrimesList.length; i++)
			for(int j = 0; j < fivePrimesList[i].length; j++)
				fivePrimes[i][j] = list2int(fivePrimesList[i][j]);
		fivePrimesList = null;
		for(int i = 0; i < hitCountsList.length; i++)
			for(int j = 0; j < hitCountsList[i].length; j++)
				hitCounts[i][j] = list2float(hitCountsList[i][j]);
		hitCountsList = null;
		System.gc();
		generateStats();
	}//end of populateArrays method
	
	public void generateStats() {
		// count readHit numbers in 1bp bins
		int max = 200;
		binCounts = new int[max+1];
		for(int i = 0; i < hitCounts.length; i++)
			for(int j = 0; j < hitCounts[i].length; j++)
				for(int k = 0; k < hitCounts[i][j].length; k++){
					int count = (int)hitCounts[i][j][k];
					count = count>max?max:count;
					binCounts[count]++;
				}
		// count readHit numbers in BINSIZE (500bp) bins
		max = 2000;		
		bin500Counts = new int[max+1];		
		for(int i = 0; i < hitCounts.length; i++){
			String chrom = id2Chrom.get(i);
			int totalLength = gen.getChromLength(chrom);
			for (int j=0;j<totalLength;j+=BINSIZE){
				int start=j+1;
				int end=j+BINSIZE;
				int count = (int)countHits(new Region(gen, chrom, start, end));
				count = count>max?max:count;
				bin500Counts[count]++;
			}
		}	
	}//end of generateStats method
	
	/**
	 * This method renormalizes the counts of this channel 
	 * by multiplying with a constant factor.
	 * @param factor
	 */
	public void normalizeCounts(double factor) {
		for(int i = 0; i < hitCounts.length; i++)
			for(int j = 0; j < hitCounts[i].length; j++)
				for(int k = 0; k < hitCounts[i][j].length; k++)
					hitCounts[i][j][k] = (float)(hitCounts[i][j][k]*factor);
//		hitCounts[i][j][k] *= factor;
		
		updateTotalHits();
	}//end of normalizeCounts method
	
	public void filterAllBases(float maxReadperBP){
		for(int i = 0; i < hitCounts.length; i++)
			for(int j = 0; j < hitCounts[i].length; j++)
				for(int k = 0; k < hitCounts[i][j].length; k++)
					if (hitCounts[i][j][k] > maxReadperBP)
						hitCounts[i][j][k] = maxReadperBP;
						
		updateTotalHits();
	}
	/*
	 * Reset bases with huge number of reads to 1
	 * Return the base positions that are reset
	 */
	public ArrayList<Pair<Point, Float>> resetHugeBases(float threshold){
		ArrayList<Pair<Point, Float>> bases = new ArrayList<Pair<Point, Float>>();
		for(int i = 0; i < hitCounts.length; i++)
			for(int j = 0; j < hitCounts[i].length; j++)
				for(int k = 0; k < hitCounts[i][j].length; k++)
					if (hitCounts[i][j][k] > threshold){
						Point p = new Point(gen, id2Chrom.get(i), fivePrimes[i][j][k]);
						bases.add(new Pair<Point, Float>(p, hitCounts[i][j][k]));
						hitCounts[i][j][k] = 1;
					}
						
		updateTotalHits();
		return bases;
	}
	/*
	 * Exclude the data from specified regions
	 * It will truncate the data arrays.
	 */
	public void excludeRegions(ArrayList<Region> regions){
		for (Region r:regions){
			int chrID = chrom2ID.get(r.getChrom());
			for(int j = 0; j < fivePrimes[chrID].length; j++){
				int[] tempStarts = fivePrimes[chrID][j];		
				if(tempStarts.length != 0) {
					int start_ind = Arrays.binarySearch(tempStarts, r.getStart());
					int end_ind   = Arrays.binarySearch(tempStarts, r.getEnd());					
					if( start_ind < 0 ) { start_ind = -start_ind - 1; }
					if( end_ind < 0 )   { end_ind   = -end_ind-1-1 ; }
					
					int length = start_ind + tempStarts.length-end_ind-1;
					int[] newData = new int[length];
					float[] tempCounts = hitCounts[chrID][j];
					float[] newCounts = new float[length];
					if (start_ind>0){
						System.arraycopy(tempStarts, 0, newData, 0, start_ind);
						System.arraycopy(tempCounts, 0, newCounts, 0, start_ind);
					}
					if (tempStarts.length-end_ind-1>0){
						System.arraycopy(tempStarts, end_ind+1, newData, start_ind, tempStarts.length-end_ind-1);
						System.arraycopy(tempCounts, end_ind+1, newCounts, start_ind, tempCounts.length-end_ind-1);
					}
					fivePrimes[chrID][j] = newData;
					hitCounts[chrID][j] = newCounts;
				}	
			}					
		}
		updateTotalHits();
		
	}
	private void updateTotalHits(){
		totalHits = 0.0;
		for(int i = 0; i < hitCounts.length; i++)
			for(int j = 0; j < hitCounts[i].length; j++)
				for(int k = 0; k < hitCounts[i][j].length; k++)
					totalHits += hitCounts[i][j][k];
	}
	private int[] list2int(List<Integer> list) {
		int[] out = new int[list.size()];
		for(int i = 0; i < out.length; i++)
			out[i] = list.get(i);
		return out;
	}
	
	private float[] list2float(List<Float> list) {
		float[] out = new float[list.size()];
		for(int i = 0; i < out.length; i++)
			out[i] = list.get(i);
		return out;
	}	
	
	//Accessors	
	/**
	 * get the total number of hits (of the all alignment/files)
	 * @return
	 */
	public double getHitCount(){
		return totalHits;
	} 
	
	public double getBaseCount(){
		return totalBases;
	} 
	
	public String getName(){
		return name;
	}
	
	public void displayStats(){
		System.out.println(name+"\tBases: "+totalBases+"\tHitCounts: "+totalHits);
	}
	
	public void printBinCounts(){
		StringBuilder sb = new StringBuilder();
		for (int i=0;i<binCounts.length;i++){
			sb.append(i+"\t"+binCounts[i]+"\n");
		}
		writeFile(name.trim()+"_1bpCount.txt", sb.toString());
	}
	public void printBin500Counts(){
		StringBuilder sb = new StringBuilder();
		for (int i=0;i<bin500Counts.length;i++){
			sb.append(i+"\t"+bin500Counts[i]+"\n");
		}
		writeFile(name.trim()+"_500bpCount.txt", sb.toString());
	}
	
	public int getMaxHitPerBP(double fraction){
		double toKeep = (1-fraction) * totalHits;
		int accumulative = 0;
		for (int i=0;i<binCounts.length;i++){
			accumulative += binCounts[i]*i;
			if (accumulative>=toKeep)
				return i;
		}
		return binCounts.length;
	}
	
	private void writeFile(String fileName, String text){
		try{
			FileWriter fw = new FileWriter(fileName, false); //new file
			fw.write(text);
			fw.close();
		} 
		catch (IOException e) {
			e.printStackTrace();
		}
	}
	
}//end of ReadCache class

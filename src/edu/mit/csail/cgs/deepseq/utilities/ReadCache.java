package edu.mit.csail.cgs.deepseq.utilities;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;

import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.StrandedBase;
import edu.mit.csail.cgs.utils.probability.NormalDistribution;
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
	
	private TreeMap<String, Integer> chrom2ID=new TreeMap<String,Integer>();	
	private HashMap<Integer,String> id2Chrom=new HashMap<Integer,String>();
	
	public ReadCache(Genome g, String name, HashMap<String, Integer> chrom2ID, HashMap<Integer,String> id2Chrom){
		totalHits=0;
		totalBases=0;
		gen=g;
		List<String> chromList = g.getChromList();
		numChroms = chromList.size();
		this.name = name;
		if (chrom2ID!=null && id2Chrom!=null){
			this.chrom2ID.putAll(chrom2ID);
			this.id2Chrom=id2Chrom;
		}
		else{
			int i=0; 
			for(String c:chromList){
				this.chrom2ID.put(c, i);
				this.id2Chrom.put(i, c);
				i++;
			}
		}
	
		//Initialize the data structures
		fivePrimes    = new int[numChroms][2][];
		hitCounts = new float[numChroms][2][];
		
		fivePrimesList = new ArrayList[numChroms][2];
		for(int i = 0; i < fivePrimesList.length; i++) { for(int j = 0; j < fivePrimesList[i].length; j++) { fivePrimesList[i][j] = new ArrayList<Integer>(); } }
		
		hitCountsList = new ArrayList[numChroms][2];
		for(int i = 0; i < hitCountsList.length; i++) { for(int j = 0; j < hitCountsList[i].length; j++) { hitCountsList[i][j] = new ArrayList<Float>(); } }
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
				while(ctrl_idx<ctrl_bases.size() && 
						ctrl_bases.get(ctrl_idx).getCoordinate()<b.getCoordinate()){
					ctrl_idx++;
				}
				if (ctrl_idx==ctrl_bases.size()){	// reach end of control reads
					bases.add(b);
					continue;
				}
					
				// if there is control reads at the same position, subtract it
				if ( ctrl_bases.get(ctrl_idx).getCoordinate()==b.getCoordinate()){
					float count = (float) (b.getCount() - ctrl_bases.get(ctrl_idx).getCount()*ratio);
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
			
            while (start_ind > 0 && tempStarts[start_ind - 1] >= r.getStart() ) {
                start_ind--;
            }
            while (end_ind < tempStarts.length && tempStarts[end_ind] <= r.getEnd()) {
                end_ind++;
            }
			for(int k = start_ind; k < end_ind; k++) {
				bases.add(new StrandedBase(strand, tempStarts[k], hitCounts[chrID][j][k]));
			}	
		}	
		return bases;
	}//end of getStrandedBases method
	
	/**
	 * Count hits in the region, both strands
	 * @param r
	 * @return
	 */	
	public float countHits(Region r) {
		return countStrandedBases(r,'+')+countStrandedBases(r,'-');
	}
	
    public float countStrandedBases(Region r, char strand) {
		String chr = r.getChrom();
		int chrID = chrom2ID.get(chr);
		int j = (strand=='+') ? 0 : 1;
		int[] tempStarts = fivePrimes[chrID][j];		
        float count = 0;
		if(tempStarts.length != 0) {
			int start_ind = Arrays.binarySearch(tempStarts, r.getStart());
			int end_ind   = Arrays.binarySearch(tempStarts, r.getEnd());
			if( start_ind < 0 ) { start_ind = -start_ind - 1; }
			if( end_ind < 0 )   { end_ind   = -end_ind - 1; }
			
            while (start_ind > 0 && tempStarts[start_ind - 1] >= r.getStart() ) {
                start_ind--;
            }
            while (end_ind < tempStarts.length && tempStarts[end_ind] <= r.getEnd()) {
                end_ind++;
            }
			for(int k = start_ind; k < end_ind; k++) {
                count += hitCounts[chrID][j][k];
            }
        }
        return count;
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
				try{
					int[][][] tmp = allStarts.get(0);		// assuming only one text file
					int[] allPositions = tmp[i][j];
					tmp[i][j]=null;
					for (int k=1;k<allStarts.size();k++){		// for each of files/replicates
						tmp = allStarts.get(k);
						allPositions = mergeOrderedList(allPositions, tmp[i][j]);
						tmp[i][j]=null;
					}
					System.gc();
					if (allPositions==null || allPositions.length==0)
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
					fivePrimesList[i][j].trimToSize();
					hitCountsList[i][j].trimToSize();
					
					// update stats
					totalBases += fivePrimesList[i][j].size();
					for (float c: hitCountsList[i][j])
						totalHits += c;
				}catch (Exception e){
					System.err.println("Error: loading chomosome "+id2Chrom.get(i)+" "+(j==0?"+":"-")+" strand.");
					e.printStackTrace(System.err);
					System.exit(-1);
				}
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
	 * Converts lists of Integers to integer arrays, deletes the lists for saving memory. <br>
	 * This is usually called after addHits() or addAllFivePrimes().
	 * All array elements are ordered in terms of the <tt>five primes of reads</tt>.
	 */
	public void populateArrays(boolean generateStats) {
		for(int i = 0; i < fivePrimesList.length; i++)
			for(int j = 0; j < fivePrimesList[i].length; j++)
				fivePrimes[i][j] = list2int(fivePrimesList[i][j]);
		fivePrimesList = null;
		for(int i = 0; i < hitCountsList.length; i++)
			for(int j = 0; j < hitCountsList[i].length; j++)
				hitCounts[i][j] = list2float(hitCountsList[i][j]);
		hitCountsList = null;
		System.gc();
		if (generateStats)
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
	/**
	 * Reset duplicate reads that pass Poisson threshold. 
	 * The Poisson lambda parameter is calculated by an Gaussian Kernel smoothing
	 * that puts more weight on nearby bases (same chrom, same strand),
	 * the position being considered is excluded for Guassian Kernel computation
	 */
	public void applyPoissonGaussianFilter(double threshold, int width, int rand_seed){
        //init the Guassian kernel prob. for smoothing the read profile of called events
		double g[] = new double[width*4+1];
		NormalDistribution gaussianDist = new NormalDistribution(0, width*width);
		for (int i=0;i<g.length;i++)
			g[i]=gaussianDist.calcProbability((double)i);
		
		DRand re = new DRand(rand_seed);
		Poisson P = new Poisson(0, re);
			
		for(int i = 0; i < hitCounts.length; i++)
			for(int j = 0; j < hitCounts[i].length; j++){
				float counts[] = hitCounts[i][j];
				int pos[] = fivePrimes[i][j];
				for(int k = 0; k < counts.length; k++){
					int posK = pos[k]; 
					double sum = 0;
					for (int x=1;x<=width*4;x++){		// at most extend out 250 idx
						if (k+x>=counts.length|| pos[k+x]-posK>width*4)
							break;
						sum += counts[k+x]*g[pos[k+x]-posK];
					}
					for (int x=1;x<=width*4;x++){		// at most extend out 250 idx
						if (k-x<0 || posK-pos[k-x]>width*4)
							break;
						sum += counts[k-x]*g[posK-pos[k-x]];
					}
					sum = sum/(1-g[0]);				// exclude this position for evaluation
					
					double countThres=0;
					P.setMean(sum);
					double pvalue=1;
					for(int b=1; pvalue>threshold; b++){
						pvalue=1-P.cdf(b);	//p-value as the tail of Poisson
						countThres=b;
					}
										
//					System.out.println(String.format("%s%d\t%.1f\t%.1f\t%.3f", (j==0?"+":"-"), posK, counts[k], countThres, sum));
					if (counts[k] > Math.max(1,countThres))
						counts[k] = (float) Math.max(1,countThres);					
				}
			}
		
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
	
	public String getStatsString(){
		return (name+"\tBases: "+totalBases+"\tHitCounts: "+totalHits);
	}
	
	public void deleteUnenrichedReadData(ArrayList<Region> enrichedRegions){
		fivePrimesList = new ArrayList[numChroms][2];
		for(int i = 0; i < fivePrimesList.length; i++) { for(int j = 0; j < fivePrimesList[i].length; j++) { fivePrimesList[i][j] = new ArrayList<Integer>(); } }
		hitCountsList = new ArrayList[numChroms][2];
		for(int i = 0; i < hitCountsList.length; i++) { for(int j = 0; j < hitCountsList[i].length; j++) { hitCountsList[i][j] = new ArrayList<Float>(); } }

		for(Region r:enrichedRegions){
			int chrID   = chrom2ID.get(r.getChrom());
			for (int strandInd=0; strandInd<=1;strandInd++){
				List<StrandedBase> bases = getStrandedBases(r, strandInd==0?'+':'-');
				for (StrandedBase b:bases){
					fivePrimesList[chrID][strandInd].add(b.getCoordinate());
					hitCountsList[chrID][strandInd].add(b.getCount());
				}
			}
		}
		fivePrimes = new int[numChroms][2][];
		hitCounts = new float[numChroms][2][];		
		
		populateArrays(false);
	}
	
	public void printBinCounts(String prefix){
		StringBuilder sb = new StringBuilder();
		for (int i=0;i<binCounts.length;i++){
			sb.append(i+"\t"+binCounts[i]+"\n");
		}
		CommonUtils.writeFile(prefix+"_"+name.trim()+"_1bpCount.txt", sb.toString());
	}
	public void printBin500Counts(String prefix){
		StringBuilder sb = new StringBuilder();
		for (int i=0;i<bin500Counts.length;i++){
			sb.append(i+"\t"+bin500Counts[i]+"\n");
		}
		CommonUtils.writeFile(prefix+"_"+name.trim()+"_500bpCount.txt", sb.toString());
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
	/** 
	 * Write Read Start Count (RSC) file
	 */
	public void writeRSC(){
		StringBuilder sb = new StringBuilder();
		for(int i = 0; i < hitCounts.length; i++){
			String chrom = id2Chrom.get(i);
			for(int j = 0; j < hitCounts[i].length; j++){
				StringBuilder sb_sub = new StringBuilder();
				char strand = j==0?'+':'-';
				int subCount = hitCounts[i][j].length;
				if (subCount>0){
					for(int k = 0; k < subCount; k++){
						sb_sub.append(String.format("%d\t%.1f\n",  fivePrimes[i][j][k], hitCounts[i][j][k]));
					}
					sb.append(String.format("%s\t%c\t%d\n",  chrom, strand, subCount));
					sb.append(sb_sub).append("\n");
				}
			}
		}
		CommonUtils.writeFile(name.trim()+".rsc", sb.toString());
	}
	/** 
	 * Write genetrack format file
	 * # genetrack tab-delimited file lists the gnomic coordinates (chrome and index) of the 5' ends of sequencing tags on the forward and reverse strands.
# chrome: the name of the chromosome
# index: the 5' ends position of sequencing tags
# forward: number of tag counts on the forward (+) strand
# reverse: number of tag counts on the reverse (-) strand
# value: the sum of tag counts on both strands
	 */
	public void writeGeneTrack(){
		StringBuilder sb = new StringBuilder();
		sb.append("#chrom\tindex\tforward\treverse\tvalue\n");
		for(int i = 0; i < hitCounts.length; i++){
			String chrom = id2Chrom.get(i);
			HashMap<Integer, Pair<Integer, Integer>> strandedCounts = new HashMap<Integer, Pair<Integer, Integer>>();
			for(int j = 0; j < hitCounts[i].length; j++){
				int subCount = hitCounts[i][j].length;
				if (subCount>0){
					for(int k = 0; k < subCount; k++){
						int pos = fivePrimes[i][j][k];
						if (j==0)
							strandedCounts.put(pos, new Pair<Integer, Integer>((int)hitCounts[i][j][k],0));
						else{
							if (!strandedCounts.containsKey(pos))
								strandedCounts.put(pos, new Pair<Integer, Integer>(0,(int)hitCounts[i][j][k]));
							else
								strandedCounts.put(pos, new Pair<Integer, Integer>(strandedCounts.get(pos).car(),(int)hitCounts[i][j][k]));
						}
					}
				}
			}
			ArrayList<Integer> coords = new ArrayList<Integer>();
			coords.addAll(strandedCounts.keySet());
			Collections.sort(coords);
			for(int pos:coords){
				Pair<Integer, Integer> pair = strandedCounts.get(pos);
				sb.append(String.format("chr%s\t%d\t%d\t%d\t%d\n", chrom, pos, pair.car(),pair.cdr(),pair.car()+pair.cdr()));
			}
			CommonUtils.appendFile(name.trim()+".genetrack", sb.toString());
			sb = new StringBuilder();
		}
		
	}
	/** 
	 * Read Read Start Count (RSC) file
	 */
	public void readRSC(String filename) throws IOException{
		BufferedReader bin = null;
        bin = new BufferedReader(new InputStreamReader(new FileInputStream(filename)));
		
		for(int i = 0; i < fivePrimes.length; i++){
			for(int j = 0; j < fivePrimes[i].length; j++){
				fivePrimes[i][j] = new int[0];
				hitCounts[i][j] = new float[0];
			}
		}

        int[] currentCoords = null;
        float[] currentCounts = null;
        int idx = -1;
        totalHits = 0;
        totalBases = 0;
        String line;
        while((line = bin.readLine()) != null) { 
            line = line.trim();
            if (line.equals(""))
            	continue;
            String[] f=line.split("\t");
            if (f.length==3){		// header
            	int i = chrom2ID.get(f[0]);
            	int j = f[1].charAt(0)=='+'?0:1;
            	int baseCount = Integer.parseInt(f[2]);
            	totalBases += baseCount;
            	currentCoords = new int[baseCount];
            	currentCounts = new float[baseCount];
            	fivePrimes[i][j] = currentCoords;
            	hitCounts[i][j] = currentCounts;
            	idx = 0;
	            continue;
            }
            if (f.length==2){		// data
            	currentCoords[idx] = Integer.parseInt(f[0]);
            	float count = Float.parseFloat(f[1]);
            	currentCounts[idx] = count;
            	totalHits += count;
            	idx++;
            }
        }			
        if (bin != null) {
            bin.close();
        }
	}
}//end of ReadCache class

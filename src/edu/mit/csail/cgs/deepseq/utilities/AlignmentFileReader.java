package edu.mit.csail.cgs.deepseq.utilities;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;

import edu.mit.csail.cgs.datasets.chipseq.ChipSeqAlignment;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.Read;
import edu.mit.csail.cgs.deepseq.ReadHit;
import edu.mit.csail.cgs.projects.readdb.ClientException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.stats.StatUtil;

/**
 * This class basically stores the hits coming from the correspoding files. <br>
 * We have made use of an unusual convention for reducing running time purposes. <br>
 * The hits are basically being represented by 3 main fields: <tt>starts, hitCounts</tt>
 * and <tt>hitIDs</tt>. <br>
 * Each of these fields are 3D arrays where:  <br>
 * - the first dimension corresponds to the chromosome that a hit belongs to (based on
 * the mapping from a chromosome as a <tt>String</tt> to an integer via the 
 * <tt>chrom2ID</tt> map).    <br>
 * - the second dimension corresponds to the strand. 0 for '+' (Watson), 1 for '-' (Crick). <br>
 * - the third dimension contains information for a hit (e.g. its start, counts, ID).
 * 
 * @author shaunmahony
 *
 */
public abstract class AlignmentFileReader {

	protected File inFile;
	protected double totalHits;
	protected double totalWeight;
	protected Genome gen;
	protected int misMatch;
	protected boolean useNonUnique;
	protected int numChroms;
		
	//Data structures for pre-loading
	
	/**
	 * Chromosomes are stored as IDs for efficiency (saving memory) 
	 */
	//protected int[] chrs=null;
	
	/**
	 * Starts of the read hits. <br>
	 * First dimension represents the corresponding chromosome ID. <br>
	 * Second dimension represents the strand. 0 for '+', 1 for '-' <br>
	 * Third dimension contains the starts of the hits
	 */
	protected int[][][] starts=null;
	
	protected ArrayList<Integer>[][] startsList = null;
	
	/**
	 * Number of hits that corresponds to the <tt>Read</tt> of each <tt>ReadHit</tt>
	 * First dimension represents the corresponding chromosome ID. <br>
	 * Second dimension represents the strand. 0 for '+', 1 for '-' <br>
	 * Third dimension contains the number of hits of the <tt>Read</tt> of each <tt>ReadHit</tt> 
	 */
	protected int[][][] hitCounts=null;
	
	protected ArrayList<Integer>[][] hitCountsList = null;
	
	/**
	 * The IDs of the read hits
	 * First dimension represents the corresponding chromosome ID. <br>
	 * Second dimension represents the strand. 0 for '+', 1 for '-' <br>
	 * Third dimension contains the IDs of the read hits
	 */
	protected int[][][] hitIDs=null;//Only necessary if the ID field in ReadHit becomes important
	
	protected ArrayList<Integer>[][] hitIDsList = null;
	
	/**
	 * Strands of the read hits
	 */
	//protected char[][] strands=null;
	
	protected HashMap<String, Integer> chrom2ID=new HashMap<String,Integer>();
	protected HashMap<Integer,String> id2Chrom=new HashMap<Integer,String>();
	protected int readLength;
	protected int currID=0;
	
	public AlignmentFileReader(File f, Genome g, int mis, boolean nonUnique, int idSeed){
		gen=g;
		totalHits=0;
		totalWeight=0;
		inFile=f;
		misMatch=mis;
		useNonUnique = nonUnique;
		currID=idSeed;
		if(gen==null)
			estimateGenome();
		List<String> chromList = gen.getChromList();
		numChroms = chromList.size();
		
		System.out.print("Loading reads from: "+f.getName()+" ... ");
		
		//Initialize the chromosome name lookup tables
		int i=0; 
		for(String c:chromList){
			chrom2ID.put(c, i);
			id2Chrom.put(i, c);
			i++;
		}
	
		//Initialize the data structures
		starts    = new int[numChroms][2][];
		hitCounts = new int[numChroms][2][];
		hitIDs    = new int[numChroms][2][];
		
		startsList    = new ArrayList[numChroms][2];
		for(i = 0; i < startsList.length; i++) { for(int j = 0; j < startsList[i].length; j++) { startsList[i][j] = new ArrayList<Integer>(); } }
		
		hitCountsList = new ArrayList[numChroms][2];
		for(i = 0; i < hitCountsList.length; i++) { for(int j = 0; j < hitCountsList[i].length; j++) { hitCountsList[i][j] = new ArrayList<Integer>(); } }
		
		hitIDsList    = new ArrayList[numChroms][2];
		for(i = 0; i < hitIDsList.length; i++) { for(int j = 0; j < hitIDsList[i].length; j++) { hitIDsList[i][j] = new ArrayList<Integer>(); } }
		
		countReads();
		System.out.println("Loaded");		
	}
	
	//Accessors
	public Genome getGenome(){return gen;}
	public void setGenome(Genome g){gen = g;}
	
	/**
	 * Loads reads that have hits in coords that pass the mismatch thresholds
	 * @param coords
	 * @return
	 */
	public List<ReadHit> loadHits(Region coords) {
		List<ReadHit> hits = new ArrayList<ReadHit>();
		String chr = coords.getChrom();
		if(chrom2ID.containsKey(chr)){
			int chrID = chrom2ID.get(chr);
			
			// j = 0 <=> strand = '+', j = 1 <=> strand = '-'
			for(int j = 0; j <= 1; j++) {
				
				//assumes tempStarts is sorted
				int[] tempStarts = starts[chrID][j];
				
				if(tempStarts.length != 0) {
					int start_ind = Arrays.binarySearch(tempStarts, coords.getStart());
					int end_ind   = Arrays.binarySearch(tempStarts, coords.getEnd());
					
					if( start_ind < 0 ) { start_ind = -start_ind - 1; }
					if( end_ind < 0 )   { end_ind   = -end_ind - 1; }
					
					start_ind = StatUtil.searchFrom(tempStarts, ">=", coords.getStart(), start_ind);
					end_ind   = StatUtil.searchFrom(tempStarts, "<=",   coords.getEnd(), end_ind);
					
					char str = j == 0 ? '+' : '-';			
					for(int k = start_ind; k <= end_ind; k++) {
						double weight = 1/(double)hitCounts[chrID][j][k];
						hits.add(new ReadHit(gen, hitIDs[chrID][j][k], chr, tempStarts[k], tempStarts[k]+readLength-1, str, weight ));
					}	
				}						
			}
		}
		return hits;
	}//end of loadHits method

	// load paired read hit 5' coordinates (sorted, can be duplicates) and counts 
	public Pair<ArrayList<Integer>,ArrayList<Float>> loadStrandedFivePrimeCounts(Region coords , char strand){
		ArrayList<Integer> hits = new ArrayList<Integer>();
		ArrayList<Float> counts = new ArrayList<Float>();
		String chr = coords.getChrom();
		if(chrom2ID.containsKey(chr)){
			int chrID = chrom2ID.get(chr);
			int j = strand == '+'?0:1;	// j = 0 <=> strand = '+', j = 1 <=> strand = '-'
			int shift = strand=='+'?0:readLength-1; 		// shift for minus strand
			int[] tempStarts = starts[chrID][j];	//assumes tempStarts is sorted	
			if(tempStarts.length != 0) {
				int start_ind = Arrays.binarySearch(tempStarts, coords.getStart()-shift);
				int end_ind   = Arrays.binarySearch(tempStarts, coords.getEnd()-shift);
				if( start_ind < 0 ) { start_ind = -start_ind - 1; }
				if( end_ind < 0 )   { end_ind   = -end_ind - 1; }
				start_ind = StatUtil.searchFrom(tempStarts, ">=", coords.getStart()-shift, start_ind);
				end_ind   = StatUtil.searchFrom(tempStarts, "<=",   coords.getEnd()-shift, end_ind);
					
				for(int k = start_ind; k <= end_ind; k++) {
					hits.add(tempStarts[k]+shift);
					counts.add((float)hitCounts[chrID][j][k]);
				}
			}
		}
		return new Pair<ArrayList<Integer>,ArrayList<Float>>(hits, counts);
	}
	/**
	 * Count total reads & weight
	 */
	protected abstract void countReads();
		
	/**
	 * Estimate a fake genome from the observed read coordinates
	 */
	protected abstract void estimateGenome();
	
	/**
	 * Gets the stranded weight of all hits (of all chromosomes) for the specified strand
	 * @param strand 
	 * @return
	 */
	protected double getStrandedTotalWeight(char strand) {
		int strandInd = strand == '+' ? 0 : 1;
		double weight = 0;
		for(int i = 0; i < hitCounts.length; i++) {
			int[] hitCountsTemp = hitCounts[i][strandInd];
			for(Integer el:hitCountsTemp)
				weight += 1/(double)el;
		}
		
		return weight;
	}//end of getStrandedTotalWeight method
	
	
	//Add hits to data structure
	protected void addHits(Read r){
		int numHits = r.getHits().size();
		for(ReadHit h : r.getHits()){
			int chrID   = chrom2ID.get(h.getChrom());
			char strand = h.getStrand();
			int strandInd = strand == '+' ? 0 : 1;
			//System.out.println(h.getChrom()+"\t"+h.getStart()+"\t"+h.getStrand());
			startsList[chrID][strandInd].add(h.getStart());
			hitIDsList[chrID][strandInd].add(h.getID());
			hitCountsList[chrID][strandInd].add(numHits);
			totalHits++;
		}
	}//end of addHits method
	
	/**
	 * Converts lists of integers to integer arrays, deletes the lists for saving
	 * memory and permutes all array elements so that they are ordered in terms of
	 * the array <tt>starts</tt>.
	 */
	protected void populateArrays() {
		
		for(int i = 0; i < startsList.length; i++)
			for(int j = 0; j < startsList[i].length; j++)
				starts[i][j] = list2int(startsList[i][j]);
		
		for(int i = 0; i < startsList.length; i++) { for(int j = 0; j < startsList[i].length; j++) { startsList[i][j].clear(); } }
		startsList = null;

		////////////////////////////////////////////////////////////////////////

		for(int i = 0; i < hitIDsList.length; i++)
			for(int j = 0; j < hitIDsList[i].length; j++)
				hitIDs[i][j] = list2int(hitIDsList[i][j]);
		
		for(int i = 0; i < hitIDsList.length; i++) { for(int j = 0; j < hitIDsList[i].length; j++) { hitIDsList[i][j].clear(); } }
		hitIDsList = null;
		
		////////////////////////////////////////////////////////////////////////

		for(int i = 0; i < hitCountsList.length; i++)
			for(int j = 0; j < hitCountsList[i].length; j++)
				hitCounts[i][j] = list2int(hitCountsList[i][j]);
		
		for(int i = 0; i < hitCountsList.length; i++) { for(int j = 0; j < hitCountsList[i].length; j++) { hitCountsList[i][j].clear(); } }
		hitCountsList = null;
		
		////////////////////////////////////////////////////////////////////////
		
		for(int i = 0; i < starts.length; i++) {  // chr
			for(int j = 0; j < starts[i].length; j++) { // strand
				int[] inds = StatUtil.findSort(starts[i][j]);
				hitIDs[i][j] = StatUtil.permute(hitIDs[i][j], inds);
				hitCounts[i][j] = StatUtil.permute(hitCounts[i][j], inds);
			}
		}
	}//end of populateArrays method

	
	/***************
	 ** ACCESSORS **
	 **************/
	
	/**
	 * Get the read length
	 */
	public int getReadLength(){return readLength;}
	
	/**
	 * get the total number of hits (of the alignment)
	 * @return
	 */
	public double getTotalHits(){if(totalHits==0){countReads();}return totalHits;}
	
	/**
	 * get the total weight of hits (of the alignment)
	 * @return
	 */
	public double getTotalWeight(){if(totalWeight==0){countReads();}return totalWeight;}
	
	/**
	 * get the ID of the alignment
	 * @return
	 */
	public int getCurrID(){return currID;}
	
	/**
	 * Count the lines in a file (this is typically, but not always, the maximum possible number of valid hits) 
	 */
	protected int countMaxHits(){
		try {
			LineNumberReader lineCounter = new LineNumberReader(new FileReader(inFile));
			String nextLine = null;
			while ((nextLine = lineCounter.readLine()) != null) {
				if (nextLine == null)
					break;
			}
			int count = lineCounter.getLineNumber();
			lineCounter.close();
			return(count);
		} catch (Exception done) {
			done.printStackTrace();
		}
		return(0);
	}//end of countMaxHits method
	
	/**
	 * get all the start positions
	 * clean the memory space
	 * @return
	 */
	public int[][][] getStarts(){
		hitIDs = null;
		hitCounts = null;
		System.gc();
		return starts;
	}//end of getStarts method
	
	/**
	 * Get sorted start positions of all reads (regardless of strand) in one chrom
	 */
	public int[] getStartCoords(String chrom){
		if(chrom2ID.containsKey(chrom)){
			int chromID = chrom2ID.get(chrom);
			int[] watsonReads = starts[chromID][0];
			int[] crickReads = starts[chromID][0];
			int[] allReads = new int[watsonReads.length+crickReads.length];
			int watson_idx=0; 
			int crick_idx=0;
			for (int i=0;i<allReads.length;i++){
				if (watsonReads[watson_idx]>crickReads[crick_idx]){
					allReads[i]=crickReads[crick_idx];
					if (crick_idx<crickReads.length-1)
						crick_idx++;
				}
				else{
					allReads[i]=watsonReads[watson_idx];
					if (watson_idx<watsonReads.length-1)
						watson_idx++;
				}
			}
			return allReads;
		}else{
			return null;
		}
	}//end of getStartCoords method
		
	private int[] list2int(List<Integer> list) {
		int[] out = new int[list.size()];
		for(int i = 0; i < out.length; i++)
			out[i] = list.get(i);
	   
		return out;
	}//end of list2int method

	public void cleanup(){
		starts=null;
		hitCounts=null;
		hitIDs=null;
	}
}//end of AlignmentFileReader class

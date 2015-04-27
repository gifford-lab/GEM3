package edu.mit.csail.cgs.deepseq.utilities;

import java.io.File;
import java.io.FileReader;
import java.io.LineNumberReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.Read;
import edu.mit.csail.cgs.deepseq.ReadHit;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public abstract class PairedAlignmentFileReader {

    protected File inFile;
    protected double totalHits;
    protected double totalWeight;
    protected Genome gen;
    protected int misMatch;
    protected boolean useNonUnique;
    protected int numChroms;
    protected HashMap<String, Integer> chrom2ID=new HashMap<String,Integer>();
    protected HashMap<Integer,String> id2Chrom=new HashMap<Integer,String>();
    
    //Data structures for pre-loading
    
    /**
     * Chromosomes are stored as IDs for efficiency (saving memory) 
     */
    //protected int[] chrs=null;
    
    /**
     * The first of the paired reads. <br>
     * First dimension represents the corresponding chromosome ID. <br>
     * Second dimension represents the strand. 0 for '+', 1 for '-' <br>
     * Third dimension contains the starts of the hits
     */
    protected int[][][] readOne=null;
    
    protected ArrayList<Integer>[][] readOneList = null;
    
    /**
     * The second of the paired reads. <br>
     * First dimension represents the corresponding chromosome ID. <br>
     * Second dimension represents the strand. 0 for '+', 1 for '-' <br>
     * Third dimension contains the ends of the hits
     */
    protected int[][][] readTwo=null;
    
    protected ArrayList<Integer>[][] readTwoList = null;
    
    /**
     * Number of hits that corresponds to the <tt>Read</tt> of each <tt>ReadHit</tt><br>
     * Should be the same between paired ends? <br>
     * First dimension represents the corresponding chromosome ID. <br>
     * Second dimension represents the strand. 0 for '+', 1 for '-' <br>
     * Third dimension contains the number of hits of the <tt>Read</tt> of each <tt>ReadHit</tt> 
     */
    protected float[][][] hitCounts=null;
    
    protected ArrayList<Float>[][] hitCountsList = null;
    
    /**
     * The IDs of the read hits<br>
     * Currently just ID of the first Read <br>
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
    
    protected int readLengthOne;
    protected int readLengthTwo;
    protected int insertLength;
    protected int currID=0;
    
    
    /**
     * This constructor use the file to initialize genome information
     * @param f file    
     * @param g genome
     */
    public PairedAlignmentFileReader(File f){
        estimateGenome(f);
    }

    
    public PairedAlignmentFileReader(File f, Genome g, int mis, boolean nonUnique, int idSeed,
            HashMap<String, Integer> chrom2ID, HashMap<Integer,String> id2Chrom){
        gen=g;
        totalHits=0;
        totalWeight=0;
        inFile=f;
        misMatch=mis;
        useNonUnique = nonUnique;
        currID=idSeed;
        numChroms = gen.getChromList().size();
        
        this.chrom2ID=chrom2ID;
        this.id2Chrom=id2Chrom;
        
        System.out.print("    Loading reads from: "+f.getName()+" ... ");
        
        //Initialize the data structures
        readOne   = new int[numChroms][2][];
        readTwo   = new int[numChroms][2][];
        hitCounts = new float[numChroms][2][];
        hitIDs    = new int[numChroms][2][];
        
        readOneList    = new ArrayList[numChroms][2];
        for(int i = 0; i < readOneList.length; i++) { for(int j = 0; j < readOneList[i].length; j++) { readOneList[i][j] = new ArrayList<Integer>(); } }

        readTwoList    = new ArrayList[numChroms][2];
        for(int i = 0; i < readTwoList.length; i++) { for(int j = 0; j < readTwoList[i].length; j++) { readTwoList[i][j] = new ArrayList<Integer>(); } }

        hitCountsList = new ArrayList[numChroms][2];
        for(int i = 0; i < hitCountsList.length; i++) { for(int j = 0; j < hitCountsList[i].length; j++) { hitCountsList[i][j] = new ArrayList<Float>(); } }
        
        hitIDsList    = new ArrayList[numChroms][2];
        for(int i = 0; i < hitIDsList.length; i++) { for(int j = 0; j < hitIDsList[i].length; j++) { hitIDsList[i][j] = new ArrayList<Integer>(); } }
        
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
                int[] tempStarts = readOne[chrID][j];
                
                if(tempStarts.length != 0) {
                    int start_ind = Arrays.binarySearch(tempStarts, coords.getStart());
                    int end_ind   = Arrays.binarySearch(tempStarts, coords.getEnd());
                    
                    if( start_ind < 0 ) { start_ind = -start_ind - 1; }
                    if( end_ind < 0 )   { end_ind   = -end_ind - 1; }
                    
                    start_ind = StatUtil.searchFrom(tempStarts, ">=", coords.getStart(), start_ind); //read closest to start site strictly greater than start
                    end_ind   = StatUtil.searchFrom(tempStarts, "<=",   coords.getEnd(), end_ind); //read closest to end strictly less than end
                    
                    char str = j == 0 ? '+' : '-';          
                    for(int k = start_ind; k <= end_ind; k++) {
                        double weight = 1/(double)hitCounts[chrID][j][k];
                        if (j==0)
                            hits.add(new ReadHit(gen, hitIDs[chrID][j][k], chr, tempStarts[k], tempStarts[k]+readLengthOne-1, str, weight ));
                        else
                            hits.add(new ReadHit(gen, hitIDs[chrID][j][k], chr, tempStarts[k]-readLengthOne+1, tempStarts[k], str, weight ));
                    }   
                }                       
            }
        }
        return hits;
    }//end of loadHits method

    // load paired read hit 5' coordinates (sorted, can be duplicates) and counts 
    public Pair<Pair<ArrayList<Integer>, ArrayList<ArrayList<Integer>>>,ArrayList<ArrayList<Float>>> loadStrandedFivePrimeCounts(Region coords , char strand){
        ArrayList<Integer> fiveHits = new ArrayList<Integer>();
        ArrayList<ArrayList<Integer>> threeHits = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<Float>> counts = new ArrayList<ArrayList<Float>>();
        String chr = coords.getChrom();
        if(chrom2ID.containsKey(chr)){
            int chrID = chrom2ID.get(chr);
            int j = strand == '+'?0:1;  // j = 0 <=> strand = '+', j = 1 <=> strand = '-'
            int[] tmp5Primes = readOne[chrID][j];    //assumes tempStarts is sorted  
            int[] tmp3Primes = readTwo[chrID][j];
            if(tmp5Primes.length != 0) {
                int start_ind = Arrays.binarySearch(tmp5Primes, coords.getStart());
                int end_ind   = Arrays.binarySearch(tmp5Primes, coords.getEnd());
                if( start_ind < 0 ) { start_ind = -start_ind - 1; }
                if( end_ind < 0 )   { end_ind   = -end_ind - 1; }
                start_ind = StatUtil.searchFrom(tmp5Primes, ">=", coords.getStart(), start_ind);
                end_ind   = StatUtil.searchFrom(tmp5Primes, "<=",   coords.getEnd(), end_ind);
                  
                //initializing the matrix
                int prevHit = tmp5Primes[start_ind];
                ArrayList<Integer> tempHits = new ArrayList<Integer>();
                tempHits.add(tmp3Primes[start_ind]);
                ArrayList<Float> tempCounts = new ArrayList<Float>();
                tempCounts.add(hitCounts[chrID][j][start_ind]);
                for(int k = start_ind+1; k <= end_ind; k++) {
                    if (prevHit != tmp5Primes[k]) {
                        fiveHits.add(prevHit);
                        threeHits.add(tempHits);//finish up row in matrix
                        counts.add(tempCounts);
                        
                        tempHits = new ArrayList<Integer>();//start new row in matrix
                        tempCounts = new ArrayList<Float>();
                        prevHit = tmp5Primes[k];
                    }
                    tempHits.add(tmp3Primes[k]);
                    tempCounts.add(hitCounts[chrID][j][k]);
                
                    //fiveHits.add(tmp5Primes[k]);
                    //counts.add(hitCounts[chrID][j][k]);
                }
                if (tempHits.size() != 0) {
                    fiveHits.add(prevHit);
                   threeHits.add(tempHits);
                   counts.add(tempCounts);
                   
                }
            }
        }
        Pair<ArrayList<Integer>, ArrayList<ArrayList<Integer>>> hitMatrix = 
                new Pair<ArrayList<Integer>, ArrayList<ArrayList<Integer>>>(fiveHits, threeHits);
        return new Pair<Pair<ArrayList<Integer>, ArrayList<ArrayList<Integer>>>,ArrayList<ArrayList<Float>>>(hitMatrix, counts);
    }
    
    /**
     * Count total reads & weight
     */
    protected abstract void countReads();
        
    /**
     * Estimate a fake genome from the observed read coordinates
     */
    protected abstract void estimateGenome(File f);
    
    /**
     * Gets the stranded weight of all hits (of all chromosomes) for the specified strand
     * @param strand 
     * @return
     */
    protected double getStrandedTotalWeight(char strand) {
        int strandInd = strand == '+' ? 0 : 1;
        double weight = 0;
        for(int i = 0; i < hitCounts.length; i++) {
            float[] hitCountsTemp = hitCounts[i][strandInd];
            for(Float el:hitCountsTemp)
                weight += 1/(double)el;
        }
        
        return weight;
    }//end of getStrandedTotalWeight method
    
    
    //Add hits to data structure
    //assumes same number of hits in each read
    protected void addHits(Read r1, Read r2){
        for(int i = 0; i < r1.getHits().size(); i++){
            ReadHit h = r1.getHits().get(i);
            ReadHit h2 = r2.getHits().get(i);
            if (!chrom2ID.containsKey(h.getChrom()))
                continue;
            int chrID   = chrom2ID.get(h.getChrom());
            char strand = h.getStrand();
            int strandInd = strand == '+' ? 0 : 1;
            //System.out.println(h.getChrom()+"\t"+h.getStart()+"\t"+h.getStrand());
            readOneList[chrID][strandInd].add(strand == '+' ?h.getStart():h.getEnd());
            readTwoList[chrID][strandInd].add(strand == '+' ?h2.getStart():h2.getEnd());
            hitIDsList[chrID][strandInd].add(h.getID());
            hitCountsList[chrID][strandInd].add((float)h.getWeight());
            totalHits++;
        }
    }//end of addHits method
    
    /**
     * Converts lists of integers to integer arrays, deletes the lists for saving
     * memory and permutes all array elements so that they are ordered in terms of
     * the array <tt>fivePrimes</tt>.
     */
    protected void populateArrays() {
        
        for(int i = 0; i < readOneList.length; i++)
            for(int j = 0; j < readOneList[i].length; j++)
                readOne[i][j] = list2int(readOneList[i][j]);
        
        for(int i = 0; i < readOneList.length; i++) { for(int j = 0; j < readOneList[i].length; j++) { readOneList[i][j].clear(); } }
        readOneList = null;

        ////////////////////////////////////////////////////////////////////////

        for(int i = 0; i < readTwoList.length; i++)
            for(int j = 0; j < readTwoList[i].length; j++)
                readTwo[i][j] = list2int(readTwoList[i][j]);
        
        for(int i = 0; i < readTwoList.length; i++) { for(int j = 0; j < readTwoList[i].length; j++) { readTwoList[i][j].clear(); } }
        readTwoList = null;
        
        ////////////////////////////////////////////////////////////////////////

        for(int i = 0; i < hitIDsList.length; i++)
            for(int j = 0; j < hitIDsList[i].length; j++)
                hitIDs[i][j] = list2int(hitIDsList[i][j]);
        
        for(int i = 0; i < hitIDsList.length; i++) { for(int j = 0; j < hitIDsList[i].length; j++) { hitIDsList[i][j].clear(); } }
        hitIDsList = null;
        
        ////////////////////////////////////////////////////////////////////////

        for(int i = 0; i < hitCountsList.length; i++)
            for(int j = 0; j < hitCountsList[i].length; j++)
                hitCounts[i][j] = list2Float(hitCountsList[i][j]);
        
        for(int i = 0; i < hitCountsList.length; i++) { for(int j = 0; j < hitCountsList[i].length; j++) { hitCountsList[i][j].clear(); } }
        hitCountsList = null;
        
        ////////////////////////////////////////////////////////////////////////
        
        for(int i = 0; i < readOne.length; i++) {  // chr
            for(int j = 0; j < readOne[i].length; j++) { // strand
                int[] inds = StatUtil.findSort(readOne[i][j]);
                readTwo[i][j] = StatUtil.permute(readTwo[i][j], inds);
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
    public int getReadLengthOne(){return readLengthOne;}
    
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
     * get all the readOne positions
     * should modify to ensure it's 5 prime?
     * clean the memory space
     * @return
     */
    public int[][][] getreadOne(){
        hitIDs = null;
        hitCounts = null;
        System.gc();
        return readOne;
    }
    

    private int[] list2int(List<Integer> list) {
        int[] out = new int[list.size()];
        for(int i = 0; i < out.length; i++)
            out[i] = list.get(i);
       
        return out;
    }//end of list2int method

    private float[] list2Float(List<Float> list) {
        float[] out = new float[list.size()];
        for(int i = 0; i < out.length; i++)
            out[i] = list.get(i);
       
        return out;
    }//end of list2int method
    
    public void cleanup(){
        readOne=null;
        readTwo=null;
        hitCounts=null;
        hitIDs=null;
    }
    
    /**
     * @param args
     */
    public static void main(String[] args) {
        // TODO Auto-generated method stub

    }

}

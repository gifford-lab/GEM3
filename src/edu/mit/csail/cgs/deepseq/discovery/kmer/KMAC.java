package edu.mit.csail.cgs.deepseq.discovery.kmer;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;
import java.util.TreeMap;
import java.util.TreeSet;

import javax.imageio.ImageIO;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.analysis.MotifInstance;
import edu.mit.csail.cgs.deepseq.analysis.MotifScan;
import edu.mit.csail.cgs.deepseq.discovery.Config;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KMAC0.KmerCluster;
import edu.mit.csail.cgs.deepseq.discovery.kmer.mtree.MTree;
import edu.mit.csail.cgs.deepseq.discovery.kmer.mtree.MTree.MTreeNode;
import edu.mit.csail.cgs.deepseq.discovery.kmer.mtree.MTree.TreeObject;
import edu.mit.csail.cgs.deepseq.features.ComponentFeature;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.utils.stats.ROC;
import edu.mit.csail.cgs.utils.stats.StatUtil;
import edu.mit.csail.cgs.utils.stats.StatUtil.DensityClusteringPoint;
import edu.mit.csail.cgs.utils.strings.multipattern.AhoCorasick;
import edu.mit.csail.cgs.utils.strings.multipattern.SearchResult;

public class KMAC {
	public final static String KMAC_VERSION = "1.0";

	public static final int RC=100000;		// extra bp add to indicate negative strand match of kmer
	private static final int UNALIGNED=9999;	// the special shift for unaligned kmer
	public static final char[] LETTERS = {'A','C','G','T'};
	public static final int MAXLETTERVAL = Math.max(Math.max(Math.max('A','C'),Math.max('T','G')),
            Math.max(Math.max('a','c'),Math.max('t','g'))) + 1;
	
	private boolean standalone = false;
	public void setStandalone(){
		standalone = true;
	}
	Config config = new Config();	

	private boolean engineInitialized =false;
	private double[] bg= new double[4];	// background frequency based on GC content
	private String outName, params;
	private boolean isDebugging = false;
	public void setIsDebugging(){isDebugging = true;}
	/** region-->index_seqsNeg for negative sequences */
	private TreeMap<Region, Integer> neg_region_map;

	private double[] profile;
	
	private ArrayList<Sequence> seqList;
	private ArrayList<Sequence> seqListNeg;
	private String[] seqs;			// DNA sequences around binding sites
	private double[] seq_weights;	// sequence hit count weighted by binding strength, corresponding to the seqs[]
	public double[] getSequenceWeights(){return seq_weights;}
	public void setSequenceWeights(double[] w){
		seq_weights=w;
	}

	private String[] seqsNeg;		// DNA sequences in negative sets
	private ArrayList<String> seqsNegList=new ArrayList<String>(); // Effective negative sets, excluding overlaps in positive sets
	public int getNegSeqCount(){return negSeqCount;}
	public int getPosSeqCount(){return posSeqCount;}
    private int posSeqCount;
    private int negSeqCount;
    private double posNegSeqRatio;
    public void setTotalSeqCount(int pos, int neg){
    	posSeqCount = pos;
    	negSeqCount = neg;
	    posNegSeqRatio = posSeqCount/(double)negSeqCount;
    }
	private int negRegionDistance;
	public String[] getPositiveSeqs(){return seqs;};
	public double get_NP_ratio(){return (double)negSeqCount/posSeqCount;}

	/** the info about sequence hit, total k-mer coverage */
	int[] posCoveredWidth;
	int[] negCoveredWidth;
	String[][] posHitStrings;		// the matched sequences of the best hit of current motif in each training sequence
	String[][] negHitStrings;
	public void setCoveredWidth(int[] posCoveredWidth, int[] negCoveredWidth){
		this.posCoveredWidth = posCoveredWidth;
		this.negCoveredWidth = negCoveredWidth;
	}
	public void setHitStrings(String[][] posHitStrings, String[][] negHitStrings){
		this.posHitStrings = posHitStrings;
		this.negHitStrings = negHitStrings;
	}
	private double[] pseudoCountRatios;
	
	private HashMap<Integer, HashMap<String, Kmer>> allKmerMap = null;
	
	/** AhoCorasick algorithm tree object for multi-pattern search<br>
	 if run in parallel, this instance variable should be move to the thread local environment<br>
	 initAhoCorasick(): Pre-processing is to build the tree with all the patterns (kmers). <br>
	 findKSMGroupHits(): individual search can be done after init*/
	private AhoCorasick treeAhoCorasick;
	/** For mapping k-mer string from AhoCorasick tree to Kmer Object<br>One string may map to multiple gapped k-mers */
	private HashMap<String, Kmer[]> str2kmers = new HashMap<String, Kmer[]>();	
	/** Paired with str2kmers HashMap, one string may have different offsets in different gapped k-mers <br>
	 * Here the offset is relative to the end of the k-mer, called EndOffset, 
	 * which is computed in initAhoCorasick() method, then used in findKSMGroupHits() to compute the KSM match position */
	private HashMap<String, int[]> str2kmerEndOffsets = new HashMap<String, int[]>();
	
	public boolean isInitialized(){ return engineInitialized;}
	
	private SequenceGenerator<Region> seqgen;
	private long tic;	

	/** motif clusters */
	ArrayList<MotifCluster> clusters = new ArrayList<MotifCluster>();
	public MotifCluster getPrimaryCluster(){
		if (clusters.size()>=1)
			return clusters.get(0);
		return null;
	}
	public ArrayList<MotifCluster> getMotifClusters(){
		return clusters;
	}	
	public void setTotalSeqCounts(int posSeqCount, int negSeqCount){
		this.posSeqCount = posSeqCount;
		this.negSeqCount = negSeqCount;
	}
	
	public void setConfig(Config config, String outName, String params){
		this.config = config;
	    this.outName = outName;
	    this.params = params;
	    Kmer.set_use_weighted_hit_count(config.use_weighted_kmer);
	}
	
	// called by standalone main() method
	public void setSequences(ArrayList<String> pos_seqs, ArrayList<String> neg_seqs, ArrayList<Double> pos_w){
		if (config.k_seqs==-1)
			config.k_seqs = pos_seqs.size();
		int seqNum = Math.min(pos_seqs.size(), config.k_seqs);
		seqs = new String[seqNum];	
		for (int i=0;i<seqNum;i++){
			String seq = pos_seqs.get(i);
			if (config.repeat_fraction<1){
				int count = 0;
				for (char c:seq.toCharArray())
					if (Character.isLowerCase(c) || c=='N')				// assuming lower case sequences are repeats
						count++;
				if (count>seq.length()*config.repeat_fraction)			// if repeat fraction in sequence is too high, skip
					continue;
				if (count>1){									// convert lower case repeat to N
					char[] chars = seq.toCharArray();
					for (int j=0;j<chars.length;j++)
						if (Character.isLowerCase(chars[j]))
							chars[j] = 'N';
					seq = new String(chars);
				}
			}
			// if repeat_fraction>=1, allow all repeats, convert to upper case
			seqs[i] = seq.toUpperCase();
		}
		
		seqNum = seqs.length;
		seq_weights = new double[seqNum];
		for (int i=0;i<seqNum;i++){
			switch (config.seq_weight_type){
				case -1:	seq_weights[i]=1/pos_w.get(i);break;
				case 0:	seq_weights[i]=1.0;break;
				case 1: seq_weights[i]=pos_w.get(i);break;
				case 2: seq_weights[i]=Math.sqrt(pos_w.get(i));break;
				case 3: seq_weights[i]=Math.log(pos_w.get(i));
					if(seq_weights[i]<=0) 
						System.err.println("Non-positive sequence weight after log():"+seq_weights[i]);
					break;
				default: System.err.println("Sequence weighting type is not defined!");System.exit(-1);
			}
		}
		StatUtil.mutate_normalize(seq_weights);
		for (int i=0;i<seqNum;i++){
			seq_weights[i] = seq_weights[i]*seqNum;	// scale weights such that average weight = 1
		}
		
		// If neg seqs are not provided, use shuffled sequences as negative sequences
		if (neg_seqs.isEmpty()){
			if (config.k_neg_dinu_shuffle){
				System.out.println("Use di-nucleotide shuffled sequences as negative sequences.\n");
//				// append all sequences, then shuffle
//				StringBuilder sb = new StringBuilder();
//				for (int i=0;i<seqNum;i++){
//					sb.append(seqs[i]);
//				}
//				for (int j=0;j<config.neg_pos_ratio;j++){
//					Random randObj = new Random(config.rand_seed+j);
//					String shuffled = SequenceUtils.dinu_shuffle(sb.toString(), randObj);
//					int start = 0;
//					for (int i=0;i<seqNum;i++){
//						seqsNegList.add(shuffled.substring(start, start+seqs[i].length()));
//						start+=seqs[i].length();
//					}
//				}
				// shuffle each sequence
				for (int j=0;j<config.neg_pos_ratio;j++){		
					Random randObj = new Random(config.rand_seed+j);
					for (int i=0;i<seqNum;i++){
						seqsNegList.add(SequenceUtils.dinu_shuffle(seqs[i], randObj));
					}
				}
			}
			else{		// single nucleotide shuffling
				System.out.println("Use shuffled sequences as negative sequences.\n");
				for (int j=0;j<config.neg_pos_ratio;j++){		
					Random randObj = new Random(config.rand_seed+j);
					for (int i=0;i<seqNum;i++){
						seqsNegList.add(SequenceUtils.shuffle(seqs[i], randObj));
					}
				}
			}
		}
		else{	// use the supplied negative sequences
			int negSeqNum = seqNum * config.neg_pos_ratio;
			HashSet<Integer> ids = new HashSet<Integer>();
			Random randObj = new Random(config.rand_seed);
			for (int i=0;i<negSeqNum;i++){
				ids.add(randObj.nextInt(neg_seqs.size()));
			}
			int makeup = (negSeqNum-ids.size())*2;
			for (int i=0;i<makeup;i++){
				ids.add(randObj.nextInt(neg_seqs.size()));
				if (ids.size()==negSeqNum)
					break;
			}
			
			for (int i:ids){
				String seq = neg_seqs.get(i);
				seqsNegList.add(seq.toUpperCase());
			}
		}
		posSeqCount = seqs.length;
	    negSeqCount = seqsNegList.size();
	    posNegSeqRatio = posSeqCount/(double)negSeqCount;
	    updateSequenceInfo();
	}
	
	private void updateSequenceInfo(){
		
		if (config.use_pos_weight){
			int seqLen = seqs[0].length();
			
			// logistic distribution to fit the spatial resolution shape, with a more heavy tail than Gaussian
			// http://en.wikipedia.org/wiki/Logistic_distribution
			// ctcf_sigma = 9.53; GABP_sigma = 15.98;
		    profile = new double[seqLen+1];
		    double sigma = 13;
		    for (int i=0; i<=seqLen/2; i++){
		    	double e = Math.exp(-i/sigma);
		    	profile[seqLen/2-i] = e/(sigma*(1+e)*(1+e));
		    	profile[seqLen/2+i] = profile[seqLen/2-i];
		    }
		    StatUtil.normalize(profile);
	//	   	System.out.println(CommonUtils.arrayToString(profile, "%.4f"));
		}
		
	    // count cg-content
		int gcCount = 0;
		int negLength = 0;
		for (String seq:seqsNegList){
			negLength += seq.length();
			for (char c:seq.toCharArray())
				if (c=='C'||c=='G')
					gcCount ++;
		}
		double gcRatio = (double)gcCount/negLength;
		System.out.println(String.format("Estimated GC content is %.2f. Set [--gc -1] to use the estimated GC content.", gcRatio));
		if (config.gc>0){
			System.out.println(String.format("Provided  GC content is %.2f.", config.gc));
			bg[0]=(1-config.gc)/2; 
	    	bg[1]=config.gc/2; 
	    	bg[2]=bg[1]; 
	    	bg[3]=bg[0];
		}
		else{
			bg[0]=(1-gcRatio)/2; 
	    	bg[1]=gcRatio/2; 
	    	bg[2]=bg[1]; 
	    	bg[3]=bg[0];
		}
		System.out.println();
	}
	/**
	 * Stand-alone KMAC constructor
	 */
	public KMAC(){
	}

	public KMAC(Genome g, boolean useCache, boolean use_db_genome, String genomePath){
//		setUseKmerWeight();

		seqgen = new SequenceGenerator<Region>();
    	if (use_db_genome)
    		seqgen.useLocalFiles(false);
		if (useCache)
			seqgen.useCache(true);	
		seqgen.setGenomePath(genomePath);
	}
	
	/** 
	 * Contruct a Kmer Engine from a list of Kmers.<br>
	 * This is used for KSM motif scanning outside of KMAC.
	 */
	public KMAC(ArrayList<Kmer> kmers, Config config){
		if (config!=null)
			this.config = config;
		if (!kmers.isEmpty()){
			initAhoCorasick(kmers);
		}
	}
	/**
	 * Set up the light weight genome cache. Only load the sequences for the specified regions.<br>
	 * At the same time, retrieve negative sequences (for only once, no caching)
	 * @param regions
	 */
	public double setupRegionCache(ArrayList<Region> cacheRegions, ArrayList<Region> negativeRegions, int negRegionDistance){
		this.negRegionDistance = negRegionDistance;
		double gcRatio=0;
		if (!seqgen.isRegionCached()){
			seqsNeg = seqgen.setupRegionCache_new(cacheRegions, negativeRegions);
			neg_region_map = new TreeMap<Region, Integer>();
			for (int i=0;i<negativeRegions.size();i++){
				neg_region_map.put(negativeRegions.get(i), i);
			}
			// count cg-content
			int gcCount = 0;
			for (String s:seqsNeg){
				for (char c:s.toCharArray())
					if (c=='C'||c=='G')
						gcCount ++;
			}
			gcRatio = (double)gcCount/seqsNeg.length/seqsNeg[0].length();
			bg[0]=0.5-gcRatio/2; 
        	bg[1]=gcRatio/2; 
        	bg[2]=bg[1]; 
        	bg[3]=bg[0];
		}
		return gcRatio;
	}

	public void updateOutPrefix(String outPrefix){
	    this.outName = outPrefix;
	}

	/**
	 * Load pos/neg test sequences based on event positions, run from GEM<br>
	 * Skip repeat masked sequences according to config.repeat_fraction, otherwise convert repeat characters into 'N'
	 * 
	 * @param events
	 * @param winSize
	 * @param winShift
	 */
	public void loadTestSequences(ArrayList<ComponentFeature> events, int winSize){
	
		int eventCount = events.size();
		ArrayList<Region> posImpactRegion = new ArrayList<Region>();			// to make sure negative region is not within negRegionDistance of positive regions.
		ArrayList<String> posSeqs = new ArrayList<String>();
		ArrayList<Double> posSeqWeights = new ArrayList<Double>();
		for(int i=0;i<eventCount;i++){
			Region posRegion = events.get(i).getPeak().expand(winSize/2);
			String seq = seqgen.execute(posRegion);
			if (config.repeat_fraction<1){
				int count = 0;
				for (char c:seq.toCharArray())
					if (Character.isLowerCase(c) || c=='N')				// assuming lower case sequences are repeats
						count++;
				if (count>seq.length()*config.repeat_fraction)			// if repeat fraction in sequence is too high, skip
					continue;
				if (count>1){									// convert lower case repeat to N
					char[] chars = seq.toCharArray();
					for (int j=0;j<chars.length;j++)
						if (Character.isLowerCase(chars[j]))
							chars[j] = 'N';
					seq = new String(chars);
				}
			}
			posSeqs.add(seq.toUpperCase());		// if repeat_fraction>=1, allow all repeats, convert to upper case
			switch (config.seq_weight_type){
			case 0:	posSeqWeights.add(10.0);break;
			case 1: posSeqWeights.add(events.get(i).getTotalEventStrength());break;
			case 2: posSeqWeights.add(Math.sqrt(events.get(i).getTotalEventStrength()));break;
			case 3: posSeqWeights.add(Math.log(events.get(i).getTotalEventStrength()));break;
			default: System.err.println("Sequence weighting type is not defined! No weighting.");posSeqWeights.add(1.0);
			}
			posImpactRegion.add(events.get(i).getPeak().expand(negRegionDistance));
		}
		seqs = new String[posSeqs.size()];	// DNA sequences around binding sites
		posSeqs.toArray(seqs);
		
		int seqNum = seqs.length;
	    seq_weights = new double[seqNum];
		for (int i=0;i<seqNum;i++)
			seq_weights[i]=posSeqWeights.get(i);
		StatUtil.mutate_normalize(seq_weights);
		for (int i=0;i<seqNum;i++)
			seq_weights[i] = seq_weights[i]*seqNum;	// scale weights with total sequence count, and total weight
		
		seqsNegList.clear();
		if (config.k_neg_dinu_shuffle){
			System.out.println("Use di-nucleotide shuffled sequences as negative sequences.\n");
			Random randObj = new Random(config.rand_seed);
			for (int i=0;i<seqs.length;i++)
				seqsNegList.add(SequenceUtils.dinu_shuffle(seqs[i], randObj));
		}

		//		else{
//			/** Negative sequences has been retrieved when setting up region caches */
//			ArrayList<Region> negRegions = new ArrayList<Region>();
//			negRegions.addAll(neg_region_map.keySet());
//			Region.filterOverlapRegions(negRegions, posImpactRegion);	// make sure negative region is not within negRegionDistance of positive regions.
//			int negCount = 0;
//			int len = winSize/2*2+1;
//			for (Region r:negRegions){
//				String seq_retrieved = seqsNeg[neg_region_map.get(r)];
//				if (seq_retrieved.length()<len)
//					continue;
//				String seq = seq_retrieved.substring(0, len);
//				if (config.repeat_fraction<1){
//					int count = 0;
//					for (char c:seq.toCharArray())
//						if (Character.isLowerCase(c) || c=='N')
//							count++;
//					if (count>seq.length()*config.repeat_fraction)			// if repeat fraction in sequence is too high, skip
//						continue;
//					if (count>1){									// convert lower case repeat to N
//						char[] chars = seq.toCharArray();
//						for (int j=0;j<chars.length;j++)
//							if (Character.isLowerCase(chars[j]))
//								chars[j] = 'N';
//						seq = new String(chars);
//					}
//				}
//				seqsNegList.add(seq.toUpperCase());		// if repeat_fraction>=1, allow repeats, convert to upper case
//				negCount++;
//				if (negCount==seqs.length)				// limit the neg region count to be same or less than positive region count
//					break;
//			}
//		}
		posSeqCount = seqs.length;
		seqsNegList.trimToSize();
	    negSeqCount = seqsNegList.size();
	    if (config.verbose>1 || config.repeat_fraction!=1)
			System.out.println(String.format("From %d events, loaded %d positive sequences, skipped %d(%.1f%%) repeat sequences\n", 
					events.size(), posSeqCount, events.size()-posSeqCount, 100-100.0*posSeqCount/events.size()));
	    
	    updateSequenceInfo();			
		    
	//		cern.jet.random.engine.RandomEngine randomEngine = new cern.jet.random.engine.MersenneTwister();
	//		ArrayList<String> negSeqList = new ArrayList<String>();
	//		for(int i=0;i<eventCount;i++){
	//			// getting negative sequences
	//			// exclude negative regions that overlap with positive regions, or exceed start of chrom
	//			// it is OK if we lose a few sequences here
	//			Region posRegion = seqCoors[i];
	//			int start = 0;
	//			double rand = randomEngine.nextDouble();
	//			if (rand>0.5)
	//				start = (int) (posRegion.getEnd()+1 + winShift*rand);
	//			else
	//				start =(int) (posRegion.getStart()-1 - winShift*(1-rand));
	//			int end = start + posRegion.getWidth()-1;			// end inclusive
	//			if (start < 0 || end >= genome.getChromLength(posRegion.getChrom()))
	//				continue;
	//			Region negRegion = new Region(genome, posRegion.getChrom(), start, end);			
	//			if (i>0 && seqCoors[i-1].overlaps(negRegion))
	//				continue;
	//			if (i<(eventCount-2) && seqCoors[i+1].overlaps(negRegion))
	//				continue;
	//			negSeqList.add(seqgen.execute(negRegion).toUpperCase());
	//		}
	//		seqsNeg = new String[negSeqList.size()];	// DNA sequences in negative sets
	//		for (int i=0; i<seqsNeg.length;i++){
	//			seqsNeg[i] = negSeqList.get(i);
	//		}
	}
	
	public void printInputSequences(String outName){
	    StringBuilder sb = new StringBuilder();
	    for (int i=0;i<seqs.length;i++)
	    	sb.append(String.format("%s\t%.2f\n", seqs[i], seq_weights[i]));
	    CommonUtils.writeFile(outName+"_pos_seqsw.txt", sb.toString());
	    sb = new StringBuilder();
	    for (int i=0;i<seqsNegList.size();i++)
	    	sb.append(String.format("%s\n", seqsNegList.get(i)));
	    CommonUtils.writeFile(outName+"_neg_seqs.txt", sb.toString());
	}
	
	/** 
	 * Run through a range of k values to discovery motifs
	 */
	public void discoverMotifs (int k_min, int k_max, int[] eventCounts){
		
		ArrayList<MotifCluster> allClusters = new ArrayList<MotifCluster>();		
		cern.jet.random.engine.RandomEngine randomEngine = new cern.jet.random.engine.MersenneTwister();
		
		if (seqs.length<500){	// relax cutoff if not many sequences
			config.pwm_noise = 0.1;
			config.kmer_hgp = -1.3;
		}
		
		/** Initialize the sequence objects */
		seqList = new ArrayList<Sequence>();
		for (int i=0;i<seqs.length;i++){
			Sequence s = new Sequence(seqs[i], i);
			seqList.add(s);
		}
		seqList.trimToSize();
		seqListNeg = new ArrayList<Sequence>();
		for (int i=0;i<seqsNegList.size();i++){
			Sequence s = new Sequence(seqsNegList.get(i), i);
			seqListNeg.add(s);
		}
		seqListNeg.trimToSize();
				
		/**
		 * For each k, generate exact k-mers and gapped kmers, density clustering, KMAC
		 */
		allKmerMap = new HashMap<Integer, HashMap<String, Kmer>>();
		this.pseudoCountRatios = new double[k_max+1];
		for (int k=k_min;k<=k_max;k++){
			if (config.verbose>1)
				System.out.println("\nmemory used = "+
					(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1048576  +"M");
			pseudoCountRatios[k] = Math.max((seqs[0].length()-k+1)*2 / Math.pow(4, k), 0.005);
			StringBuilder sb = new StringBuilder();
			System.out.println("\n----------------------------------------------------------\nRunning k="+k+" ...\n");
			Pair<ArrayList<Kmer>, ArrayList<Kmer>> pair = selectEnrichedKmers(k);
			ArrayList<Kmer> allSignificantKmers = pair.car();
			ArrayList<Kmer> kmers = pair.cdr();
		
			if (kmers.isEmpty()){
				System.out.println("\nNo enriched k-mer!");
				continue;
			}
		
//			COMMENT block start, to SKIP kmac, only get k-mers

			kmers.trimToSize();
			Collections.sort(kmers);
			
			// setup the matrix form of gapped k-mer representation, for distance calculation
			for (Kmer km:kmers){
				km.setMatrix();
			}
			
			System.out.println("\n------------------------- k = "+ k +" ----------------------------\n");
			System.out.println("Total number of k-mers: n = "+kmers.size());

			ArrayList<Kmer> centerKmers = new ArrayList<Kmer> ();
			if (config.mtree==-1){
				// printout m-tree performance information
				System.gc();
				long mem = (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1048576;
				System.out.println("n^2 / 2 = "+kmers.size()*(kmers.size()-1)/2 + "\tmem="+mem+"M");
				for (int d=3;d<=4;d++){
					for (int c=8;c<=8;c+=2){
//					for (int c=8;c<=20;c+=2){
//					for (int c=5;c<=55;c+=10){
						System.gc();
						long tic=System.currentTimeMillis();
						numDistCalcuation = 0;
						System.out.print("d="+d+"\t c="+c+"\t");
						MTree dataPoints = MTree.constructTree(kmers, c);
						System.out.print("n_tree="+numDistCalcuation+"\t");
						centerKmers = densityClusteringWithMTree(kmers, dataPoints, d);
						System.out.print("\ttime="+CommonUtils.timeElapsed(tic));
						System.out.println("\tmem="+
						((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1048576 - mem) +"M");					
					}
					System.out.println();
				}
				
				// printout distance matrix performance information
				System.gc();
				long tic=System.currentTimeMillis();
				float[][] distanceMatrix = computeWeightedDistanceMatrix2(kmers, config.print_dist_matrix);
				centerKmers = densityClusteringWithDistMatrix(kmers, distanceMatrix, config.dc>0?config.dc:k/3.0);
				System.out.print("Distance Matrix: time="+CommonUtils.timeElapsed(tic));
				System.out.println("\tmem="+
				((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1048576 - mem) +"M");					
				distanceMatrix = null;
				continue;
			}
			else if (config.mtree!=0 || kmers.size()>10000){		// Explicitly setting capacity or better for mtree
				long tic=System.currentTimeMillis();
				int capacity = 0;
				if (config.mtree==0)
					capacity = kmers.size()/10000 + 8;
				else
					capacity = config.mtree;
				if (config.verbose>1)
					System.out.println("Construct m-tree for k-mers, node capacity = "+capacity);		
				
				MTree dataPoints = MTree.constructTree(kmers, capacity);
				centerKmers = densityClusteringWithMTree(kmers, dataPoints, config.dc>0?config.dc:k/3.0);
				dataPoints = null;
				System.gc();
				if (config.verbose>1)
					System.out.println("Density clustering: " + CommonUtils.timeElapsed(tic));
			}
			else{
				long tic=System.currentTimeMillis();
				if (config.verbose>1)
					System.out.print("Computing k-mer pairwise distance matrix ...");
				
				float[][] distanceMatrix = computeWeightedDistanceMatrix2(kmers, config.print_dist_matrix);
				
				if (config.verbose>1)
					System.out.println(" OK: "+CommonUtils.timeElapsed(tic));
				
				centerKmers = densityClusteringWithDistMatrix(kmers, distanceMatrix, config.dc>0?config.dc:k/3.0);
				distanceMatrix = null;
				System.gc();
				if (config.verbose>1)
					System.out.println("Density clustering: " + CommonUtils.timeElapsed(tic));
			}
			
			// use all the significant k-mers to get center-kmer neighbors
			ArrayList<ArrayList<Kmer>> neighbourList = new ArrayList<ArrayList<Kmer>>();  // centerKmers and neighbourList are matched lists
			double cutoff = config.kmer_deviation_factor*k;	// maximum kmer distance to be considered as neighbors
			
			int numKmer = 5000;
	        for (int j=0;j<centerKmers.size();j++){	
	        	Kmer seedKmer = centerKmers.get(j);
				
				ArrayList<Kmer> tmp = new ArrayList<Kmer>();
				for (Kmer km: allSignificantKmers){		// get all the k-mers
					if (KMAC.editDistance(seedKmer, km) <= cutoff)
						tmp.add(km);
				}
				if (tmp.size()<numKmer){
					neighbourList.add(tmp);					
				}
				else{
					ArrayList<Kmer> neighbours = new ArrayList<Kmer>();
					for (Kmer km: kmers){		// get all the selected k-mers first
						if (KMAC.editDistance(seedKmer, km) <= cutoff)
							neighbours.add(km);
					}
					System.out.print(neighbours.size()+" --> ");
					
					// randomly sample some allSignificantKmers, upto total 5000 kmers
					for (int r=0;r<tmp.size();r++){
						int idx = (int)(randomEngine.nextDouble()*tmp.size());
						Kmer km = tmp.get(idx);
						if (!neighbours.contains(km)){
							neighbours.add(km);
							if (neighbours.size()>=numKmer)
								break;
						}
					}
					System.out.println(neighbours.size());
					neighbourList.add(neighbours);
				}
	        }
	        
			// clear matrix, it is only used for calculating distance
			for (Kmer km:allSignificantKmers){
				km.clearMatrix();
			}
			kmers = null;
			allSignificantKmers = null;
			pair = null;
			System.gc();
			System.out.println();
			
//			numKmerToTry=0; // skip KMAC step, used only for testing
			
	        ArrayList<MotifCluster> tmp = new ArrayList<MotifCluster>();
	        for (int j=0;j<centerKmers.size();j++){	

	        	Kmer seedKmer = centerKmers.get(j);
				ArrayList<Kmer> neighbours = neighbourList.get(j);
				neighbourList.set(j, null);	// clean up
				System.gc();
				if (config.verbose>1)
	    			System.out.println("\n------------------------------------------------------------");
				System.out.println("Aligning and clustering k-mers with seed "+seedKmer.kmerString+",   \t#"+j);
	    		if (config.verbose>1)
	    			System.out.println("\nmemory used = "+
						(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1048576  +"M");		

	        	MotifCluster c = KmerMotifAlignmentClustering(seqList, neighbours, seedKmer, k);

	        	if (c!=null && c.wm!=null){
	        		tmp.add(c);
		        	if (config.kg_hit_adjust_type==2)
		        		c.setCoveredWidth(posCoveredWidth, negCoveredWidth);	// save a copy in the cluster
		        	if (config.kg_hit_adjust_type==1)
		        		c.setHitStrings(posHitStrings, negHitStrings);
	        	}
	        	if (tmp.size()==config.k_top)
	        		break;
	        }
			sortMotifClusters(tmp, true);
			
//			
//			COMMENT block to SKIP kmac, only get k-mers, for testing density clustering.
			
			// print all the motifs for k, before merging
			sb.append("\n");
			printMotifClusters(tmp, sb);
			if (config.verbose>1)
				System.out.println(sb.toString());
			
			if (config.verbose>1)
				System.out.println("\nmemory used = "+
					(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1048576  +"M");		

			if (tmp.size()>1){
				System.out.println("\nMerging redundant motifs ...");
				boolean[][] checked = new boolean[tmp.size()][tmp.size()];	// a table indicating whether a pair has been checked
				
				tic = System.currentTimeMillis();
				mergeOverlapPwmMotifs (tmp, seqList, checked, config.use_ksm, 0);
				sb.append("After merging:\n");
				printMotifClusters(tmp, sb);
			}
			
			// print motifs for k, after merging
			System.out.println("\n------------------------- k = "+ k +" ----------------------------");
			System.out.println(sb.toString());
			
			allClusters.addAll(tmp);
		} // for each k
		
//		allK_allPatterns.clear(); allK_allPatterns = null;
		allKmerMap.clear();	allKmerMap = null;
		clusters = allClusters;
		allClusters = null;
		
		/** merge motifs from all K values */
		if (clusters.size()>1){
			
			sortMotifClusters(clusters, true);
			
			StringBuilder sb_all = new StringBuilder("\n------------------------- "+ new File(outName).getName() +", k = ALL -------------------------\n");
			if (config.verbose>1)
				System.out.println("\n---------------------------------------------\nMotifs from all k values:\n");
			printMotifClusters(clusters, sb_all);
			System.out.print(sb_all.toString());
			
			if (config.verbose>1)
				System.out.println("\nmemory used = "+
					(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1048576  +"M");		
			tic = System.currentTimeMillis();
			boolean[][] checked = new boolean[clusters.size()][clusters.size()];	// a table indicating whether a pair has been checked
			mergeOverlapPwmMotifs (clusters, seqList, checked, config.use_ksm, 0);
			
			if (config.verbose>1)
				System.out.println("\nFinished merging motifs.\n");
			sb_all.append("After merging:\n");
			printMotifClusters(clusters, sb_all);
			System.out.print(sb_all.toString());
		}
		
		if (config.verbose>1)
			System.out.println("\nmemory used = "+
				(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1048576  +"M");		

		/** Output final motifs, set binding positions, etc */
    	System.out.println(CommonUtils.timeElapsed(tic)+": Finalizing "+ clusters.size() +" motifs ...\n");
		// remove clusters with low hit count
		ArrayList<MotifCluster> toRemove = new ArrayList<MotifCluster>();
		for (int i=0;i<clusters.size();i++){
			MotifCluster c = clusters.get(i);
			double hitRatio = (double)c.pwmThreshold.posHit / posSeqCount;
			if (i>=10 && hitRatio<config.motif_hit_factor_report || hitRatio<config.motif_hit_factor)
					toRemove.add(c);
			//TODO: add a new config.motif_significance, or use fold change to remove motif here
//			if (config.evaluate_by_ksm && c.ksmThreshold.motif_significance>config.hgp)
//				toRemove.add(c);
//			if (!config.evaluate_by_ksm && c.pwmThreshold.motif_significance>config.hgp)
//				toRemove.add(c);
		}
		clusters.removeAll(toRemove);
		
		sortMotifClusters(clusters, true);

		for (int i=0;i<clusters.size();i++){
			MotifCluster cluster = clusters.get(i);
			
			/** use all aligned sequences to find expected binding sites, set kmer offset */
	    	// average all the binding positions to decide the expected binding position
			initAhoCorasick(cluster.alignedKmers);
			if (config.kg_hit_adjust_type==1){
				this.posHitStrings = cluster.posHitStrings.clone();
				this.negHitStrings = cluster.negHitStrings.clone();		
			}
			else if (config.kg_hit_adjust_type==2){
				this.posCoveredWidth = cluster.posCoveredWidth.clone();
				this.negCoveredWidth = cluster.negCoveredWidth.clone();
			}

			HashMap<Integer, KmerGroup> seq2kg = alignByKSM(seqList, cluster.alignedKmers, cluster);	// get seq2kg map for building KSM logo
	    	int leftmost = Integer.MAX_VALUE;
	    	int total_aligned_seqs = 0;
	    	for (Sequence s : seqList){
				if (s.pos==UNALIGNED)
					continue;
				if (s.pos < leftmost )
					leftmost = s.pos;		
				total_aligned_seqs++;
			}
	    	cluster.total_aligned_seqs = total_aligned_seqs;
	    	double[] bs = new double[total_aligned_seqs];
	    	int count = 0;
	    	int midPos=seqList.get(0).seq.length()/2;		// assume all the seqs are of the same length
			StringBuilder sb = new StringBuilder();
			for (Sequence s : seqList){
				if (s.pos==UNALIGNED)
					continue;
				String seq = s.getAlignedSeq();
				if (config.print_aligned_seqs)
					sb.append(String.format("%d\t%d\t%s\t%s%s\n", s.id, s.pos, s.isOriginalOrientation?"F":"R", CommonUtils.padding(-leftmost+s.pos, '.'), seq));
				bs[count]=midPos+s.pos;
				count++;
			}
			
			// make KSM logo
			// sort the seqs by the strongest Kmer, then KGscore
			ArrayList<Seq> allSeqList = new ArrayList<Seq>();
			ArrayList<Seq> seqSortList = new ArrayList<Seq>();
			HashMap<Kmer, ArrayList<Sequence>> km2seqs = new HashMap<Kmer, ArrayList<Sequence>>();
			for (Sequence s: seqList){
				Seq seq = new Seq();
				seq.id = s.id;
				allSeqList.add(seq);			
				if (s.pos!=UNALIGNED){
					for (Kmer km: seq2kg.get(s.id).kmers){
						if (!km2seqs.containsKey(km)){
							km2seqs.put(km, new ArrayList<Sequence>());
						}
						km2seqs.get(km).add(s);
					}					
					seq.kgScore = seq2kg.get(s.id).getScore();
					seqSortList.add(seq);
				}
			}
			seqSortList.trimToSize();
			allSeqList.trimToSize();
			
//			if (i==10)
//				System.out.println("m"+i);

			int totalHitCount = 0;
			int sortId = 0;
			ArrayList<Kmer> bestKmers = new ArrayList<Kmer>();
			ArrayList<Integer> kmHitCounts = new ArrayList<Integer>();
			while(totalHitCount<seqSortList.size()){
				int bestHit=0;
				Kmer bestKm = null;
				for (Kmer km:km2seqs.keySet()){
					if (bestHit<km2seqs.get(km).size()){
						bestHit=km2seqs.get(km).size();
						bestKm=km;
					}
				}
				totalHitCount += bestHit;
				if (totalHitCount>= seqSortList.size())
					break;
				bestKmers.add(bestKm);
				kmHitCounts.add(bestHit);
				ArrayList<Sequence> seqs = km2seqs.get(bestKm);
				km2seqs.remove(bestKm);
				for (Sequence s:seqs)
					if (allSeqList.get(s.id).kmerSortId > sortId)
						allSeqList.get(s.id).kmerSortId = sortId;
				for (Kmer km:km2seqs.keySet())
					km2seqs.get(km).removeAll(seqs);
//				System.out.println((bestKm.isSeedOrientation?bestKm.kmerString:bestKm.kmerRC) +"\t"+bestHit+"\t"+sortIdSetCount);
				sortId++;
			}
			Collections.sort(seqSortList);

			int pictHeight = 1200;		// KSM logo image height
			int minCount = seqSortList.size()/pictHeight*15;
			int lastCountIdx = kmHitCounts.size();
			int accumCount = 0;
			for (int ik = 0; ik<kmHitCounts.size(); ik++){
				accumCount += kmHitCounts.get(ik);
				if (kmHitCounts.get(ik)<minCount){
					lastCountIdx = ik;
					break;
				}
			}
			
			// print
			String[] ss = new String[seqSortList.size()];
			int seqAlignmentLength = cluster.k*3;
			for (int id=0; id<seqSortList.size(); id++){
				if (id==accumCount){
					ss[id]=CommonUtils.padding(seqAlignmentLength, 'N');
					continue;
				}
				Sequence s = seqList.get(seqSortList.get(id).id);
				String seq = s.getAlignedSeq();
				if (config.print_aligned_seqs)
					sb.append(String.format("%d\t%d\t%s\t%s%s\n", s.id, s.pos, s.isOriginalOrientation?"F":"R", CommonUtils.padding(-leftmost+s.pos, '.'), seq));
				int start = -s.pos - cluster.k;
				int end = start + seqAlignmentLength;
				int padding = 0;
				int endpadding = 0;
				if (start<0){
					padding = -start;
					start = 0;
				}
				if (end>=seq.length()){
					endpadding = end - seq.length();
					end = seq.length();
				}
				StringBuilder sb_logo = new StringBuilder();
				sb_logo.append(CommonUtils.padding(padding, 'N')).append(seq.substring(start, end))
					.append(CommonUtils.padding(endpadding, 'N'));
				ss[id] = sb_logo.toString();
			}
			CommonUtils.visualizeSequences(ss, 10, 1, new File(outName+".m"+i+".sequenceHits.png"));
			
			// KSM logo
			int currentSeqId = 0;
			ArrayList<WeightMatrix> wms = new ArrayList<WeightMatrix>();
			for (int ik = 0; ik<lastCountIdx; ik++){
//				if (currentSeqId+kmHitCounts.get(ik)>=ss.length)
//					break;
				float[][] pfm = new float[seqAlignmentLength][MAXLETTERVAL];  // positions x letters
				for (int p=0;p<pfm.length;p++){
					for (char base:LETTERS)			// 0 count can cause log(0), set pseudo-count 0.375 to every pos, every base
						pfm[p][base] = 0.375f; 		//http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2490743/
				} 
		    	for (int is=currentSeqId;is<currentSeqId+kmHitCounts.get(ik);is++){
					for (int p=0;p<pfm.length;p++){
		    			char base = ss[is].charAt(p);
		    			pfm[p][base] += 1;
		    		}
		    	}
				currentSeqId += kmHitCounts.get(ik);

				WeightMatrix wm = new WeightMatrix(pfm);
				if (config.ksm_logo_text)
					wm.setNameVerType("m"+i, "kmer"+ik, kmHitCounts.get(ik)+"/"+ss.length);
				else
					wm.setNameVerType(null,null,null);
				wms.add(wm);
			}
			CommonUtils.printKSMMotifLogo(wms, kmHitCounts.subList(0, lastCountIdx), new File(outName+".m"+i+".KSM.png"), pictHeight, 800/seqAlignmentLength);

			// median BS position relative to seed k-mer start 
			if (bs.length==0){
		    	if (config.verbose>1){
		    		System.out.println("!!! No binding site match !!!");
        			System.out.println(String.format("%s: #%d KSM %.2f\thit %d+/%d- seqs\tkAUC=%.1f\t%s", CommonUtils.timeElapsed(tic), i,
        					cluster.ksmThreshold.motif_cutoff, cluster.ksmThreshold.posHit, cluster.ksmThreshold.negHit, cluster.ksmThreshold.motif_significance, cluster.seedKmer.kmerString));
        			System.out.println(String.format("%s: #%d PWM %.2f/%.2f\thit %d+/%d- seqs\tpAUC=%.1f\t%s", CommonUtils.timeElapsed(tic), i,
        					cluster.pwmThreshold.motif_cutoff, cluster.wm.getMaxScore(), cluster.pwmThreshold.posHit, cluster.pwmThreshold.negHit, cluster.pwmThreshold.motif_significance, WeightMatrix.getMaxLetters(cluster.wm)));
		    	}
		    	continue;
			}
				
			cluster.pos_BS_seed=(int)Math.ceil(StatUtil.median(bs));		
			if (config.print_aligned_seqs)		// Note: the motif id of this seqs_aligned.txt may not be the same as final motif id.
				CommonUtils.writeFile(outName+"_"+i+"_seqs_aligned.txt", sb.toString());

			for (Kmer km: cluster.alignedKmers){			// set k-mer offset
				km.setKmerStartOffset(km.shift-cluster.pos_BS_seed);
			}			
		}	// for each motif found

		// print KSM and motif hits
		for (MotifCluster cluster : clusters){
			if (config.kg_hit_adjust_type==2){
				// wrap int into string[], a hack
				String[][] pos = new String[cluster.posCoveredWidth.length][];
				for (int i=0;i<pos.length;i++)
					pos[i] = new String[] {String.valueOf(cluster.posCoveredWidth[i])};
				String[][] neg = new String[cluster.negCoveredWidth.length][];
				for (int i=0;i<neg.length;i++)
					neg[i] = new String[] {String.valueOf(cluster.negCoveredWidth[i])};
				GappedKmer.printKSM(cluster.alignedKmers, pos, neg, seq_weights, cluster.k, 0, posSeqCount, negSeqCount, 
					cluster.ksmThreshold.motif_cutoff, outName+".m"+cluster.clusterId, false, true, false);
			}
			else // config.kg_hit_adjust_type ==1 or ==0 (null)
				GappedKmer.printKSM(cluster.alignedKmers, cluster.posHitStrings, cluster.negHitStrings, seq_weights, cluster.k, 0, posSeqCount, negSeqCount, 
						cluster.ksmThreshold.motif_cutoff, outName+".m"+cluster.clusterId, false, true, false);
		}
		if (config.print_motif_hits){		// PWM motif hits
			ArrayList<WeightMatrix> pwms=new ArrayList<WeightMatrix>();
			ArrayList<Double> thresholds=new ArrayList<Double>();
			for (MotifCluster c : clusters){
				if (c.wm!=null){
					pwms.add(c.wm);
					thresholds.add(c.pwmThreshold.motif_cutoff);
				}
			}
			ArrayList<MotifInstance> instances = MotifScan.getPWMInstances(seqs, pwms, thresholds);
			StringBuilder sb = new StringBuilder("#Motif\tSeqID\tMatch\tSeqPos\tScore\n");
		    for (int i=0;i<instances.size();i++){
		    	MotifInstance mi = instances.get(i);
		    	sb.append(mi.motifID).append("\t").append(mi.seqID).append("\t").append(mi.matchSeq).append("\t")
		    	.append(mi.strand=='+'?mi.position:-mi.position).append("\t").append(String.format("%.2f", mi.score)).append("\n");
		    }
		    CommonUtils.writeFile(outName+".PWM.motifInstances.txt", sb.toString());
		}
		
		/** final outputs */
		if (clusters.isEmpty()){
			System.out.println("\n----------------------------------------------\nNone of the k values form an enriched PWM, stop here!\n");
			File f = new File(outName);
			String name = f.getName();
			StringBuffer html = new StringBuffer("<style type='text/css'>/* <![CDATA[ */ table, td{border-color: #600;border-style: solid;} table{border-width: 0 0 1px 1px; border-spacing: 0;border-collapse: collapse;} td{margin: 0;padding: 4px;border-width: 1px 1px 0 0;} /* ]]> */</style>");
			html.append("<script language='javascript' type='text/javascript'><!--\nfunction popitup(url) {	newwindow=window.open(url,'name','height=75,width=400');	if (window.focus) {newwindow.focus()}	return false;}// --></script>");
			html.append("<table><th bgcolor='#A8CFFF'><font size='5'>");
			html.append(name).append("</font></th>");
			html.append("<tr><td valign='top' width='500'><br>");
			if (!this.standalone && eventCounts!=null){
				html.append("<a href='"+name+".GEM_events.txt'>Significant Events</a>&nbsp;&nbsp;: "+eventCounts[0]);
				html.append("<br><a href='"+name+".GEM_insignificant.txt'>Insignificant Events</a>: "+eventCounts[1]);
				html.append("<br><a href='"+name+".GEM_filtered.txt'>Filtered Events</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;: "+eventCounts[2]);
			}
			html.append("<p><p>Motif can not be found!<p>");
			html.append("</td></tr></table>");
			CommonUtils.writeFile(outName+".results.htm", html.toString());
			return;
		}
		
		// print PWM spatial distribtution
		computeMotifDistanceDistribution(outName);
		outputClusters(eventCounts);
	}
	
	
	/** 
	 * Generate ungapped k-mers from the positive sequences, upto a relax enrichment cutoff
	 * */
	public HashMap<String, Kmer> generateKmers(int k){
		tic = System.currentTimeMillis();
		
		double relaxFactor = 0.6;
		double relaxed_fold = Math.max(config.k_fold*relaxFactor, 1);
		double relaxed_hgp = Math.min(relaxFactor*config.kmer_hgp, -1.3);
		
		// compute the smallest PosCount needed to be significant, with negCount=0
		int minSigBaseKmerCount;
		for (minSigBaseKmerCount=3;minSigBaseKmerCount<posSeqCount;minSigBaseKmerCount++){
			double hgp = computeHGP(minSigBaseKmerCount, 0);
			if (hgp<relaxed_hgp){
				break;
			}
		}
		// expected count of kmer = total possible unique occurences of kmer in sequence / total possible kmer sequence permutation
		// num_seq * num_kmer_per_seq / total_num_kmer_for_k, 
		// sequence: scan for reverse compliment, k-mer: merge reverse compliment, therefore cancel out
		int expectedBaseKmerCount = Math.max(minSigBaseKmerCount, (int) Math.round(seqs.length*(seqs[0].length()-k+1) / Math.pow(4, k)));

		
		/*******************************************************************
		 * Scan the sequences to generate kmers
		 ******************************************************************/
		HashMap<String, HashSet<Integer>> kmerstr2seqs = new HashMap<String, HashSet<Integer>>();
		for (int seqId=0;seqId<posSeqCount;seqId++){
			// only count one direction, k-mer and its RC will be merged with all hits being added together
			String seq = seqs[seqId];
			int numPos = seq.length()-k+1;		// substring() is end exclusive
			HashSet<String> uniqueKmers = new HashSet<String>();			// only count repeated kmer once in a sequence
			for (int i=0;i<numPos;i++){
				String kstring = seq.substring(i, i+k);
				if (kstring.contains("N"))									// ignore 'N', converted from repeat when loading the sequences
					continue;
				uniqueKmers.add(kstring);
			}
			for (String s: uniqueKmers){
				if (!kmerstr2seqs.containsKey(s)){
					 kmerstr2seqs.put(s, new HashSet<Integer>());
				}
				kmerstr2seqs.get(s).add(seqId);
			}
		}

		// Merge kmer and its reverse compliment(RC)
		// Remove low-count base k-mers
		ArrayList<String> kmerStrings = new ArrayList<String>();
		kmerStrings.addAll(kmerstr2seqs.keySet());
		
		// create kmers from its and RC's counts
		for (String key:kmerStrings){
			if (!kmerstr2seqs.containsKey(key))		// this kmer has been removed, represented by RC
				continue;
			// consolidate kmer and its reverseComplment kmer, remove if count is not high enough
			String key_rc = SequenceUtils.reverseComplement(key);				
			if (!key_rc.equals(key)){	// if it is not reverse compliment itself
				if (kmerstr2seqs.containsKey(key_rc)){
					int kCount = kmerstr2seqs.get(key).size();
					int rcCount = kmerstr2seqs.get(key_rc).size();
					// if k-mer hit is less than expected count, remove
					if (kCount+rcCount < expectedBaseKmerCount){
						kmerstr2seqs.remove(key);	
						kmerstr2seqs.remove(key_rc);	
					}
					else{
						String winner = kCount>=rcCount?key:key_rc;
						String loser = kCount>=rcCount?key_rc:key;
						kmerstr2seqs.get(winner).addAll(kmerstr2seqs.get(loser));	// winner takes all
						kmerstr2seqs.remove(loser);					// remove the loser kmer because it is represented by its RC
					}
				}
				else if (kmerstr2seqs.get(key).size()<expectedBaseKmerCount){
					kmerstr2seqs.remove(key);						
				}				
			}
		}
		if (config.verbose > 1)
			System.out.println(String.format("k=%d, mapped %d base k-mers, min_base_kmer_Hit=%d, %s", k,
				kmerstr2seqs.keySet().size(), expectedBaseKmerCount, CommonUtils.timeElapsed(tic)));

		
		/*****************************************************************
		 * Select significantly over-representative kmers 
		 * Count kmer hits in the negative sequences, then compute hgp
		 *****************************************************************/
		
		tic = System.currentTimeMillis();
		//Aho-Corasick for searching Kmers in negative sequences
		//ahocorasick_java-1.1.tar.gz is an implementation of Aho-Corasick automata for Java. BSD license.
		//from <http://hkn.eecs.berkeley.edu/~dyoo/java/index.html> 
		AhoCorasick tmp = new AhoCorasick();
		/** Base K-mer String to k-mer map */
		HashMap<String, Kmer> bkMap = new HashMap<String, Kmer>();
		for (String s:kmerstr2seqs.keySet()){	
			Kmer kmer = new Kmer(s, kmerstr2seqs.get(s), seq_weights);		// create k-mers
			bkMap.put(kmer.kmerString,kmer);
			tmp.add(s.getBytes(), s);
		}
		kmerstr2seqs=null;	// clean up
		System.gc();

		tmp.prepare();
		
		// count hits in the negative sequences
		HashMap<String, HashSet<Integer>> kmerstr2negSeqs = new HashMap<String, HashSet<Integer>>();
		for (int negSeqId=0; negSeqId<negSeqCount;negSeqId++){
			String seq = seqsNegList.get(negSeqId);
			HashSet<Object> kmerHits = new HashSet<Object>();	// to ensure each sequence is only counted once for each kmer
			Iterator searcher = tmp.search(seq.getBytes());
			while (searcher.hasNext()) {
				SearchResult result = (SearchResult) searcher.next();
				kmerHits.addAll(result.getOutputs());
			}
			String seq_rc = SequenceUtils.reverseComplement(seq);
			searcher = tmp.search(seq_rc.getBytes());
			while (searcher.hasNext()) {
				SearchResult result = (SearchResult) searcher.next();
				kmerHits.addAll(result.getOutputs());
			}
			for (Object o: kmerHits){
				String ks = (String) o;
				if (!kmerstr2negSeqs.containsKey(ks))					
					kmerstr2negSeqs.put(ks, new HashSet<Integer>());
				kmerstr2negSeqs.get(ks).add(negSeqId);
			}
		}
		
		// remove k-mers that do not pass the relaxed hgp cutoff
		HashSet<String> toRemove = new HashSet<String>();
		for (Kmer kmer:bkMap.values()){
			HashSet<Integer> neghits = kmerstr2negSeqs.get(kmer.kmerString);
			if (neghits==null)
				neghits = new HashSet<Integer>();
			if (kmer.getPosHitCount() > neghits.size() / get_NP_ratio() * relaxed_fold){
				double hgp = computeHGP(kmer.getPosHitCount(), neghits.size());
				if (hgp < relaxed_hgp){		// the hgp and fold change here is slightly relaxed for base-kmers
					kmer.setNegHits(neghits);	
					kmer.setHgp(hgp);
				}
				else
					toRemove.add(kmer.kmerString);	
			}
			else{
				toRemove.add(kmer.kmerString);		
			}
		}
		for (String ks:toRemove)
			bkMap.remove(ks);
		
		allKmerMap.put(k, bkMap);
		if (config.verbose > 1)
			System.out.println(String.format("k=%d, relaxed_hgp=%.2f, total exact kmer=%d, %s", 
				k, relaxed_hgp, bkMap.size(), CommonUtils.timeElapsed(tic)));
		return bkMap;
	}
	
	private ArrayList<char[]> getPatterns(int k, int gap){
		/***************************************************************************************
		 * Compute pair-wise k-mer mismatch distance, prepare for making gapped k-mers
		 ***************************************************************************************/
		int kFull = k+gap;
		
		ArrayList<Kmer> kmers = new ArrayList<Kmer>();
		kmers.addAll(allKmerMap.get(kFull).values());
		kmers.trimToSize();
		Collections.sort(kmers);

		ArrayList<char[]> patterns = new ArrayList<char[]>();	// gap patterns	
		ArrayList<Integer> mismatchPos = new ArrayList<Integer>();
		int mismatchCount = 0;
		String s1=null;
		String s2=null;
		for (int i=0;i<kmers.size();i++){
			for (int j=i+1;j<kmers.size();j++){
				s1 = kmers.get(i).kmerString;
				s2 = kmers.get(j).kmerString;
				if (s1.charAt(0)==s2.charAt(0) && s1.charAt(kFull-1)==s2.charAt(kFull-1)){
					for (int id=1;id<kFull-1;id++){	// exclude edge gaps	// s1 and s2 are of same length
						if (s1.charAt(id)!=s2.charAt(id)){
							mismatchCount++;
							mismatchPos.add(id);
						}
						if (mismatchCount>gap)			// early stop: only care mismatch numGap
							break;
					}
					if (mismatchCount==gap){	
						char[] pattern = s1.toCharArray();
						for (int id:mismatchPos)
							pattern[id]='N';
						patterns.add(pattern);
//						System.out.println(s1+" "+s2+" "+String.valueOf(pattern));
					}
					mismatchPos.clear();
					mismatchCount=0;
				}
				// RC  (the following code is duplicated, don't want to make a function to avoid overhead, b/c this will loop NxN.)
				s2 = kmers.get(j).kmerRC;
				if (s1.charAt(0)==s2.charAt(0) && s1.charAt(kFull-1)==s2.charAt(kFull-1)){	// edge must be the same, otherwise the pattern is the same as ungapped k-mer
					for (int id=1;id<s1.length()-1;id++){	// exclude edge gaps	// s1 and s2 are of same length
						if (s1.charAt(id)!=s2.charAt(id)){
							mismatchCount++;
							mismatchPos.add(id);
						}
						if (mismatchCount>gap)			// early stop: only care mismatch numGap
							break;
					}
					if (mismatchCount==gap){	
						char[] pattern = s1.toCharArray();
						for (int id:mismatchPos)
							pattern[id]='N';
						patterns.add(pattern);
//						System.out.println(s1+" "+s2+" "+String.valueOf(pattern));
					}
					mismatchPos.clear();
					mismatchCount=0;
				}
			}
		}
		patterns.trimToSize();
//		System.out.print(String.format("%d+%d(%d) ", k, gap, patterns.size())); 
		
		return patterns;
	}

	/**
	 * Assemble exact and gapped k-mers for k using pre-defined gapped wildcard patterns
	 */
	private Pair<ArrayList<Kmer>, ArrayList<Kmer>> selectEnrichedKmers(int k){ 
		ArrayList<Kmer> allSignificantKmers = new ArrayList<Kmer>();	// object to save all the significant k-mers
		ArrayList<Kmer> kmers = new ArrayList<Kmer>();
		// exact k-mers
		HashMap<String, Kmer> ungappedKmerMap = null;
		
		if (allKmerMap.containsKey(k))
			ungappedKmerMap = allKmerMap.get(k);
		else
			ungappedKmerMap = generateKmers(k);
		for (Kmer kmer: ungappedKmerMap.values()){
			if (kmer.getHgp() <= config.kmer_hgp){		// get those pass HGP cutoff
				kmers.add(kmer);	
			}
		}
		allSignificantKmers.addAll(kmers);
		if (config.print_all_kmers){
			ArrayList<Kmer> kmersAll = new ArrayList<Kmer>();
			kmersAll.addAll(ungappedKmerMap.values());
			Collections.sort(kmersAll);
			GappedKmer.printKSM(kmersAll, null, null, null, k, 0, posSeqCount, negSeqCount, 0, outName+"_all_w"+seqs[0].length(), true, false, true);
		}
		allKmerMap.remove(k);
		ungappedKmerMap = null;		// remove the k-mers with length k, assuming selectEnrichedKmers() is called with increasing k, they will not be used any more
		Collections.sort(kmers, new Comparator<Kmer>(){
            public int compare(Kmer o1, Kmer o2) {
                return o1.compareByHGP(o2);
            }
        });
		System.out.println(String.format("k=%d, exact kmers=%d", k, kmers.size()));
		
		/**
		 * Construct gapped k-mers from the significant (relaxed fold change and hgp) base k-mers
		 */
		for (int numGap=1;numGap<=config.gap;numGap++){
			HashMap<String, GappedKmer> gkMap = new HashMap<String, GappedKmer>();
			int kFull = k+numGap;
			
			HashMap<String,Kmer> bkMap = null;
			if (allKmerMap.containsKey(kFull))
				bkMap =	allKmerMap.get(kFull);
			else{
				bkMap = generateKmers(kFull);
				if (k<config.k_max)
					allKmerMap.put(kFull, bkMap);
			}
			ArrayList<char[]> allPatterns = getPatterns(k, numGap);
			
			// prepare the variants
			int numVariants = 1;
			for (int j=0;j<numGap;j++){
				numVariants *= LETTERS.length;
			}
			char[][] variants = new char[numVariants][numGap];
			for (int j=0;j<numVariants;j++){
				char[] v = new char[numGap];
				int div = j;
				for (int l=0;l<v.length;l++){
					v[l] = LETTERS[div % LETTERS.length];
					div /= LETTERS.length;
				}
				variants[j]=v;
			}
			
			for (char[] p : allPatterns){
				// check whether this gapped kmer has been made
				String gkStr = String.valueOf(p);
				if (gkMap.containsKey(gkStr)||gkMap.containsKey(SequenceUtils.reverseComplement(gkStr)))
					continue;	// this gapped kmer has been made, skip to next pattern
				
				ArrayList<Integer> gapPos = new ArrayList<Integer>();
				for (int i=0;i<p.length;i++)
					if (p[i]=='N')
						gapPos.add(i);
				
				GappedKmer gk = new GappedKmer(gkStr);

				for (int j=0;j<numVariants;j++){
					char[] v = variants[j];
					for (int l=0;l<v.length;l++){
						p[gapPos.get(l)] = v[l];
					}
					String m = String.valueOf(p);
					if (bkMap.containsKey(m))
						gk.addBaseKmer(bkMap.get(m), true);
					else{
						String mrc = SequenceUtils.reverseComplement(m);
						if (bkMap.containsKey(mrc)){
							gk.addBaseKmer(bkMap.get(mrc), false); 
						}
					}
				}
				
				//TODO: think!!! ignore if a singleton longer kmer, b/c it will be covered by a basic k-mer
				if (gk.getBaseKmers().size()>1){	
					gk.mergePosHits(seq_weights);
					gkMap.put(gk.kmerString, gk);
				}
				else{
					System.err.println(gk.toString());
				}
			}// loop all patterns
			allPatterns = null;
			
			ArrayList<GappedKmer> gks = new ArrayList<GappedKmer>();
			HashSet<Kmer> baseKmers = new HashSet<Kmer>();
			for (String key:gkMap.keySet()){	
				GappedKmer gk = gkMap.get(key);
				gks.add(gk);
				baseKmers.addAll(gk.getBaseKmers());
			}
			System.out.println(String.format("k=%d+%d, gapped kmer=%d, base kmer=%d", 
					k, numGap, gks.size(), baseKmers.size()));
			
			
			/** 
			 *  score the gapped kmers, hypergeometric p-value, select significant k-mers
			 */
			ArrayList<Kmer> results = new ArrayList<Kmer>();
			for (GappedKmer gk:gks){
				gk.mergeNegHits();		// aggregate base-kmer negative hits to the gapped kmer
				if (gk.getPosHitCount() >= gk.getNegHitCount()/get_NP_ratio() * config.k_fold ){
					// optimize the base k-mers to get best hgp
					ArrayList<Kmer> baseKmerList = new ArrayList<Kmer>();
					baseKmerList.addAll(gk.getBaseKmers());
					if (config.optimize_base_kmers)
						optimizeKSM(baseKmerList, pseudoCountRatios[k]);
					HashSet<Kmer> toremove = new HashSet<Kmer>();
					for (Kmer baseKm:gk.getBaseKmers())
						if (!baseKmerList.contains(baseKm))
							toremove.add(baseKm);
					for (Kmer km: toremove)
						gk.removeBasekmers(km);

					if (gk.getBaseKmers().size()<=1)		// skip if only has 1 subkmer
						continue;
					gk.mergePosHits(seq_weights);
					gk.mergeNegHits();
					gk.setHgp(computeHGP(gk.getPosHitCount(), gk.getNegHitCount()));
					if (gk.getHgp() <= config.kmer_hgp){		// Gapped k-mers passing cutoff are significant
						results.add(gk);
					}
				}
			}
			results.trimToSize();
			Collections.sort(results, new Comparator<Kmer>(){
	            public int compare(Kmer o1, Kmer o2) {
	                return o1.compareByHGP(o2);
	            }
	        });	
			// prepare matrix for distance calculation, will be clear in discoverMotif() method
			for (Kmer km:results){
				km.setMatrix();
			}
			allSignificantKmers.addAll(results);
		
			/** 
			 *  reduce the k-mer set by removing similar gapped k-mers
			 */
			if (config.cluster_gapped){
				ArrayList<Kmer> results_final = new ArrayList<Kmer>();
				double dist_cutoff_to_reduce = 0;
				int numTop = config.max_gkmer*5;
				if (results.size()>numTop){
					if (config.verbose>1)
						System.out.println(String.format("Reduce gapped-kmer, n=%d, hgp1x=%.1f, hgp5x=%.1f", 
								results.size(), results.get(config.max_gkmer).hgp_lg10, results.get(numTop).hgp_lg10));
	
					for (int i=0;i<numTop;i++)
						results_final.add(results.get(i));
					results = results_final;
				}
				if (results.size()>config.max_gkmer){
					
	//				// direct reduce
	//				if (results.size()>1000){
	//					long tic=System.currentTimeMillis();
	//					dist_cutoff_to_reduce = 1.5;
	//					for (Kmer km: results)
	//						km.setMatrix();
	//					int[] idx = reduceGappedKmerSet(results, dist_cutoff_to_reduce);
	//					for (int i=0;i<results.size();i++){
	//						if (idx[i]==1)
	//							results_final.add(results.get(i));
	//					}
	//					if (config.verbose>1)
	//					System.out.println("Reduce gapped-kmer set, dist="+dist_cutoff_to_reduce+", from " + 
	//							idx.length + " to " + results_final.size() + ", " + CommonUtils.timeElapsed(tic));
	//				}
					
					//  otherwise too heavy load for mtree range search
	
					// reduce by distance, limit cutoff to 1.5
//					while (results.size()>config.max_gkmer && !(dist_cutoff_to_reduce>=1.5)){	
//						long tic=System.currentTimeMillis();
//						dist_cutoff_to_reduce += 0.5;
//						// the results list has been sorted by hgp
//						MTree dataPoints = MTree.constructTree(results, 10);
//						int[] idx = reduceGappedKmerSetMtree(dataPoints, dist_cutoff_to_reduce);
//						results_final = new ArrayList<Kmer>();
//						for (int i=0;i<results.size();i++){
//							if (idx[i]==1)
//								results_final.add(results.get(i));
//						}
//						if (config.verbose>1)
//							System.out.println("Reduce gapped-kmer, dist="+dist_cutoff_to_reduce+", from " + 
//									idx.length + " to " + results_final.size() + ", " + CommonUtils.timeElapsed(tic));
//						if (results_final.size()<config.max_gkmer/2){		// if reduce too much, not reduce, take the top k-mers
//							results_final = new ArrayList<Kmer>();
//							int final_count = Math.min(results.size(), config.max_gkmer);
//							for (int i=0;i<final_count;i++){
//								results_final.add(results.get(i));
//							}
//							break;
//						}
//						else	// prepare for next round of reduction
//							results = results_final;
//					}
					
					// reduce by similar hit counts
					double sharedHitsRatio = 0.9;
					while (results.size()>config.max_gkmer && !(sharedHitsRatio<=0.5)){	
						long tic=System.currentTimeMillis();
						sharedHitsRatio -= 0.1;
						// the results list has been sorted by hgp
						int[] idx = reduceGappedKmerSetByHits(results, sharedHitsRatio);  // cutoff = shared / min(iHits, jHits)
						results_final = new ArrayList<Kmer>();
						for (int i=0;i<results.size();i++){
							if (idx[i]==1)
								results_final.add(results.get(i));
						}
						if (config.verbose>1)
							System.out.println(String.format("Reduce gapped-kmers with shared ratio >= %.1f, from %d to %d, %s", 
									sharedHitsRatio, idx.length, results_final.size(), CommonUtils.timeElapsed(tic)));
						if (results_final.size()<config.max_gkmer/2){		// if reduce too much, not reduce, take the top k-mers
							results_final = new ArrayList<Kmer>();
							int final_count = Math.min(results.size(), config.max_gkmer);
							for (int i=0;i<final_count;i++){
								results_final.add(results.get(i));
							}
							break;
						}
						else	// prepare for next round of reduction
							results = results_final;
					}					
					if (results_final.size()>config.max_gkmer){		// if still too many, take top k-mers
						results = results_final;
						results_final = new ArrayList<Kmer>();
						int final_count = Math.min(results.size(), config.max_gkmer);
						for (int i=0;i<final_count;i++){
							results_final.add(results.get(i));
						}
					}
				}
				else		// no need to reduce
					results_final = results;
				
	//			System.out.println(String.format("k=%d+%d, selected %d gapped k-mers (hgp=%.1f) from %d+/%d- sequences, %s", 
	//					kOrginal, numGap, results_final.size(), results_final.size()==0?0:results_final.get(results_final.size()-1).hgp_lg10, posSeqCount, negSeqCount, CommonUtils.timeElapsed(tic)));
				System.out.println(String.format("k=%d+%d, selected %d gapped k-mers (hgp=%.1f), %s", 
						k, numGap, results_final.size(), results_final.size()==0?0:results_final.get(results_final.size()-1).hgp_lg10, CommonUtils.timeElapsed(tic)));
				
				kmers.addAll(results_final);
			}
			// print all gapped k-mers 
			if (config.print_all_kmers){
				ArrayList<Kmer> kms = new ArrayList<Kmer>();
				for (GappedKmer gk:gks){
					kms.add(gk);			
				}
				Collections.sort(kms);
				GappedKmer.printKSM(kms, null, null, null, k, numGap, posSeqCount, negSeqCount, 0, outName+"_all_w"+seqs[0].length(), true, false, true);
			}
			
		}// loop all gap numbers
		
		return new Pair<ArrayList<Kmer>, ArrayList<Kmer>>(allSignificantKmers, kmers);
	}

	// assuming the 
	private static int[] reduceGappedKmerSetMtree(MTree mtreeDataPoints, double distanceCutoff){
		// keys = kmer indices (0 to dataPoints.getSize() - 1)
		// values = kmers in range
		int[] idx = new int[mtreeDataPoints.getSize()];
		HashMap<Integer, ArrayList<Kmer>> rangeResults = new HashMap<Integer, ArrayList<Kmer>>();
		ArrayList<TreeObject> traversal = MTree.traverse(mtreeDataPoints.getRoot());

		// for every k-mer, find k-mers within distance, register the pair bi-directionally
		for (TreeObject o : traversal) {
			ArrayList<Kmer> rangeResult = mtreeDataPoints.rangeSearch(o.getData(), distanceCutoff);
			for (Kmer kmer : rangeResult) {
				if (rangeResults.get(o.getData().getIndex()) == null) 
					rangeResults.put(o.getData().getIndex(), new ArrayList<Kmer>());
				rangeResults.get(o.getData().getIndex()).add(kmer);

				if (rangeResults.get(kmer.getIndex()) == null) 
					rangeResults.put(kmer.getIndex(), new ArrayList<Kmer>());
				rangeResults.get(kmer.getIndex()).add(o.getData());
			}
			// remove o object from the tree
			MTreeNode container = o.getContainer();
			container.getObjects().remove(o.getIndex());
			if (container.getObjects().size() == 0) {
				if (container.getParent() != null) {
					container.getParent().setChild(null);
				}
			}
		}

		// remove all the k-mers that were marked as in-range with other k-mers
		// idx array is indexed by k-mer id (i.e. treeIndex)
		// the m-tree objects and index array has been sorted before this method being called
		// so the following will keep the high ranking k-mers and remove the weaker k-mers
		for (int i = 0; i < idx.length; i++) {
			if (idx[i] == -1)
				continue;
			ArrayList<Kmer> kmers = rangeResults.get(i);
			for (Kmer km: kmers)
				if (idx[km.getIndex()]!=-1)
					idx[km.getIndex()]=-1;
			idx[i] = 1;
		}
		
		return idx;
	}

	private static int[] reduceGappedKmerSetByHits(ArrayList<Kmer> kmers, double orCutoff){
		// keys = kmer indices (0 to dataPoints.getSize() - 1)
		// values = kmers in range
		int[] idx = new int[kmers.size()];		// 1: selected;  -1: to be removed
		// idx array is indexed by k-mer id (i.e. treeIndex)
		// the kmers list has been sorted before this method being called
		// so the following will keep the high ranking k-mers and remove the weaker k-mers
		for (int i = 0; i < idx.length-1; i++) {
			if (idx[i] == -1)
				continue;
			Kmer ki = kmers.get(i);
			int iHits = ki.getPosBits().cardinality();
			idx[i] = 1;
			int count=0;
//			System.out.println(ki.toShortString());
			for (int j = i+1; j < idx.length; j++){	// this is O( n^2 ) looping, but we have removed some indexes
				if (idx[j] == -1)
					continue;
				Kmer kj = kmers.get(j);
				BitSet cp = (BitSet) ki.getPosBits().clone();	// clone, b/c BitSet.add() will modify BitSet obj
				int jHits = kj.getPosBits().cardinality();
				cp.and(kj.getPosBits());		
				int sharedHits = cp.cardinality();
				// kmer i hits as positive
//				double or = StatUtil.odds_ratio(iHits, posSeqCount-iHits, sharedHits, jHits-sharedHits, 10, 10);
				double or  =  ((double)sharedHits )/ Math.min(iHits, jHits);
//				System.out.print(String.format("%s\t%d\t%d\t%.2f\t%s\t", ki.getKmerStrRC(), jHits, sharedHits, or, kj.getKmerStrRC()));
				if (or>orCutoff){
					idx[j]=-1;
					count++;
				}
			} // j
//			System.out.print(" "+count);
		} // i
		return idx;
	}
	
	private static int[] reduceGappedKmerSet(ArrayList<Kmer> kmers, double distanceCutoff){
		// keys = kmer indices (0 to dataPoints.getSize() - 1)
		// values = kmers in range
		int[] idx = new int[kmers.size()];
		// idx array is indexed by k-mer id (i.e. treeIndex)
		// the kmers list has been sorted before this method being called
		// so the following will keep the high ranking k-mers and remove the weaker k-mers
		for (int i = 0; i < idx.length-1; i++) {
			if (idx[i] == -1)
				continue;
			Kmer kmi = kmers.get(i);
			for (int j = i+1; j < idx.length; j++)		// this is O( n^2 ) looping, but we have removed some indexes
				if (idx[j]!=-1 && KMAC.editDistance(kmi, kmers.get(j))<=distanceCutoff)
					idx[j]=-1;
			idx[i] = 1;
		}
		
		return idx;
	}
	
	private double[][] computeDistanceMatrix2(ArrayList<Kmer> kmers, boolean print_dist_matrix, int cutoff) {
		int kmerCount = kmers.size();
		double[][] distanceMatrix = new double[kmerCount][kmerCount];
		for (int i = 0; i < kmerCount; i++) {
			System.out.println(i);
			for (int j = 0; j <= i; j++) {
				double d = MTree.testDist(kmers.get(i).kmerString, kmers.get(j).kmerString);
				distanceMatrix[i][j] = d;
				distanceMatrix[j][i] = d;
			}
		}
		return distanceMatrix;
	}
	
	private double[][] computeDistanceMatrix(ArrayList<Kmer> kmers, boolean print_dist_matrix, int cutoff) {
		int kmerCount = kmers.size();
		double[][]distanceMatrix = new double[kmerCount][kmerCount];
		for (int i=0;i<kmerCount;i++){
			String ks_i = kmers.get(i).kmerString;
			for (int j=0;j<=i;j++){
//				String ks_j = kmers.get(j).getKmerString();
				distanceMatrix[i][j]=CommonUtils.strMinDistanceWithCutoff(ks_i, kmers.get(j).kmerString, cutoff);
//				distanceMatrix[i][j]=Math.min(CommonUtils.strMinDistance(ks_i, ks_j), CommonUtils.strMinDistance(krc_i, ks_j));  // much slower
			}
		}
		for (int i=0;i<kmerCount;i++){
			for (int j=i+1;j<kmerCount;j++){
				distanceMatrix[i][j]=distanceMatrix[j][i];
			}
		}
//		System.out.println("computeDistanceMatrix, "+CommonUtils.timeElapsed(tic));

		if (print_dist_matrix){
	        StringBuilder output = new StringBuilder();
	        for (int j=0;j<kmerCount;j++){
	        	output.append(String.format("%s\t",kmers.get(j).kmerString));
	        }
	        CommonUtils.replaceEnd(output, '\n');
	        for (int j=0;j<kmerCount;j++){
	        	output.append(String.format("%d\t",kmers.get(j).getPosHitCount()));
	        }
	        CommonUtils.replaceEnd(output, '\n');
	        for (int j=0;j<kmerCount;j++){
		        for (int i=0;i<distanceMatrix[j].length-1;i++){
		        	output.append(String.format("%.1f\t",distanceMatrix[j][i]));
		        }
		        output.append(String.format("%.1f",distanceMatrix[j][kmerCount-1])).append("\n");
	        }
	        CommonUtils.writeFile(outName+".distance_matrix.txt", output.toString());
	        
	        // Gephi format
	        // http://gephi.github.io/users/supported-graph-formats/csv-format/, Matrix
	        output = new StringBuilder();
	        // header line
	        output.append(";");
	        for (int j=0;j<kmerCount;j++){
	        	output.append(String.format("%s;",kmers.get(j).kmerString));
	        }
	        CommonUtils.replaceEnd(output, '\n');
	        // data lines
	        for (int j=0;j<kmerCount;j++){
	        	output.append(kmers.get(j).kmerString+";");
		        for (int i=0;i<distanceMatrix[j].length;i++){
		        	output.append(String.format("%.1f;",distanceMatrix[j][i]));
		        }
		        CommonUtils.replaceEnd(output, '\n');
	        }
	        CommonUtils.writeFile(outName+".distance_matrix.csv", output.toString());
		}
		return distanceMatrix;
	}
	
	/**
	 * Compute Weighted Distance Matrix<br>
	 * Old version, for gapped kmers, the distance is the weighted average distance between base kmers of the gapped kmers<br>
	 * The weight is the net hit count of the kmers
	 * @param kmers
	 * @param print_dist_matrix
	 * @return
	 */
	private double[][] computeWeightedDistanceMatrix(ArrayList<Kmer> kmers, boolean print_dist_matrix, int cutoff, int k) {
		long tic=System.currentTimeMillis();
		if (config.verbose>1)
			System.out.print("Computing k-mer distance matrix ... ");
		
		/** basicKmers include the exact kmers and the base-kmers of the gapped kmers */
		HashSet<Kmer> basicKmerSet = new HashSet<Kmer>();
		for (Kmer km: kmers){
			km.addBaseKmersToSet(basicKmerSet);
	    }
		ArrayList<Kmer> basicKmers = new ArrayList<Kmer>();
		basicKmers.addAll(basicKmerSet);
		HashMap<Kmer,Integer> basicKmerMap = new HashMap<Kmer,Integer>();
		for (int i=0;i<basicKmers.size();i++){
			basicKmerMap.put(basicKmers.get(i), i);
		}
		double[][] bDist = computeDistanceMatrix(basicKmers, false, cutoff);
		
		int kmerCount = kmers.size();
		double[][]distanceMatrix = new double[kmerCount][kmerCount];
		for (int i=0;i<kmerCount;i++){
			Kmer km = kmers.get(i);
			ArrayList<Kmer> sks = new ArrayList<Kmer>();
			if (km instanceof GappedKmer)
				sks.addAll( ((GappedKmer)km).getBaseKmers() );
			else
				sks.add(km);
					
			for (int j=0;j<=i;j++){
				Kmer kmj = kmers.get(j);
				ArrayList<Kmer> skj = new ArrayList<Kmer>();
				if (kmj instanceof GappedKmer)
					skj.addAll( ((GappedKmer)kmj).getBaseKmers() );
				else
					skj.add(kmj);
				
				double distSum = 0;
				double weightSum = 0;
				
				for (Kmer k1:sks)
					for (Kmer k2: skj){
						double w2 = k1.getNetHitCount(posNegSeqRatio)*k2.getNetHitCount(posNegSeqRatio);
						weightSum += w2;
						distSum += w2 * bDist[basicKmerMap.get(k1)][basicKmerMap.get(k2)];
					}
				distanceMatrix[i][j]=distSum/weightSum;
			}
		}
		for (int i=0;i<kmerCount;i++){
			for (int j=i+1;j<kmerCount;j++){
				distanceMatrix[i][j]=distanceMatrix[j][i];
			}
		}
		for (int i=0;i<kmerCount;i++)
			distanceMatrix[i][i]=0;
		
//		System.out.println("computeDistanceMatrix, "+CommonUtils.timeElapsed(tic));
	
		if (print_dist_matrix){
	        StringBuilder output = new StringBuilder();
	        for (int j=0;j<kmerCount;j++){
	        	output.append(String.format("%s\t",kmers.get(j).kmerString));
	        }
	        CommonUtils.replaceEnd(output, '\n');
	        for (int j=0;j<kmerCount;j++){
	        	output.append(String.format("%d\t",kmers.get(j).getNetHitCount(posNegSeqRatio)));
	        }
	        CommonUtils.replaceEnd(output, '\n');
	        for (int j=0;j<kmerCount;j++){
		        for (int i=0;i<distanceMatrix[j].length-1;i++){
		        	output.append(String.format("%.1f\t",distanceMatrix[j][i]));
		        }
		        output.append(String.format("%.1f",distanceMatrix[j][kmerCount-1])).append("\n");
	        }
	        CommonUtils.writeFile(outName+".weighted_distance_matrix.txt", output.toString());
	        
	        // Gephi format
	        // http://gephi.github.io/users/supported-graph-formats/csv-format/, Matrix
	        output = new StringBuilder();
	        output.append("Id,Label,Count\n");
	        for (int j=0;j<kmerCount;j++){
	        	output.append(String.format("%d,%s,%d\n",j,kmers.get(j).kmerString,kmers.get(j).getNetHitCount(posNegSeqRatio)));
	        }
	        CommonUtils.writeFile(String.format("%s.k%d.dist.gephi_nodes.csv", outName, k), output.toString());
	        
	        output = new StringBuilder();
	        double dist_cutoff = 1.5;
	        output.append("Target,Source,Weight,Type,Dist\n");
	        for (int j=0;j<kmerCount;j++){
		        for (int i=0;i<distanceMatrix[j].length;i++){
		        	if (i!=j && distanceMatrix[j][i] <= dist_cutoff)
		        	output.append(String.format("%d,%d,1,Undirected,%.1f\n",j,i,distanceMatrix[j][i]));
		        }
	        }
	        CommonUtils.writeFile(String.format("%s.k%d.dist%.1f.gephi_edges.csv", outName, k, dist_cutoff), output.toString());
		}
		if (config.verbose>1)
			System.out.println(" " + CommonUtils.timeElapsed(tic));
		
		return distanceMatrix;
	}
	
	/**
	 * Compute Weighted Distance Matrix<br>
	 * Gapped k-mer distance as Levenshtien distance of matrix form
	 */
	public float[][] computeWeightedDistanceMatrix2(ArrayList<Kmer> kmers, boolean print_dist_matrix) {
		int kmerCount = kmers.size();
		float[][]distanceMatrix = new float[kmerCount][kmerCount];
		for (int i=0;i<kmerCount;i++){
			Kmer kmi = kmers.get(i);
			distanceMatrix[i][i] = 0;
			for (int j = 0; j < i; j++) {
				Kmer kmj = kmers.get(j);
				float dist = KMAC.editDistance(kmi, kmj);
				distanceMatrix[i][j] = dist;
				distanceMatrix[j][i] = dist;
			}
		}
		if (print_dist_matrix){
	        StringBuilder output = new StringBuilder();
	        for (int j=0;j<kmerCount;j++){
	        	output.append(String.format("%s\t",kmers.get(j).kmerString));
	        }
	        CommonUtils.replaceEnd(output, '\n');
	        for (int j=0;j<kmerCount;j++){
	        	output.append(String.format("%d\t",kmers.get(j).getNetHitCount(posNegSeqRatio)));
	        }
	        CommonUtils.replaceEnd(output, '\n');
	        for (int j=0;j<kmerCount;j++){
		        for (int i=0;i<distanceMatrix[j].length-1;i++){
		        	output.append(String.format("%.1f\t",distanceMatrix[j][i]));
		        }
		        output.append(String.format("%.1f",distanceMatrix[j][kmerCount-1])).append("\n");
	        }
	        CommonUtils.writeFile(outName+".weighted_distance_matrix.txt", output.toString());
		}
		return distanceMatrix;
	}
	
	/**
	 * This method computes the edit distance between 2 k-mers<br>
	 * It is similar to Levenshtein distance, but not considering internal insertion/deletion.<br>
	 * It takes the min of distances from forward and reverse compliment orientation.
	 */
	public static int numDistCalcuation = 0;
	public static float editDistance(Kmer k1, Kmer k2) {
		if (k1.kmerString.length() > k2.kmerString.length()) {
			return KMAC.editDistance(k2, k1);
		}
		// k1 is the shorter kmer (or they are equal length)
		float[][] m1 = k1.getMatrix();
		float[][] m2 = k2.getMatrix();
		float[][] m2RC = k2.getMatrixRC();
		float forwardDistance = KMAC.editDistanceByMatrix(m1, m2, m1.length+m2.length);
		float best = KMAC.editDistanceByMatrix(m1, m2RC, forwardDistance);
//		double forwardDistance2 = KMAC1.editDistanceByMatrix(m1, m2);
//		double best2 = Math.min(KMAC1.editDistanceByMatrix(m1, m2RC), forwardDistance2);
//		assert (best==best2);
		numDistCalcuation++;

		return best;
	}
	
	/**
	 * This method computes the edit distance using matrix representation for gapped k-mers<br>
	 * It is similar to Levenshtein distance, but not considering internal insertion/deletion.
	 */
	public static float editDistanceByMatrix(float[][] m1, float[][] m2) {
		int dist = m1.length + m2.length; // final distance
		for (int i = 0; i < m1.length + m2.length - 1; i++) {
			// we slide m1 along m2, where i is the index at m2 that m1's max index is aligned with
			// for example, for i = 0; position m1.length - 1 in m1 is aligned with position 0 in m2
			int slideDist = 0; // this is the total distance for this given slide of m1 along m2
			
			// for each case, we first add the hanging ends as part of the edit distance
			// next, we compare the overlapping indices to find the dist
			if (i < m1.length - 1) {
				slideDist += m1.length + m2.length - 2 * i - 2;
				for (int j = 0; j < i + 1; j++) {
					int compare1 = m1.length - 1 - i + j;
					int compare2 = j;
					float compareDist = 0;
					for (int k = 0; k < 4; k++) {
						compareDist += Math.abs(m1[compare1][k] - m2[compare2][k]);
					}
					compareDist /= 2;
					slideDist += compareDist;
				}
			}
			else if (i > m2.length - 1) {
				slideDist += 2 * i + 2 - m1.length - m2.length;
				for (int j = 0; j < m1.length + m2.length - 1 - i; j++) {
					int compare1 = j;
					int compare2 = i + 1 - m1.length + j;
					float compareDist = 0;
					for (int k = 0; k < 4; k++) {
						compareDist += Math.abs(m1[compare1][k] - m2[compare2][k]);
					}
					compareDist /= 2;
					slideDist += compareDist;
				}
			}
			else {
				slideDist += m2.length - m1.length;
				for (int j = 0; j < m1.length; j++) {
					int compare1 = j;
					int compare2 = i - m1.length + 1 + j;
					float compareDist = 0;
					for (int k = 0; k < 4; k++) {
						compareDist += Math.abs(m1[compare1][k] - m2[compare2][k]);
					}
					compareDist /= 2;
					slideDist += compareDist;
				}
			}
			dist = Math.min(slideDist, dist);
		}
		
		return dist;
	}
	/**
	 * This method computes the edit distance using matrix representation for gapped k-mers<br>
	 * It is similar to Levenshtein distance, but not considering internal insertion/deletion.<br>
	 * It stops if the distance is larger than the cutoff value.
	 */
	public static float editDistanceByMatrix(float[][] m1, float[][] m2, float cutoff) {
		float dist = m1.length + m2.length; // final distance
		
		// fill the sliding indexs such that the sliding start from maximum overlap of m1 and m2
		int[] idxs = new int [m1.length + m2.length-2];
		int mid = idxs.length/2;
		for (int i = 0; i < mid; i++){
			idxs[i*2] = mid+i;
			if (i*2+1<idxs.length)
				idxs[i*2+1] = mid-1-i;
		}
eachSliding:for (int it = 0; it < idxs.length; it++) {
			int i = idxs[it];
			// we slide m1 along m2, where i is the index at m2 that m1's max index is aligned with
			// for example, for i = 0; position m1.length - 1 in m1 is aligned with position 0 in m2
			float slideDist = 0; // this is the total distance for this given slide of m1 along m2
			
			// for each case, we first add the hanging ends as part of the edit distance
			// next, we compare the overlapping indices to find the dist
			if (i < m1.length - 1) {	// if m1 left is hanging
				slideDist += m1.length + m2.length - 2 * i - 2;
				if (slideDist>=cutoff){
					dist = Math.min(cutoff, dist);
					break eachSliding;
				}
				for (int j = 0; j < i + 1; j++) {
					int compare1 = m1.length - 1 - i + j;
					int compare2 = j;
					float compareDist = 0;
					for (int k = 0; k < 4; k++) {
						compareDist += Math.abs(m1[compare1][k] - m2[compare2][k]);
					}
					compareDist /= 2;
					slideDist += compareDist;
					if (slideDist>=cutoff){
						dist = Math.min(cutoff, dist);
						continue eachSliding;
					}
				}
			}
			else if (i > m2.length - 1) {	// if m1 right is hanging (i.e. slide passed m2)
				slideDist += 2 * i + 2 - m1.length - m2.length;
				if (slideDist>=cutoff){
					dist = Math.min(cutoff, dist);
					break eachSliding;
				}
				for (int j = 0; j < m1.length + m2.length - 1 - i; j++) {
					int compare1 = j;
					int compare2 = i + 1 - m1.length + j;
					float compareDist = 0;
					for (int k = 0; k < 4; k++) {
						compareDist += Math.abs(m1[compare1][k] - m2[compare2][k]);
					}
					compareDist /= 2;
					slideDist += compareDist;
					if (slideDist>=cutoff){
						dist = Math.min(cutoff, dist);
						continue eachSliding;
					}
				}
			}
			else {	// all m1 positions are aligned with m2 positions
				slideDist += m2.length - m1.length;
				if (slideDist>=cutoff){
					dist = Math.min(cutoff, dist);
					break eachSliding;
				}
				for (int j = 0; j < m1.length; j++) {
					int compare1 = j;
					int compare2 = i - m1.length + 1 + j;
					float compareDist = 0;
					for (int k = 0; k < 4; k++) {
						compareDist += Math.abs(m1[compare1][k] - m2[compare2][k]);
					}
					compareDist /= 2;
					slideDist += compareDist;
					if (slideDist>=cutoff){
						dist = Math.min(cutoff, dist);
						continue eachSliding;
					}
				}
			}
			dist = Math.min(slideDist, dist);
			cutoff = dist;
		}
		return dist;	
	}	
	
	public static double ycDistance(Kmer k1, Kmer k2) {
		int numOperations = 1;
		ArrayList<Kmer> sk1 = new ArrayList<Kmer>();
		ArrayList<Kmer> sk2 = new ArrayList<Kmer>();
		if (k1 instanceof GappedKmer) {
			sk1.addAll(((GappedKmer)k1).getBaseKmers());
			numOperations *= ((GappedKmer)k1).getBaseKmers().size();
		}
		else {
			sk1.add(k1);
		}
		if (k2 instanceof GappedKmer) {
			sk2.addAll(((GappedKmer)k2).getBaseKmers());
			numOperations *= ((GappedKmer)k2).getBaseKmers().size();
		}
		else {
			sk2.add(k2);
		}
		
		double distSum = 0;
		double weightSum = 0;
		
		for (Kmer s1: sk1) {
			for (Kmer s2: sk2) {
				double w2 = k1.getPosHitCount()*k2.getPosHitCount();
				weightSum += w2;
				distSum += w2 * MTree.testDist(s1.kmerString, s2.kmerString);
				//distSum += w2 * CommonUtils.strMinDistanceWithCutoff(s1.getKmerString(), s2.getKmerString(), cutoff);
			}
		}
		//System.out.println(numOperations);
		return distSum/weightSum;
	}
	
	public static double ycDistanceWC(Kmer k1, Kmer k2, double posNegRatio, int cutoff) {
		// this one has cutoff
		int numOperations = 1;
		ArrayList<Kmer> sk1 = new ArrayList<Kmer>();
		ArrayList<Kmer> sk2 = new ArrayList<Kmer>();
		if (k1 instanceof GappedKmer) {
			sk1.addAll(((GappedKmer)k1).getBaseKmers());
			numOperations *= ((GappedKmer)k1).getBaseKmers().size();
		}
		else {
			sk1.add(k1);
		}
		if (k2 instanceof GappedKmer) {
			sk2.addAll(((GappedKmer)k2).getBaseKmers());
			numOperations *= ((GappedKmer)k2).getBaseKmers().size();
		}
		else {
			sk2.add(k2);
		}
		
		double distSum = 0;
		double weightSum = 0;
		
		for (Kmer s1: sk1) {
			for (Kmer s2: sk2) {
				double w2 = k1.getNetHitCount(posNegRatio)*k2.getNetHitCount(posNegRatio);
				weightSum += w2;
				distSum += w2 * MTree.testDist(s1.kmerString, s2.kmerString);
				//distSum += w2 * CommonUtils.strMinDistanceWithCutoff(s1.getKmerString(), s2.getKmerString(), cutoff);
			}
		}
		//System.out.println(numOperations);
		return distSum/weightSum;
	}
	
	private ArrayList<Kmer> densityClusteringWithDistMatrix(ArrayList<Kmer> kmers, float[][] distanceMatrix, double distance_cutoff) {
		if (config.verbose>1)
			System.out.println(String.format("Density clustering k-mers, critical distance (dc) = %.1f", distance_cutoff));		
		long tic = System.currentTimeMillis();
		int k = kmers.get(0).getK();
		
		ArrayList<BitSet> posHitList = new ArrayList<BitSet>();
		ArrayList<BitSet> negHitList = new ArrayList<BitSet>();
		for (Kmer km:kmers){
			posHitList.add(km.posBits);
			negHitList.add(km.negBits);
        }
		
		ArrayList<DensityClusteringPoint> centers = hitWeightedDensityClustering(distanceMatrix, 
				posHitList, negHitList, seq_weights, posNegSeqRatio, distance_cutoff, config.k_top*3,
				config.refine_centerKmers, config.use_self_density);
		ArrayList<Kmer> results = new ArrayList<Kmer>();
		boolean topKmerIsNotCenter = true;
		for (DensityClusteringPoint p:centers){
			results.add(kmers.get(p.id));
			if (p.id == 0)
				topKmerIsNotCenter = false;
		}
		if (topKmerIsNotCenter)	{	// add the strongest k-mer if it has not been selected as a cluster center
			results.add(kmers.get(0));
//			System.err.println("Add top kmer!!");
		}
		results.trimToSize();

//		System.out.println(String.format("cluster_num=%d, kmer_distance<=%d, delta>=%d", centers.size(), dc, delta));
		System.out.println("\nTop "+config.k_top+" cluster center k-mers:");
		System.out.println(Kmer.toShortHeader(k)+"\tId\tDensity\tDelta\tGamma\tCluster_size");
		int displayCount = Math.min(config.k_top, centers.size());
		for (int i=0;i<displayCount;i++){
			DensityClusteringPoint p = centers.get(i);
			System.out.println(String.format("%s    \t%d\t%.1f\t%.1f\t%.1f\t%d",
					results.get(i).toShortString(), p.id, p.density, p.delta, p.gamma, p.members.size()));
		}
		if (topKmerIsNotCenter)
			System.out.println(results.get(results.size()-1).toShortString());
		if (config.verbose>1)
			System.out.println(CommonUtils.timeElapsed(tic));
		
		return results;
	}
	
	private ArrayList<Kmer> densityClusteringWithMTree(ArrayList<Kmer> kmers, MTree dataPoints, double distance_cutoff) {
		if (config.verbose>1)
			System.out.println(String.format("Density clustering k-mers, critical distance (dc) = %.1f", distance_cutoff));		

		long tic = System.currentTimeMillis();
		int k = kmers.get(0).getK();
		
		ArrayList<BitSet> posHitList = new ArrayList<BitSet>();
		ArrayList<BitSet> negHitList = new ArrayList<BitSet>();
		for (Kmer km:kmers){
			posHitList.add(km.posBits);
			negHitList.add(km.negBits);
        }

//		if (config.verbose>1)
//			System.out.println("distance_cutoff="+distance_cutoff);
		ArrayList<DensityClusteringPoint> centers = KMAC.hitWeightedDensityClusteringMTree(dataPoints, config.mtree, 
				posHitList, negHitList, seq_weights, posNegSeqRatio, distance_cutoff, config.k_top*3,
				config.refine_centerKmers, config.use_self_density);
		ArrayList<Kmer> results = new ArrayList<Kmer>();
		//TODO: for each cluster, select the strongest kmer that is far away from other stronger cluster center kmers
		for (DensityClusteringPoint p:centers){
			// find the best representative kmer
			Kmer km = kmers.get(p.id);
			results.add(km);
		}
		results.trimToSize();

		if (config.mtree!=-1){
			System.out.println("\nTop "+config.k_top+" cluster center k-mers:");
			System.out.println(Kmer.toShortHeader(k)+"\tId\tDensity\tDelta\tGamma\tCluster_size");
			int displayCount = Math.min(config.k_top, centers.size());
			for (int i=0;i<displayCount;i++){
				DensityClusteringPoint p = centers.get(i);
				System.out.println(String.format("%s    \t%d\t%.1f\t%.1f\t%.1f\t%d",
						results.get(i).toShortString(), p.id, p.density, p.delta, p.gamma, p.members.size()));
			}
			if (config.verbose>1)
				System.out.println(CommonUtils.timeElapsed(tic));
		}
		return results;
	}	
	/**
	 * This implements the clustering method introduced by "Clustering by fast search and find of density peaks, Science. 2014 Jun 27;344(6191):"<br>
	 * It is extended to incorporate weights for the data points. <br> 
	 * Currently it is not used by KMAC
	 * @param distanceMatrix
	 * @param weights
	 * @param distanceCutoff kmer distance cutoff, kmers with equal or less distance are consider neighbors when computing local density
	 * @param deltaCutoff delta value cutoff, kmers with equal or higher delta values are used for selecting cluster centers
	 * @return
	 */
	public static ArrayList<DensityClusteringPoint> weightedDensityClustering(float[][] distanceMatrix, double[] weights, int distanceCutoff, int deltaCutoff){
		ArrayList<DensityClusteringPoint> data = new ArrayList<DensityClusteringPoint>();
		StatUtil util = new StatUtil();
		// compute local density for each point
		for (int i=0;i<distanceMatrix.length;i++){
			DensityClusteringPoint p = util.new DensityClusteringPoint();
			p.id = i;
			float[] distanceVector = distanceMatrix[i];
			for (int j=0;j<distanceVector.length;j++){
				if (distanceVector[j] <= distanceCutoff)
					p.density+=weights[j];				
			}
			data.add(p);
		}
		data.trimToSize();		
		Collections.sort(data);
		
		// compute the shortest distance to the potential centers
		// for the top point
		data.get(0).delta = StatUtil.getMax(distanceMatrix[data.get(0).id]);
		if (data.get(0).delta <= deltaCutoff)		// in rare situation, at least return the top point
			data.get(0).delta = deltaCutoff + 1;
		data.get(0).gamma = data.get(0).delta * data.get(0).density;
		
		// for the rest of points
		for (int i=1;i<data.size();i++){
			float min = Float.MAX_VALUE;
			int id = data.get(i).id;
//			if (id==284)
//				id=id;
			for (int j=0;j<i;j++){		// for the points have higher (or equal) density than point i
//				if (data.get(i).density<data.get(j).density 	// need this comparison because of tied ranking
				if (min>distanceMatrix[id][data.get(j).id])
					min = distanceMatrix[id][data.get(j).id];
			}
			data.get(i).delta = min;
			data.get(i).gamma = data.get(i).delta * data.get(i).density;
		}
		
//		Collections.sort(data, new Comparator<DensityClusteringPoint>(){
//            public int compare(DensityClusteringPoint o1, DensityClusteringPoint o2) {
//                return o1.compareToById(o2);
//            }
//        });
//		for (int i=0;i<50;i++){
//			DensityClusteringPoint p = data.get(i);
//			System.out.println(String.format("#%d\t%.1f\t%.1f\t%.1f\t%.1f",p.id, weights[p.id], p.density, p.delta, p.gamma));
//		}
		
		// sort results by gamma, excluding points with delta smaller than the cutoff value
		ArrayList<DensityClusteringPoint> data_small_delta = new ArrayList<DensityClusteringPoint>();
		for (DensityClusteringPoint p: data){
			if (p.delta<deltaCutoff)
				data_small_delta.add(p);
		}
		data.removeAll(data_small_delta);
		Collections.sort(data, new Comparator<DensityClusteringPoint>(){
            public int compare(DensityClusteringPoint o1, DensityClusteringPoint o2) {
                return o1.compareToByGamma(o2);
            }
        });
		data.trimToSize();
		
		return data;
	}
	
	/**
	 * This implements the clustering method introduced by "Clustering by fast search and find of density peaks, Science. 2014 Jun 27;344(6191):"<br>
	 * It is extended to incorporate weights for the data points. <br>
	 * This method use net (positive-negative) hit as weight (e.g. sequence hit by the k-mers)
	 * @param distanceMatrix
	 * @param weights
	 * @param distanceCutoff kmer distance cutoff, kmers with equal or less distance are consider neighbors when computing local density
	 * @param deltaCutoff delta value cutoff, kmers with equal or higher delta values are used for selecting cluster centers
	 * @return Return points with high delta values, then points that are far from the high delta points
	 */
	public static ArrayList<DensityClusteringPoint> hitWeightedDensityClustering(float[][] distanceMatrix, 
			ArrayList<BitSet> posHitList, ArrayList<BitSet> negHitList, double[] seq_weights, 
			double posNegSeqRatio, double distanceCutoff, int numCenters, boolean refine_centers, boolean use_self_density){
		ArrayList<DensityClusteringPoint> data = new ArrayList<DensityClusteringPoint>();
		StatUtil util = new StatUtil();
		// compute local density for each point
		for (int i=0;i<distanceMatrix.length;i++){
			DensityClusteringPoint p = util.new DensityClusteringPoint();
			p.id = i;
			// self_density: individual k-mer hit count
			double self_density = CommonUtils.calcWeightedHitCount(posHitList.get(i),seq_weights) - negHitList.get(i).cardinality()*posNegSeqRatio;
			float[] distanceVector = distanceMatrix[i];
			// sum up to get total hit count of this point and its neighbors
			BitSet b_pos = new BitSet();
			BitSet b_neg = new BitSet();
			for (int j=0;j<distanceVector.length;j++){
				if (distanceVector[j] <= distanceCutoff){
					b_pos.or(posHitList.get(j));
					b_neg.or(negHitList.get(j));
				}
			}
			if (use_self_density)
				// group hit count * self_density, to down-weight weak kmers from being selected as center
				p.densitySxN = (float)Math.sqrt((CommonUtils.calcWeightedHitCount(b_pos,seq_weights) - b_neg.cardinality()*posNegSeqRatio) * self_density);
			else
				// group hit count 
				p.densitySxN = (float)(CommonUtils.calcWeightedHitCount(b_pos,seq_weights) - b_neg.cardinality()*posNegSeqRatio);
			if (Double.isNaN(p.densitySxN))
				continue;
			p.density = p.densitySxN;
//			p.density = b_pos.cardinality()-b_neg.cardinality()*posNegSeqRatio;
			data.add(p);
		}
		data.trimToSize();	
		Collections.sort(data);
		
		// compute the shortest distance to the potential centers
		// for the top point
		data.get(0).delta = StatUtil.getMax(distanceMatrix[data.get(0).id]);
		data.get(0).gamma = data.get(0).delta * data.get(0).density;
		data.get(0).delta_id = 0;
		data.get(0).members.add(data.get(0));
//		data.get(0).memberIds.add(data.get(0).id);
		
		// for the rest of points
		for (int i=1;i<data.size();i++){
			float min = Float.MAX_VALUE;
			DensityClusteringPoint point_i = data.get(i);
			point_i.members.add(point_i);
			int id = point_i.id;
//			if (id==91)
//				id=id;
			// find the nearest stronger point of i, point i to it
			for (int j=0;j<i;j++){		// for the points j have higher (or equal) density than point i
				if (distanceMatrix[id][data.get(j).id] < min){
					min = distanceMatrix[id][data.get(j).id];
					point_i.delta_id = data.get(j).id;
				}
			}
			for (int j=0;j<i;j++){		// mark point i as the members of all stronger points with small distance
				if (distanceMatrix[id][data.get(j).id] <= distanceCutoff){
					data.get(j).members.add(point_i);
//					data.get(j).memberIds.add(id);
				}
			}
			point_i.delta = min;
			point_i.gamma = point_i.delta * point_i.density;	
		}
		for (int i=0;i<data.size();i++){
			if (data.get(i).delta > distanceCutoff)
				data.get(i).delta_id = data.get(i).id;
		}

		Collections.sort(data, new Comparator<DensityClusteringPoint>(){
            public int compare(DensityClusteringPoint o1, DensityClusteringPoint o2) {
                return o1.compareToByGamma(o2);
            }
        });

		double center_min_distance = distanceCutoff * 1.5;
		ArrayList<DensityClusteringPoint> results = new ArrayList<DensityClusteringPoint>();
		if (refine_centers){
			while(!data.isEmpty()){
				DensityClusteringPoint p = data.get(0);
				data.remove(0);
				for (DensityClusteringPoint r:results){
					if (distanceMatrix[p.id][r.id]<center_min_distance){
						p = null;;
						break;
					}
				}
				if (p != null){
					results.add(p);
					if (results.size()>=numCenters)
						break;
				}
			}
		}
		else{
			for (int i=0;i<numCenters;i++)
				results.add(data.get(i));
		}

		return results;
	}	
		
	/** 
	 * Density clustering with m-tree data structre<br>
	 * posHitList, negHitList, and seq_weights are assumed to be in the same indices as m-tree indices of the k-mers
	 * @param mtreeDataPoints
	 * @param posHitList
	 * @param negHitList
	 * @param seq_weights
	 * @param posNegSeqRatio
	 * @param distanceCutoff
	 * @return
	 */
	public static ArrayList<DensityClusteringPoint> hitWeightedDensityClusteringMTree(MTree mtreeDataPoints, int mtree,
			ArrayList<BitSet> posHitList, ArrayList<BitSet> negHitList, double[] seq_weights, 
			double posNegSeqRatio, double distanceCutoff, int numCenters, boolean refine_centers, boolean use_self_density){
		// keys = kmer indices (0 to dataPoints.getSize() - 1)
		// values = kmers in range
		HashMap<Integer, ArrayList<Kmer>> rangeResults = new HashMap<Integer, ArrayList<Kmer>>();
		ArrayList<DensityClusteringPoint> data = new ArrayList<DensityClusteringPoint>();
		StatUtil util = new StatUtil();
		int tmp = numDistCalcuation;
		ArrayList<TreeObject> traversal = MTree.traverse(mtreeDataPoints.getRoot());
		for (TreeObject o : traversal) {
			ArrayList<Kmer> rangeResult = mtreeDataPoints.rangeSearch(o.getData(), distanceCutoff);
			for (Kmer kmer : rangeResult) {
				if (rangeResults.get(o.getData().getIndex()) != null) {
					rangeResults.get(o.getData().getIndex()).add(kmer);
				}
				else {
					rangeResults.put(o.getData().getIndex(), new ArrayList<Kmer>());
					rangeResults.get(o.getData().getIndex()).add(kmer);
				}
				if (rangeResults.get(kmer.getIndex()) != null) {
					rangeResults.get(kmer.getIndex()).add(o.getData());
				}
				else {
					rangeResults.put(kmer.getIndex(), new ArrayList<Kmer>());
					rangeResults.get(kmer.getIndex()).add(o.getData());
				}
			}
			MTreeNode container = o.getContainer();
			container.getObjects().remove(o.getIndex());
			if (container.getObjects().size() == 0) {
				if (container.getParent() != null) {
					container.getParent().setChild(null);
				}
			}
		}
//		TODO: printout m-tree performance information
		if (mtree==-1){
			System.out.print("n_range="+(numDistCalcuation-tmp)+"\t");
			tmp = numDistCalcuation;
		}
		
		// compute density using the range search results
		// dataPoints arrayList is indexed by k-mer id (i.e. treeIndex)
		for (int i = 0; i < mtreeDataPoints.getSize(); i++) {
			DensityClusteringPoint p = util.new DensityClusteringPoint();
			p.id = i;
			// self_density: individual k-mer hit count
			double self_density = CommonUtils.calcWeightedHitCount(posHitList.get(i),seq_weights) - negHitList.get(i).cardinality()*posNegSeqRatio;
			ArrayList<Kmer> inRange = rangeResults.get(i);
			// sum up to get total hit count of this point and its neighbors
			BitSet b_pos = new BitSet();
			BitSet b_neg = new BitSet();
			for (Kmer kmer: inRange) {
				b_pos.or(posHitList.get(kmer.getIndex()));
				b_neg.or(negHitList.get(kmer.getIndex()));
			}
			if (use_self_density)
				// group hit count * self_density, to down-weight weak kmers from being selected as center
				p.densitySxN = (float)Math.sqrt((CommonUtils.calcWeightedHitCount(b_pos,seq_weights) - b_neg.cardinality()*posNegSeqRatio) * self_density);
			else
				// group hit count 
				p.densitySxN = (float)(CommonUtils.calcWeightedHitCount(b_pos,seq_weights) - b_neg.cardinality()*posNegSeqRatio);
			if (Double.isNaN(p.densitySxN))
				continue;
			p.density = p.densitySxN;
//				p.density = b_pos.cardinality()-b_neg.cardinality()*posNegSeqRatio;
			data.add(p);
//			System.out.println(i);
		}
		data.trimToSize();	
		Collections.sort(data);		// by default, sort by density

		int[] kmerIdx2dataIdx = new int[mtreeDataPoints.getSize()];	// indexed by kmer index, to find the index of data list.
		for (int i=0;i<kmerIdx2dataIdx.length;i++)
			kmerIdx2dataIdx[i] = -1;		// init
		for (int i=0;i<data.size();i++)
			kmerIdx2dataIdx[data.get(i).id] = i;
		
		/** compute the shortest distance to the potential centers */
		
		// for the 2nd to all the rest of points
		float maxDist = -1;
		for (int i=1;i<data.size();i++){
			DensityClusteringPoint point_i = data.get(i);
			point_i.members.add(point_i);
//			point_i.memberIds.add(i);
			int id = point_i.id;
			Kmer kmer = mtreeDataPoints.getData().get(id);
			double density_i = point_i.density;

//			if(kmer.getKmerString().equals("GGGGAGGGG"))
//				i=i;

			// find the nearest stronger point (have higher or equal density than point i), point i to it
			ArrayList<Kmer> inRange = rangeResults.get(id);
			ArrayList<Kmer> stronger = new ArrayList<Kmer>();
			for (Kmer km: inRange){
				int idx = kmerIdx2dataIdx[km.getIndex()];
				if (idx==-1 || idx==i)
					continue;
				if (data.get(idx).density>=density_i)
					stronger.add(km);
			}
			
			if (stronger.isEmpty()){	// all the inRange neighbors are not stronger, do the exhaustive search
				float min = Float.MAX_VALUE;
				for (int j=0;j<i;j++){		// for the points j have higher (or equal) density than point i
					int jd = data.get(j).id;
					float ijDistance = KMAC.editDistance(kmer, mtreeDataPoints.getData().get(jd));
					if (ijDistance < min) {
						min = ijDistance;
//						data.get(i).delta_id = data.get(j).id;
					}
					if (ijDistance <= distanceCutoff) {
						data.get(j).members.add(point_i);
//						data.get(j).memberIds.add(i);
					}
				}
				point_i.delta = min;
			}
			else{
				float min = Float.MAX_VALUE;
				for (Kmer km:stronger){		// for each stronger point
					int j = kmerIdx2dataIdx[km.getIndex()];
					float ijDistance = KMAC.editDistance(kmer, mtreeDataPoints.getData().get(km.getIndex()));
					if (ijDistance < min) {
						min = ijDistance;
//						data.get(i).delta_id = data.get(j).id;
					}
					if (ijDistance <= distanceCutoff) {
						data.get(j).members.add(point_i);
//						data.get(j).memberIds.add(i);
					}
				}
				point_i.delta = min;
			}
			point_i.gamma = point_i.delta * point_i.density;
			if (point_i.delta > maxDist)
				maxDist = point_i.delta;
		}
//		for (int i=0;i<data.size();i++){
//			if (data.get(i).delta > distanceCutoff)
//				data.get(i).delta_id = data.get(i).id;
//		}
		
		// for the top point, get max distance
		data.get(0).delta = maxDist;
		data.get(0).gamma = data.get(0).delta * data.get(0).density;
//		data.get(0).delta_id = 0;
		data.get(0).members.add(data.get(0));
//		data.get(0).memberIds.add(data.get(0).id);
		
		if (mtree==-1){
			System.out.print("n_delta="+(numDistCalcuation-tmp));
			System.out.print("\tn_total="+numDistCalcuation);
		}
		Collections.sort(data, new Comparator<DensityClusteringPoint>(){
            public int compare(DensityClusteringPoint o1, DensityClusteringPoint o2) {
                return o1.compareToByGamma(o2);
            }
        });
		double center_min_distance = distanceCutoff * 1.5;
		ArrayList<DensityClusteringPoint> results = new ArrayList<DensityClusteringPoint>();
		if (refine_centers){
			while(!data.isEmpty()){
				DensityClusteringPoint p = data.get(0);
				data.remove(0);
				for (DensityClusteringPoint r:results){
					if (KMAC.editDistance(mtreeDataPoints.getData().get(r.id), mtreeDataPoints.getData().get(p.id))<center_min_distance){
						p = null;;
						break;
					}
				}
				if (p != null){
					results.add(p);
					if (results.size()>=numCenters)
						break;
				}
			}
		}
		else{
			for (int i=0;i<numCenters;i++)
				results.add(data.get(i));
		}

		return results;
	}

/**
 * This is the main method for KMAC v1 (density clustering, with gappedKmer) motif discovery<br>
 * It discovers one single motif starting from the seed k-mer and return. It is called with different seed k-mers for multiple motifs for all k values.<br>
 * A few important points:<br>
 * 1. Definition of k-mer positions<br>
 * 		The start of the seed k-mer is the global reference point. All the other positions can be represented as relatve to seed position.<br>
 * 		Shift: kmer_seed = start_of_kmer - start_of_seed. The start of a k-mer relative to the start of the seed k-mer when they are aligned<br>
 * 		StartOffset: kmer_bs = kmer_seed - seed_bs, The start of a k-mer relative to the binding site (BS) position (esitmated by averaging the sequence midpoints)<br>
 * 					During KMAC motif finding, startOffset is set as kmer.Shift in extractKSM(). After the motif is found, it is set to be relative to BS.<br>
 * 2. Definition of the sequence alignment position<br>
 * 		seq_seed = -(seed_seq) = -(seed_kmer - kmer_seq) = kmer_seed + kmer_seq. The start of the sequence relative to the start of the seed k-mer. This is done by linking the position of the matched k-mer group (or simply just a k-mer) to the seed k-mer.<br>
 * 3. PWM position:<br>
 * 		PWM position is the start position of the PWM. Cluster.pos_pwm_seed = pwm - seed, for converting pwm pos to seed position.
 * 4. Orientation of k-mers and sequences<br>
 * 		Each k-mer is merged with its reverse-compliment k-mer. The k-mer string stores one orientation. In KMAC1, we do NOT flip the k-mer string. Instead, the Kmer object keep track of the orientation as seedOrientation (whether the k-mer is in the orientation consistent with the seed k-mer orientation, i.e. they are extracted from the sequence alignment built from the seed k-mer).<br>
 * 		The baseKmers hashTable of a gapped k-mer, stores the info of whether the base-kmers are the same orientationo as the gapped-kmer. Because a base-kmer could potentially be the base-kmers of gapped k-mers that has conflicting orientation. Thus needing to store orientation info in each gappedKmer-baseKmer pair.<br>
 * 		The orientation of a sequence is stored in the Sequence object, to note whether it is in the forward orientation as the original input sequence.
 */
	public MotifCluster KmerMotifAlignmentClustering (ArrayList<Sequence> seqList, ArrayList<Kmer> kmers, Kmer seed, int k){
//		if (seed.getKmerString().equals("AGCAGCCAG"))
//			k+=0;
		tic = System.currentTimeMillis();
		int seed_range = k;
		MotifCluster cluster = new MotifCluster();
		
		if (kmers.size()==0)
			return null;

		if (config.verbose>1)
			System.out.println(CommonUtils.timeElapsed(tic)+": kmer num = "+kmers.size());

		// make a copy so that the modified k-mer information does not affect other motifClusters
		cluster.inputKmers = Kmer.deepCloneKmerList(kmers, seed, seq_weights);
		cluster.seedKmer = cluster.inputKmers.get(0);
		cluster.k = k;
		kmers = cluster.inputKmers;
		seed = cluster.seedKmer;
		for (Kmer km: kmers){
			km.setShift(0);
			km.setKmerStartOffset(0);
			km.setSeedOrientation(true);
		}
//		if (seed.kmerString.equals("CNTTGTNNT"))
//			k=k+0;
		
		/** init kmerSet with seed family of seed kmer, by adding mismatch k-mers, order by #mm */
		ArrayList<Kmer> seedFamily = findSeedFamily(kmers, seed.kmerString);
		seed.setAlignString(seed.kmerString);
		seedFamily.add(seed);
		
		Collections.sort(seedFamily);
		ArrayList<Kmer> tmp = new ArrayList<Kmer>();
		seed.setMatrix();
		for (Kmer km: seedFamily){
			km.setMatrix();
			if (editDistance(seed, km) <= 2.5)
				tmp.add(km);
//			System.out.println(String.format("%.2f\t%s", editDistance(seed, km), km.toString2()));
		}
		for (Kmer km: seedFamily)
			km.clearMatrix();
		seedFamily = tmp;
		
		// init AC search engine and init coveredWidth = k
		if (config.kg_hit_adjust_type==2){
			if (posCoveredWidth==null || posCoveredWidth.length!=seqList.size())
				posCoveredWidth = new int[seqList.size()];
			for (int i=0;i<posCoveredWidth.length;i++)
				posCoveredWidth[i] = k;
			if (negCoveredWidth==null || negCoveredWidth.length!=seqListNeg.size())
				negCoveredWidth = new int[seqListNeg.size()];
			for (int i=0;i<negCoveredWidth.length;i++)
				negCoveredWidth[i] = k;
		}
		else if (config.kg_hit_adjust_type==1){
			if (posHitStrings==null || posHitStrings.length!=seqList.size())
				posHitStrings = new String[seqList.size()][];
			for (int i=0;i<posHitStrings.length;i++)
				posHitStrings[i] = new String[]{""};
			if (negHitStrings==null || negHitStrings.length!=seqListNeg.size())
				negHitStrings = new String[seqListNeg.size()][];
			for (int i=0;i<negHitStrings.length;i++)
				negHitStrings[i] = new String[]{""};
		}
		
		initAhoCorasick(seedFamily);
		
		KmerGroup kg = config.use_weighted_kmer ? 
				new KmerGroup(posCoveredWidth, negCoveredWidth, seedFamily, 0, seq_weights) : 
				new KmerGroup(posCoveredWidth, negCoveredWidth, seedFamily, 0, null);
		if (config.verbose>1) {
//			System.out.println(CommonUtils.timeElapsed(tic)+": evaluateKsmROC(seedFamily)");
			Pair<double[],double[] > pair = scoreKsmSequences (seqList, seqListNeg, seedFamily);
			double kAUC = evaluateMotifROC(pair.car(), pair.cdr(), config.fpr).motif_significance;
			System.out.println(CommonUtils.timeElapsed(tic)+String.format(": Seed family %d kmers, motif kAUC = %.1f", seedFamily.size(), kAUC));
		}
		
		/** get sequence alignment using seedFamily kmer positions */
//		System.out.println(getKmerClusterAlignmentString(seedFamily, seedFamily.size()));
    	cluster.alignedKmers = seedFamily;
    	MotifThreshold thresh = new MotifThreshold();
    	if (config.use_odds_ratio)
    		thresh.motif_cutoff = 1;
    	else
    		thresh.motif_cutoff = 3;
    	thresh.posHit = kg.getGroupHitCount();
    	thresh.negHit = kg.getGroupNegHitCount();
    	thresh.motif_significance = 0;			// TODOTODO
    	cluster.ksmThreshold = thresh;
		int hitCount = alignByKSM (seqList, seedFamily, cluster).size();
		if (hitCount<seqs.length*config.motif_hit_factor){
			if (config.verbose>1)
	    		System.out.println(String.format("%s: Seed and family kmers match too few (%d) sequences, stop here.", CommonUtils.timeElapsed(tic), hitCount));
			return null;
		}

    	seedFamily = null;
    	
    	// build first PWM
    	if (config.verbose>1 && config.pwm_noise!=0)
			System.out.println(CommonUtils.timeElapsed(tic)+ ": PWM noise = " + config.pwm_noise);
    	NewPWM newPWM = buildPWM(seqList, cluster, config.pwm_noise, tic, true);	
		
    	if (newPWM!=null){
			newPWM.updateClusterPwmInfo(cluster);
			int hit_count = alignByPWM(seqList, cluster, false);
	    	if (config.verbose>1)
	    		System.out.println(CommonUtils.timeElapsed(tic)+": PWM "+WeightMatrix.getMaxLetters(cluster.wm)
	    				+" align " + hit_count+" sequences.");
	    	   	
			/** Iteratively build PWM and KSM */
	    	int flag = iteratePWMKSM (cluster, seqList, seed_range, config.use_ksm);
	    	if (flag==0)
	    		return cluster;
	    	else
	    		return null;
		}
    	else
    		return null;
	}
	
	/** 
 * Recursively merge overlapped motif clusters using PWMs<br>
 * Assuming the clusters are sorted by cluster id (0-based)<br>
 * This method consider the motif hit distance, thus avoiding merging dis-similar motifs that co-bind many sequences<br>
 * checked matrix in indexed by the clusterid, thus is ok to sort the clusters without changing clusterid 
 */
@SuppressWarnings({ "unchecked", "rawtypes" })
private void mergeOverlapPwmMotifs (ArrayList<MotifCluster> clusters, ArrayList<Sequence> seqList, boolean[][] checked, 
		boolean useKSM, int recursion){	
	
	if (recursion>10)
		return;
	if (config.verbose>1)
		System.out.println("\n"+CommonUtils.timeElapsed(tic)+": Merging redundant motifs, iteration " + recursion);
	boolean isChanged = false;
	int maxClusterId=0;
	for (MotifCluster c:clusters)
		if (c.clusterId > maxClusterId)
			maxClusterId = c.clusterId;
	
	// hits matrix contains the sequence-motif hit, indexed by clusterid value
	ArrayList[][] hits = new ArrayList[seqs.length][maxClusterId+1];
	for (int j=0;j<clusters.size();j++){
		MotifCluster c = clusters.get(j);
		if (c.wm==null)
			continue;
		WeightMatrixScorer scorer = new WeightMatrixScorer(c.wm);
		for (int i=0;i<seqs.length;i++){
			hits[i][c.clusterId]=CommonUtils.getAllPWMHit(seqs[i], c.wm.length(), scorer, c.pwmThreshold.motif_cutoff);
		}
	}
			
	int seqLen = seqs[0].length();
	ArrayList<MotifCluster> toRemove = new ArrayList<MotifCluster>();
	for (int m=0;m<clusters.size();m++){
		for (int j=m+1;j<clusters.size();j++){
			MotifCluster cluster1 = clusters.get(m);
			MotifCluster cluster2 = clusters.get(j);
			// one of the cluster may be marked as toRemoved, but if the other cluster has just merged, 
			// then all its checked[][] entries will be set as false, including the remove-marked cluster
			// thus need to also check the alignedKmers object reference
			if(checked[cluster1.clusterId][cluster2.clusterId] || cluster1.alignedKmers==null || cluster2.alignedKmers==null)
				continue;
			// multiply k_ratio ( >1 for larger k value) to prefer longer k-mer
			double k_ratio = cluster1.k>cluster2.k ? config.k_ratio : cluster1.k<cluster2.k? 1/config.k_ratio : 1;
			if (config.evaluate_by_ksm){
				if ((cluster1.pwmThreshold.motif_significance+cluster1.ksmThreshold.motif_significance)*k_ratio
						< cluster2.pwmThreshold.motif_significance+cluster2.ksmThreshold.motif_significance){
	//				if (cluster1.pwmThreshold.motif_significance+cluster1.ksmThreshold.motif_significance 
	//						< cluster2.pwmThreshold.motif_significance+cluster2.ksmThreshold.motif_significance){
					cluster1 = clusters.get(j);
					cluster2 = clusters.get(m);
				}
			}
			else{
				if ((cluster1.ksmThreshold.motif_significance)*k_ratio
						< cluster2.ksmThreshold.motif_significance){
					cluster1 = clusters.get(j);
					cluster2 = clusters.get(m);
				}
			}
			if (cluster1.wm==null||cluster2.wm==null)
				continue;
			
			int range = seqLen - cluster1.wm.length()/2 - cluster2.wm.length()/2 + 10;  // add 10 for preventing arrayOutofBound
			int[] same = new int[range*2+1];
			int[] diff = new int[range*2+1];
			for (int i=0;i<seqs.length;i++){
				ArrayList<Integer> hitm = hits[i][cluster1.clusterId];
				ArrayList<Integer> hitj = hits[i][cluster2.clusterId];
				if (hitm.isEmpty()||hitj.isEmpty())
					continue;

				for (int pm:hitm){
					for (int pj:hitj){
						if ((pm>=0&&pj>=0) || (pm<0&&pj<0))
							same[pj-pm+range]++;
						else
							diff[-pj-pm+range]++;			// -pj to get the coord on the same strand as pm
					}
				}
			}
			
			// minOverlap limit the motif hit distance to be short, avoiding merging dis-similar motifs that co-bind many sequences
			int minOverlap = Math.min(cluster1.wm.length(), cluster2.wm.length())/2;
			int minHitCount = Math.min(cluster1.pwmThreshold.posHit, cluster2.pwmThreshold.posHit);
			int maxCount = 0;		// maximum co-hit sequence counts at one distance between two motif hits
			int maxDist = 0;		// the motif distance that gives the max count
			boolean isRC = false;
			for (int p=-minOverlap;p<=minOverlap;p++){
				if (same[p+range]>maxCount){
					maxCount = same[p+range];
					maxDist = p;
				}
			}
			for (int p=-minOverlap;p<=minOverlap;p++){
				if (diff[p+range]>maxCount){
					maxCount = diff[p+range];
					isRC = true;
					maxDist = p;
				}
			}
			
			if (maxCount > minHitCount*config.pwm_hit_overlap_ratio){				// if there is large enough overlap, try to merge 2 clusters
				if (config.verbose>1) {
					System.out.println("\n-------------------------------------------------------------");
		    		System.out.println(String.format("%s: Motif pair have %d overlap hits, motif hit distance = %d%s:\n#%d PWM %s\t(hit=%d, pAUC=%.1f, kAUC=%.1f)\n#%d PWM %s\t(hit=%d, pAUC=%.1f, kAUC=%.1f)\n", 
		    				CommonUtils.timeElapsed(tic), maxCount, maxDist, isRC?"rc":"", 
		    				cluster1.clusterId, WeightMatrix.getMaxLetters(cluster1.wm), cluster1.pwmThreshold.posHit, cluster1.pwmThreshold.motif_significance, cluster1.ksmThreshold.motif_significance,
		    				cluster2.clusterId, WeightMatrix.getMaxLetters(cluster2.wm), cluster2.pwmThreshold.posHit, cluster2.pwmThreshold.motif_significance, cluster2.ksmThreshold.motif_significance));
				}
				MotifCluster newCluster = cluster1.clone(true);
				for (Sequence s:seqList)
					s.resetAlignment();
				int count_pwm_aligned = alignByPWM(seqList, newCluster, false);	// align with first PWM
		    	if (config.verbose>1)
		    		System.out.println(CommonUtils.timeElapsed(tic)+": PWM "+WeightMatrix.getMaxLetters(newCluster.wm)
		    				+" align " + count_pwm_aligned+" sequences.");

				// align additional sequences matching second PWM
				WeightMatrix wm = isRC?WeightMatrix.reverseComplement(cluster2.wm):cluster2.wm;
		        WeightMatrixScorer scorer = new WeightMatrixScorer(wm);		
				count_pwm_aligned=0;
				for (Sequence s:seqList){
		    	  String seq = s.getAlignedSeq();			// PWM to scan unaligned sequence, and align if pass threshold
		    	  if (s.pos!=UNALIGNED )
		    		  continue;
		    	      	  
		          WeightMatrixScoreProfile profiler = scorer.execute(seq);
		          double maxSeqScore = Double.NEGATIVE_INFINITY;
		          int maxScoringShift = 0;
		          char maxScoringStrand = '+';
		          for (int p=0;p<profiler.length();p++){
		        	  double score = profiler.getHigherScore(p);
		        	  if (maxSeqScore<score || (maxSeqScore==score && maxScoringStrand=='-')){	// equal score, prefer on '+' strand
		        		  maxSeqScore = score;
		        		  maxScoringShift = p;
		        		  maxScoringStrand = profiler.getHigherScoreStrand(p);
		        	  }
		          }
		          // if a sequence pass the motif score, align with PWM hit
		          if (maxSeqScore >= cluster2.pwmThreshold.motif_cutoff){
					if (maxScoringStrand =='-'){
						maxScoringShift = seqLen-maxScoringShift-wm.length();
						s.RC();
						// i.e.  (seq.length()-1)-maxScoringShift-(wm.length()-1);
					}
					s.pos = cluster1.pos_pwm_seed-(maxScoringShift+cluster2.wm.length()/2-maxDist-cluster1.wm.length()/2);
					count_pwm_aligned ++;
		          }
		          else
		        	  s.pos = UNALIGNED;
		        }	// each sequence
				if (config.verbose>1)
		    		System.out.println(CommonUtils.timeElapsed(tic)+": PWM "+WeightMatrix.getMaxLetters(wm)+" align additional "+count_pwm_aligned+" sequences.");
				
				// if cluster1 contains all cluster2 hits, remove cluster2 and skip the rest
				if (count_pwm_aligned==0){
					cluster2.wm = null;
					toRemove.add(cluster2);
					cluster2.cleanup();
					for (MotifCluster c: clusters){
						checked[c.clusterId][cluster2.clusterId]=true;
						checked[cluster2.clusterId][c.clusterId]=true;
					}
		    		System.out.println(CommonUtils.timeElapsed(tic)+": PWM "+WeightMatrix.getMaxLetters(wm)+" has been removed.");
					continue;
				}
					
		    	// build PWM
				NewPWM newPWM = buildPWM(seqList, newCluster, 0, tic, false);
				if (newPWM!=null){
					newPWM.updateClusterPwmInfo(newCluster);				// update to use this PWM motif for alignment
					newCluster.ksmThreshold.motif_significance = 0; 		// set it to low to allow a new start
					newCluster.pwmThreshold.motif_significance = 0; 		// set it to low to allow a new start
					count_pwm_aligned = alignByPWM(seqList, newCluster, false);
			    	if (config.verbose>1)
			    		System.out.println(CommonUtils.timeElapsed(tic)+": PWM "+WeightMatrix.getMaxLetters(newCluster.wm)
			    				+" align " + count_pwm_aligned+" sequences.");
				}
				
				// Iteratively improve PWM and sequences alignment
				if (config.kg_hit_adjust_type==1){
					this.posHitStrings = newCluster.posHitStrings.clone();
					this.negHitStrings = newCluster.negHitStrings.clone();		
				}
				else if (config.kg_hit_adjust_type==2){
					this.posCoveredWidth = newCluster.posCoveredWidth.clone();
					this.negCoveredWidth = newCluster.negCoveredWidth.clone();
				}
				int flag = iteratePWMKSM (newCluster, seqList, newCluster.k, useKSM);
				
				// merge if the new PWM is more enriched TODO
//				System.out.println(String.format("merPWM=%.1f, merKSM=%.1f -- c1_PWM=%.1f, c1_KSM=%.1f", newCluster.pwmThreshold.motif_significance, newCluster.ksmThreshold.motif_significance,
//						cluster1.pwmThreshold.motif_significance, cluster1.ksmThreshold.motif_significance));
				boolean isNewClusterBetter = false;
				if (config.evaluate_by_ksm){
					isNewClusterBetter = newCluster.ksmThreshold.motif_significance > cluster1.ksmThreshold.motif_significance;
				}
				else{
					isNewClusterBetter = newCluster.pwmThreshold.motif_significance+newCluster.ksmThreshold.motif_significance > cluster1.pwmThreshold.motif_significance+cluster1.ksmThreshold.motif_significance;					
				}
				if (flag==0 && isNewClusterBetter){	
					if (config.kg_hit_adjust_type==1)
						newCluster.setHitStrings(posHitStrings, negHitStrings);
					else if (config.kg_hit_adjust_type==2)
						newCluster.setCoveredWidth(posCoveredWidth, negCoveredWidth);
					clusters.set(m, newCluster);
					isChanged = true;
					for (int d=0;d<checked.length;d++){
						checked[d][cluster1.clusterId]=false;
						checked[cluster1.clusterId][d]=false;
					}
					cluster1 = newCluster;
					newCluster = null;
					if (config.verbose>1)
			    		System.out.println(String.format("\n%s: Motifs #%d and #%d merge to %s, pAUC=%.1f, kAUC=%.1f", 
			    				CommonUtils.timeElapsed(tic), cluster1.clusterId, cluster2.clusterId,
			    				WeightMatrix.getMaxLetters(cluster1.wm), 
			    				cluster1.pwmThreshold.motif_significance, cluster1.ksmThreshold.motif_significance));	
				}
				else{	
					checked[cluster1.clusterId][cluster2.clusterId]=true;
					checked[cluster2.clusterId][cluster1.clusterId]=true;
					if (config.verbose>1)
			    		System.out.println(String.format("%s: Merged PWM is not more enriched, do not merge", CommonUtils.timeElapsed(tic)));
				}
				// no matter merge or not, remove cluster 2
				toRemove.add(cluster2);
				cluster2.cleanup();
				System.out.println(String.format("%s: Remove motif #%d.", 
		    			CommonUtils.timeElapsed(tic), cluster2.clusterId));
				for (MotifCluster c: clusters){
					checked[c.clusterId][cluster2.clusterId]=true;
					checked[cluster2.clusterId][c.clusterId]=true;
				}
			}		
			else{					// if overlap is not big enough
				checked[cluster1.clusterId][cluster2.clusterId] = true;
				checked[cluster2.clusterId][cluster1.clusterId] = true;
			}
		}
	}
	clusters.removeAll(toRemove);
	System.gc();
	
	// if merged, recursively merge, otherwise return
	if (isChanged){		
		// checked matrix is indexed by the clusterid, thus is ok to sort the clusters without changing clusterid 
		sortMotifClusters(clusters, false);			
		mergeOverlapPwmMotifs (clusters, seqList, checked, useKSM, recursion+1);
	}
}

	private void outputClusters(int[] eventCounts){
		// output cluster information, PFM, and PWM
		File f = new File(outName);
		String name = f.getName();
		f=null;
			
		System.out.println("\n------------------------- "+ name +" final motifs ----------------------------\n");

		// output kmer alignments
		StringBuilder alignedKmer_sb = new StringBuilder();
		for (MotifCluster c:clusters){
			// print aligned k-mers of this cluster
			ArrayList<Kmer> alignedKmers = c.alignedKmers;
	    	if (!alignedKmers.isEmpty())
	    		alignedKmer_sb.append("Motif m"+c.clusterId+", n="+alignedKmers.size()+"\n");
	    	// sort for output
//			Collections.sort(alignedKmers, new Comparator<Kmer>(){
//			    public int compare(Kmer o1, Kmer o2) {
//			    	return o1.compareByHGP(o2);
//			    }
//			});	
	    	alignedKmer_sb.append(getKmerClusterAlignmentString(alignedKmers, alignedKmers.size()));
		}
		CommonUtils.writeFile(outName+".KSM_Alignements.txt", alignedKmer_sb.toString());
		alignedKmer_sb = null;

		// output PWM info
		System.out.println();
		StringBuilder pfm_sb = new StringBuilder();		// default TRANSFAC/STAMP format
		StringBuilder pfm_jasper_sb = new StringBuilder();		// JASPAR format
		StringBuilder pfm_meme_sb = new StringBuilder();		// MEME format
		StringBuilder pfm_homer_sb = new StringBuilder();		// HOMER format
		for (MotifCluster c:clusters){
			if (c.wm==null)
				continue;
     		WeightMatrix wm = c.wm;
			if (config.evaluate_by_ksm && c.ksmThreshold!=null){
				if (c.ksmThreshold.posHit>c.total_aligned_seqs)
	    			c.total_aligned_seqs = c.ksmThreshold.posHit;
			}
			else{
				if (c.pwmThreshold.posHit>c.total_aligned_seqs)
	    			c.total_aligned_seqs = c.pwmThreshold.posHit;
			}
    		System.out.println(String.format("%s motif #%d", name, c.clusterId));
    		if (config.use_ksm && c.ksmThreshold!=null){
    			System.out.println(String.format("\nKSM top k-mer: %s, total %d k-mers", c.seedKmer.getKmerStrRC(), c.alignedKmers.size()));    			
    			System.out.println(String.format("KSM threshold: %.2f, \t\thit=%d+/%d-, kAUC=%.1f\n", c.ksmThreshold.motif_cutoff, c.ksmThreshold.posHit, c.ksmThreshold.negHit, c.ksmThreshold.motif_significance));
    		}
			int pos = c.pos_BS_seed-c.pos_pwm_seed;
    		if (pos>=0)
    			System.out.println(CommonUtils.padding(pos, ' ')+"|\n"+ WeightMatrix.printMatrixLetters(wm));
    		else
    			System.out.println(WeightMatrix.printMatrixLetters(wm).trim());
    		System.out.println(String.format("PWM threshold: %.2f/%.2f, \thit=%d+/%d-, pAUC=%.1f", c.pwmThreshold.motif_cutoff, c.wm.getMaxScore(), c.pwmThreshold.posHit, c.pwmThreshold.negHit, c.pwmThreshold.motif_significance));
    		pfm_sb.append(CommonUtils.makeTRANSFAC (c.pfm, c.pwmThreshold.posHit, 
    				String.format("DE %s_m%d %.2f %d k%d_c%d", name, c.clusterId, c.pwmThreshold.motif_cutoff, pos, c.k, c.pwmThreshold.posHit)));
			if (config.outputMEME)
				pfm_meme_sb.append(CommonUtils.makeMEME (c.pfm, c.pwmThreshold.posHit, 
						String.format("%s_m%d_p%d_k%d_c%d", name, c.clusterId, pos, c.k, c.pwmThreshold.posHit)));
			if (config.outputJASPAR)
				pfm_jasper_sb.append(CommonUtils.makeJASPAR (c.pfm, c.pwmThreshold.posHit, 
						String.format("%s_m%d_p%d_k%d_c%d", name, c.clusterId, pos, c.k, c.pwmThreshold.posHit)));
			if (config.outputHOMER)				
				pfm_homer_sb.append(CommonUtils.makeHOMER (c.pfm, c.pwmThreshold.posHit, 
						String.format("%s_m%d_p%d_k%d_c%d", name, c.clusterId, pos, c.k, c.pwmThreshold.posHit)));
			System.out.println("--------------------------------------------------------------\n");

			
			// paint motif logo
			c.wm.setNameVerType(name, "m"+c.clusterId, null);
			CommonUtils.printMotifLogo(c.wm, new File(outName+".m"+c.clusterId+".PWM.png"), 75);
			
			WeightMatrix wm_rc = WeightMatrix.reverseComplement(wm);
			wm_rc.setNameVerType(name, "m"+c.clusterId, "rc");
			CommonUtils.printMotifLogo(wm_rc, new File(outName+".m"+c.clusterId+".PWM_rc.png"), 75);
		}
		
		System.out.println();	// to break from the motif output list
		
		CommonUtils.writeFile(outName+".all.PFM.txt", pfm_sb.toString());
		if (config.outputMEME)
			CommonUtils.writeFile(outName+".all.PFM_MEME.txt", pfm_meme_sb.toString());
		if (config.outputJASPAR)
			CommonUtils.writeFile(outName+".all.PFM_JASPAR.txt", pfm_jasper_sb.toString());
		if (config.outputHOMER)
			CommonUtils.writeFile(outName+".all.PFM_HOMER.txt", pfm_homer_sb.toString());

		// output HTML report
		boolean isGEM = !this.standalone && eventCounts!=null;
		StringBuffer html = new StringBuffer("<style type='text/css'>/* <![CDATA[ */ table.table {border-width: 0 0 1px 1px; border-spacing: 0;border-collapse: collapse;} th,td{border-color: #600;border-style: solid; margin: 0;padding: 4px;border-width: 1px 1px 1px 1px;} table.noborder, th.noborder, td.noborder{border: 0px solid black;}/* ]]> */</style>");
		html.append("<script language='javascript' type='text/javascript'><!--\nfunction popitup(url) {	newwindow=window.open(url,'name','height=75,width=400');	if (window.focus) {newwindow.focus()}	return false;}// --></script>");
		// Name of the analysis
		html.append("\n<center><table class=\"noborder\"><th bgcolor='#A8CFFF' colspan=3 class=\"noborder\"><font size='5'>");
		html.append(isGEM?"GEM/":"").append("KMAC results: ");
		html.append(name).append("</font></th>");
		// Links to result files
		html.append("\n<tr>");
		if (isGEM){
			html.append("<td valign='top' class=\"noborder\">");
			html.append("<b>Binding Event calling</b>:<p>");
			html.append("<a href='"+name+".GEM_events.txt'>Significant Events</a>&nbsp;&nbsp;: "+eventCounts[0]);
			if (eventCounts[1]!=0)
				html.append("<br><a href='"+name+".GEM_insignificant.txt'>Insignificant Events</a>: "+eventCounts[1]);
			else
				html.append("<br>Insignificant Events: 0");
			if (eventCounts[2]!=0)
				html.append("<br><a href='"+name+".GEM_filtered.txt'>Filtered Events</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;: "+eventCounts[2]);
			else
				html.append("<br>Filtered Events: 0");
			html.append("</td>\n<td valign='top' class=\"noborder\">Read distribution<br><img src='"+name.substring(0,name.length()-2)+".all.Read_Distributions.png' width='350'><hr>");
			html.append("</td>\n<td valign='top' class=\"noborder\">");
		}
		else{
			html.append("\n<td colspan=3 class=\"noborder\">");
		}
		html.append("<b>Motif Discovery</b>:<p>");
		html.append(String.format("<p>Total number of sequences: %d+/%d-", posSeqCount, negSeqCount));
		html.append("<br><ul><li>KSM files:");
		for (MotifCluster c:clusters)
			html.append(" <a href='"+name+".m"+c.clusterId+".KSM.txt'>m"+c.clusterId+"</a>");
		html.append("<li><a href='"+name+".KSM_Alignements.txt'>KSM alignment file</a>");
		html.append("<li><a href='"+name+".all.PFM.txt'>All motif PFMs</a></ul></td></tr>");
		html.append("<table></center>\n");
		
		html.append("<table border=1 class=\"table\">");
		html.append("<th>Motif ID</th><th>High scoring k-mers of the KSM motif</th><th>KSM motif</th><th>Aligned bound sequences</th>");
		html.append("<th>PWM and spatial distribution</th>");
		
    	for (int j=0;j<clusters.size();j++){
    		MotifCluster c = clusters.get(j);
    		html.append("\n<tr><td align='center'>m"+c.clusterId+"</td><td>");
    		html.append("<table border=1 class=\"table\"><th>K-mer</th><th>Offset</th><th>Pos Hit</th><th>Neg Hit</th><th>HGP</th>");
	    	int leftmost_km = Integer.MAX_VALUE;
	    	ArrayList<Kmer> outputKmers = new ArrayList<Kmer>();		
    		// clone kmers, needed to set clusterId
    		ArrayList<Kmer> tmp = new ArrayList<Kmer>();
    		ArrayList<Kmer> alignedKmers = c.alignedKmers;		// each kmer in diff cluster has been clone, would not overwrite
    		skipKmer: for (int i=0;i<alignedKmers.size();i++){
	    		Kmer km = alignedKmers.get(i);
	    		if (km instanceof GappedKmer){
	    			km.setMatrix();
	    			for (Kmer kmer: tmp){
	    				kmer.setMatrix();
	    				if (editDistance(km, kmer)<2)
	    					continue skipKmer;	// skip if this kmer is a gapped k-mer and similar to other output kmers
	    			}
	    		}
	    		km.setClusterId(j);
				if (km.getKmerStartOffset()<leftmost_km)
					leftmost_km = km.getKmerStartOffset();
				tmp.add(km.clone());	// clone kmer for output only
				if (tmp.size()>=10)
					break;
			}
	    	outputKmers.addAll(tmp);
	    	
			for (Kmer km: outputKmers){
				html.append("<tr><td>");
				html.append("<b><font size='4' face='Courier New'>");
				String kmString = km.isSeedOrientation()?km.kmerString:km.getKmerRC();
				char[] kmStr = kmString.toCharArray();
				html.append(CommonUtils.padding(-leftmost_km+km.getKmerStartOffset(), '-'));
				for (char b:kmStr){
					switch(b){
					case 'A': html.append("<font color='green'>A</font>");break;
					case 'C': html.append("<font color='blue'>C</font>");break;
					case 'G': html.append("<font color='orange'>G</font>");break;
					case 'T': html.append("<font color='red'>T</font>");break;
					case 'N': html.append("<font color='grey'>N</font>");break;
					}
				}
				html.append("</font></b></td>");
				html.append(String.format("<td>%d</td><td>%d</td><td>%d</td><td>%.1f</td></tr>", 
						km.getKmerStartOffset(), km.getPosHitCount(), km.getNegHitCount(), km.getHgp()));
			}
			html.append("</table>");
			String prefix = name+".m"+c.clusterId;
			html.append("</td><td valign='top'><a href='"+prefix+".KSM.txt'>"+"<img src='"+prefix+".KSM.png"+"' height='300' width='200'></a>");
    		html.append(String.format("<br>%.2f, hit=%d+/%d-, auc=%.1f</td>", 
    				c.ksmThreshold.motif_cutoff, c.ksmThreshold.posHit, c.ksmThreshold.negHit, c.ksmThreshold.motif_significance));
			html.append("<td valign='top' align='center'><img src='"+prefix+".sequenceHits.png"+"' height='300' width='"+c.k*30+"'></td>");
			if (c.wm==null){
				html.append("<td</td></tr>");
				continue;
			}
    		html.append("<td align='center' valign='bottom'><img src='"+prefix+".PWM.png"+"'>");
    		html.append("<br><img src='"+prefix+".PWM_rc.png"+"'><br>");
    		prefix = name+".Spatial_dist.m0_m"+c.clusterId;
    		html.append("<br><a href='"+prefix+".txt'>"+"<img src='"+prefix+".png"+"' width='220'></a>");
    		html.append(String.format("<p>%.2f/%.2f, hit=%d+/%d-, auc=%.1f<br>", 
    				c.pwmThreshold.motif_cutoff, c.wm.getMaxScore(), c.pwmThreshold.posHit, c.pwmThreshold.negHit, c.pwmThreshold.motif_significance));
    		html.append("</td></tr>");
    	}
		html.append("</table>");
		html.append("</td></tr></table>");
		html.append("<p>"+params+"<p>");
		CommonUtils.writeFile(outName+".results.htm", html.toString());
		
//		for (int i=0;i<Math.min(clusters.size(), 20);i++){
//			ArrayList<Kmer> clusterKmers = clusters.get(i).alignedKmers;
//			MotifThreshold t = this.estimateClusterKgsThreshold(clusterKmers);
//			if (t!=null)
//				System.out.println(String.format("%d\t%.2f\t%d\t%d\t%.1f", i, t.score, t.posHit, t.negHit, t.hgp ));
//		}
		
	}
	
//	
//	/** Index k-mers and pos/neg sequences, remove un-enriched k-mers <br>
// * 	record hit sequence id in the k-mers
// * 	keep track of the k-mer positions in the sequences so that we can use one to find the other<br>
// * 	This will set the sequences as unaligned.
// * */
//private static void indexKmerSequences(ArrayList<Kmer> kmers, double[]seq_weights, ArrayList<Sequence> seqList, 
//		ArrayList<Sequence> seqListNeg, double kmer_hgp){
//
//	/* Initialization, setup sequence list, and update kmers */
//	/** basicKmers include the exact kmers and the base-kmers of the gapped kmers */
//	HashSet<Kmer> basicKmers = new HashSet<Kmer>();
//	for (Kmer km: kmers){
//		km.addBaseKmersToSet(basicKmers);
//    }
//	int totalPosCount = seqList.size();
//	int totalNegCount = seqListNeg.size();
//
//	HashMap <Kmer, HashSet<Integer>> kmer2seq = mapBasicKmerSequences(basicKmers, seqList);
//	HashMap <Kmer, HashSet<Integer>> kmer2seqNeg = mapBasicKmerSequences(basicKmers, seqListNeg);
//	
//	// update kmerCount, and hgp()
//	HashSet<Kmer> unenriched = new HashSet<Kmer>();
//	for (Kmer km:basicKmers){
//		if (kmer2seq.containsKey(km))
//			km.setPosHits(kmer2seq.get(km), seq_weights);
//		else
//			km.setPosHits(new HashSet<Integer>(), seq_weights);
//		if (kmer2seqNeg.containsKey(km))
//			km.setNegHits(kmer2seqNeg.get(km));
//		else
//			km.setNegHits(new HashSet<Integer>());
//	}
//	for (Kmer km:kmers){
//		if (km instanceof GappedKmer){	
//			((GappedKmer)km).update(seq_weights);	// the pos hits of the base-kmers has been updated
//			if (km.getPosHitCount()==0)
//				continue;
//			km.setHgp(computeHGP(totalPosCount, totalNegCount, km.getPosHitCount(), km.getNegHitCount()));	// neg hit count is not change
//			if (km.getHgp()>kmer_hgp)
//				unenriched.add(km);
//		}
//		else{		// for exact kmers
//			if (kmer2seq.containsKey(km)){
//				if (km.getPosHitCount()==0)
//					continue;
//				km.setHgp(computeHGP(totalPosCount, totalNegCount, km.getPosHitCount(), km.getNegHitCount()));	
//				if (km.getHgp()>kmer_hgp)
//					unenriched.add(km);
//			}
//			else{
//				unenriched.add(km);
//			}
//		}
//	}
//	kmers.removeAll(unenriched);
//	
//	// update seqList with GK hits
//	HashSet<Kmer> basicKmersUsed = new HashSet<Kmer>();
//	for (Kmer km: kmers){
//		km.addBaseKmersToSet(basicKmersUsed);
//    }
//	basicKmers.removeAll(basicKmersUsed); // now basicKmers contains all the un-used k-mers
//	for (Sequence s:seqList){
//		s.removeAllKmers(basicKmers);
//	}
//}
//	/** Index k-mers for a set of sequences <br>
//	 * 	record hit sequence id in the k-mers<br>
//	 * 	keep track of the k-mer positions in the sequences so that we can use one to find the other<br>
//	 * 	Sequences are modified to index the kmer match positions, and are set as unaligned.<br>
//	 *  @return returns the mapping of kmer--[sequence ids]
//	 * */
//	private static HashMap <Kmer, HashSet<Integer>> mapBasicKmerSequences(HashSet<Kmer> kmers, ArrayList<Sequence> seqList){
//		/** Next, k-mer search and index is done with basicKmers */
//		// build the kmer search tree
//		AhoCorasick oks = new AhoCorasick();
//		for (Kmer km: kmers){
//			oks.add(km.getKmerString().getBytes(), km);
//	    }
//		oks.prepare();
//		
//		// index kmer->seq id, seq->(kmer and position)
//		// it start with the basicKmers, then the base-kmers will be replaced by the WC-kmers
//		HashMap <Kmer, HashSet<Integer>> kmer2seq = new HashMap <Kmer, HashSet<Integer>>();
//		ArrayList<HashMap<Kmer, ArrayList<Integer>>> fPosArray = new ArrayList<HashMap<Kmer, ArrayList<Integer>>>();
//		ArrayList<HashMap<Kmer, ArrayList<Integer>>> rPosArray = new ArrayList<HashMap<Kmer, ArrayList<Integer>>>();
//
//		for (Sequence s:seqList){
//			/** forward strand (as the original orientation of input sequence) matches */
//			HashMap<Kmer, ArrayList<Integer>> fPos = new HashMap<Kmer, ArrayList<Integer>>();		// forward
//			/** reverse strand (as the original orientation of input sequence) matches */
//			HashMap<Kmer, ArrayList<Integer>> rPos = new HashMap<Kmer, ArrayList<Integer>>();		// reverse
//			
//			String seq = s.seq;						// get the sequence from original strand
//			s.clearKmerPosIndex();
//			s.resetAlignment();
//			HashSet<Kmer> results = findMatchedKmers (seq, oks);
//			if (!results.isEmpty()){
//				for (Kmer km: results){		
//					if (!kmer2seq.containsKey(km)){
//						kmer2seq.put(km, new HashSet<Integer>());
//					}
//					kmer2seq.get(km).add(s.id);
//
//					// This is DIFFERENT from alignSequencesByCoOccurence(ArrayList<Kmer>)
//					// forward strand
//					ArrayList<Integer> pos = StringUtils.findAllOccurences(seq, km.getKmerString());
//					ArrayList<Integer> posRC = StringUtils.findAllOccurences(seq, km.getKmerRC());
////						if(pos.isEmpty()&&posRC.isEmpty())
////							System.err.println("pos.isEmpty()&&posRC.isEmpty()");
//					if (!fPos.containsKey(km))
//						fPos.put(km, new ArrayList<Integer>());
//					for (int p:pos){
//						fPos.get(km).add(p);
//					}
//					for (int p:posRC){
//						fPos.get(km).add(p+RC);					// match kmer RC
//					}
//					// reverse strand
//					pos = StringUtils.findAllOccurences(s.rc, km.getKmerString());
//					posRC = StringUtils.findAllOccurences(s.rc, km.getKmerRC());
//					if (!rPos.containsKey(km))
//						rPos.put(km, new ArrayList<Integer>());
//					for (int p:pos){
//						rPos.get(km).add(p);
//					}
//					for (int p:posRC){
//						rPos.get(km).add(p+RC);					// match kmer RC
//					}
//				} // each matched k-mer
//			}
//			fPosArray.add(fPos);
//			rPosArray.add(rPos);
//		}
//		fPosArray.trimToSize();
//		rPosArray.trimToSize();
//		
//		// for each sequence, replace the base-kmer hits with its WC kmers 
//		// sub-kmer to WC-kmer mapping is many to many, need to collect all the positions for each WK, from all base-kmers, 
//		// then add WK-mers and remove base-kmers.
//		HashMap<Kmer, HashSet<Integer>> wkPosMap = new HashMap<Kmer, HashSet<Integer>>();
//		HashSet<Kmer> replacedBaseKmers = new HashSet<Kmer>();
//		for (int i=0;i<seqList.size();i++){
//			Sequence s = seqList.get(i);
//			HashMap<Kmer, ArrayList<Integer>> fPos = fPosArray.get(i);
//			HashMap<Kmer, ArrayList<Integer>> rPos = rPosArray.get(i);
//			for (Kmer km: fPos.keySet()){		
//				if (km.getGappedKmers()!=null){	// only process base-kmers, ignore exact kmers
//					for (GappedKmer wk: km.getGappedKmers()){
//						if (! wkPosMap.containsKey(wk)){
//							wkPosMap.put(wk, new HashSet<Integer>());
//						}
//						Set<Kmer> bks = wk.getBaseKmers();
//						if (!bks.contains(km)){
//							System.err.println("Inconsistent wk-sk");
//							System.err.println(km.toShortString()+"\t"+km.hashCode()+"\n");
//							for (Kmer sk:bks)
//								System.err.println(sk.toShortString()+"\t"+sk.hashCode());
//						}
//						if (wk.getBaseKmerOrientation(km)){		// same orientation, directly copy
//							wkPosMap.get(wk).addAll(fPos.get(km));
//						}
//						else{		// diff orientation, reverse strand for the match kmer position
//							for (int pos: fPos.get(km)){
//								if (pos>RC/2)
//									wkPosMap.get(wk).add(pos-RC);
//								else
//									wkPosMap.get(wk).add(pos+RC);
//							}
//						}
//					}
//					replacedBaseKmers.add(km);
//				}
//			}
////			k=k+0;
//			for (Kmer rk:replacedBaseKmers)
//				fPos.remove(rk);
//			replacedBaseKmers.clear(); 	// clear for each sequence
//			for (Kmer wk: wkPosMap.keySet()){
//				ArrayList<Integer> newPos = new ArrayList<Integer>();
//				newPos.addAll(wkPosMap.get(wk));
//				fPos.put(wk, newPos);
//			}
//			wkPosMap.clear(); 			// clear for each sequence
//			
//			// same for the reverse sequence
//			for (Kmer km: rPos.keySet()){		
//				if (km.getGappedKmers()!=null){	// only process base-kmers, ignore exact kmers
//					for (GappedKmer wk: km.getGappedKmers()){
//						if (! wkPosMap.containsKey(wk)){
//							wkPosMap.put(wk, new HashSet<Integer>());
//						}
//						if (wk.getBaseKmerOrientation(km)){		// same orientation, directly copy
//							wkPosMap.get(wk).addAll(rPos.get(km));
//						}
//						else{		// diff orientation, reverse strand for the match kmer position
//							for (int pos: rPos.get(km)){
//								if (pos>RC/2)
//									wkPosMap.get(wk).add(pos-RC);
//								else
//									wkPosMap.get(wk).add(pos+RC);
//							}
//						}
//					}
//					replacedBaseKmers.add(km);
//				}
//			}
//			for (Kmer rk:replacedBaseKmers)
//				rPos.remove(rk);
//			replacedBaseKmers.clear(); 	// clear for each sequence
//			
//			for (Kmer wk: wkPosMap.keySet()){
//				ArrayList<Integer> newPos = new ArrayList<Integer>();
//				newPos.addAll(wkPosMap.get(wk));
//				rPos.put(wk, newPos);
//			}
//			wkPosMap.clear(); 			// clear for each sequence
//			
//			// copy the fPos and rPos info into arrays to save space
//			for (Kmer km: fPos.keySet()){
//				s.fPos.put(km, posList2PairArray(fPos.get(km)));
//			}				
//			for (Kmer km: rPos.keySet()){
//				s.rPos.put(km, posList2PairArray(rPos.get(km)));
//			}
//		} // each sequence
//		return kmer2seq;
//	}
//	
	private String getKmerClusterAlignmentString (ArrayList<Kmer> alignedKmers, int topNum){
		StringBuilder sb = new StringBuilder();
		Collections.sort(alignedKmers);		    	// sort by positive hit count
		int leftmost_km = Integer.MAX_VALUE;
		topNum = Math.min(alignedKmers.size(), topNum);
		for (int i=0;i<topNum;i++){
			Kmer km = alignedKmers.get(i);
			if (km.getKmerStartOffset()<leftmost_km)
				leftmost_km = km.getKmerStartOffset();
		}
		for (int i=0;i<topNum;i++){
			Kmer km = alignedKmers.get(i);
			sb.append(km.getKmerStartOffset()+"\t"
					+CommonUtils.padding(-leftmost_km+km.getKmerStartOffset(), '.')
					+(km.isSeedOrientation?km.kmerString:km.getKmerRC())+"\t"
					+km.getPosHitCount()+"\t"+km.getNegHitCount()+"\t"
					+String.format("%.1f", km.getHgp())+"\t"+km.getAlignString()+"\n");
		}
		return sb.toString();
	}
	
	private static Pair<int[],int[]> posList2PairArray(ArrayList<Integer> ps){
		ArrayList<Integer> ps_rc = new ArrayList<Integer>();
		for (int p:ps)
			if (p>RC/2)
				ps_rc.add(p);
		ps.removeAll(ps_rc);
		int k[] = new int[ps.size()];
		for (int i=0;i<ps.size();i++)
			k[i]=ps.get(i);
		int krc[] = new int[ps_rc.size()];
		for (int i=0;i<ps_rc.size();i++)
			krc[i]=ps_rc.get(i);
		return new Pair<int[],int[]>(k, krc);
	}
	
//	private void alignSequencesUsingSeedFamily(ArrayList<Sequence> seqList, ArrayList<Kmer> kmers, Kmer seed){
//		ArrayList<Kmer> seedFamily = getSeedKmerFamily(kmers, seed);	// from all kmers
//		
//		/** align sequences using kmer positions */
//    	for (Sequence s : seqList){						// use all sequences
//    		s.reset();
//    		for (Kmer km:seedFamily){
//				int seed_seq = s.getSeq().indexOf(km.getAlignString());
//				if (seed_seq<0){
//					s.RC();
//					seed_seq = s.getSeq().indexOf(km.getAlignString());
//					if (seed_seq<0)
//						continue;
//				}
////				if (km.getKmerString().equals("CCACGCG")||km.getKmerRC().equals("CCACGCG"))
////					km.getK();
//				s.pos = -seed_seq;
//				break;				// seq is aligned, do not try with weaker k-mers
//    		}
//		}
//	}

	/**
	 * Get all the kmers that has 1 mismatch to the kmerStr<br>
	 * and set alignString and its shift for the kmers, the shift is wrt the input kmerStr orientation
	 * @param kmers
	 * @param kmerStr the length should be the same as kmers
	 * @return
	 */
	private ArrayList<Kmer> findSeedFamily(ArrayList<Kmer> kmers, String kmerStr) {
		ArrayList<Kmer> family = new ArrayList<Kmer>();
		if (kmers.isEmpty())
			return family;
		ArrayList<Kmer> kmerListCopy = Kmer.copyKmerList(kmers);
		ArrayList<Kmer> toRemove = new ArrayList<Kmer>();
		int mm = 1;
//			int mm = config.dc;
//			if (k>=8)
//				mm = config.dc+1;
		// progressively allow more mismatch, this will give the more similar kmer priorty for sequence matching
		for (int i=1;i<=mm;i++){
			for (Kmer kmer: kmerListCopy){
		    	if (CommonUtils.strMinDistance(kmerStr, kmer.kmerString)==i){
		    		Pair<Integer,Integer> p = CommonUtils.strMinDistanceAndShift(kmerStr, kmer.kmerString);
		    		kmer.setAlignString(kmer.kmerString);
		    		kmer.setShift(p.cdr());
		    		kmer.setKmerStartOffset(p.cdr());
		    		kmer.setSeedOrientation(true);
		    		family.add(kmer);
			    	toRemove.add(kmer);
		    	}
		    	else if (CommonUtils.strMinDistance(kmerStr, kmer.getKmerRC())==i){	
		    		// do not RC kmer, set RC string as the alignString for alignment, 
		    		// shift is relative to the input kmerStr orientation
		    		Pair<Integer,Integer> p = CommonUtils.strMinDistanceAndShift(kmerStr, kmer.getKmerRC());
		    		kmer.setAlignString(kmer.getKmerRC());
		    		kmer.setShift(p.cdr());
		    		kmer.setKmerStartOffset(p.cdr());
		    		kmer.setSeedOrientation(false);
		    		family.add(kmer);
			    	toRemove.add(kmer);
		    	}
			}
			kmerListCopy.removeAll(toRemove);
			toRemove.clear();
		}

		return family;
	}
	
	/**
	 * alignSequencesUsingKmers will reset the sequences then align them
	 * @param kmers
	 * @param cluster
	 */
			
	private HashMap<Integer, KmerGroup> alignByKSM (ArrayList<Sequence> seqList, ArrayList<Kmer> kmers, MotifCluster cluster){
		
		BitSet bitSeqWithKmer = new BitSet();		// set all the sequences that contain at least one k-mer
		for (Kmer km:kmers)
			bitSeqWithKmer.or(km.posBits);
		
		HashMap<Integer, KmerGroup> kgs = new HashMap<Integer, KmerGroup>();
		KmerGroup bestKG = null;
		for (Sequence s : seqList){
    		s.resetAlignment();	
			if (!bitSeqWithKmer.get(s.id))
				continue;
    		KmerGroup[] matches = findKsmGroupHits(s.seq, s.rc);
    		if (matches==null)
    			continue;
    		else
    			bestKG = matches[0];
			if (bestKG.kg_score>=cluster.ksmThreshold.motif_cutoff){
				if (bestKG.bs>RC/2){		// match on reverse strand
					s.RC();
					s.pos = -(bestKG.bs-RC); // seq_seed = - seed_seq
				}
				else
					s.pos = -bestKG.bs;		// seq_seed = - seed_seq
			}
			kgs.put(s.id, bestKG);
		}
//		if (config.verbose>1)
//			System.out.println(CommonUtils.timeElapsed(tic)+ ": KSM align "+ ksm_hit_count +" sequences.");
		return kgs;
	}
	
	/** align seqList using a PWM<br>
	 * sequences that have a PWM hit are aligned according to the hit positions<br>
	 * if the sequence has been aligned with KSM, only align those that don't pass the relaxed PWM cutoff
	 */
	private int alignByPWM(ArrayList<Sequence> seqList, MotifCluster cluster, boolean relaxAlignedSequences){
			    	
		WeightMatrix wm = cluster.wm;
        WeightMatrixScorer scorer = new WeightMatrixScorer(wm);		
        double motif_cutoff = cluster.pwmThreshold.motif_cutoff;
		int count_pwm_aligned=0;
		int wm_len = wm.length();
		for (Sequence s:seqList){
			String seq = s.getAlignedSeq();	
			int seqLen = seq.length();
			if (seqLen<wm_len)
				continue;
			
//			boolean ksmAlignOK = false;
			if (s.pos!=UNALIGNED && relaxAlignedSequences){
				int idx = cluster.pos_pwm_seed-s.pos;	// pwm_seq = pwm_seed - seq_seed;
				if (idx+wm_len<seqLen && idx>=0){
					String match = seq.substring(idx, idx+wm_len);
					WeightMatrixScoreProfile profiler = scorer.execute(match);
					// The KSM aligned position can be matched by PWM with a relax cutoff, keep it the same
					if (profiler.getHigherScore(0) >= motif_cutoff * 0.5)
						continue;
				}
			}
    	  
			// if the KSM match is not good, use PWM to scan and align using the PWM if pass threshold    	  
			WeightMatrixScoreProfile profiler = scorer.execute(seq);
			double maxSeqScore = Double.NEGATIVE_INFINITY;
			int maxScoringShift = 0;
			char maxScoringStrand = '+';
			for (int j=0;j<profiler.length();j++){
				double score = profiler.getHigherScore(j);
				if (maxSeqScore<score || (maxSeqScore==score && maxScoringStrand=='-')){	// equal score, prefer on '+' strand
					maxSeqScore = score;
					maxScoringShift = j;
					maxScoringStrand = profiler.getHigherScoreStrand(j);
				}
			}
			// if a sequence pass the motif score, align with PWM hit
			if (maxSeqScore >= motif_cutoff){
				if (maxScoringStrand =='-'){
					maxScoringShift = seqLen-maxScoringShift-wm.length(); // (seq.length()-1)-maxScoringShift-(wm.length()-1);
					s.RC();
				}
				s.pos = cluster.pos_pwm_seed-maxScoringShift;
				count_pwm_aligned ++;
			}
//			else{		// not match by PWM
//				if (!ksmAlignOK)			// if not (match by KSM, and by relaxed PWM)
//					s.pos = UNALIGNED;
//			}
		}	// each sequence

		return count_pwm_aligned;
	}
	
//	/** 
//	 * Return KmerGroup matches in the indexed sequence<br>
//	 * This implementation is rely on pre-index_seq_kmers.
//	 * Using the input kmers set as the KSM (passed in KSM with a threshold??)
//	 * @param s
//	 * @param kmers
//	 * @return return KmerGroup object, return null if not match is found
//	 */
//	private KmerGroup[] findIndexedKsmGroupHits(Sequence s, HashSet<Kmer> kmers){
//		
//		// negPositionPadding is added to to seed_seq for indexing the arrayList because seed_seq may be slightly negative.
//
//		System.currentTimeMillis();
//		for (Kmer km: s.fPos.keySet()){
//			if(!kmers.contains(km))
//				continue;
//			// the string contains this k-mer
//			if (km.isSeedOrientation()){
//				for (int km_seq: s.fPos.get(km).car()){ // string match kmer
//					int seed_seq = km_seq-km.getShift();		// seed_seq = km_seq - km_seed 
//					forward.get(seed_seq+negPositionPadding).add(km);
//				}
//				for (int km_seq: s.rPos.get(km).car()){// RC string match kmer, use reverse list
//					int seed_seq = km_seq-km.getShift();		// seed_seq = km_seq - km_seed 
//					reverse.get(seed_seq+negPositionPadding).add(km);
//				}
//			}
//			else{
//				for (int km_seq: s.fPos.get(km).cdr()){ // string match kmer RC (seed orientation)
//					km_seq -= RC;		// remove RC label, set orientation
//					int seed_seq = km_seq-km.getShift();		// seed_seq = km_seq - km_seed 
//					forward.get(seed_seq+negPositionPadding).add(km);
//				}
//				for (int km_seq: s.rPos.get(km).cdr()){  // RC string match kmer RC (seed orientation), use reverse list
//					// no change, use RC_kmer as RC_string
//					int seed_seq = km_seq-km.getShift();		// seed_seq = km_seq - km_seed 
//					reverse.get(seed_seq-RC+negPositionPadding).add(km);
//				}
//			}
//		}
//		// count position hits in both strand
//		int forwardCount=0,reverseCount=0;
//		for (ArrayList<Kmer>kms:forward)
//			if (!kms.isEmpty())
//				forwardCount++;
//		for (ArrayList<Kmer>kms:reverse)
//			if (!kms.isEmpty())
//				reverseCount++;
//		
//		if (forwardCount+reverseCount==0)
//			return null;
////		System.out.println("Index k-mer positions "+CommonUtils.timeElapsed(t));
//		
//		// The most time consuming step is after this:  optimizeKSM() takes 60ms, computeHGP() takes 8ms; other steps less than 1ms.
//		KmerGroup[] matches = new KmerGroup[forwardCount+reverseCount];
//		int idx = 0;
//		for (int p=0;p<forward.size();p++){
//			if (!forward.get(p).isEmpty()){
////				System.out.println("for (int p=0;p<forward.size();p++) "+CommonUtils.timeElapsed(t));
//				ArrayList<Kmer> kms = forward.get(p);
//				if (config.optimize_KG_kmers && kms.size()>1)
//					optimizeKSM(kms);
////				System.out.println("optimizeKSM "+CommonUtils.timeElapsed(t));
//				KmerGroup kg = config.use_weighted_kmer ? new KmerGroup(kms, p-negPositionPadding, seq_weights) : new KmerGroup(kms, p-negPositionPadding);
//				matches[idx]=kg;
////				System.out.println("Made KG "+CommonUtils.timeElapsed(t));
//				kg.setScore(computeSiteSignificanceScore(kg.getGroupHitCount(), kg.getGroupNegHitCount()));
//				idx++;
////				System.out.println("computeHGP "+CommonUtils.timeElapsed(t));
//				kms.clear();
//			}
//		}	
////		System.out.println(CommonUtils.timeElapsed(t));
//		for (int p=0;p<reverse.size();p++){
//			if (!reverse.get(p).isEmpty()){
//				ArrayList<Kmer> kms = reverse.get(p);
//				if (config.optimize_KG_kmers && kms.size()>1)
//					optimizeKSM(kms);
//				KmerGroup kg = config.use_weighted_kmer ? new KmerGroup(kms, p-negPositionPadding+RC, seq_weights) : new KmerGroup(kms, p-negPositionPadding+RC);
//				matches[idx]=kg;
////				System.out.println("Made KG "+CommonUtils.timeElapsed(t));
//				kg.setScore(-computeHGP(kg.getGroupHitCount(), kg.getGroupNegHitCount()));
//				idx++;
//				kms.clear();
//			}
//		}
////		System.out.println("Done "+CommonUtils.timeElapsed(t));
//
//		return matches;
//	}

//	/** Iteratively build PWM and KSM <br>
//	 *  The seqList should have been aligned so that a new PWM can be built. <br>
//	 *  The seqList should have been indexed with kmers. <br>
//	 *  Need to have use_KSM flag here because of PWM only refinement */
//	private void iteratePWMandKSM (MotifCluster cluster, ArrayList<Sequence> seqList, int seed_range, boolean use_KSM){	
//    	while(true){	
//			NewPWM newPWM = buildPWM(seqList, cluster, 0, tic, true);	
//			if (newPWM==null)
//				return;
//			
//			// test if we want to accept the new PWM
//			if (config.evaluate_by_both){
//				if  (newPWM.hgp >= cluster.pwmThreshold.motif_significance)//TODOTODO
//					return;		// previous pwm+ksm is more enriched, stop here
//			}
//			else{
//				if( (!config.evaluate_by_ksm) && (newPWM.hgp >= cluster.pwmThreshold.motif_significance))//TODOTODO
//					return;		// previous pwm is more enriched, stop here
//			}
//			newPWM.updateClusterPwmInfo(cluster);
//			int pwm_hit_count = alignByPWM(seqList, cluster, config.evaluate_by_ksm||config.evaluate_by_both);
//	    	if (config.verbose>1)
//	    		System.out.println(CommonUtils.timeElapsed(tic)+": PWM "+WeightMatrix.getMaxLetters(cluster.wm)+" align "
//	    				+ (config.evaluate_by_ksm?"additional ":"") + pwm_hit_count+" sequences.");
//
//			if (use_KSM){
//				NewKSM newKSM = extractKSM (seqList, seed_range, null);
//				if (newKSM==null ||newKSM.threshold==null)
//					return;
//				
//				if (config.evaluate_by_both){
//					if  (newKSM.threshold.motif_significance < cluster.ksmThreshold.motif_significance)
//						return;		// previous pwm+ksm is more enriched, stop here
//				}
//				else{
//					if (config.evaluate_by_ksm && (newKSM.threshold.motif_significance < cluster.ksmThreshold.motif_significance))
//						return;
//				}
//				cluster.alignedKmers = newKSM.kmers;
//				cluster.ksmThreshold = newKSM.threshold;
//
//				alignByKSM(seqList, cluster.alignedKmers, cluster);
//			}
//    	}  // Iteratively improving the HGP of the motif
//	}
	/** Iteratively build both PWM and KSM for each set of aligned sequences <br>
	 *  The seqList should have been aligned so that a new PWM can be built. <br>
	 *  The seqList should have been indexed with kmers. <br>
	 *  Need to have use_KSM flag here because of PWM only refinement */
	private int iteratePWMKSM (MotifCluster cluster, ArrayList<Sequence> seqList, int seed_range, boolean use_KSM){	
		ArrayList<Kmer> inputBackup = null;
		ArrayList<Kmer> ksmBackup = null;
		Kmer seedBackup = null;
		while(true){
    		// build PWM
			NewPWM newPWM = buildPWM(seqList, cluster, 0, tic, true);
			if (newPWM==null)
				return -1;
			
			// build KSM
			NewKSM newKSM=null;
			if (use_KSM){
				// update kmer Engine with all input k-mers for extracting KSM, do not use base kmer for matching
				inputBackup = cluster.inputKmers;
				ksmBackup = cluster.alignedKmers;
				seedBackup = cluster.seedKmer;
				cluster.inputKmers = Kmer.deepCloneKmerList(cluster.inputKmers, cluster.seedKmer, seq_weights);
				cluster.seedKmer = cluster.inputKmers.get(0);
				initAhoCorasick(cluster.inputKmers, false, false);		// for extractKSM(), set both to false
				newKSM = extractKSM (seqList, seed_range, pseudoCountRatios[cluster.k]);
				if (newKSM==null || newKSM.threshold==null)
					return -1;
			}
			
			// test if we want to accept the new PWM and KSM
//			System.out.println(String.format("newPWM=%.1f, newKSM=%.1f -- oldPWM=%.1f, oldKSM=%.1f", newPWM.motif_significance, newKSM.threshold.motif_significance,
//					cluster.pwmThreshold.motif_significance, cluster.ksmThreshold.motif_significance));
			if  (newPWM.motif_significance+newKSM.threshold.motif_significance <= 
					cluster.pwmThreshold.motif_significance+cluster.ksmThreshold.motif_significance + config.kmac_iteration_delta){
				// restore k-mers
				cluster.inputKmers = inputBackup;
				cluster.alignedKmers = ksmBackup;
				cluster.seedKmer = seedBackup;
				return 0;		// previous pwm+ksm is more enriched, stop here
			}

			newPWM.updateClusterPwmInfo(cluster);
			cluster.alignedKmers = newKSM.kmers;
			cluster.ksmThreshold = newKSM.threshold;
			
			if (use_KSM){
				initAhoCorasick(cluster.alignedKmers);
				int hit_count = alignByKSM(seqList, cluster.alignedKmers, cluster).size();
				if (config.verbose>1)
		    		System.out.println(CommonUtils.timeElapsed(tic)+": KSM "+cluster.seedKmer.kmerString+" align "
		    				 + hit_count+" sequences.");
			}
			int pwm_hit_count = alignByPWM(seqList, cluster, config.pwm_align_new_only);
	    	if (config.verbose>1)
	    		System.out.println(CommonUtils.timeElapsed(tic)+": PWM "+WeightMatrix.getMaxLetters(cluster.wm)+" align "
	    				+ pwm_hit_count+" sequences.");
    	}// Iteratively improve the combined PWM and KSM motif significance
	}
	
	/** Iteratively build both PWM and KSM for each set of aligned sequences <br>
	 *  The seqList should have been aligned so that a new PWM can be built. <br>
	 *  The seqList should have been indexed with kmers. <br>
	 *  Need to have use_KSM flag here because of PWM only refinement */
//	private void iteratePWMKSM2 (MotifCluster cluster, ArrayList<Sequence> seqList, int seed_range, boolean use_KSM){	
//		int hit_count;
//    	while(true){
//    		// build PWM
//			NewPWM newPWM = buildPWM(seqList, cluster, 0, tic, true);
//			if (newPWM==null)
//				return;
//			
//			// build KSM
//			NewKSM newKSM=null;
//			if (use_KSM){
//				newKSM = extractKSM (seqList, seed_range, null);
//				if (newKSM==null ||newKSM.threshold==null)
//					return;
//			}
//			
//			// test if we want to accept the new PWM and KSM
//			if  (newPWM.motif_significance+newKSM.threshold.motif_significance <= cluster.pwmThreshold.motif_significance+cluster.ksmThreshold.motif_significance)
//					return;		// previous pwm+ksm is more enriched, stop here
//
//			newPWM.updateClusterPwmInfo(cluster);
//			cluster.alignedKmers = newKSM.kmers;
//			cluster.ksmThreshold = newKSM.threshold;
//
//			if (use_KSM)
//				hit_count = alignByKSM(seqList, cluster.alignedKmers, cluster);
//			
//    		// build PWM
//			newPWM = buildPWM(seqList, cluster, 0, tic, true);
//			if (newPWM==null)
//				return;
//			
//			// build KSM
//			if (use_KSM){
//				newKSM = extractKSM (seqList, seed_range, null);
//				if (newKSM==null ||newKSM.threshold==null)
//					return;
//			}
//			
//			// test if we want to accept the new PWM and KSM
//			if  (newPWM.motif_significance+newKSM.threshold.motif_significance <= cluster.pwmThreshold.motif_significance+cluster.ksmThreshold.motif_significance)
//					return;		// previous pwm+ksm is more enriched, stop here
//
//			newPWM.updateClusterPwmInfo(cluster);
//			cluster.alignedKmers = newKSM.kmers;
//			cluster.ksmThreshold = newKSM.threshold;
//			
//			hit_count = alignByPWM(seqList, cluster, false);
//	    	if (config.verbose>1)
//	    		System.out.println(CommonUtils.timeElapsed(tic)+": PWM "+WeightMatrix.getMaxLetters(cluster.wm)
//	    				+" align " + hit_count+" sequences.");
//    	}// Iteratively improve the combine PWM and KSM motif significance
//	}
	
	/**
	 * Compute the distance distributions between primary and secondary motifs
	 */
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public void computeMotifDistanceDistribution (String name){		
		System.out.println("\nCompute motif distance distribution ...");
		
		ArrayList[][] hits = new ArrayList[seqs.length][clusters.size()];
		for (int j=0;j<clusters.size();j++){
			MotifCluster c = clusters.get(j);
//			if (c.wm!=null)
//				System.err.println("Cluster "+j+" PWM length="+c.wm.length());
			if (c.wm!=null){
				WeightMatrixScorer scorer = new WeightMatrixScorer(c.wm);
				for (int i=0;i<seqs.length;i++){
					hits[i][j]=CommonUtils.getAllPWMHit(seqs[i], c.wm.length(), scorer, c.pwmThreshold.motif_cutoff);
				}
			}
			else{
				for (int i=0;i<seqs.length;i++){
					hits[i][j]=new ArrayList<Integer>();
				}
			}
				
		}
		int seqLen = seqs[0].length();
		for (int m=0;m<1;m++){
			if (clusters.get(m).wm==null)
				continue;
			for (int j=0;j<clusters.size();j++){
				if (clusters.get(j).wm==null)
					continue;
				int range = seqLen - clusters.get(m).wm.length()/2 - clusters.get(j).wm.length()/2 + 1;	// add 1 to correct for ceiling effect
				int[] same = new int[range*2+1];
				int[] diff = new int[range*2+1];
				for (int i=0;i<seqs.length;i++){
					ArrayList<Integer> hitm = hits[i][m];
					ArrayList<Integer> hitj = hits[i][j];
					if (hitm.isEmpty()||hitj.isEmpty())
						continue;
					if (m==j){		//self comparison
						for (int a=0;a<hitm.size();a++){
							int pm = hitm.get(a);
							for (int b=a;b<hitm.size();b++){
								int pj = hitm.get(b);
								if ((pm>=0&&pj>=0) || (pm<0&&pj<0))
									same[pj-pm+range]++;
								else
									diff[-pj-pm+range]++;			// -pj to get the coord on the same strand as pm
							}
						}
					}
					else{
						for (int pm:hitm){
							for (int pj:hitj){
								if ((pm>=0&&pj>=0) || (pm<0&&pj<0))
									same[pj-pm+range]++;
								else
									diff[-pj-pm+range]++;			// -pj to get the coord on the same strand as pm
							}
						}
					}
				}
				StringBuilder sb = new StringBuilder();
				int x[]=new int[range*2+1];
				for (int i=0;i<same.length;i++){
					x[i]=i-range;
					sb.append(String.format("%d\t%d\t%d\n", x[i], same[i], diff[i]));
				}
				String fileSuffix = name+".Spatial_dist.m"+clusters.get(m).clusterId+"_m"+clusters.get(j).clusterId;
				plotMotifDistanceDistribution(x, same, diff, fileSuffix+".png");
				
				CommonUtils.writeFile(fileSuffix+".txt", sb.toString());
			}
		}
	}
	
	public void plotMotifDistanceDistribution(int[]x, int[]same, int[]diff, String filename){
		File f = new File(filename);
		int w = 660;
		int h = 300;
		int margin= 30;
		int x_frame = 0;
		int y_frame = 0;
		
		System.setProperty("java.awt.headless", "true");
	    BufferedImage im = new BufferedImage(w, h, BufferedImage.TYPE_INT_ARGB);
	    Graphics g = im.getGraphics();
	    Graphics2D g2 = (Graphics2D)g;
	    g2.setRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
	    g2.setColor(Color.white);
	    g2.fillRect(0, 0, w, h);	
	    
	    int total = 0;
	    int max_same = 0;
	    for (int s:same){
	    	if(s>max_same)
	    		max_same=s;
	    	total+=s;
	    }
	    int max_diff = 0;
	    for (int s:diff){
	    	if(s>max_diff)
	    		max_diff=s;
	    	total+=s;
	    }
	    
	    float x_ratio = 2*60.0f/(w-margin*2-x_frame);
	    float y_ratio = (max_same+max_diff)/(float)(h-margin*3-y_frame);
	    int y_zero = Math.round((float)max_same*(h-margin*3-y_frame)/(max_same+max_diff))+margin;
	    int x_y_axis = x_frame+margin;
	    int x_zero = x_frame+(w-x_frame)/2;
	    int y_bottom = h-margin-y_frame;
	    
	    // frame and tick
	    g2.setColor(Color.gray);
	    g2.drawLine(x_y_axis, y_bottom, w-margin, y_bottom);					// x-frame
//	    g2.drawLine(x_y_axis, margin, x_y_axis, y_bottom);						// y-frame    
	    g2.drawLine(x_y_axis, y_zero, w-margin, y_zero);						// x-axis    
	    int font = 20;
	    g.setFont(new Font("Arial",Font.PLAIN,font));
	    for (int p=-3;p<=3;p++){
	    	int x_coor= x_zero+Math.round(p*20/x_ratio);
	    	g2.drawLine(x_coor, y_bottom-4, x_coor, y_bottom);	// tick  
	    	g2.drawString(p==0?" 0":""+p*20, x_zero+Math.round(p*20/x_ratio)-font/2, y_bottom+font);			// tick label
	    }
	    
	    // message
	    g2.setColor(Color.black);
	    if (max_same>max_diff)
	    	g2.drawString("Total:"+total, x_y_axis+font, y_zero-font*3);
	    else
	    	g2.drawString("Total:"+total, x_y_axis+font, y_zero+font*3);
	    
	    // plot the data
	    int diameter = 8;
	    for (int i=0;i<x.length;i++){
	    	int x_coor= x_zero+Math.round(x[i]/x_ratio);
	    	if (same[i]!=0){
		    	g2.setColor(Color.blue);
		    	g2.drawLine(x_coor, y_zero-Math.round(same[i]/y_ratio), x_coor, y_zero);	// same bar  
		    	g2.drawOval(x_coor-diameter/2, y_zero-Math.round(same[i]/y_ratio)-diameter/2, diameter, diameter);
		    	if (same[i]==max_same){
		    		g2.drawString(String.format("%d,%d", x[i], same[i]), x_coor-font, y_zero-Math.round(same[i]/y_ratio)-font/2);
		    	}
	    	}
	    	if (diff[i]!=0){
		    	g2.setColor(Color.red);
		    	g2.drawLine(x_coor, y_zero, x_coor, y_zero+Math.round(diff[i]/y_ratio));	// diff bar  
		    	g2.drawOval(x_coor-diameter/2, y_zero+Math.round(diff[i]/y_ratio)-diameter/2, diameter, diameter);
		    	if (diff[i]==max_diff){
		    		g2.drawString(String.format("%d,%d", x[i], diff[i]), x_coor-font, y_zero+Math.round(diff[i]/y_ratio)+font+5);
		    	}
	    	}
	    }
	    
	    try{
	    	ImageIO.write(im, "png", f);
	    }
	    catch(IOException e){
	    	System.err.println("Error in printing file "+filename);
	    }
	}
	
	private NewPWM buildPWM(ArrayList<Sequence> seqList, MotifCluster cluster, double noiseRatio, long tic, boolean onlyBetter){	    	
		ArrayList<String> alignedSeqs = new ArrayList<String>();
		ArrayList<Double> weights = new ArrayList<Double>();
		int k = cluster.seedKmer.kmerString.length();
		
		for (Sequence seq:seqList){
			if (seq.pos==UNALIGNED)
				continue;
			// align window = k*2+1		
			int start = k/2-k-seq.pos;		// k/2 is mid_seed (middle of seed k-mer), so start = mid_seed -k - seq_seed = mid_seq -k
			int end = start+2*k+1;
			String startPadding = "";
			String endPadding = "";
			if (start<0){
				startPadding = CommonUtils.padding(-start, "N");
				start=0;
			}
			int seqLen = seq.getAlignedSeq().length();
			if (end>seqLen){
				endPadding = CommonUtils.padding(end-seqLen, "N");
				end=seqLen;
			}
			if (start>=end)
				continue;
 			String s = startPadding+seq.getAlignedSeq().substring(start, end)+endPadding;
 			alignedSeqs.add(s);
 			
 			// determine the weight for the sequence
 			double weight = 1.0;
			if (config.use_pwm_binding_strength_weight)
				weight *= seq_weights[seq.id];
			//TODO: now we apply this to all possible motifs (b.c. don't know which is the primary with KMAC2Gap, is it OK?
			if (config.use_pos_weight){			
				int prof_pos = k/2-seq.pos;
				if (prof_pos<0)
					prof_pos = 0;
				else if (prof_pos>seqLen-1)
					prof_pos = seqLen-1;
				weight *=profile[prof_pos];
			}
			weights.add(weight);
    	}
		alignedSeqs.trimToSize();
		
		if (alignedSeqs.size()<seqs.length*config.motif_hit_factor){
			if (config.verbose>1)
	    		System.out.println(String.format("%s: Seq_Count %d is too few, stop building PWM.", CommonUtils.timeElapsed(tic), alignedSeqs.size()));
			return null;
		}
		
		return buildPWMfromAlignedSequences(alignedSeqs, weights, cluster, noiseRatio, onlyBetter);
    }
	
//	public NewPWM buildPWMfromHits(ArrayList<Sequence> seqList, MotifCluster cluster, Iterator<PWMHit> hits){
//		ArrayList<String> alignedSeqs = new ArrayList<String>();
//		ArrayList<Double> weights = new ArrayList<Double>();
//		while(hits.hasNext()){
//			PWMHit hit = hits.next();
//			String seq = seqList.get(hit.seqId).getSeqStrand(true);
//			String s = seq.substring(hit.start, hit.end+1);
//			if (!hit.isForward)
//				s = SequenceUtils.reverseComplement(s);
//			alignedSeqs.add(s);
//			weights.add(hit.weight*hit.responsibility);
//		}
//		// make PWM from aligned sequence segament
//		return buildPWMfromAlignedSequences(alignedSeqs, weights, cluster, 0, true);	
//	}
//	

	private NewPWM buildPWMfromAlignedSequences(ArrayList<String> alignedSeqs, ArrayList<Double> weights, 
		MotifCluster cluster, double noiseRatio, boolean onlyBetter){
		if (alignedSeqs.isEmpty())
			return null;
		int seqLength = alignedSeqs.get(0).length();
		double[][] pfm = new double[seqLength][MAXLETTERVAL];  // positions x letters
    	if (config.verbose>1)
    		System.out.println(String.format("%s: %d seqs to build PWM.", CommonUtils.timeElapsed(tic), alignedSeqs.size()));
		// count base frequencies
    	double meanWeight = 0;
    	for (double w:weights)
    		meanWeight+=w;
    	meanWeight /= weights.size();
		for (int p=0;p<pfm.length;p++){
			for (char base:LETTERS)			// 0 count can cause log(0), set pseudo-count 0.375 to every pos, every base
				pfm[p][base]=0.375*meanWeight; 		//http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2490743/
		} 
    	for (int i=0;i<alignedSeqs.size();i++){
			for (int p=0;p<pfm.length;p++){
    			char base = alignedSeqs.get(i).charAt(p);
    			pfm[p][base] += weights.get(i);
    		}
    	}
    	
		double[][] pwm = pfm.clone();
		for (int p=0;p<pfm.length;p++){
    		double countN = pfm[p]['N'];
    		if (noiseRatio!=0){
	    		double countLetters = 0;
	    		for (int b=0;b<LETTERS.length;b++){
	    			countLetters += pfm[p][LETTERS[b]];	
	    		}
	    		countLetters+=countN;								// total weighted base count
	    		countN += countLetters * noiseRatio;				// add noise to N count
    		}
    		if (countN!=0){
    			for (int b=0;b<LETTERS.length;b++){
        			char base = LETTERS[b];						// add the fraction of 'N' according to bg dist
	    			pfm[p][base] += countN*bg[b];
	    		}   
    		}
			pwm[p]=pfm[p].clone();
		}
    	// make the PWM
    	// normalize, compare to background, and log2
    	double[] ic = new double[pwm.length];						// information content
    	for (int p=0;p<pwm.length;p++){						// for each position
    		double sum=0;
    		for (char base:LETTERS){						// do not count 'N'
    			sum += pwm[p][base];
    		}
    		for (int b=0;b<LETTERS.length;b++){
    			char base = LETTERS[b];
    			double f = pwm[p][base]/sum;						// normalize freq
    			pwm[p][base] = Math.log(f/bg[b])/Math.log(2.0);		//log base 2
    			ic[p] += f*pwm[p][base];
    		}
    	}
    	
    	// // trim low ic ends (simple method)
    	int leftIdx=ic.length-1;
    	for (int p=0;p<ic.length;p++){
    		if (ic[p]>=config.ic_trim){
    			leftIdx = p;
    			break;
    		}
    	}
    	int rightIdx=0;
    	for (int p=ic.length-1;p>=0;p--){
    		if (ic[p]>=config.ic_trim){    			
    			rightIdx=p;
    			break;
    		}
    	}
    	
    	// special treatment for 'N': set it to lowest score
		if (rightIdx-leftIdx+1>3){		// pwm is long enough
	    	for(int p=leftIdx;p<=rightIdx;p++){
	    		double lowest = 2;
	    		for (char base:LETTERS){
	    			if (lowest>pwm[p][base])
	    				lowest = pwm[p][base];
	    		}
	    		pwm[p]['N']=lowest;
	    	}
		}
		else {
			if (config.verbose>1)
				System.out.println("PWM is too short, W="+(rightIdx-leftIdx+1));
			return null;
		}
		
		/* try all pwm length with the most IC-rich columns, find the best PWM */
		int[] left, right;
		int minLength = cluster.k-2;
		if (cluster.k<=7)
			minLength = cluster.k-1;
		int seedLength = cluster.seedKmer.kmerString.length();
		int seed_start = seqLength/2 - seedLength/2;
		int seed_end = seed_start + seedLength-1;

		if (rightIdx-leftIdx+1>minLength){		// length > k-1
			if (config.bestIC_PWM_trim){
				left=new int[rightIdx-leftIdx+1-minLength];
				right=new int[rightIdx-leftIdx+1-minLength];
				for (int i=0;i<left.length;i++){
					int bestLeft = -1;
					double bestSumIC = 0;
					for(int p=leftIdx;p<=rightIdx-minLength-i;p++){
						int end = minLength+i+p;
						if (ic[p]<config.ic_trim || ic[end]<config.ic_trim)			// if the ends have low ic, skip
							continue;
						double sumIC=0; 
						for (int j=p;j<=end;j++)
							sumIC += ic[j];
						if(sumIC<1*(end-p+1))							// average IC >= 1
							continue;
						if (bestSumIC<sumIC){
							bestSumIC=sumIC;
							bestLeft = p;
						}
					}
					left[i]=bestLeft;
					right[i]=bestLeft+minLength+i;
				}
			}
			else{	// centered on the seed_middle position
				ArrayList<Integer> startList = new ArrayList<Integer>();
				ArrayList<Integer> endList = new ArrayList<Integer>();
				int width = Math.min(seed_start-leftIdx, rightIdx-seed_end);
				if (width>=-1){
					for (int i=-1;i<=width;i++){
						startList.add(seed_start-i);endList.add(seed_end+i);
						if (seed_start-i-1<leftIdx || seed_end+i+1>rightIdx)
							continue;
						if (ic[seed_start-i-1]>ic[seed_end+i+1]){
							startList.add(seed_start-i-1);endList.add(seed_end+i);
						}
						else{
							startList.add(seed_start-i);endList.add(seed_end+i+1);
						}
					}
					left=new int[startList.size()];
					right=new int[startList.size()];
					for (int i=0;i<startList.size();i++){
						left[i] = startList.get(i);
						right[i] = endList.get(i);
					}
				}
				else{				// if high IC bases have too little overlap with seed bases
					left=new int[1];
					right=new int[1];
					left[0]=Math.max(seed_start-1, leftIdx);
					right[0]=Math.min(rightIdx, seed_end+1);
				}
			}
		}
		else{				// if it is not very long
			left=new int[1];
			right=new int[1];
			left[0]=leftIdx;
			right[0]=rightIdx;
		}
		
		// limit PWM length to k and k+1
		ArrayList<Integer> newLefts = new ArrayList<Integer>();
		ArrayList<Integer> newRights = new ArrayList<Integer>();
		for (int i=0;i<left.length;i++){
			// limit PWM length to k and k+1
//			if ( (config.k_PWM_trim && (right[i]-left[i]+1==cluster.k ||right[i]-left[i]+1==cluster.k+1)) ||
//					!config.k_PWM_trim){
			// limit PWM length to seed k-mer length
			if ( (config.k_PWM_trim && (right[i]-left[i]+1==cluster.seedKmer.k)) || !config.k_PWM_trim){
				newLefts.add(left[i]);
				newRights.add(right[i]);
			}
		}
		if (newLefts.isEmpty()){	// if not match length requirement, just take the longest PWM
			newLefts.add(left[left.length-1]);
			newRights.add(right[right.length-1]);
		}
    	double best_significance = -0.001;
    	WeightMatrix bestWM = null;
    	int bestLeft=0;
    	int bestRight=0;    
		for (int i=0;i<newLefts.size();i++){
			int l = newLefts.get(i);
			int r = newRights.get(i);
	    	float[][] matrix = new float[r-l+1][MAXLETTERVAL];   
	    	for(int p=l;p<=r;p++){
	    		for (char base:LETTERS){							// ignore 'N' count
	    			matrix[p-l][base]=(float) pwm[p][base];
	    		}
	    	}
	    	WeightMatrix wm = new WeightMatrix(matrix);
   	
	    	// Check the AUROC of new PWM
	    	double motif_significance = evaluatePwmROC(wm, config.fpr).motif_significance;
    		if (config.verbose>1)
    			if (motif_significance==0)			// TODOTODO
    				System.out.println(String.format("%s: PWM %s is not enriched", CommonUtils.timeElapsed(tic), WeightMatrix.getMaxLetters(wm)));
    		if (motif_significance > best_significance){
    			bestWM = wm;
    			best_significance = motif_significance;
    			bestLeft=l;
    			bestRight=r;
    		}
		}
		if (bestWM==null){
			if (config.verbose>1)
				System.out.println(CommonUtils.timeElapsed(tic)+": None of PWM is enriched.");
			return null;
		}
		
    	MotifThreshold estimate = optimizePwmThreshold(bestWM, "", bestWM.getMaxScore()*0.5, bestWM.getMaxScore()*0.7, bestWM.length());
		if (estimate==null) {
			if (config.verbose>1)
				System.out.println(CommonUtils.timeElapsed(tic)+": None of PWM is enriched.");
			return null;
		}
		else{
			if (config.verbose>1)
				System.out.println(String.format("%s: PWM %s\thit %d+/%d- seqs\tpAUC=%.1f\t%.2f/%.2f", CommonUtils.timeElapsed(tic), 
					WeightMatrix.getMaxLetters(bestWM), estimate.posHit, estimate.negHit, best_significance, estimate.motif_cutoff, bestWM.getMaxScore()));
		}

    	// normailze and store pfm
    	float[][] pfm_trim = new float[bestRight-bestLeft+1][MAXLETTERVAL];   
    	for(int p=bestLeft;p<=bestRight;p++){
    		for (char base:LETTERS){
    			pfm_trim[p-bestLeft][base]=(float) pfm[p][base];
    		}
    	}
		for (int p=0;p<pfm_trim.length;p++){
			float sum = 0;
			for (char base : LETTERS)
				sum += pfm_trim[p][base];
			for (char base : LETTERS)
				pfm_trim[p][base] /= sum;
		}

		NewPWM newPWM = new NewPWM();
		newPWM.wm = bestWM;
		newPWM.pwmGoodQuality = (estimate.posHit>seqs.length*config.motif_hit_factor);
		newPWM.threshold = estimate.motif_cutoff;
		newPWM.motif_significance = best_significance;
		newPWM.pwmPosHitCount = estimate.posHit;
		newPWM.pwmNegHitCount = estimate.negHit;		
		newPWM.pfm = pfm_trim;
		newPWM.pos_pwm_seed = bestLeft - seed_start;		// pwm_seed = pwm_start-seed_start

    	return newPWM;
	}
	
	private class NewPWM{
		double threshold = 0;
		/** motif_significance is defined as partial AUC of ROC curve */
    	double motif_significance = -0.1;
    	WeightMatrix wm = null;
    	boolean pwmGoodQuality = false;
    	int pwmPosHitCount = 0;
    	int pwmNegHitCount = 0;
    	float[][] pfm = null;
    	int pos_pwm_seed = 0;   
    	
    	private void updateClusterPwmInfo(MotifCluster cluster){
    		if (wm==null)
    			return;
    		cluster.wm = wm;
        	cluster.pwmGoodQuality = pwmGoodQuality;
        	cluster.pwmThreshold.motif_cutoff = threshold;
        	cluster.pwmThreshold.motif_significance = motif_significance;
        	cluster.pwmThreshold.posHit = pwmPosHitCount;
        	cluster.pwmThreshold.negHit = pwmNegHitCount;
        	cluster.pfm = pfm;
        	cluster.pos_pwm_seed = pos_pwm_seed;
    	}
	}
	

	
//	private void alignSequencesUsingSeedFamily(ArrayList<Sequence> seqList, ArrayList<Kmer> kmers, Kmer seed){
//		ArrayList<Kmer> seedFamily = getSeedKmerFamily(kmers, seed);	// from all kmers
//		
//		/** align sequences using kmer positions */
//    	for (Sequence s : seqList){						// use all sequences
//    		s.reset();
//    		for (Kmer km:seedFamily){
//				int seed_seq = s.getSeq().indexOf(km.getAlignString());
//				if (seed_seq<0){
//					s.RC();
//					seed_seq = s.getSeq().indexOf(km.getAlignString());
//					if (seed_seq<0)
//						continue;
//				}
////				if (km.getKmerString().equals("CCACGCG")||km.getKmerRC().equals("CCACGCG"))
////					km.getK();
//				s.pos = -seed_seq;
//				break;				// seq is aligned, do not try with weaker k-mers
//    		}
//		}
//	}


	/** Create a NewKSM to store kmer set and the threshold */
	private class NewKSM{
		private ArrayList<Kmer> kmers=null;
		private MotifThreshold threshold = null;
		private NewKSM(KMAC kmac, ArrayList<Kmer> kmers, double pseudoCountRatio){
			this.kmers = kmers;
			initAhoCorasick(kmers);
			kmac.initKgHitSequences(seqList, seqListNeg, kmers);
			Pair<double[],double[] > pair = scoreKsmSequences (seqList, seqListNeg, kmers);
			threshold = optimizeThreshold(pair.car(), pair.cdr(), 0, Double.POSITIVE_INFINITY, pseudoCountRatio);
			threshold.motif_significance = evaluateMotifROC(pair.car(), pair.cdr(), config.fpr).motif_significance;
			if (threshold==null)
				return;
			if (config.verbose>1)
				System.out.println(String.format("%s: KSM cutoff %.2f\thit %d+/%d- seqs\tkAUC=%.1f", 
					CommonUtils.timeElapsed(tic), threshold.motif_cutoff, threshold.posHit, threshold.negHit, threshold.motif_significance));
		}
	}
	/**
	 * Get aligned kmers within seed_range using the aligned sequences<br>
	 * Assuming the kEngine has been initialized with appropriate k-mers
	 * @param seqList
	 * @param seed_range
	 * @param excludes
	 * @return
	 */
	private NewKSM extractKSM (ArrayList<Sequence> seqList, int seed_range, double pseudoCountRatio){

		/** kmer2pos: record all the hit positions (in reference to the seed position) of all k-mers in the alignment */
		HashMap<Kmer, ArrayList<Integer>> kmer2pos_seed = new HashMap<Kmer, ArrayList<Integer>>();
		for (Sequence s:seqList){
			if (s.pos != UNALIGNED){		// aligned seqs
				// find all the k-mers in the sequence, and all of their positions
				String seq = s.getAlignedSeq();		// aligned seqs are in seed orientation
				Iterator searcher = treeAhoCorasick.search(seq.getBytes());
				while (searcher.hasNext()) {
					SearchResult sr = (SearchResult) searcher.next();
					for (Object m: sr.getOutputs()){
						for(Kmer km :str2kmers.get(m)){				// these have to be gapped k-mer
//							if (km.getKmerStr().equals("AGCCTNCCTCC"))
//								seed_range=seed_range+0;
							int km_seq = sr.getLastIndex() - km.k + 1;	// minus k to get the k-mer position (i.e. first base of the k-mer)
							// km_seed = km_seq + seq_seed
							int km_seed = km_seq+s.pos;
							if (km_seed<-seed_range || km_seed>seed_range)
								continue;
							if (!kmer2pos_seed.containsKey(km))
								kmer2pos_seed.put(km, new ArrayList<Integer>());
							kmer2pos_seed.get(km).add(km_seed);	
						}
					}
				}
				// the reverse compliment of alignment orientation
				String seq_rc = s.getAlignedSeqRC();
				searcher = treeAhoCorasick.search(seq_rc.getBytes());
				while (searcher.hasNext()) {
					SearchResult sr = (SearchResult) searcher.next();
					for (Object m: sr.getOutputs()){
						for(Kmer km :str2kmers.get(m)){	
							if (km.kmerString.equals("AGCCTNCCTCC"))
								seed_range=seed_range+0;
							int km_seq = seq_rc.length() - sr.getLastIndex() + 1;	// the start position of kmer_rc in seq
							int km_seed = km_seq+s.pos;  // km_seed = km_seq + seq_seed
							if (km_seed<-seed_range || km_seed>seed_range)
								continue;
							if (!kmer2pos_seed.containsKey(km))
								kmer2pos_seed.put(km, new ArrayList<Integer>());
							kmer2pos_seed.get(km).add(km_seed + RC);	// +RC: "found on RC strand"
						}
					}
				}
			}
		}
		// comment out the method using indexed k-mer
//		for (Sequence s:seqList){
//			if (s.pos != UNALIGNED){		// aligned seqs
//				HashMap<Kmer, Pair<int[],int[]>> kmerPos_seq = s.getKmerPos();	// get from the aligned strand orientation
//				for (Kmer km:kmerPos_seq.keySet()){
//					if (!kmer2pos_seed.containsKey(km))
//						kmer2pos_seed.put(km, new ArrayList<Integer>());
//					Pair<int[],int[]> hits_seq = kmerPos_seq.get(km);
//					if (hits_seq==null)
//						continue;
//					for (int km_seq:hits_seq.cdr()){	// if it is kmerRC match on the aligned strand
//						int km_seed = (km_seq-RC)+s.pos;	// km_seed = km_seq + seq_seed
//						if (km_seed>=-seed_range && km_seed<=seed_range)
//							kmer2pos_seed.get(km).add(km_seed+RC);		// RC to label that kmer is in opposite orientation of seed kmer
//					}
//					for (int km_seq:hits_seq.car()){
//						int pos = km_seq+s.pos;
//						if (pos>=-seed_range && pos<=seed_range)
//							kmer2pos_seed.get(km).add(pos);	
//					}
//				}
//			}
//		}

    	/** find k-mers that are consistently aligned, set kmer consensus position */
		ArrayList<Kmer> alignedKmers = new ArrayList<Kmer>();
		for (Kmer km:kmer2pos_seed.keySet()){
//			if (km.getKmerStr().equals("AGCCTNCCTCC"))
//				seed_range=seed_range+0;
			ArrayList<Integer> posKmer = kmer2pos_seed.get(km);		// all km_seed positions of this kmer
			// The kmer hit in the 2*k region should be at least 1/3 of total hit (why??)
			if (posKmer==null || posKmer.size() < km.getPosHitCount()*config.kmer_inRange_fraction){			
				km.setAlignString("Too few hit in range "+posKmer.size());
//				if (config.verbose>1)
//					System.out.println(String.format("%s\t%d\t%d", km.getKmerRCStr(), km.getPosHitCount(), posKmer.size()));
				continue;
			}	
			// find the most frequent kmerPos
			Pair<int[], int[]> sorted = StatUtil.sortByOccurences(posKmer);
			int counts[] = sorted.cdr();
			if (counts.length<1)
				continue;
			int posSorted[] = sorted.car();
			int maxCount = counts[counts.length-1];
			// The number of consistently aligned k-mer should be more than some fraction of the k-mer hits
			// posKmer.size() only count in seqs in the alignment, getPosHitCount() count all seqs
			if (maxCount < Math.min(posKmer.size(),km.getPosHitCount()) * config.kmer_consistent_fraction){
				km.setAlignString("Low consistent hit count "+maxCount);
//				if (config.verbose>1)
//					System.out.println(String.format("%s\t%d\t%d\t%d", km.getKmerRCStr(), km.getPosHitCount(), posKmer.size(), maxCount));
				continue;
			}
			
			// pass the test, now set the k-mer shift, i.e. position relative to the seed kmer
			ArrayList<Integer> maxPos = new ArrayList<Integer>();
			for (int i=counts.length-1;i>=0;i--){
				if (counts[i]==maxCount){
					int p = posSorted[i];
					maxPos.add(p);
				}	
				else		// do not need to count for non-max elements
					break;
			}
			int shift=0;
			if (maxPos.size()>1){		// if tie with 1+ positions, get the one closest to seed kmer
				int min = Integer.MAX_VALUE;
				int minIdx = 0;
				for (int i=0;i<maxPos.size();i++){
					int pos = maxPos.get(i);
					if (pos>RC/2)
						pos -= RC;
					int distance = Math.abs(pos);
					if (distance<min){
						minIdx = i;
						min = distance;
					}
				}
				shift=maxPos.get(minIdx);
			}
			else{
				shift=maxPos.get(0);
			}
			if (shift>RC/2){		// if kmer is in opposite orientation of seed kmer
				km.setSeedOrientation(false);
				shift -= RC;
			}
			else
				km.setSeedOrientation(true);
			assert (!(shift<-seed_range || shift>seed_range));
			km.setShift(shift);			
			km.setAlignString(maxCount+"/"+posKmer.size());
			km.setKmerStartOffset(km.getShift());
			alignedKmers.add(km);
		}	
		
		kmer2pos_seed = null;
//		System.out.println(String.format("%s: Extracted %d k-mers.", CommonUtils.timeElapsed(tic), alignedKmers.size()));

		if (alignedKmers.isEmpty())
			return null;
		Collections.sort(alignedKmers);
		if (config.optimize_kmer_set){
			int tmp = alignedKmers.size();
			optimizeKSM(alignedKmers, pseudoCountRatio);
			if (config.verbose>1)
				System.out.println(String.format("%s: Extract new KSM, optimize %d to %d k-mers.", CommonUtils.timeElapsed(tic), tmp, alignedKmers.size()));
		}
		else{
			if (config.verbose>1)
				System.out.println(String.format("%s: Extract new KSM, %d k-mers.", CommonUtils.timeElapsed(tic), alignedKmers.size()));
		}
		alignedKmers.trimToSize();
		
//		System.out.println(getKmerClusterAlignmentString(alignedKmers, 20));
		
		return new NewKSM(this, alignedKmers, pseudoCountRatio);
	}

	/**
	 * Optimize the HGP of the kmer set by removing non-essential k-mers<br>
	 * optimizeKSM() should not change the state of the individual k-mers, it will only remove k-mers from the kmers list.
	 * @param alignedKmers
	 */
	private void optimizeKSM(ArrayList<Kmer> kmers, double pseudoCountRatio){
		ArrayList<Kmer> kmerCopy = Kmer.copyKmerList(kmers);
		if (kmers.size()<=1)
			return;

		Collections.sort(kmers, new Comparator<Kmer>(){
		    public int compare(Kmer o1, Kmer o2) {
		    	return o1.compareByHGP(o2);
		    }
		});	
		Collections.reverse(kmers);		// reverse sort so that weaker hit are remove before stronger hit
		boolean changed=true;
		
		// mapping from sequence id to kmers
		HashMap<Integer, HashSet<Kmer>> seq2kmers = new HashMap<Integer, HashSet<Kmer>>();
		for (Kmer km: kmers){
			BitSet bitset = km.getPosBits();
			for (int i = bitset.nextSetBit(0); i >= 0; i = bitset.nextSetBit(i+1)) {
				if (!seq2kmers.containsKey(i))
					seq2kmers.put(i, new HashSet<Kmer>());
				seq2kmers.get(i).add(km);
	 		}
		}
		if (isDebugging)
			System.err.println();
		HashMap<Integer, HashSet<Kmer>> seq2kmers_neg = new HashMap<Integer, HashSet<Kmer>>();
		for (Kmer km: kmers){
			BitSet bitset = km.getNegBits();
			if (isDebugging)
				System.err.println(km.toShortString()+" "+bitset.toString());
			for (int i = bitset.nextSetBit(0); i >= 0; i = bitset.nextSetBit(i+1)) {
				if (!seq2kmers_neg.containsKey(i))
					seq2kmers_neg.put(i, new HashSet<Kmer>());
				seq2kmers_neg.get(i).add(km);
	 		}
		}
		
		
		ArrayList<Kmer> kmers_toRemove = new ArrayList<Kmer>();
		while(changed){
			changed = false;

			// If a k-mer is removed, the hit count is changed only when those hit sequences have only this 1 kmer hit.
			// otherwise, this sequence will still consider a hit because of the other k-mers.
			// Thus, the reduction in hit count is only those single-kmer sequence hits.

			kmers.removeAll(kmers_toRemove);
			for (Kmer km: kmers_toRemove){
				BitSet bitset = km.getPosBits();
				for (int i = bitset.nextSetBit(0); i >= 0; i = bitset.nextSetBit(i+1)) {
					if (seq2kmers.get(i)!=null)
						seq2kmers.get(i).remove(km);
		 		}
				bitset = km.getNegBits();
				for (int i = bitset.nextSetBit(0); i >= 0; i = bitset.nextSetBit(i+1)) {
					if (seq2kmers_neg.get(i)!=null)
						seq2kmers_neg.get(i).remove(km);
		 		}
			}
			kmers_toRemove.clear();
			if (isDebugging)
				System.err.println("kmers.size() "+kmers.size());
			float posHitCount = 0;	// weighted count
			for (int id: seq2kmers.keySet())
				if (config.use_weighted_kmer)
					posHitCount += seq_weights[id];
				else
					posHitCount ++;
			int negHitCount = seq2kmers_neg.size();
			double score_all = computeMotifSignificanceScore(Math.round(posHitCount), negHitCount, pseudoCountRatio);

			for (int i=0;i<kmers.size();i++){
				
				Kmer km = kmers.get(i);
				float count_with_single_kmer = 0;			// weighted count
				int count_with_single_kmer_neg = 0;
				BitSet hits_to_remove = new BitSet();
				BitSet hits_to_remove_neg = new BitSet();
				BitSet hits = km.getPosBits();
				BitSet hits_neg = km.getNegBits();
				// for all the sequences that has a hit for this k-mer km
				for (int id = hits.nextSetBit(0); id >= 0; id = hits.nextSetBit(id+1)) {
					if (seq2kmers.get(id).size()==1){
						if (config.use_weighted_kmer)
							count_with_single_kmer += seq_weights[id];
						else
							count_with_single_kmer++;
						hits_to_remove.set(id);
					}
				}
				for (int id = hits_neg.nextSetBit(0); id >= 0; id = hits_neg.nextSetBit(id+1)) {
					if (seq2kmers_neg.get(id).size()==1){
						count_with_single_kmer_neg++;
						hits_to_remove_neg.set(id);
					}
				}
				// the score can only be improved if we remove the k-mer with stronger effect on negative sequences.
				if (count_with_single_kmer_neg>0){
					int pos_new = Math.round(posHitCount - count_with_single_kmer);
					int neg_new = negHitCount - count_with_single_kmer_neg;
					double score_remove_this_km = computeMotifSignificanceScore( pos_new<=0?0:pos_new, neg_new<=0?0:neg_new, pseudoCountRatio);
					// test whether removing this k-mer will improve enrichment significance hgp
					if (score_remove_this_km>score_all){
						// remove this k-mer 
						if (isDebugging)
							System.err.println(String.format("%s: p%.1f n%d p-%.1f n-%d score: %.1f-->%.1f", 
							km.toShortString(), posHitCount, negHitCount, count_with_single_kmer, count_with_single_kmer_neg,
							score_all, score_remove_this_km));

						changed = true;
						kmers_toRemove.add(km);
						// remove the sequence hit containing only this kmer
						for (int id = hits_to_remove.nextSetBit(0); id >= 0; id = hits_to_remove.nextSetBit(id+1)) 
							seq2kmers.remove(id);
						for (int id = hits_to_remove_neg.nextSetBit(0); id >= 0; id = hits_to_remove_neg.nextSetBit(id+1)) 
							seq2kmers_neg.remove(id);		
						if (isDebugging){
							System.err.println(hits_to_remove_neg+" "+hits_to_remove_neg.toString());
							System.err.println(km.toShortString()+" "+seq2kmers_neg.keySet().toString());
						}
						break;	// break the for loop, continue with while loop because changed = true;
					}
				}
			}
		}
		kmers.trimToSize();
		if (kmers.isEmpty()){
			for (Kmer km:kmerCopy)
				System.out.println(km.toShortString());
		}
	}
//
//	/**
//	 * Optimize the kmer set by removing non-essential k-mers to improve the overall HGP
//	 * @param alignedKmers
//	 */
//	private void optimizeKSM_HGP(ArrayList<Kmer> kmers){
//		ArrayList<Kmer> kmerCopy = Kmer.copyKmerList(kmers);
//		if (kmers.size()<=1)
//			return;
//
//		Collections.sort(kmers, new Comparator<Kmer>(){
//		    public int compare(Kmer o1, Kmer o2) {
//		    	return o1.compareByHGP(o2);
//		    }
//		});	
//		Collections.reverse(kmers);		// reverse so that weaker hit are remove before stronger hit
//		boolean changed=true;
//		
//		// mapping from sequence id to kmers
//		HashMap<Integer, HashSet<Kmer>> seq2kmers = new HashMap<Integer, HashSet<Kmer>>();
//		for (Kmer km: kmers){
//			BitSet bitset = km.getPosBits();
//			for (int i = bitset.nextSetBit(0); i >= 0; i = bitset.nextSetBit(i+1)) {
//				if (!seq2kmers.containsKey(i))
//					seq2kmers.put(i, new HashSet<Kmer>());
//				seq2kmers.get(i).add(km);
//	 		}
//		}
//		HashMap<Integer, HashSet<Kmer>> seq2kmers_neg = new HashMap<Integer, HashSet<Kmer>>();
//		for (Kmer km: kmers){
//			BitSet bitset = km.getNegBits();
//			for (int i = bitset.nextSetBit(0); i >= 0; i = bitset.nextSetBit(i+1)) {
//				if (!seq2kmers_neg.containsKey(i))
//					seq2kmers_neg.put(i, new HashSet<Kmer>());
//				seq2kmers_neg.get(i).add(km);
//	 		}
//		}
//		
//		float posHitCount = 0;	// weighted count
//		for (int id: seq2kmers.keySet())
//			if (config.use_weighted_kmer)
//				posHitCount += seq_weights[id];
//			else
//				posHitCount ++;
//		int negHitCount = seq2kmers_neg.size();
//		double hgp_all = computeHGP(Math.round(posHitCount), negHitCount);
//		
//		while(changed){
//			changed = false;
//			ArrayList<Kmer> kmers_toRemove = new ArrayList<Kmer>();
//
//			for (int i=0;i<kmers.size();i++){
//				
//				Kmer km = kmers.get(i);
//				float count_with_single_kmer = 0;			// weighted count
//				int count_with_single_kmer_neg = 0;
//				HashSet<Integer> hits_to_remove = new HashSet<Integer>();
//				HashSet<Integer> hits_to_remove_neg = new HashSet<Integer>();
//				HashSet<Integer> hits = km.getPosHits();
//				HashSet<Integer> hits_neg = km.getNegHits();
//				for (int id:hits){
//					if (seq2kmers.get(id).size()==1){
//						if (config.use_weighted_kmer)
//							count_with_single_kmer += seq_weights[id];
//						else
//							count_with_single_kmer++;
//						hits_to_remove.add(id);
//					}
//				}
//				for (int id:hits_neg){
//					if (seq2kmers_neg.get(id).size()==1){
//						count_with_single_kmer_neg++;
//						hits_to_remove_neg.add(id);
//					}
//				}
//				if (count_with_single_kmer_neg>0){
//					double hgp_remove_this = computeHGP( Math.round(posHitCount-count_with_single_kmer), negHitCount-count_with_single_kmer_neg);
//					// test whether removing this k-mer will improve enrichment significance hgp
//					if (hgp_remove_this<hgp_all){
////						System.err.println(String.format("%s: p%d n%d p-%d n-%d hgp: %.1f-->%.1f", 
////								km.toString(), posHitCount, negHitCount, count_with_single_kmer, count_with_single_kmer_neg,
////								hgp_all, hgp_remove_this));
//						hgp_all = hgp_remove_this;
//						// remove this k-mer 
//						changed = true;
//						kmers_toRemove.add(km);
//						// remove the sequence hit containing only this kmer
//						for (int h: hits_to_remove)
//							seq2kmers.remove(h);
//						for (int h: hits_to_remove_neg)
//							seq2kmers_neg.remove(h);						
//					}
//				}
//			}
//			kmers.removeAll(kmers_toRemove);
//		}
//		kmers.trimToSize();
//		if (kmers.isEmpty()){
//			for (Kmer km:kmerCopy)
//				System.out.println(km.toShortString());
//		}
//	}
//
//    
    class PWMHit implements Comparable<PWMHit>{
    	int clusterId;
    	WeightMatrix wm;
    	int seqId;
    	boolean isForward;
    	int start;
    	int end;
    	double score;
    	/** emission probability (prob of hit seq given the pwm) */
//    	double eProb;		
    	/** responsibility : the fraction of this data (hit group) is explained by this PWM hit*/
    	double responsibility=1;
    	String str;
    	double weight=1;	// the weight of hit when computing PWM, Weight = EventStrength * PositionalLogitProb. 
    	
		public int compareTo(PWMHit h) {					// descending score
			if(score<h.score){return(1);}
			else if(score>h.score){return(-1);}
			else return(0);
		}
		public int compareByPosition(PWMHit h) {					// ascending position
			if(start<h.start){return(-1);}
			else if(start>h.start){return(1);}
			else {
				if(end<h.end){return(-1);}
				else if(end>h.end){return(1);}
				else return(0);
			}
		}
		boolean overlaps(PWMHit h){
			int minHalf = Math.min(end-start+1, h.end-h.start+1)/2;
			if (start+minHalf>h.end || h.start+minHalf>end)
				return false;
			else
				return true;
		}
		boolean overlaps(int start, int end){
			if (this.start>end || start>this.end)
				return false;
			else
				return true;
		}
		public String toString(){
			return String.format("%d:%s==>%s%d:%d-%d, r%.2f", clusterId, WeightMatrix.getMaxLetters(wm), isForward?"":"-", seqId, start, end, responsibility);
		}
    }
    // sort by kmerSortId, then by kgScore, for KSM LOGO
    private class Seq implements Comparable<Seq>{
		int id;				// original input id
		int kmerSortId=99999;
		double kgScore;
		public int compareTo(Seq s){
			int diff = kmerSortId<s.kmerSortId ? -1 : ((kmerSortId==s.kmerSortId) ? (kgScore>s.kgScore?-1:(kgScore==s.kgScore?0:1)) : 1);
			return diff;
		}
		public String toString(){
			return String.format("%d\t%.3f\t%d", kmerSortId, kgScore, id);
		}
    }  
    
    /**
     * 
     * @author Yuchun
     *
     */
    private class Sequence{
		int id;				// original input id
		String seq;			// original input sequence
		String rc;			// reverse compliment of the original input sequence
		/** seq_seed: the position of sequence relative to the seed kmer position, the squence should be RC() to match the seed orientation */
		int pos=UNALIGNED;
		/** is the sequence in the original orientation as the input? */
		private boolean isOriginalOrientation = true;
		
//		/** forward strand (as the original orientation of input sequence) matches, <br>
//		 * Pair<int[],int[]> stores the positions for kmer and kmerRC, respectively*/
//		HashMap<Kmer, Pair<int[],int[]>> fPos = new HashMap<Kmer, Pair<int[],int[]>>();		// forward
//		/** reverse strand (as the original orientation of input sequence) matches */
//		HashMap<Kmer, Pair<int[],int[]>> rPos = new HashMap<Kmer, Pair<int[],int[]>>();		// reverse
		
		private Sequence(String seq, int id){
			this.id = id;
			this.seq = seq;
			this.rc = SequenceUtils.reverseComplement(seq);
		}
		/** reverse the sequence orientation to match the seed orientation */
		private void RC(){
			isOriginalOrientation = !isOriginalOrientation;
		}
		/**
		 * Get sequence string in the aligned orientation
		 */
		private String getAlignedSeq(){
			return isOriginalOrientation?seq:rc;
		}
		/**
		 * Get sequence string in the RC of aligned orientation
		 */
		private String getAlignedSeqRC(){
			return isOriginalOrientation?rc:seq;
		}
		private void resetAlignment(){
			pos = UNALIGNED;
			isOriginalOrientation = true;
		}
		
		public String toString(){
			return String.format("%d\t%s\t%d\t%s", id, isOriginalOrientation?"F":"R", pos, isOriginalOrientation?seq:rc);
		}
	}

	
	public class MotifCluster implements Comparable<MotifCluster>{
		public int pos_pwm_seed;
		public int pos_BS_seed;
		public MotifThreshold ksmThreshold = new MotifThreshold();
		public MotifThreshold pwmThreshold = new MotifThreshold();
		public WeightMatrix wm;
		
		int clusterId;
		int k;
		Kmer seedKmer;
		float[][] pfm;
		boolean pwmGoodQuality = false;
		ArrayList<Kmer> alignedKmers;			// The K-mer set motif, a set of aligned k-mers
		ArrayList<Kmer> inputKmers;				// The whole set of input K-mers
		int total_aligned_seqs;
		double pi;
		/** the info about sequence hit, total k-mer coverage */
		int[] posCoveredWidth;
		int[] negCoveredWidth;
		String[][] posHitStrings;		// the covered sequence of the best hit of current motif
		String[][] negHitStrings;
		/** save a copy in the cluster, clone()*/
		public void setCoveredWidth(int[] posCoveredWidth, int[] negCoveredWidth){
			this.posCoveredWidth = posCoveredWidth.clone();
			this.negCoveredWidth = negCoveredWidth.clone();
		}
		/** save a copy in the cluster, clone()*/
		public void setHitStrings(String[][] posHitStrings, String[][] negHitStrings){
			this.posHitStrings = posHitStrings.clone();
			this.negHitStrings = negHitStrings.clone();
		}		
		MotifCluster(){}	// empty constructor;

		protected MotifCluster clone(boolean cloneKmers){
			MotifCluster cluster = new MotifCluster();
			cluster.clusterId = clusterId;
			if (pfm!=null){
				cluster.pfm = pfm.clone();
				cluster.wm = wm;
			}
			if (ksmThreshold!=null)
				cluster.ksmThreshold = (MotifThreshold)ksmThreshold.clone();
			if (pwmThreshold!=null)
				cluster.pwmThreshold = (MotifThreshold)pwmThreshold.clone();
			cluster.pwmGoodQuality = pwmGoodQuality;
			cluster.pos_pwm_seed = pos_pwm_seed;
			cluster.pos_BS_seed = pos_BS_seed;
			if (cloneKmers){
				cluster.inputKmers = (ArrayList<Kmer>)Kmer.deepCloneKmerList(inputKmers, seedKmer, seq_weights);
				cluster.seedKmer = cluster.inputKmers.get(0);
				cluster.alignedKmers = alignedKmers;
			}
			else{
				cluster.seedKmer = seedKmer;
				cluster.inputKmers = inputKmers;		// this is read only
				cluster.alignedKmers = alignedKmers;
			}
			cluster.k = k;
			if (posHitStrings!=null)
				cluster.setHitStrings(posHitStrings, negHitStrings);
			if (posCoveredWidth!=null)
				cluster.setCoveredWidth(posCoveredWidth, negCoveredWidth);
			
			return cluster;
		}
		
		private void cleanup(){
			pfm=null;
			inputKmers=null;
			alignedKmers=null;
			posHitStrings=null;
			negHitStrings=null;
			posCoveredWidth=null;
			negCoveredWidth=null;
		}
		
		public ArrayList<Kmer> getAlignedKmers(){
			return alignedKmers;
		}
		/**
		 * Sort by PWM HGP, then by KSM HGP
		 */
		public int compareTo(MotifCluster c) {					// descending pwmThreshold.motif_significance
			if(pwmThreshold.motif_significance > c.pwmThreshold.motif_significance){return(-1);}
			else if(pwmThreshold.motif_significance < c.pwmThreshold.motif_significance){return(+1);}
			else // if same pwmHGP
				if(ksmThreshold.motif_significance > c.ksmThreshold.motif_significance){return(-1);}
				else if(ksmThreshold.motif_significance < c.ksmThreshold.motif_significance){return(+1);}
				else return(0);
		}
		
		/**
		 * Sort by KSM Significance, then by PWM Significance
		 */
		public int compareToByKsmSignificance(MotifCluster c) {					// descending ksmThresholdHGP
			if(ksmThreshold.motif_significance > c.ksmThreshold.motif_significance){return(-1);}
			else if(ksmThreshold.motif_significance < c.ksmThreshold.motif_significance){return(+1);}
			else 	// if same ksm Significance
				if(pwmThreshold.motif_significance > c.pwmThreshold.motif_significance){return(-1);}
				else if(pwmThreshold.motif_significance < c.pwmThreshold.motif_significance){return(+1);}
				else return(0);
		}
		
		/**
		 * Sort by PWM+KSM significance, then by KSM significance
		 */
		public int compareToByKsmPwmSignificance(MotifCluster c) {					// descending pwmThreshold.motif_significance
			if(pwmThreshold.motif_significance+ksmThreshold.motif_significance > c.pwmThreshold.motif_significance+c.ksmThreshold.motif_significance)
				return(-1);
			else
				if(pwmThreshold.motif_significance+ksmThreshold.motif_significance < c.pwmThreshold.motif_significance+c.ksmThreshold.motif_significance)
					return(+1);
				else // if same 
				if(ksmThreshold.motif_significance > c.ksmThreshold.motif_significance){return(-1);}
				else if(ksmThreshold.motif_significance < c.ksmThreshold.motif_significance){return(+1);}
				else return(0);
		}

		/**
		 * This method is only used for selecting K when the hgp of the clusters are very close (<5%)<br>
		 * It first select for k that gives pwm width most close to k, then select for pwm that have a highest pwm score per position
		 * @param c
		 * @return
		 */
		public int compareForSelectingK(MotifCluster c) {					
			int diff = Math.abs(wm.length()+1-seedKmer.getK());				// expect k to be 1 base longer than pwm length	
			double pwmAvg = wm.getMaxScore()/wm.length();
			int diffC = Math.abs(c.wm.length()+1-c.seedKmer.getK());	
			double pwmAvgC = c.wm.getMaxScore()/c.wm.length();
			if(diff==diffC){
				if(pwmAvg<pwmAvgC){return(1);}
				else if(pwmAvg>pwmAvgC){return(-1);}
				else return(0);
			}
			else if (diff<diffC)
				return (-1);
			else
				return (1);
		}
		
		public String toString(){
			return String.format("Motif %d: %s, pAUC=%.1f\n%s, kAUC= %.1f\n", clusterId, wm!=null?WeightMatrix.getMaxLetters(wm):"----", pwmThreshold.motif_significance,
					seedKmer.getKmerStrRC(), ksmThreshold.motif_significance);
		}
	}
	/**
	 * Sort motif clusters and set cluster id
	 * @param motifs
	 */
	private void sortMotifClusters(ArrayList<MotifCluster> motifs, boolean resetClusterId){
		// sort clusters, set clusterid
		if (config.evaluate_by_ksm){
			Collections.sort(motifs, new Comparator<MotifCluster>() {
	            public int compare(MotifCluster o1, MotifCluster o2) {
	                return o1.compareToByKsmSignificance(o2);
	            }
	        });
		}
		else{
			Collections.sort(motifs, new Comparator<MotifCluster>() {
	            public int compare(MotifCluster o1, MotifCluster o2) {
	                return o1.compareToByKsmPwmSignificance(o2);
	            }
	        });
		}
		if (resetClusterId){
			for (int j=0;j<motifs.size();j++){
				motifs.get(j).clusterId = j;
			}
		}
	}
	
	private void printMotifClusters(ArrayList<MotifCluster> motifs, StringBuilder sb){
		for (MotifCluster c:motifs){
//			if (config.evaluate_by_ksm){
//				sb.append(String.format("k=%d\tthresh=%.2f\thit=%d\thgp=1e%.1f\tTopKmer= %s\n", c.k, c.ksmThreshold.kg_score, 
//						c.ksmThreshold.posHit, c.ksmThreshold.motif_hgp, c.seedKmer.getKmerString()));
//			}
//			else if (c.wm!=null){
//				sb.append(String.format("k=%d\tthresh=%.2f\thit=%d\thgp=1e%.1f\tW=%d\tPWM= %s.\n", c.k, c.pwmThreshold,
//						c.pwmThreshold.posHit, c.pwmThreshold.motif_significance, c.wm.length(), WeightMatrix.getMaxLetters(c.wm)));
//			}
			sb.append(String.format("#%d\tk=%d\tKSM= %s \t%.2f, %d, kAUC=%.1f", c.clusterId, c.k, c.seedKmer.kmerString, 
					c.ksmThreshold.motif_cutoff, c.ksmThreshold.posHit, c.ksmThreshold.motif_significance));
			if (c.wm!=null)
				sb.append(String.format("\tPWM= %s\t%.2f, %d, pAUC=%.1f\n", 
						WeightMatrix.getMaxLetters(c.wm), c.pwmThreshold.motif_cutoff, c.pwmThreshold.posHit, c.pwmThreshold.motif_significance));
			else
				sb.append("\tNo significant PWM\n");
		}
		sb.append("\n");
	}
	
	/**
	 * Compute motif significance score [ OR or -log10(hgp) ] from total hits <br>More significant, higher score
	 */	
	public double computeMotifSignificanceScore(double posHitCount, double negHitCount, double pseudoCountRatio){
		if (config.use_odds_ratio)
			return StatUtil.odds_ratio(posSeqCount, negSeqCount, posHitCount, negHitCount, posSeqCount*pseudoCountRatio, negSeqCount*pseudoCountRatio);
		else
			return -computeHGP((int)posHitCount, (int)negHitCount);
	}
	
	/**
	 * Compute matched site (KmerGroup) significance score [ OR or -log10(hgp) ] based on config.use_odds_ratio <br>More significant, higher score
	 */	
	public double computeSiteSignificanceScore(double posHitCount, double negHitCount){
		if (config.use_odds_ratio)
			return StatUtil.odds_ratio(posSeqCount, negSeqCount, posHitCount, negHitCount, 10, 10*config.neg_pos_ratio);
		else
			return -computeHGP((int)posHitCount, (int)negHitCount);
	}

//	/**
//	 * Compute motif significance score [ odds ratio ]<br>More significant ratio, higher score
//	 */	
//	public double computeMotifSignificanceScore(int posHitCount, int negHitCount){
//		return StatUtil.odds_ratio(posSeqCount, negSeqCount, posHitCount, negHitCount, 3);
//	}
//	
	/**
	 * Compute hgp (log10) using the positive/negative sequences<br>
	 * More negative hgp, more significant p-value
	 */	
	public double computeHGP(int posHitCount, int negHitCount){
		if (posHitCount==0)
			return 0;
		return computeHGP(posSeqCount, negSeqCount, posHitCount, negHitCount);
	}
	/**
	 * Compute hgp (log10) using the positive/negative sequences<br>
	 * More negative hgp, more significant p-value
	 */
	public static double computeHGP(int posSeq, int negSeq, int posHit, int negHit){
		int allHit = posHit + negHit;
		int allSeq = posSeq + negSeq;
		if (posHit<negHit){		// select smaller x for hyperGeometricCDF_cache(), to reduce # of x sum operations
			double hgcdf = StatUtil.hyperGeometricCDF_cache(posHit, allSeq, allHit, posSeq);
			if (hgcdf>0.99)
				return computeHGP_TINY(posSeq, negSeq, posHit, negHit);
			else
				return Math.log(1-hgcdf);
		}
		else{	// flip the problem, compute cdf of negative count, CDF for negative hit do not include negHit
			double hgcdf=0;
			if (negHit==0)
				hgcdf = StatUtil.hyperGeometricCDF_cache(0, allSeq, allHit+2+1, negSeq);		// add 1 negHit, 2 posHit as pseudo count
			else
				hgcdf = StatUtil.hyperGeometricCDF_cache(negHit-1, allSeq, allHit, negSeq);
			if (hgcdf==0||hgcdf<=Double.MIN_VALUE)
				return computeHGP_TINY(posSeq, negSeq, posHit, negHit);
			else
				return Math.log10(hgcdf);
		}
	}
	/**
	 * Compute hgp using the positive/negative sequences, high precision approximation
	 * Only use for very small p-value (<MIN_VALUE, 2^-1074)
	 */
	public static double computeHGP_TINY(int posSeq, int negSeq, int posHit, int negHit){
		int allHit = posHit + negHit;
		int allSeq = posSeq + negSeq;
		// flip the problem, compute cdf of negative count
		double hgcdf_log10=0;
		if (negHit==0)
			hgcdf_log10 = StatUtil.log10_hyperGeometricCDF_cache_appr(0, allSeq, allHit+2+1, negSeq); // add 1 negHit, 2 posHit as pseudo count
		else
			hgcdf_log10 = StatUtil.log10_hyperGeometricCDF_cache_appr(negHit-1, allSeq, allHit, negSeq);
		return hgcdf_log10;
	}
	
	/**
	 * Optimize the threshold of a PWM (larger than startingScore) using the positive/negative sequences<br>
	 * Approximate grid search to find best HGP, to reduce run time
	 */
	private MotifThreshold optimizePwmThreshold(WeightMatrix wm, String outName, double startingScore, double endingScore, double pseudoCountRatio){
		double[] posSeqScores = new double[posSeqCount];
		double[] negSeqScores = new double[negSeqCount];
		for (int i=0;i<posSeqCount;i++){
			posSeqScores[i]=WeightMatrixScorer.getMaxSeqScore(wm, seqs[i], config.strand_type==1);
		}
		for (int i=0;i<negSeqCount;i++){
			negSeqScores[i]=WeightMatrixScorer.getMaxSeqScore(wm, seqsNegList.get(i), config.strand_type==1);
		}
		MotifThreshold score = optimizeThreshold(posSeqScores, negSeqScores, startingScore, endingScore, pseudoCountRatio);
		return score;
	}

	private MotifThreshold evaluatePwmROC(WeightMatrix wm, double falsePositiveRate){
		double[] posSeqScores = new double[posSeqCount];
		double[] negSeqScores = new double[negSeqCount];
		for (int i=0;i<posSeqCount;i++){
			posSeqScores[i]=WeightMatrixScorer.getMaxSeqScore(wm, seqs[i], config.strand_type==1);
		}
		for (int i=0;i<negSeqCount;i++){
			negSeqScores[i]=WeightMatrixScorer.getMaxSeqScore(wm, seqsNegList.get(i), config.strand_type==1);
		}
		return evaluateMotifROC(posSeqScores, negSeqScores, falsePositiveRate);
	}

	/**
	 * Compute the partial AUROC given the scores for positive and negative data points<br>
	 * Do not determine the optimal cutoff score
	 * @param posSeqScores
	 * @param negSeqScores
	 * @param falsePositiveRate
	 * @return
	 */
	private MotifThreshold evaluateMotifROC(double[] posSeqScores, double[] negSeqScores, double falsePositiveRate){
		ROC roc = new ROC(posSeqScores, negSeqScores);
		MotifThreshold score = new MotifThreshold();
		score.motif_significance = roc.partialAUC(falsePositiveRate)/falsePositiveRate*100;
		return score;
	}
		
	/**
	 * compute motif scores<br>
	 * @param idxs	The indices of elements to compute
	 * @param poshits
	 * @param neghits
	 * @param motifScores
	 * @return
	 */
	private Pair<Double, Integer> findBestScore(ArrayList<Integer> idxs, double[] poshits, double[] neghits, double[] motifScores, double pseudoCountRatio){
		for (int i:idxs){
			if (i==0 && poshits[0]==posSeqCount)	//why?
    			motifScores[0]=0;
			motifScores[i]= computeMotifSignificanceScore(poshits[i], neghits[i], pseudoCountRatio);
		}

		Pair<Double, TreeSet<Integer>> maxScore = StatUtil.findMax(motifScores);
		int maxIdx = maxScore.cdr().last();
		
		return new Pair<Double, Integer>(maxScore.car(),maxIdx);
	}
	
	public class MotifThreshold{
		/** motif score cutoff */
		public double motif_cutoff;
		public int posHit;
		public int negHit;
		/** motif_significance is defined as partial AUC of ROC curve */
		public double motif_significance=0;	
		public MotifThreshold clone(){
			MotifThreshold thresh = new MotifThreshold();
			thresh.motif_cutoff = this.motif_cutoff;
			thresh.motif_significance = this.motif_significance;
			thresh.posHit = this.posHit;
			thresh.negHit = this.negHit;
			return thresh;			
		}
	}
	
	/**
	 * Grid search the best score cutoff given the positive scores and negative scores<br>
	 * This is applicable to both PWM and KSM
	 * @param posSeqScores	motif scanning score of positive sequences, should be in the same order
	 * @param negSeqScores motif scanning score of negative sequences
	 * @return
	 */
	private MotifThreshold optimizeThreshold(double[] posSeqScores, double[] negSeqScores, double startingScore, double endingScore, double pseudoCountRatio){
		int[] posIdx = StatUtil.findSort(posSeqScores);		
		Arrays.sort(negSeqScores);
		
		// find the threshold motif score
		TreeSet<Double> posScoreUnique = new TreeSet<Double>();
		for (double s:posSeqScores)
			posScoreUnique.add(s);
		Double[] posScores_u = new Double[posScoreUnique.size()];
		posScoreUnique.toArray(posScores_u);
		double[] poshits = new double[posScoreUnique.size()];		// #pos_hit given the corresponding score cutoff
		double[] neghits = new double[posScoreUnique.size()];		// #neg_hit given the corresponding score cutoff
		double[] motifScores = new double[posScoreUnique.size()];
		
		int startIdx = Arrays.binarySearch(posScores_u, startingScore);
		if( startIdx < 0 ) { 
			startIdx = -startIdx-1; 	//insert point
		}
		int endIdx = Arrays.binarySearch(posScores_u, endingScore);
		if( endIdx < 0 ) { 
			endIdx = -endIdx-1; 		//insert point
		}
		if (endIdx-1<startIdx)		// if all scores are not in the staring-ending range, just take the nearest of starting score
			endIdx = startIdx+1;
		for (int i=endIdx-1;i>=startIdx;i--){
			double key = posScores_u[i];
			int index = CommonUtils.findKey(posSeqScores, key);
			if (config.use_weighted_kmer){
				double weightedHit = 0;
				for (int s=index; s<posSeqScores.length; s++)
					weightedHit += seq_weights[ posIdx[s] ];
				poshits[i] = weightedHit;
			}
			else
				poshits[i] = posSeqScores.length-index;
			index = CommonUtils.findKey(negSeqScores, key);
			neghits[i] = negSeqScores.length-index;
			if (neghits[i]>negSeqScores.length*config.fpr){		// do not consider cutoff that gives 0.1 false positive rate
				startIdx = i+1;
				break;
			}

		}

		ArrayList<Integer> idxs = new ArrayList<Integer>();			// the score ids to compute HGP
		for (int i=endIdx-1;i>=startIdx;i--){
			if (poshits[i]*1.0 >= neghits[i]*1.5*posSeqScores.length/negSeqScores.length){	// posHit should be at least 2 fold
				idxs.add(i);
//				System.out.println(String.format("%d\t%.0f\t%.0f\t%.2f", i, poshits[i], neghits[i], posScores_u[i]));
			}
		}
		if (idxs.isEmpty())
			idxs.add(endIdx-1);
		
		Pair<Double, Integer> best;
		
		if (idxs.size()>100 && config.use_grid_search){
		
			// coarse search
			int gridStep = (int)Math.ceil(Math.sqrt((double)idxs.size()/2));
			ArrayList<Integer> idxCoarse = new ArrayList<Integer>();	
			for (int i=0;i<idxs.size();i+=gridStep){
				idxCoarse.add(idxs.get(i));
			}
			if (idxCoarse.get(idxCoarse.size()-1)!=idxs.get(idxs.size()-1))
				idxCoarse.add(idxs.get(idxs.size()-1));
			
			best = findBestScore(idxCoarse, poshits, neghits, motifScores, pseudoCountRatio);
			
			// finer resolution search
			int bestIdx = idxs.indexOf(best.cdr());
			int start = Math.max(bestIdx-gridStep+1, 0);
			int end = Math.min(bestIdx+gridStep-1,idxs.size()-1) ;
			ArrayList<Integer> idxFine = new ArrayList<Integer>();	
			for (int i=start;i<=end;i++){
				idxFine.add(idxs.get(i));
			}
			
			best = findBestScore(idxFine, poshits, neghits, motifScores, pseudoCountRatio);
		}
		else
			best = findBestScore(idxs, poshits, neghits, motifScores, pseudoCountRatio);
		
		MotifThreshold score = new MotifThreshold();
		score.motif_cutoff = posScores_u[best.cdr()];
		score.motif_significance = best.car();
		score.posHit = (int) Math.round(poshits[best.cdr()]);
		score.negHit = (int) Math.round(neghits[best.cdr()]);
		
		return score;		
	}
	/**
	 * Compute the KSM scores for pos/neg sequences<br>
	 * The KSM k-mers are assumed to have been loaded into the Engine
	 * @returns the KSM scores for pos/neg sequences.
	 */
	private Pair<double[],double[] > scoreKsmSequences (ArrayList<Sequence> seqList, ArrayList<Sequence> seqListNeg, ArrayList<Kmer> kmers){
//		if (config.verbose>1)
//			System.out.println(CommonUtils.timeElapsed(tic)+ ": scoreKsmSequences, start.");

		BitSet bitSeqWithKmer = new BitSet();
		BitSet bitSeqWithKmerNeg = new BitSet();
		for (Kmer km:kmers){
			bitSeqWithKmer.or(km.posBits);
			bitSeqWithKmerNeg.or(km.negBits);
		}
		
		double[] posSeqScores = new double[seqList.size()];
		double[] negSeqScores = new double[seqListNeg.size()];
		KmerGroup[] kgs=null;
		for (int i=0;i<seqList.size();i++){
			Sequence s = seqList.get(i);
			if (!bitSeqWithKmer.get(s.id))
				continue;
			kgs = findKsmGroupHits(s.seq, s.rc);				// both sequence orientation will be scan
			if (kgs==null)
				posSeqScores[i]=0;
			else
				posSeqScores[i]=kgs[0].getScore();
		}

		for (int i=0;i<seqListNeg.size();i++){
			Sequence s = seqListNeg.get(i);
			if (!bitSeqWithKmerNeg.get(s.id))
				continue;
			kgs = findKsmGroupHits(s.seq, s.rc);				// both sequence orientation will be scan
			if (kgs.length==0)
				negSeqScores[i]=0;
			else
				negSeqScores[i]=kgs[0].getScore();
		}
		return new Pair<double[],double[]>(posSeqScores,negSeqScores);
	}

	/**
	 * Init the KmerGroup hit covered sequence and width for pos/neg sequences<br>
	 * The KSM k-mers are assumed to have been loaded into the Engine
	 * @returns the KSM scores for pos/neg sequences.
	 */
	private void initKgHitSequences (ArrayList<Sequence> seqList, ArrayList<Sequence> seqListNeg, ArrayList<Kmer> kmers){
		BitSet bitSeqWithKmer = new BitSet();
		BitSet bitSeqWithKmerNeg = new BitSet();
		for (Kmer km:kmers){
			bitSeqWithKmer.or(km.posBits);
			bitSeqWithKmerNeg.or(km.negBits);
		}
		
		if (config.kg_hit_adjust_type==2){
			int[] posCoveredWidth = new int[this.posCoveredWidth.length];
			KmerGroup[] kgs=null;
			for (int i=0;i<seqList.size();i++){
				Sequence s = seqList.get(i);
				if (!bitSeqWithKmer.get(s.id)){
					posCoveredWidth[i]=0;
					continue;
				}
				kgs = findKsmGroupHits(s.seq, s.rc);				// both sequence orientations will be scanned
				if (kgs==null)
					posCoveredWidth[i]=0;
				else
					posCoveredWidth[i]=kgs[0].getCoveredWidth();
			}
			this.posCoveredWidth = posCoveredWidth;

			int[] negCoveredWidth = new int[this.negCoveredWidth.length];
			for (int i=0;i<seqListNeg.size();i++){
				Sequence s = seqListNeg.get(i);
				if (!bitSeqWithKmerNeg.get(s.id)){
					negCoveredWidth[i]=0;
					continue;
				}
				kgs = findKsmGroupHits(s.seq, s.rc);				// both sequence orientations will be scanned
				if (kgs.length==0)
					negCoveredWidth[i]=0;
				else
					negCoveredWidth[i]=kgs[0].getCoveredWidth();
			}
			this.negCoveredWidth = negCoveredWidth;
		}
		else if (config.kg_hit_adjust_type==1){
			String[][] posHitStrings = new String[posSeqCount][];
			KmerGroup[] kgs=null;
			for (int i=0;i<seqList.size();i++){
				Sequence s = seqList.get(i);
				if (!bitSeqWithKmer.get(s.id)){
					posHitStrings[i]=null;
					continue;
				}
				kgs = findKsmGroupHits(s.seq, s.rc);				
				if (kgs==null)
					posHitStrings[i]=null;
				else{
					posHitStrings[i]=new String[kgs.length];
					for (int j=0;j<kgs.length;j++)
						posHitStrings[i][j] = kgs[j].getCoveredSequence();
				}
			}

			String[][] negHitStrings = new String[negSeqCount][];
			for (int i=0;i<seqListNeg.size();i++){
				Sequence s = seqListNeg.get(i);
				if (!bitSeqWithKmerNeg.get(s.id)){
					negHitStrings[i]=null;
					continue;
				}
				kgs = findKsmGroupHits(s.seq, s.rc);				
				if (kgs.length==0)
					negHitStrings[i]=null;
				else{
					negHitStrings[i]=new String[kgs.length];
					for (int j=0;j<kgs.length;j++)
						negHitStrings[i][j] = kgs[j].getCoveredSequence();
				}
			}
			// need to update after the whole process, because findKsmGroupHits() uses the HitStrings
			this.posHitStrings = posHitStrings;
			this.negHitStrings = negHitStrings;
		}
	}

	
	/**
	 * A wrapper for updateEngine(kmers, boolean, boolean). Used before KSM scanning, or alignByKSM().<br>
	 * It set areKmersAligned to be true, set use_base_kmers using config setting.
	 * @param kmers
	 */
	public void initAhoCorasick(ArrayList<Kmer> kmers){
		initAhoCorasick(kmers, config.match_base_kmer, true);
	}
	
	/** Prepare the AhoCorasic search Engine, assuming the kmers are unique<br>
	 * Also setup str2kmers and str2kmerOffsets hash table to retrieve K-mer object once k-mer strings are matched
	 * @param kmers List of kmers (with kmerString, sequence hit count)
	 * @param use_base_kmers true: use base kmers for KG scoring; false: use whole gapped k-mer for scoring
	 * @param areKmersAligned true: set string orientation as seed, for alignKSM(); false: set String orientation as kmer, for extractKSM()
	 */
	public void initAhoCorasick(ArrayList<Kmer> kmers, boolean use_base_kmers, boolean areKmersAligned){
		if (kmers.isEmpty()){
			engineInitialized = false;
			return;
		}		
		
		//Init Aho-Corasick (AC) algorithm for searching multiple Kmers in sequences
		//ahocorasick_java-1.1.tar.gz is an implementation of Aho-Corasick automata for Java. BSD license.
		//from <http://hkn.eecs.berkeley.edu/~dyoo/java/index.html> 
		treeAhoCorasick = new AhoCorasick();
		// str2Kmers use a list/array to store kmers because a base k-mer string may represent different gapped k-mers, 
		// with different offsets
		HashMap<String, ArrayList<Kmer>> tmpStr2kmers = new HashMap<String, ArrayList<Kmer>>();	
		HashMap<String, ArrayList<Integer>> tmpStr2kmerEndOffsets = new HashMap<String, ArrayList<Integer>>();	
//		int a=0;
		for (Kmer km: kmers){
//			if (km.getKmerStr().equals("ATGGT"))	// for debugging
//				a = kmers.size();
			if (km instanceof GappedKmer){
				GappedKmer gk = (GappedKmer)km;
				if (use_base_kmers){
					for (Kmer bk: gk.getBaseKmers()){
						String kmerStr;
						Kmer baseKmer = null;
						if (areKmersAligned){
							if ((gk.isSeedOrientation() && gk.getBaseKmerOrientation(bk)) || (!gk.isSeedOrientation() && !gk.getBaseKmerOrientation(bk))){
								baseKmer = bk;
							}
							else{
								baseKmer = bk.clone();		// clone the base k-mer so that it can be modified without affecting other gKmers
								baseKmer.setKmerString(bk.kmerRC);
								baseKmer.setKmerStartOffset(gk.kmerStartOffset);
								baseKmer.setShift(gk.kmerStartOffset);
							}
							kmerStr = baseKmer.kmerString;
						}
						else
							kmerStr = gk.getBaseKmerOrientation(bk)?bk.kmerString:bk.kmerRC;	// for extractKSM, use gkmer orientation
						
						if (!tmpStr2kmers.containsKey(kmerStr)){
							tmpStr2kmers.put(kmerStr, new ArrayList<Kmer>());	
							tmpStr2kmerEndOffsets.put(kmerStr, new ArrayList<Integer>());
						}
						tmpStr2kmers.get(kmerStr).add(baseKmer);			// use base-kmers for scoring
						// base kmer has same length as gkmer, thus same offset
						// Here when initializaing, convert startOffset to endOffset
						// b/c AC search returns KmEnd_seq position, add the kmer length to get seed_kmEnd
						// - kmStart_seed - kmEnd_kmStart (k) --> seed_kmEnd; 
						// When scanning KSM, compute match position seed_seq: seed_kmEnd + kmEnd_seq --> seed_seq
						tmpStr2kmerEndOffsets.get(kmerStr).add(-gk.kmerStartOffset-gk.k);
					}
				}
				else{	// use whole gapped kmer for matching KG
					for (Kmer bk: gk.getBaseKmers()){
						String kmerStr;
						if (areKmersAligned)
							kmerStr = gk.isSeedOrientation()?(gk.getBaseKmerOrientation(bk)?bk.kmerString:bk.kmerRC)
									: (!gk.getBaseKmerOrientation(bk)?bk.kmerString:bk.kmerRC);
						else
							kmerStr = gk.getBaseKmerOrientation(bk)?bk.kmerString:bk.kmerRC;
						if (!tmpStr2kmers.containsKey(kmerStr)){
							tmpStr2kmers.put(kmerStr, new ArrayList<Kmer>());	
							tmpStr2kmerEndOffsets.put(kmerStr, new ArrayList<Integer>());
						}
						tmpStr2kmers.get(kmerStr).add(km);			// use gapped kmers for scoring
						tmpStr2kmerEndOffsets.get(kmerStr).add(-gk.kmerStartOffset-gk.k);	 // base kmer has same length as gkmer, thus same offset
					}
				}
			}
			else{		// if exact kmer
				String kmerStr = null;
				if (areKmersAligned)
					kmerStr = km.isSeedOrientation()?km.kmerString:km.kmerRC;
				else
					kmerStr = km.kmerString;
				if (!tmpStr2kmers.containsKey(kmerStr)){
					tmpStr2kmers.put(kmerStr, new ArrayList<Kmer>());	
					tmpStr2kmerEndOffsets.put(kmerStr, new ArrayList<Integer>());
				}
				tmpStr2kmers.get(kmerStr).add(km);	
				// Convert startOffset to endOffset, b/c AC search returns end+1 position, add the kmer length with kmer offset, minus
				// - start_seed - end_start (k) --> seed_end; seed_end + end_seq --> seed_seq
				tmpStr2kmerEndOffsets.get(kmerStr).add(-km.kmerStartOffset-km.k);
			}			
	    }
		str2kmers.clear();
		str2kmerEndOffsets.clear();
		// convert ArrayList to array, for more efficient access, b/c these data structure are used for every k-mer match
		for (String s: tmpStr2kmers.keySet()){
			ArrayList<Kmer> ks = tmpStr2kmers.get(s);
			ArrayList<Integer> os = tmpStr2kmerEndOffsets.get(s);
			// if use base kmers, the same base k-mer may match multiple gapped k-mer pattern, thus reduce redundancy here
			if (use_base_kmers && ks.size()>1){
				HashSet<Integer> uniq_offsets = new HashSet<Integer>();
				ArrayList<Integer> ids = new ArrayList<Integer>();
				for (int i=0;i<ks.size();i++){
					int o = os.get(i);
					if (uniq_offsets.contains(o))
						ids.add(i);
					else
						uniq_offsets.add(o);
				}
				Collections.sort(ids);
				Collections.reverse(ids);
				for (int id:ids){
					ks.remove(id);
					os.remove(id);
				}
			}
			Kmer[] kms = new Kmer[ks.size()];
			for (int i=0;i<kms.length;i++)
				kms[i] = ks.get(i);
			int[] ofs = new int[os.size()];
			for (int i=0;i<ofs.length;i++)
				ofs[i] = os.get(i);
			str2kmers.put(s, kms);
			str2kmerEndOffsets.put(s, ofs);
		}
		for (String s: str2kmers.keySet())
			treeAhoCorasick.add(s.getBytes(), s);
	    treeAhoCorasick.prepare();
	    engineInitialized = true;
	}		
	
	/** 
	 * Search all k-mers (loaded in the AhoCorasick tree) in the sequence, both strands<br>
	 * Do not return k-mer positions
	 * @param seq sequence string to search k-mers
	 * @return a set of kmers found
	 */
	public static HashSet<Kmer> findMatchedKmers (String seq, AhoCorasick tree){
		HashSet<Object> kmerFound = new HashSet<Object>();	// each kmer is only used 
		//Search for all kmers in the sequences using Aho-Corasick algorithms (initialized)
		//ahocorasick_java-1.1.tar.gz is an implementation of Aho-Corasick automata for Java. BSD license.
		//from <http://hkn.eecs.berkeley.edu/~dyoo/java/index.html> 
	
		Iterator searcher = tree.search(seq.getBytes());
		while (searcher.hasNext()) {
			SearchResult result = (SearchResult) searcher.next();
			kmerFound.addAll(result.getOutputs());
		}
		// the reverse compliment
		String seq_rc = SequenceUtils.reverseComplement(seq);
		searcher = tree.search(seq_rc.getBytes());
		while (searcher.hasNext()) {
			SearchResult result = (SearchResult) searcher.next();
			kmerFound.addAll(result.getOutputs());
		}
		
		HashSet<Kmer> result = new HashSet<Kmer>();
		for (Object km: kmerFound)
			result.add((Kmer)km);
		return result;
	}
	/** 
	 * Search all k-mers in the sequence<br>
	 * This is used for GEM peak calling
	 * @param seq sequence string to search k-mers
	 * @return an array of KmerGroups, with match positions, significance, etc:<br>
	 * Each k-mer group maps to a binding position in the sequence<br>
	 * Note: matches on negative strand are combined with matches on positive strand at the same position
	 */
	public KmerGroup[] findUnstrandedKsmGroupHits (String seq){
		seq = seq.toUpperCase();
		HashMap<Integer, ArrayList<Kmer>> result = new HashMap<Integer, ArrayList<Kmer>> ();
		
		//Search for all kmers in the sequences using Aho-Corasick algorithms (initialized)
		//ahocorasick_java-1.1.tar.gz is an implementation of Aho-Corasick automata for Java. BSD license.
		//from <http://hkn.eecs.berkeley.edu/~dyoo/java/index.html> 
		Iterator searcher = treeAhoCorasick.search(seq.getBytes());
		while (searcher.hasNext()) {
			SearchResult sr = (SearchResult) searcher.next();
			for (Object s: sr.getOutputs()){
				Kmer[] kmers = str2kmers.get(s);	
				// AC search returns end+1 position, end_seq; endOffset is seed_end;
				// thus   seed_end + end_seq --> seed_seq
				int[] kmerOffsets = str2kmerEndOffsets.get(s);	
				for(int i=0;i<kmers.length;i++){	
					int x = sr.getLastIndex() + kmerOffsets[i];	
					if (!result.containsKey(x))
						result.put(x, new ArrayList<Kmer>());
					result.get(x).add(kmers[i]);	
				}
			}
		}
		// the reverse compliment
		String seq_rc = SequenceUtils.reverseComplement(seq);
		searcher = treeAhoCorasick.search(seq_rc.getBytes());
		while (searcher.hasNext()) {
			SearchResult sr = (SearchResult) searcher.next();
			for (Object s: sr.getOutputs()){
				Kmer[] kmers = str2kmers.get(s);	
				int[] kmerOffsets = str2kmerEndOffsets.get(s);	
				for(int i=0;i<kmers.length;i++){	
					int x = sr.getLastIndex() + kmerOffsets[i];		// match position in seq_rc
					x = seq.length()-1-x;		// convert to position in Seq, discarding the strand info
					if (!result.containsKey(x))
						result.put(x, new ArrayList<Kmer>());
					result.get(x).add(kmers[i]);	
				}
			}
		}
		
		KmerGroup[] matches = new KmerGroup[result.keySet().size()];
		int idx = 0;
		for (int p:result.keySet()){
			ArrayList<Kmer> kmers = result.get(p);
			int leftShift = 999;			// left most shift position
			int longest = 0;				// length from the seed position
			for (Kmer km:kmers){
				if (km.getKmerStartOffset() < leftShift)
					leftShift = km.getKmerStartOffset();
				if (km.k+km.getKmerStartOffset() > longest)
					longest = km.k+km.getKmerStartOffset();
			}
			// get the matched sequence
			String s = null;
			int pSeq = 0;
			if (p<RC/2){
				pSeq = p;
				s = (String) seq.subSequence(pSeq+leftShift, pSeq+longest);
			}
			else{
				pSeq = p-RC;
				s = (String) seq_rc.subSequence(pSeq+leftShift, pSeq+longest);
			}
			KmerGroup kg = config.kg_hit_adjust_type==2?
					config.use_weighted_kmer ? new KmerGroup(posCoveredWidth, negCoveredWidth, kmers, p, seq_weights) : new KmerGroup(posCoveredWidth, negCoveredWidth, kmers, p, null) :
					config.use_weighted_kmer ? new KmerGroup(posHitStrings, negHitStrings, kmers, p, s, seq_weights) : new KmerGroup(posHitStrings, negHitStrings, kmers, p, s, null);
			matches[idx]=kg;
			kg.setScore(computeSiteSignificanceScore(kg.getGroupHitCount(), kg.getGroupNegHitCount()));
			idx++;
		}
		return matches;
	}
	
	/** 
	 * Report all strand-specific KSM group hits in both orientations of the sequence<br>
	 * Assuming the initAhoCorasick(kmers) method had been called, i.e. "treeAhoCorasick kmer search tree" instance variable has been constructed.
	 * @param seq sequence string to search k-mers
	 * @param seq_rc sequence rc to search k-mers
	 * @return an array of KmerGroups, sorted by KG scores:<br>
	 * Each k-mer group maps to a binding position (using kmer.startOffset, relative to bs ) in the sequence<br>
	 * Note: the return value is different from query(), here the match on RC strand is labeled (pos+RC)
	 */
	public KmerGroup[] findKsmGroupHits (String seq, String seq_rc){
//		seq = seq.toUpperCase();
		//Search for all kmers in the sequences using Aho-Corasick algorithms (initialized)
		//ahocorasick_java-1.1.tar.gz is an implementation of Aho-Corasick automata for Java. BSD license.
		//from <http://hkn.eecs.berkeley.edu/~dyoo/java/index.html> 

		// matches on negative strand are combined with matches on positive strand
		HashMap<Integer, ArrayList<Kmer>> result = new HashMap<Integer, ArrayList<Kmer>> ();
		
		Iterator searcher = treeAhoCorasick.search(seq.getBytes());
		while (searcher.hasNext()) {
			SearchResult sr = (SearchResult) searcher.next();
			for (Object s: sr.getOutputs()){
				Kmer[] kmers = str2kmers.get(s);	
				// AC search returns end+1 position, end_seq; endOffset is seed_end;
				// thus   seed_end + end_seq --> seed_seq
				int[] kmerOffsets = str2kmerEndOffsets.get(s);	
				for(int i=0;i<kmers.length;i++){	
					int x = sr.getLastIndex() + kmerOffsets[i];	// get the motif position
					if (!result.containsKey(x))
						result.put(x, new ArrayList<Kmer>());
					result.get(x).add(kmers[i]);	
				}
			}
		}
		// the reverse compliment
		searcher = treeAhoCorasick.search(seq_rc.getBytes());
		while (searcher.hasNext()) {
			SearchResult sr = (SearchResult) searcher.next();
			for (Object s: sr.getOutputs()){
				Kmer[] kmers = str2kmers.get(s);	
				int[] kmerOffsets = str2kmerEndOffsets.get(s);	
				for(int i=0;i<kmers.length;i++){	
					int x = sr.getLastIndex() + kmerOffsets[i] + RC;	// get the motif position, +RC: "found on RC"
					if (!result.containsKey(x))
						result.put(x, new ArrayList<Kmer>());
					result.get(x).add(kmers[i]);	
				}
			}
		}
		if (result.isEmpty())
			return null;
		
		KmerGroup[] matches = new KmerGroup[result.keySet().size()];
		int idx = 0;
		for (int p:result.keySet()){
			ArrayList<Kmer> kmers = result.get(p);
			int leftShift = 999;			// left most shift position
			int longest = 0;				// length from the seed position
			for (Kmer km:kmers){
				if (km.getKmerStartOffset() < leftShift)
					leftShift = km.getKmerStartOffset();
				if (km.k+km.getKmerStartOffset() > longest)
					longest = km.k+km.getKmerStartOffset();
			}
			// get the matched sequence
			String s = null;
			int pSeq = 0;
			if (p<RC/2){
				pSeq = p;
				s = (String) seq.subSequence(pSeq+leftShift, pSeq+longest);
			}
			else{
				pSeq = p-RC;
				s = (String) seq_rc.subSequence(pSeq+leftShift, pSeq+longest);
			}
			KmerGroup kg = posCoveredWidth!=null?
					config.use_weighted_kmer ? new KmerGroup(posCoveredWidth, negCoveredWidth, kmers, p, seq_weights) : new KmerGroup(posCoveredWidth, negCoveredWidth, kmers, p, null) :
					config.use_weighted_kmer ? new KmerGroup(posHitStrings, negHitStrings, kmers, p, s, seq_weights) : new KmerGroup(posHitStrings, negHitStrings, kmers, p, s, null);
			matches[idx]=kg;
			kg.setScore(computeSiteSignificanceScore(kg.getGroupHitCount(), kg.getGroupNegHitCount()));
			idx++;
		}
		
		Arrays.sort(matches);		// sort by descending kgScore		
		return matches;
	}
	
	public String getSequenceUppercase(Region r){
		return seqgen.execute(r).toUpperCase();
	}

	
	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException{
		long tic = System.currentTimeMillis();
		ArrayList<String> pos_seqs = new ArrayList<String>();
		ArrayList<Double> seq_w = new ArrayList<Double>();

        Config config = new Config();
        try{
			config.parseArgs(args);   
		}
		catch (Exception e){
			e.printStackTrace();
    		System.exit(-1);
		}  
		
		String out_prefix = Args.parseString(args, "out_name", null);
		File outFolder = new File(out_prefix+"_outputs");
		outFolder.mkdir();
		out_prefix = new File(outFolder, out_prefix).getAbsolutePath();
		
		// read input fasta sequence, the fasta header line may optionally has a weight of the sequence
		String pos_file = Args.parseString(args, "pos_seq", null);
		String neg_file = Args.parseString(args, "neg_seq", null);
		if (pos_file==null || (config.k==-1&&config.k_min==-1)){
			System.err.println("Example: KMAC --pos_seq c-Myc_Crawford_HeLa-S3_61bp_GEM.fasta [--neg_seq c-Myc_Crawford_HeLa-S3_61bp_GEM_neg.fasta] --k_min 5 --k_max 8 --out_name cMyc_cMyc --seed CACGTG");
			System.exit(-1);
		}
		
        System.out.println("\nKMAC (version "+KMAC_VERSION+")");

        StringBuilder sb = new StringBuilder();
		sb.append("\nOptions:\n");
		for (String arg:args){
			if (arg.trim().indexOf(" ")!=-1)
				sb.append("\"").append(arg).append("\" ");
			else
				sb.append(arg).append(" ");
		}
		String params = sb.toString();
		System.out.println(params+"\n");
		
		String format = Args.parseString(args, "format", "fasta");
		ArrayList<String> strs = CommonUtils.readTextFile(pos_file);
        String[]f = null;
//        System.out.println("Loading positive sequences ...");
		for (String line: strs){
			if (format.equals("fasta")){
	            if (line.startsWith(">")){
	        		f = line.split("\\s+");
	        		if (f.length>1){
	        			try{
	        				seq_w.add(Double.parseDouble(f[1]));
	        			}catch(NumberFormatException nfe){
	        				seq_w.add(10.0);
	        			}
	        		}
		            else
		            	seq_w.add(10.0);
	        	}
	        	else{
	        		if (config.k_win==-1 || line.length()<=config.k_win){
	        			pos_seqs.add(line);
	        		}
	        		else{
		        		int left = line.length()/2-config.k_win/2;
		        		if (left<0)
		        			continue;
		        		pos_seqs.add(line.toUpperCase().substring(left, left+config.k_win));
	        		}
	        	}
			}
			else{		// simple format: seq<TAB or SPACE>weight
				f = line.split("\\s+");
        		if (f.length>1){
	            	seq_w.add(Double.parseDouble(f[1]));
	            	pos_seqs.add(f[0].toUpperCase());
	            }
			}
		}		
		
		ArrayList<String> neg_seqs = new ArrayList<String>();
		if (neg_file!=null){
//	        System.out.println("Loading negative sequences ...");
			strs = CommonUtils.readTextFile(neg_file);
			for (String line: strs){
				if (format.equals("fasta")){
					if (!line.startsWith(">")){
						if (config.k_win==-1 || line.length()<=config.k_win){
							neg_seqs.add(line);
		        		}
		        		else{
			        		int left = line.length()/2-config.k_win/2;
			        		if (left<0)
			        			continue;
			        		neg_seqs.add(line.toUpperCase().substring(left, left+config.k_win));
		        		}
	        		}
				}
				else{
		            neg_seqs.add(line.substring(0,config.k_win).toUpperCase());
				}
			}
		}
        
		// run motif discovery
		KMAC kmac = new KMAC();
        kmac.setStandalone();
        kmac.setConfig(config, out_prefix, params);
        System.out.println(String.format("Loaded %d input positive sequences.", pos_seqs.size()));
        kmac.setSequences(pos_seqs, neg_seqs, seq_w);
        pos_seqs.clear();
        neg_seqs.clear();
        System.out.println(String.format("Use top %d center sequences (%dbp) and %d negative sequences to find motif ...", 
        		kmac.seqs.length, config.k_win, kmac.seqsNegList.size()));
        kmac.discoverMotifs(config.k_min, config.k_max, null);
        if (config.verbose>1){
			System.out.println(StatUtil.cacheAccessCount);
			System.out.println(StatUtil.getCacheSize());
        }
        System.out.println("Done: "+CommonUtils.timeElapsed(tic));
	}
		
} 


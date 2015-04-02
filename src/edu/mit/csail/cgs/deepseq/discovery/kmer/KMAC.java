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
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import javax.imageio.ImageIO;

import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.SetTools;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.utils.strings.StringUtils;
import edu.mit.csail.cgs.utils.strings.multipattern.*;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.features.ComponentFeature;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScoreProfile.PwmMatch;
import edu.mit.csail.cgs.utils.stats.StatUtil;
import edu.mit.csail.cgs.deepseq.discovery.Config;

public class KMAC {
	private final int RC=100000;		// extra bp add to indicate negative strand match of kmer
	private final int UNALIGNED=9999;	// the special shift for unaligned kmer
	public static final char[] LETTERS = {'A','C','G','T'};
	public final int MAXLETTERVAL = Math.max(Math.max(Math.max('A','C'),Math.max('T','G')),
            Math.max(Math.max('a','c'),Math.max('t','g'))) + 1;
	
	private boolean standalone = false;
	public void setStandalone(){
		standalone = true;
	}
	Config config = new Config();
	
	private int verbose;
	private Genome genome;
	private boolean engineInitialized =false;
	private int k;
	private int minHitCount = 2;
	private int numPos;
	private double[] bg= new double[4];	// background frequency based on GC content
	private double ic_trim = 0.4;
	private String outName;
	private boolean use_PWM_MM = false;
	private boolean use_smart_mm = false;	
	
	private double seedOverrideScoreDifference=1.1;
	private Kmer primarySeed = null;
	
	private int k_win;
	private double[] profile;
	private boolean isMasked;
	private String[] seqs;			// DNA sequences around binding sites
	private double[] seq_weights;	// sequence hit count weighted by binding strength, corresponding to the seqs[]
	public double[] getSequenceWeights(){return seq_weights;}
	public void setSequenceWeights(double[] w){seq_weights=w;}
	private double totalWeight;
	private String[] seqsNeg;		// DNA sequences in negative sets for ALL event regions, but only those in seqsNegList are used for motif discovery
	private ArrayList<String> seqsNegList=new ArrayList<String>(); // Effective negative sets, excluding overlaps in positive sets
	public int getNegSeqCount(){return negSeqCount;}
	public int getPosSeqCount(){return posSeqCount;}
    private int posSeqCount;
    private int negSeqCount;
    public void setTotalSeqCount(int pos, int neg){
    	posSeqCount = pos;
    	negSeqCount = neg;
    }
	private int negRegionDistance;
	/** region-->index_seqsNeg for negative sequences */
	private TreeMap<Region, Integer> neg_region_map;
	public String[] getPositiveSeqs(){return seqs;};
	public double get_NP_ratio(){return (double)negSeqCount/posSeqCount;}
	
//	private ArrayList<Kmer> allKmers = new ArrayList<Kmer>();	// all the kmers in the sequences
//	public ArrayList<Kmer> getAllKmers() {
//		return allKmers;
//	}
	private HashMap<String, Kmer> str2kmer = new HashMap<String, Kmer>();
	
	// AhoCorasick algorithm for multi-pattern search
	// Pre-processing is to build the tree with all the patters (kmers)
	// Then each individual search can be done in scan()
	private AhoCorasick tree;
	private AhoCorasick tree_negatives;
	
	public boolean isInitialized(){ return engineInitialized;}
	
	private SequenceGenerator<Region> seqgen;
	private long tic;	

	/** motif clusters */
	ArrayList<KmerCluster> clusters = new ArrayList<KmerCluster>();
	public KmerCluster getPrimaryCluster(){
		if (clusters.size()>=1)
			return clusters.get(0);
		return null;
	}
	public ArrayList<KmerCluster> getMotifClusters(){
		return clusters;
	}
	
	public KMAC(Config config, String outPrefix){
		setConfig(config, outPrefix);
	}
	private void setConfig(Config config, String outPrefix){
		this.config = config;
	    this.outName = outPrefix;
	    this.verbose = config.verbose;
	    Kmer.set_use_weighted_hit_count(config.use_weighted_kmer);
	}

	public void updateOutPrefix(String outPrefix){
	    this.outName = outPrefix;
	}
	
	public void setTotalSeqCounts(int posSeqCount, int negSeqCount){
		this.posSeqCount = posSeqCount;
		this.negSeqCount = negSeqCount;
	}
	
	
	// called by standalone main() method
	public void setSequences(ArrayList<String> pos_seqs, ArrayList<String> neg_seqs, ArrayList<Double> pos_w){
		int seqNum = Math.min(pos_seqs.size(), config.k_seqs);
		seqs = new String[seqNum];	
		for (int i=0;i<seqNum;i++)
			seqs[i] = pos_seqs.get(i);
		
		seq_weights = new double[seqNum];
		totalWeight=0;
		for (int i=0;i<seqNum;i++){
			switch (config.seq_weight_type){
			case 0:	seq_weights[i]=1;break;
			case 1: seq_weights[i]=pos_w.get(i);break;
			case 2: seq_weights[i]=Math.sqrt(pos_w.get(i));break;
			case 3: seq_weights[i]=Math.log(pos_w.get(i));break;
			default: System.err.println("Sequence weighting type is not defined!");System.exit(-1);
			}
			totalWeight += seq_weights[i];
		}
		for (int i=0;i<seq_weights.length;i++){
			seq_weights[i] = seq_weights[i]*seqs.length/totalWeight;	// scale weights with total sequence count, and total weight
		}
		if (config.use_weighted_kmer)
			Kmer.set_seq_weights(seq_weights);
		
		if (neg_seqs.isEmpty()){
			if (config.k_neg_dinu_shuffle){
				System.out.println("Use di-nucleotide shuffled sequences as negative sequences.\n");
				Random randObj = new Random(config.rand_seed);
				for (int i=0;i<seqNum;i++)
					seqsNegList.add(SequenceUtils.dinu_shuffle(seqs[i], randObj));
			}
			else{	// single nucleotide shuffling
				System.out.println("Use shuffled sequences as negative sequences.\n");
				Random randObj = new Random(config.rand_seed);
				for (int i=0;i<seqNum;i++)
					seqsNegList.add(SequenceUtils.shuffle(seqs[i], randObj));
			}
		}
		else{
			if (neg_seqs.size()>seqNum)
				for (int i=0;i<seqNum;i++)
					seqsNegList.add(neg_seqs.get(i));
			else
				seqsNegList = neg_seqs;
		}
		posSeqCount = seqs.length;
	    negSeqCount = seqsNegList.size();
	    updateSequenceInfo();
	}
	private void resetProfile(){
		for (int i=0; i<profile.length; i++)
	    	profile[i] = 1;
	}
	
	private void updateSequenceInfo(){
		k_win = seqs[0].length();
		isMasked = false;
		
		// logistic distribution to fit the spatial resolution shape, with a more heavy tail than Gaussian
		// http://en.wikipedia.org/wiki/Logistic_distribution
		// ctcf_sigma = 9.53; GABP_sigma = 15.98;
	    profile = new double[k_win];
	    double sigma = 13;
	    for (int i=0; i<=k_win/2; i++){
	    	double e = Math.exp(-i/sigma);
	    	profile[k_win/2-i] = e/(sigma*(1+e)*(1+e));
	    	profile[k_win/2+i] = profile[k_win/2-i];
	    }
	    StatUtil.normalize(profile);
//	   	System.out.println(CommonUtils.arrayToString(profile, "%.4f"));
	   	
	    // count cg-content
		int gcCount = 0;
		for (String seq:seqsNegList){
			for (char c:seq.toCharArray())
				if (c=='C'||c=='G')
					gcCount ++;
		}
		double gcRatio = (double)gcCount/negSeqCount/k_win;
		bg[0]=0.5-gcRatio/2; 
    	bg[1]=gcRatio/2; 
    	bg[2]=bg[1]; 
    	bg[3]=bg[0];
	}
	
	public KMAC(Genome g, boolean useCache, boolean use_db_genome, String genomePath, Config config, String outPrefix){
		setConfig(config, outPrefix);
		genome = g;
		seqgen = new SequenceGenerator<Region>();
    	if (use_db_genome)
    		seqgen.useLocalFiles(false);
		if (useCache)
			seqgen.useCache(true);	
		seqgen.setGenomePath(genomePath);
	}
	
	/* 
	 * Contruct a Kmer Engine from a list of Kmers
	 */
	public KMAC(ArrayList<Kmer> kmers, Config config, String outPrefix){
		setConfig(config, outName);
		if (!kmers.isEmpty()){
			if (outPrefix!=null)
				updateEngine(kmers, outPrefix);
			else
				updateEngine(kmers);
			k=kmers.get(0).getK();
		}
	}
	
	/* 
	 * Contruct a KMAC object from a list of KSM motif, used in GEM for setting the motif positional prior 
	 */
	public KMAC(KmerSet kmerset, Config config, String outPrefix){
		setConfig(config, outPrefix);
		setTotalSeqCount(kmerset.posSeqCount, kmerset.negSeqCount);
		Set<Integer> clusterIds = kmerset.getClusterIds();
		for (int id: clusterIds){
			KmerCluster cluster = new KmerCluster();
			cluster.alignedKmers = kmerset.getKmers(id);
			cluster.ksmThreshold.score = kmerset.ksmThreshold;
			clusters.add(cluster);
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
	/**
		 * Load pos/neg test sequences based on event positions<br>
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
				if (config.strand_type==1 && events.get(i).getStrand()=='-')
					seq = SequenceUtils.reverseComplement(seq);
				posSeqs.add(seq.toUpperCase());		// if repeat_fraction>=1, allow all repeats, convert to upper case
				switch (config.seq_weight_type){
				case 0:	posSeqWeights.add(1.0);break;
				case 1: posSeqWeights.add(events.get(i).getTotalEventStrength());break;
				case 2: posSeqWeights.add(Math.sqrt(events.get(i).getTotalEventStrength()));break;
				case 3: posSeqWeights.add(Math.log(events.get(i).getTotalEventStrength()));break;
				default: System.err.println("Sequence weighting type is not defined! No weighting.");posSeqWeights.add(1.0);
				}
				posImpactRegion.add(events.get(i).getPeak().expand(negRegionDistance));
			}
			seqs = new String[posSeqs.size()];	// DNA sequences around binding sites
			posSeqs.toArray(seqs);
		    seq_weights = new double[posSeqs.size()];

			totalWeight=0;
			for (int i=0;i<seq_weights.length;i++){
				seq_weights[i]=posSeqWeights.get(i);
				totalWeight += seq_weights[i];
			}
			for (int i=0;i<seq_weights.length;i++){
				seq_weights[i] = seq_weights[i]*seqs.length/totalWeight;	// scale weights with total sequence count, and total weight
			}
			if (config.use_weighted_kmer)
				Kmer.set_seq_weights(seq_weights);
			
			seqsNegList.clear();
			if (config.k_neg_dinu_shuffle){
				System.out.println("Use di-nucleotide shuffled sequences as negative sequences.\n");
				Random randObj = new Random(config.rand_seed);
				for (int i=0;i<seqs.length;i++)
					seqsNegList.add(SequenceUtils.dinu_shuffle(seqs[i], randObj));
			}
			else{
				/** Negative sequences has been retrieved when setting up region caches */
				ArrayList<Region> negRegions = new ArrayList<Region>();
				negRegions.addAll(neg_region_map.keySet());
				Region.filterOverlapRegions(negRegions, posImpactRegion);	// make sure negative region is not within negRegionDistance of positive regions.
				int negCount = 0;
				int len = winSize/2*2+1;
				for (Region r:negRegions){
					String seq_retrieved = seqsNeg[neg_region_map.get(r)];
					if (seq_retrieved.length()<len)
						continue;
					String seq = seq_retrieved.substring(0, len);
					if (config.repeat_fraction<1){
						int count = 0;
						for (char c:seq.toCharArray())
							if (Character.isLowerCase(c) || c=='N')
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
					seqsNegList.add(seq.toUpperCase());		// if repeat_fraction>=1, allow repeats, convert to upper case
					negCount++;
					if (negCount==seqs.length)				// limit the neg region count to be same or less than positive region count
						break;
				}
			}
			posSeqCount = seqs.length;
			seqsNegList.trimToSize();
		    negSeqCount = seqsNegList.size();
		    if (verbose>1 || config.repeat_fraction!=1)
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
//	/*
//	 * Find significant Kmers that have high HyperGeometric p-value
//	 * in sequence around the binding events 
//	 * and build the kmer AhoCorasick engine
//	 */
//	public void buildEngine(int k, ArrayList<ComponentFeature> events, int winSize, double hgp, double k_fold, String outPrefix, boolean print_kmer_hits){
//		ArrayList<Kmer> kms = selectEnrichedKmers(k, events, winSize, hgp, k_fold, outPrefix);
//		updateEngine(kms, outPrefix);
//	}
//	
//	public ArrayList<Kmer> selectEnrichedKmers(int k, ArrayList<ComponentFeature> events, int winSize, double hgp, double k_fold, String outName){
//		this.hgp = hgp;
//		this.outName = outName;
//		this.k_fold = k_fold;
//
//		// collect pos/neg test sequences based on event positions
//		loadTestSequences(events, winSize);
//		
//		return selectEnrichedKmers(k);
//	}
	
	/** 
	 * Select the value of k <br>
	 * that forms cluster with the best HGP
	 */
	public int selectK(int k_min, int k_max, int[] eventCounts){
		if (k_min==k_max)
			return k_min;
		
		String[] pos_seq_backup = seqs.clone();
		String[] neg_seqs_backup = new String[seqsNegList.size()];
		seqsNegList.toArray(neg_seqs_backup);
		
		// compare different values of k to select most enriched k value
		int bestK = 0;
		double bestAllHGP = 0;
		ArrayList<KmerCluster> kClusters = new ArrayList<KmerCluster>();
		StringBuilder sb = new StringBuilder("\n------------------- "+ new File(outName).getName() +" ----------------------\n");
		for (int i=0;i<k_max-k_min+1;i++){
			int k = i+k_min;
			System.out.println("\n----------------------------------------------\nTrying k="+k+" ...\n");
			ArrayList<Kmer> kmers = selectEnrichedKmers(k);
			KmerMotifAlignmentClustering(kmers, 2, false, null,"");
			double bestclusterHGP = 0;
			KmerCluster bestCluster=null;
			for (KmerCluster c:clusters){
				if (bestclusterHGP>c.pwmThresholdHGP){
					bestclusterHGP=c.pwmThresholdHGP;
					bestCluster = c;
				}
			}
			if (bestCluster!=null){
				sb.append(String.format("k=%d\thit=%d+/%d-\thgp=1e%.1f\tW=%d\tPWM=%s.\n", k, bestCluster.pwmPosHitCount, bestCluster.pwmNegHitCount, 
						bestCluster.pwmThresholdHGP, bestCluster.wm.length(), WeightMatrix.getMaxLetters(bestCluster.wm)));
				kClusters.add(bestCluster);
				
				if (bestAllHGP>bestCluster.pwmThresholdHGP)
					bestAllHGP=bestCluster.pwmThresholdHGP;
			}
			else
				sb.append(String.format("k=%d\tcannot form a PWM.\n", k));
			
			// reload non-masked sequences
			seqs = pos_seq_backup.clone();
			seqsNegList.clear();
			for (String s:neg_seqs_backup)
				seqsNegList.add(s);
		}
		
		boolean failToFindK = false;
		KmerCluster bestCluster=null;
		if (kClusters.isEmpty()){
			failToFindK = true;
		}
		else{
			// find the k value with the best HGP
			int bestHitCount = 0;
			for (KmerCluster c : kClusters){
				if (c.pwmThresholdHGP==bestAllHGP){
					bestHitCount = c.pwmPosHitCount;
					bestCluster = c;
					break;
				}
			}
			if (bestCluster==null)
				failToFindK = true;			// failed to find motif for all k-values
		}
		
		if (failToFindK){
			System.out.println("\n----------------------------------------------\nNone of the k values form an enriched PWM, stop here!\n");
			File f = new File(outName);
			String name = f.getName();
			StringBuffer html = new StringBuffer("<style type='text/css'>/* <![CDATA[ */ table, td{border-color: #600;border-style: solid;} table{border-width: 0 0 1px 1px; border-spacing: 0;border-collapse: collapse;} td{margin: 0;padding: 4px;border-width: 1px 1px 0 0;} /* ]]> */</style>");
			html.append("<script language='javascript' type='text/javascript'><!--\nfunction popitup(url) {	newwindow=window.open(url,'name','height=75,width=400');	if (window.focus) {newwindow.focus()}	return false;}// --></script>");
			html.append("<table><th bgcolor='#A8CFFF'><font size='5'>");
			html.append(name).append("</font></th>");
			html.append("<tr><td valign='top' width='500'><br>");
			if (!this.standalone && eventCounts!=null){
				html.append("<a href='"+name+"_GEM_events.txt'>Significant Events</a>&nbsp;&nbsp;: "+eventCounts[0]);
				html.append("<br><a href='"+name+"_GEM_insignificant.txt'>Insignificant Events</a>: "+eventCounts[1]);
				html.append("<br><a href='"+name+"_GEM_filtered.txt'>Filtered Events</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;: "+eventCounts[2]);
			}
			html.append("<p><p>Motif can not be found!<p>");
			html.append("</td></tr></table>");
			CommonUtils.writeFile(outName+"_result.htm", html.toString());
			return 0;
		}
		
		bestK = bestCluster.seedKmer.getK();
		System.out.print(sb.toString());
		System.out.println(String.format("\nSelected k=%d\thit=%d\thgp=1e%.1f.\n----------------------------------------------\n", 
				bestK, bestCluster.pwmPosHitCount, bestCluster.pwmThresholdHGP));
//		
//		// check if there is discontinuity in widths 
//		ArrayList<Integer> widths = new ArrayList<Integer>();
//		for (KmerCluster c : kClusters){
//			widths.add(c.wm.length());
//		}
//		Collections.sort(widths);
//		ArrayList<Integer> continuous = new ArrayList<Integer>();
//		continuous.add(widths.get(0));
//		for (int i=0;i<widths.size()-1;i++){
//			if (widths.get(i+1)-widths.get(i)>1)	// a gap in width, could match monomer or dimer, such as Oct4_Sox2
//				break;
//			else
//				continuous.add(widths.get(i+1));
//		}
//		widths.clear();
//		widths.addAll(continuous);
//		// find the most freq pwm width(s)
//		Pair<int[], int[]> sorted = StatUtil.sortByOccurences(widths);
//		int maxCount = sorted.cdr()[sorted.cdr().length-1];
//		HashSet<Integer> consensusWidths = new HashSet<Integer>();
//		for (int i=0;i<sorted.cdr().length;i++){
//			if (sorted.cdr()[i]==maxCount)
//				consensusWidths.add(sorted.car()[i]);
//		}
//
//		// find the bestHGP with consensus width
//		bestHGP = 0;
//		for (KmerCluster c : kClusters){
//			if (consensusWidths.contains(c.wm.length()) && c.pwmThresholdHGP<bestHGP)
//				bestHGP = c.pwmThresholdHGP;
//		}
//		ArrayList<KmerCluster> goodClusters = new ArrayList<KmerCluster>();			// clusters with hgp >= 95% of bestHGP
//		for (KmerCluster c:kClusters){
//			if (consensusWidths.contains(c.wm.length()) && c.pwmThresholdHGP<=bestHGP*0.90)
//				goodClusters.add(c);
//		}
//		Collections.sort(goodClusters, new Comparator<KmerCluster>(){
//		    public int compare(KmerCluster o1, KmerCluster o2) {
//		    	return o1.compareForSelectingK(o2);
//		    }
//		});	
//		bestK = goodClusters.get(0).seedKmer.getK();
//		
//		System.out.print(sb.toString());
//		System.out.println(String.format("\nSelected k=%d\tW=%d\tbestHGP=%.1f.\n----------------------------------------------\n", 
//				bestK, goodClusters.get(0).wm.length(),goodClusters.get(0).pwmThresholdHGP));
		return bestK;
	}
	
	/** 
	 * Select the value of k <br>
	 * using the seed family hgp
	 */
	public int selectK_byTopKmer(int k_min, int k_max, int[] eventCounts){
		if (k_min==k_max)
			return k_min;
		
		String[] pos_seq_backup = seqs.clone();
		String[] neg_seqs_backup = new String[seqsNegList.size()];
		seqsNegList.toArray(neg_seqs_backup);
		
		// compare different values of k to select most enriched k value
		int bestK = 0;
		double bestAllHGP = 0;
		ArrayList<KmerCluster> kClusters = new ArrayList<KmerCluster>();
		StringBuilder sb = new StringBuilder("\n------------------- "+ new File(outName).getName() +" ----------------------\n");
		for (int i=0;i<k_max-k_min+1;i++){
			int k = i+k_min;
			System.out.println("\n----------------------------------------------\nTrying k="+k+" ...\n");
			ArrayList<Kmer> kmers = selectEnrichedKmers(k);
			KmerMotifAlignmentClustering(kmers, 2, true, null, "");	// seed kmer only
			double bestclusterHGP = 0;
			KmerCluster bestCluster=null;
			for (KmerCluster c:clusters){
				if (bestclusterHGP>c.seedKmer.familyHgp){		// use the seed family hgp to select k
					bestclusterHGP=c.seedKmer.familyHgp;
					bestCluster = c;
				}
			}
			if (bestCluster!=null){
				sb.append(String.format("k=%d\thgp=1e%.1f\tseed=%s.\n", k, bestCluster.seedKmer.familyHgp, bestCluster.seedKmer.getKmerString()));
				kClusters.add(bestCluster);
				
				if (bestAllHGP>bestCluster.seedKmer.familyHgp)
					bestAllHGP=bestCluster.seedKmer.familyHgp;
			}
			else
				sb.append(String.format("k=%d\tcan not form a seed family.\n", k));
			
			// reload non-masked sequences
			seqs = pos_seq_backup.clone();
			seqsNegList.clear();
			for (String s:neg_seqs_backup)
				seqsNegList.add(s);
		}
		if (kClusters.isEmpty()){
			System.out.println("\n----------------------------------------------\nNone of the k values form an enriched seed family, stop here!\n");
			File f = new File(outName);
			String name = f.getName();
			StringBuffer html = new StringBuffer("<style type='text/css'>/* <![CDATA[ */ table, td{border-color: #600;border-style: solid;} table{border-width: 0 0 1px 1px; border-spacing: 0;border-collapse: collapse;} td{margin: 0;padding: 4px;border-width: 1px 1px 0 0;} /* ]]> */</style>");
			html.append("<script language='javascript' type='text/javascript'><!--\nfunction popitup(url) {	newwindow=window.open(url,'name','height=75,width=400');	if (window.focus) {newwindow.focus()}	return false;}// --></script>");
			html.append("<table><th bgcolor='#A8CFFF'><font size='5'>");
			html.append(name).append("</font></th>");
			html.append("<tr><td valign='top' width='500'><br>");
			if (!this.standalone && eventCounts!=null){
				html.append("<a href='"+name+"_GEM_events.txt'>Significant Events</a>&nbsp;&nbsp;: "+eventCounts[0]);
				html.append("<br><a href='"+name+"_GEM_insignificant.txt'>Insignificant Events</a>: "+eventCounts[1]);
				html.append("<br><a href='"+name+"_GEM_filtered.txt'>Filtered Events</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;: "+eventCounts[2]);
			}
			html.append("<p><p>Motif can not be found!<p>");
			html.append("</td></tr></table>");
			CommonUtils.writeFile(outName+"_result.htm", html.toString());
			
			return 0;
		}
		
		
		KmerCluster bestCluster=null;
		for (KmerCluster c : kClusters){
			if (bestAllHGP==c.seedKmer.familyHgp){
				bestCluster = c;
				break;
			}
		}
		bestK = bestCluster.seedKmer.getK();
		System.out.print(sb.toString());
		System.out.println(String.format("\nSelected k=%d\thit=%d\thgp=1e%.1f.\n----------------------------------------------\n", 
				bestK, bestCluster.pwmPosHitCount, bestCluster.seedKmer.familyHgp));
		return bestK;
	}
	
	
	/** 
	 * Select the value of k by the coverage of enriched k-mers<br>
	 * Count the exact base of coverage from every sequence
	 * @param k_min
	 * @param k_max
	 * @return
	 */
	public int selectKbyCoverage(int k_min, int k_max, ArrayList<ComponentFeature> events, int winSize){
		ArrayList<Integer> widths = new ArrayList<Integer>(); 
		// compare different values of k to select most enriched k value
		int bestK = 0;
		int bestCoverage = 0;
		for (int i=0;i<k_max-k_min+1;i++){

			if (isMasked && events!=null)
				loadTestSequences(events, winSize);
			
			int k = i+k_min;
			ArrayList<Kmer> kmers = selectEnrichedKmers(k);
			
			/* Initialization, setup sequence list, and update kmers */
			// build the kmer search tree
			AhoCorasick oks = new AhoCorasick();
			for (Kmer km: kmers){
				oks.add(km.getKmerString().getBytes(), km);
		    }
			oks.prepare();
			
			int posCoverage = 0;
			int negCoverage = 0;
			for (String seq:seqs){
				int[] chars = new int[seq.length()];
				HashSet<Kmer> results = queryTree (seq, oks, config.strand_type==1);
				if (!results.isEmpty()){
					for (Kmer km: results){		
						ArrayList<Integer> pos = StringUtils.findAllOccurences(seq, km.getKmerString());
						for (int p:pos){
							for (int j=p;j<p+k;j++)
								chars[j] = 1;
						}
						pos = StringUtils.findAllOccurences(seq, km.getKmerRC());
						for (int p:pos){
							for (int j=p;j<p+k;j++)
								chars[j] = 1;
						}
					}
				}
				for (int c:chars)
					posCoverage+=c;
			}
			for (String seq:seqsNegList){
				int[] chars = new int[seq.length()];
				HashSet<Kmer> results = queryTree (seq, oks, config.strand_type==1);
				if (!results.isEmpty()){
					for (Kmer km: results){		
						ArrayList<Integer> pos = StringUtils.findAllOccurences(seq, km.getKmerString());
						for (int p:pos){
							for (int j=p;j<p+k;j++)
								chars[j] = 1;
						}
						pos = StringUtils.findAllOccurences(seq, km.getKmerRC());
						for (int p:pos){
							for (int j=p;j<p+k;j++)
								chars[j] = 1;
						}
					}
				}
				for (int c:chars)
					negCoverage+=c;
			}
			
			int coverage = posCoverage-negCoverage;

			System.out.println(String.format("k=%d, \tcoverage=%d.", k, coverage));
			if (bestCoverage<coverage){
				bestCoverage=coverage;
				bestK = k;
			}
		}

		System.out.println(String.format("\n------------------------\nSelected k=%d\tcoverage=%d.", bestK, bestCoverage));
		return bestK;
	}
	/** 
	 * Select the value of k by the coverage of enriched k-mers<br>
	 * approximate the coverage by sum(kmerHitCount * k)
	 * @param k_min
	 * @param k_max
	 * @return
	 */
	public int selectK_old(int k_min, int k_max){
		// compare different values of k to select most enriched k value			
		ArrayList<Integer> widths = new ArrayList<Integer>(); 
		// compare different values of k to select most enriched k value
		int bestK = 0;
		int bestKxKmerCount = 0;
		for (int i=0;i<k_max-k_min+1;i++){
			int k = i+k_min;
			ArrayList<Kmer> kmers = selectEnrichedKmers(k);
			int kmerCount = 0;
			for (Kmer km:kmers)				// kmerCount as the total difference between positive and negative hit count
				kmerCount += km.getPosHitCount()-km.getNegHitCount();
			System.out.println(String.format("k=%d, \tcoverage=%d.", k, k*kmerCount));
			if (bestKxKmerCount<k * kmerCount){
				bestKxKmerCount=k * kmerCount;
				bestK = k;
			}
		}

		System.out.println(String.format("\n------------------------\nSelected k=%d\tcoverage=%d.", bestK, bestKxKmerCount));
		return bestK;
	}
	/** 
	 * Index k-mers from the positive sequences, select enriched k-mers
	 * */
	public ArrayList<Kmer> selectEnrichedKmers(int k){
		this.k = k;
		// expected count of kmer = total possible unique occurences of kmer in sequence / total possible kmer sequence permutation
		tic = System.currentTimeMillis();
		numPos = k_win-k+1;

		// only derive k-mers from forward sequence, will merge kmer and its RC next if consider both strands
		HashMap<String, HashSet<Integer>> kmerstr2seqs = new HashMap<String, HashSet<Integer>>();
		for (int seqId=0;seqId<posSeqCount;seqId++){
			String seq = seqs[seqId];
			HashSet<String> uniqueKmers = new HashSet<String>();			// only count repeated kmer once in a sequence
			for (int i=0;i<numPos;i++){
				if ((i+k)>seq.length()) // endIndex of substring is exclusive
					break;
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
		
		ArrayList<Kmer> kms = new ArrayList<Kmer>();
		ArrayList<String> kmerStrings = new ArrayList<String>();
		kmerStrings.addAll(kmerstr2seqs.keySet());
		
		// Merge kmer and its reverse compliment (RC)	
		if (config.strand_type != 1){
			// create kmers from its and RC's counts
			for (String key:kmerStrings){
				if (!kmerstr2seqs.containsKey(key))		// this kmer has been removed, represented by RC
					continue;
				// consolidate kmer and its reverseComplment kmer
				String key_rc = SequenceUtils.reverseComplement(key);				
				if (!key_rc.equals(key)){	// if it is not reverse compliment itself
					if (kmerstr2seqs.containsKey(key_rc)){
						int kCount = kmerstr2seqs.get(key).size();
						int rcCount = kmerstr2seqs.get(key_rc).size();
						String winner = kCount>=rcCount?key:key_rc;
						String loser = kCount>=rcCount?key_rc:key;
						kmerstr2seqs.get(winner).addAll(kmerstr2seqs.get(loser));	// winner take all
						kmerstr2seqs.remove(loser);					// remove the loser kmer because it is represented by its RC
					}
				}
			}
		}
		System.out.println("k="+k+", mapped "+kmerstr2seqs.keySet().size()+" k-mers, "+CommonUtils.timeElapsed(tic));

		// compute the smallest PosCount needed to be significant, even with negCount=0
		int smallestPosCount;
		for (smallestPosCount=minHitCount;smallestPosCount<posSeqCount;smallestPosCount++){
			double hgp = computeHGP(posSeqCount, negSeqCount, smallestPosCount, 0);
			if (hgp<config.kmer_hgp){
				break;
			}
		}
		int expectedCount = (int) Math.max( smallestPosCount, Math.round(seqs.length*2*(seqs[0].length()-k+1) / Math.pow(4, k)));
		if (config.strand_type == 1){
			expectedCount /= 2;
		}
		// the purpose of expectedCount is to limit kmers to run HGP test, if total kmer number is low, it can be relaxed a little
		if (kmerstr2seqs.keySet().size()<10000){	
			expectedCount = Math.min(smallestPosCount, expectedCount);
		}
		ArrayList<String> kstrs = new ArrayList<String>();
		for (String key:kmerstr2seqs.keySet()){	
			if (kmerstr2seqs.get(key).size()< expectedCount)
				continue;	// skip low count kmers 
			kstrs.add(key);
		}		
		System.out.println("Expected kmer hit count="+expectedCount + ", kmer numbers="+kstrs.size());

		// create the kmer object
		for (String s:kstrs){	
			Kmer kmer = new Kmer(s, kmerstr2seqs.get(s));
			kms.add(kmer);
		}
		kms.trimToSize();
		Collections.sort(kms);
		kmerstr2seqs=null;	// clean up
		System.gc();
		
		/**
		 * Select significantly over-representative kmers 
		 * Search the kmer counts in the negative sequences, then compare to positive counts
		 */
		tic = System.currentTimeMillis();
		//Aho-Corasick for searching Kmers in negative sequences
		//ahocorasick_java-1.1.tar.gz is an implementation of Aho-Corasick automata for Java. BSD license.
		//from <http://hkn.eecs.berkeley.edu/~dyoo/java/index.html> 
		AhoCorasick tmp = new AhoCorasick();
		for (String s: kstrs){
			tmp.add(s.getBytes(), s);
	    }
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
			if (config.strand_type != 1){	// if not want single strand only 
				String seq_rc = SequenceUtils.reverseComplement(seq);
				searcher = tmp.search(seq_rc.getBytes());
				while (searcher.hasNext()) {
					SearchResult result = (SearchResult) searcher.next();
					kmerHits.addAll(result.getOutputs());
				}
			}
			for (Object o: kmerHits){
				String kmer = (String) o;
				if (!kmerstr2negSeqs.containsKey(kmer))					
					kmerstr2negSeqs.put(kmer, new HashSet<Integer>());
				kmerstr2negSeqs.get(kmer).add(negSeqId);
			}
		}
		
		// score the kmers, hypergeometric p-value
		ArrayList<Kmer> highHgpKmers = new ArrayList<Kmer>();
		for (Kmer kmer:kms){
			if (kmerstr2negSeqs.containsKey(kmer.getKmerString())){
				kmer.setNegHits(kmerstr2negSeqs.get(kmer.getKmerString()));				
			}
			if (kmer.getPosHitCount() < kmer.getNegHitCount()/get_NP_ratio() * config.k_fold ){
				highHgpKmers.add(kmer);	
				continue;
			}
			if (config.use_weighted_kmer){
//				kmer.setWeightedPosHitCount();		//TODO: check if weight is used
				kmer.setHgp(computeHGP(posSeqCount, negSeqCount, kmer.getWeightedHitCount(), kmer.getNegHitCount()));
			}
			else{
				kmer.setHgp(computeHGP(posSeqCount, negSeqCount, kmer.getPosHitCount(), kmer.getNegHitCount()));
			}
			if (kmer.getHgp()>config.kmer_hgp){
				highHgpKmers.add(kmer);		
				continue;
			}
		}

		Collections.sort(kms);
		// print high p-value k-mers, or low fold enrichment k-mers, although remove them for further learning
		if (config.print_all_kmers)
			Kmer.printKmers(kms, posSeqCount, negSeqCount, 0, outName+"_all_w"+seqs[0].length(), true, false, true);
		kms.removeAll(highHgpKmers);
		kms.trimToSize();
		System.out.println(String.format("k=%d, selected %d k-mers from %d+/%d- sequences, %s", k, kms.size(), posSeqCount, negSeqCount, CommonUtils.timeElapsed(tic)));
		
		return kms;
	}

	/**
	 * This is the main method for KMAC motif discovery
	 */
	public ArrayList<Kmer> KmerMotifAlignmentClustering (ArrayList<Kmer> kmers_in, int topCluster, boolean only_seed_kmer, int[] eventCounts, String gemMsg){
		int seed_range = k;
		String[] pos_seq_backup = seqs.clone();
		String[] neg_seq_backup = new String[seqsNegList.size()];
		seqsNegList.toArray(neg_seq_backup);
		
		tic = System.currentTimeMillis();
		if (kmers_in.size()==0)
			return kmers_in;

		// clone to modify locally
		ArrayList<Kmer> kmers = new ArrayList<Kmer>();
		for (Kmer km:kmers_in)
			kmers.add(km.clone());
		kmers.trimToSize();
		
		ArrayList<Sequence> seqList = new ArrayList<Sequence>();
		for (int i=0;i<seqs.length;i++){
			Sequence s = new Sequence(seqs[i], i);
			seqList.add(s);
		}
		seqList.trimToSize();
		
    	clusters.clear();
    	int clusterID = 0;
    	boolean quick_restart = false;
    	boolean primarySeed_is_immutable = false;
    	// Keep iterating to find all KmerSet motifs upto #max_cluster
		while (!kmers.isEmpty() && clusterID<=config.max_cluster){
			
			if (topCluster!=-1){			// only generate a few clusters to select optimal K
				int pwmCount = 0;
				for (KmerCluster c:clusters)
					if (c.wm!=null && c.pwmGoodQuality && c.total_aligned_seqs>=seqs.length*config.motif_hit_factor_report)
						pwmCount++;
				if (pwmCount>=topCluster){
					primarySeed = null;		// do not record primary seed here
					return null;
				}
				if (clusterID>10){
					primarySeed = null;		// do not record primary seed here
					return null;
				}
			}
			
			/** Initialization of new cluster and the remaining kmers */
			KmerCluster cluster = new KmerCluster();
			clusters.add(cluster);
			
			if (!quick_restart){
				indexKmerSequences(kmers, seqList);
				
				if (kmers.isEmpty()){
					clusters.remove(clusterID);
					break;
				}
	
//				Collections.sort(kmers, new Comparator<Kmer>(){
//				    public int compare(Kmer o1, Kmer o2) {
//				    	return o1.compareByHGP(o2);
//				    }
//				});	
				quick_restart = false;
			}
			
			if (verbose>1)
				System.out.println("------------------------------------------------\n"+CommonUtils.timeElapsed(tic)+
						": Aligning cluster #"+clusterID+",   n="+kmers.size());
			
			/** get seed kmer */
			Kmer seed=null;
			boolean isSeedInherited=false;
			if (clusterID==0){
				cluster.clusterId = 0;
				if (config.seed!=null){
					String rc = SequenceUtils.reverseComplement(config.seed);
					for (Kmer km:kmers){
						if (km.getKmerString().equals(config.seed)|| km.getKmerString().equals(rc)){
							if (km.getKmerString().equals(rc))
								km.setKmerString(config.seed);
							seed = km;
							primarySeed_is_immutable = true;		
							System.out.println("\nUse specified seed k-mer: "+seed.toShortString());
							break;
						}
					}
					if (!primarySeed_is_immutable){	// Use specified seed k-mer does not match with enriched k-mers
						seed = selectBestKmer(kmers);
						System.out.println("\nSpecified seed k-mer "+config.seed+" is not found, pick top k-mer");
					}
				}
				else{			// no user specified starting k-mer
					if (config.allow_seed_inheritance && primarySeed!=null){
						String bestStr = primarySeed.getKmerString();
						String bestRc = primarySeed.getKmerRC();
						for (Kmer km:kmers){
							if (km.getKmerString().equals(bestStr)||
									km.getKmerString().equals(bestRc)){
								seed = km;
								isSeedInherited = true;
								break;
							}
						}
						if (seed==null)	{			// not found
							seed = selectBestKmer(kmers);
							if (verbose>1)
								System.out.println(CommonUtils.timeElapsed(tic)+
										": Previous seed is not found, pick top k-mer: "+seed.toShortString());
						}
						else
							if (verbose>1)
								System.out.println(CommonUtils.timeElapsed(tic)+
										": Use seed k-mer from previous round: "+primarySeed.getKmerString()+"/"+primarySeed.getKmerRC());
					}
					else{	// no previously saved seed k-mer
						seed = selectBestKmer(kmers);
						if (verbose>1)
							System.out.println(CommonUtils.timeElapsed(tic)+
									": Start with new seed, pick top k-mer: "+seed.toShortString());
						primarySeed = seed;
					}
				}
			}
			else
				seed = selectBestKmer(kmers);
			
			// print out top kmer information
			if (clusterID==0){
				System.out.println("\nTop 5 k-mers");
				System.out.println(Kmer.toShortHeader(k));
				for (int i=0;i<Math.min(5,kmers.size());i++){
					System.out.println(kmers.get(i).toShortString()+((verbose<=1)||(isSeedInherited||kmers.get(i).familyHgp==0)?"":String.format("\t%.1f",kmers.get(i).familyHgp)));
				}
				System.out.println("Seed k-mer:\n"+seed.toShortString()+((verbose<=1)||(isSeedInherited||seed.familyHgp==0)?"":String.format("\t%.1f",seed.familyHgp))+"\n");
			}
			else
				if (verbose>1)
					System.out.println(CommonUtils.timeElapsed(tic)+": Seed k-mer: "+seed.toShortString());
			
			cluster.seedKmer = seed;
			
			System.out.println("Building k-mer cluster "+clusterID+" ...");
			
			/** init kmerSet with seed family of seed kmer, by adding mismatch k-mers, order by #mm */
			ArrayList<Kmer> seedFamily = new ArrayList<Kmer>();
			seedFamily.add(seed);
			kmers.remove(seed);
			seed.setShift(0);
			seed.setAlignString(seed.getKmerString());
			if (config.use_seed_family){
				seedFamily.addAll(getMMKmers(kmers, cluster.seedKmer.getKmerString(), 0));
				KmerGroup kg = config.use_weighted_kmer ? new KmerGroup(seedFamily, 0, seq_weights, posSeqCount, negSeqCount) : new KmerGroup(seedFamily, 0, posSeqCount, negSeqCount);
				cluster.seedKmer.familyHgp = computeHGP(kg.getGroupHitCount(), kg.getGroupNegHitCount());
				if (verbose>1)
					System.out.println(CommonUtils.timeElapsed(tic)+": Seed family hgp = "+cluster.seedKmer.familyHgp);
			}
			if (only_seed_kmer)
				return kmers;
			
			/** align sequences using seedFamily kmer positions */
	    	for (Sequence s : seqList){
	    		for (Kmer km:seedFamily){
					int seed_seq = s.getSeq().indexOf(km.getAlignString());
					if (config.strand_type != 1){	// if not only for single-strand
						if (seed_seq<0){
							s.RC();
							seed_seq = s.getSeq().indexOf(km.getAlignString());
						}
					}
					if (seed_seq<0)
						continue;
					s.pos = -seed_seq;
					break;				// seq is aligned, do not try with weaker k-mers
	    		}
			}
	    	seedFamily = null;
	    	
	    	// build first PWM with some noise
	    	if (verbose>1 && config.pwm_noise!=0)
				System.out.println(CommonUtils.timeElapsed(tic)+ ": PWM noise = " + config.pwm_noise);
			if ( buildPWM(seqList, cluster, config.pwm_noise, tic, true) > -1){    	
				alignSequencesUsingPWM(seqList, cluster);
				if (use_PWM_MM){
					alignSequencesUsingPWMmm(seqList, cluster);
				}
			}
	    	if (config.use_ksm){
		    	cluster.alignedKmers = getAlignedKmers (seqList, seed_range, new ArrayList<Kmer>());
				if (cluster.alignedKmers.isEmpty()){
					kmers.remove(seed);
					quick_restart = true;
					clusterID++;
					for (Sequence s : seqList)
						s.reset();
					continue;
				}
				alignSequencesUsingKSM(seqList, cluster);
	    	}
    	
			/** Iteratively build PWM and align sequences */
			improvePWM (cluster, seqList, seed_range, config.use_ksm, use_PWM_MM);
	    	
			// compare pwm Hgp to primary cluster Hgp, so that the primary cluster will have the best Hgp
			if (cluster.wm!=null){
				double current_cluster_hgp = cluster.pwmThresholdHGP;
				double primary_cluster_hgp = clusters.get(0).pwmThresholdHGP;
				if (config.evaluate_by_ksm){
					current_cluster_hgp = cluster.ksmThreshold.hgp;
					primary_cluster_hgp = clusters.get(0).ksmThreshold.hgp;
				}
				if (config.allow_seed_reset && clusterID!=0 && 
						current_cluster_hgp<primary_cluster_hgp*seedOverrideScoreDifference && 
						!primarySeed_is_immutable){		// this pwm is better, and primary seed is allow to change
					// reset sequences, kmers, to start over with this new bestSeed
					seqs = pos_seq_backup.clone();
					seqsNegList.clear();
					for (String s:neg_seq_backup)
						seqsNegList.add(s);

					seqList.clear();
					for (int i=0;i<seqs.length;i++){
						Sequence s = new Sequence(seqs[i], i);
						seqList.add(s);
					}
					seqList.trimToSize();
					
					kmers.clear();
					for (Kmer km:kmers_in)
						kmers.add(km.clone());
					kmers.trimToSize();
					
			    	clusters.clear();
					clusterID = 0;
					primarySeed = seed;
					if (verbose>1)
						System.out.println("**** Secondary motif is more enriched than primary motif, start over with seed="+seed.getKmerString());
					primarySeed_is_immutable=true;					// marked as "reset", only once, avoid potential infinite loop
					continue;										// start over with new seed
				}
			}
				
			// get Aligned Kmers from PWM alignement. if no PWM, it is the seed_family alignment
	    	if (cluster.wm!=null){
	    		alignSequencesUsingPWM(seqList, cluster);
	    		if (cluster.pwmPosHitCount>cluster.total_aligned_seqs)
	    			cluster.total_aligned_seqs = cluster.pwmPosHitCount;
	    	}
	    	ArrayList<Kmer> alignedKmers = getAlignedKmers (seqList, seed_range, new ArrayList<Kmer>());
	    	
			int shift_remove = 0;
			switch(config.kmer_remove_mode){
			case 0: shift_remove = 0; break;
			case 1: shift_remove = 1; break;
			case 2: shift_remove = k/2; break;
			case 3: shift_remove = k; break;
			case 4: shift_remove = 1; break;		//TODO
			}
			
			for (Kmer km: alignedKmers){	
				int shift = km.getShift();
				if (shift>RC/2)
					shift-=RC;
				if (Math.abs(shift)<shift_remove)	
					kmers.remove(km);
			}	    	
			kmers.remove(seed);
			
			/** mask aligned sequences */
			float k_mask_f=1;
			if (cluster.wm!=null){				
		        if (cluster.pwmGoodQuality){			// if PWM quality is not too bad, mask
			        WeightMatrixScorer scorer = new WeightMatrixScorer(cluster.wm);
			        int left = Math.round(cluster.wm.length()/2); // mask center base
			        int right = left+1;
			        if (!config.k_mask_1base){		// mask whole PWM
			        	left = Math.round(cluster.wm.length()*(1-k_mask_f)/2);
				        right = Math.round(cluster.wm.length()*(1-(1-k_mask_f)/2));
			        }
			        for (Sequence s : seqList){
						if (s.pos==UNALIGNED)
							continue;
						String seq = s.getSeq();
						boolean found = false;
						while(true){				// mask all possible matches
							Pair<Integer, Double> hit = CommonUtils.scanPWM(seq, cluster.wm.length(), scorer);
							if (hit.cdr()<cluster.pwmThreshold)
								break;
							found = true;
							int start = Math.abs(hit.car());// if match on rc strand, the start coordinate is the same, as implemented in WeightMatrixScorer.score()
							if (seq.substring(start+left, start+right).contains("N")){
								found = false;
								break;				// this site has been masked, avoid infinite loop because the masked site still pass PWM threshold
							}
							// replace with N
							seq=seq.substring(0, start+left)
									.concat(CommonUtils.padding(right-left, 'N'))
									.concat(seq.substring(start+right, seq.length()));
						}
						if (found){
							s.setSeq(seq);
							seqs[s.id]=s.getSeqStrand(true);			// Mask the original sequence also
						}
					}
			        // mask also negative sequences
			        for (int i=0; i<seqsNegList.size();i++){
			        	String seq = seqsNegList.get(i);
						boolean found = false;
						while(true){				// mask all possible matches
							Pair<Integer, Double> hit = CommonUtils.scanPWM(seq, cluster.wm.length(), scorer);
							if (hit.cdr()<cluster.pwmThreshold)
								break;
							found = true;
							int start = Math.abs(hit.car());// if match on rc strand, the start coordinate is the same, as implemented in WeightMatrixScorer.score()
							if (seq.substring(start+left, start+right).contains("N")){
								found = false;
								break;				// this site has been masked, avoid infinite loop because the masked site still pass PWM threshold
							}
							// replace with N
							seq=seq.substring(0, start+left)
									.concat(CommonUtils.padding(right-left, 'N'))
									.concat(seq.substring(start+right, seq.length()));
						}
						if (found){
							seqsNegList.set(i, seq);
						}
					}
			        isMasked = true;
		        }
			}		// mask with PWM

			// clone alignedKmers to avoid being modify in the next clustering iteration
			ArrayList<Kmer> clones = new ArrayList<Kmer>();
			for (Kmer km:cluster.alignedKmers)
				clones.add(km.clone());
			clones.trimToSize();
			cluster.alignedKmers = clones;
			
			clusterID++;
			quick_restart = false;
		} // Loop for each cluster
		
		if (topCluster!=-1){		// only generate a few clusters to select optimal K, should not reach here
			primarySeed = null;		// do not record primary seed here
			return null;
		}
		
		// reload non-masked sequences and kmers
		seqs = pos_seq_backup.clone();
		seqsNegList.clear();
		for (String s:neg_seq_backup)
			seqsNegList.add(s);
		for (int i=0;i<seqs.length;i++)
			seqList.set(i, new Sequence(seqs[i], i));
		
		kmers.clear();
		for (Kmer km:kmers_in)
			kmers.add(km.clone());
		kmers.trimToSize();
		
		indexKmerSequences(kmers, seqList);
		
		// refine PWMs using un-masked sequences, and merge similar PWMs
		if (config.refine_pwm){
	    	// re-build PWM with un-masked sequences
			if (!standalone){
				if (verbose>1)
					System.out.println("\nRefining motifs with un-masked sequences ...");
				for (int i=0;i<clusters.size();i++){		// do not need this for primary cluster
					KmerCluster cluster = clusters.get(i);
					if (cluster.wm==null)
						continue;
					if (verbose>1)
						System.out.println(String.format("\n%s: Refining %s hgp=1e-%.1f", CommonUtils.timeElapsed(tic), 
								WeightMatrix.getMaxLetters(cluster.wm), cluster.pwmThresholdHGP));
					for (double dec=0;dec<0.3;dec+=0.1){		// refine using more relax threshold to try more sequences
						while(true){
							HashMap<Integer, PWMHit> hits = findAllPWMHits (seqList, cluster, config.wm_factor-dec); 
							if (buildPWMfromHits(seqList, cluster, hits.values().iterator())<=-1)	// But the selection wm_factor IS NOT CHANGED
								break;
						}
					}
				}
			}

			ArrayList<KmerCluster> badClusters = new ArrayList<KmerCluster>();
			for (KmerCluster c:clusters){
	    		if (c.wm==null || (!c.pwmGoodQuality))
	    			badClusters.add(c);
			}
			clusters.removeAll(badClusters);
			
			// sort secondary clusters
			ArrayList<KmerCluster> secondaryClusters = new ArrayList<KmerCluster>();
			for (int i=1;i<clusters.size();i++){
				secondaryClusters.add(clusters.get(i));
			}
			clusters.removeAll(secondaryClusters);
			Collections.sort(secondaryClusters);
			clusters.addAll(secondaryClusters);
			for (int i=0;i<clusters.size();i++){
				clusters.get(i).clusterId = i;
			}	
			
			if (verbose>1)
				System.out.println("\nMerge overlapping motif clusters ...\n");
			boolean[][] checked = new boolean[clusters.size()][clusters.size()];	// whether a pair has been checked
			mergeOverlapClusters (outName, seqList, seed_range, config.use_ksm, use_PWM_MM, checked, 0);
		}	// refine PWMs
		
		/** post processing */
		// remove clusters if it did not form PWM, or PWM hit is too few
		ArrayList<KmerCluster> badClusters = new ArrayList<KmerCluster>();
		for (KmerCluster c:clusters){
    		if (c.wm==null || (!c.pwmGoodQuality) || c.total_aligned_seqs<seqs.length*config.motif_hit_factor_report)
    			badClusters.add(c);
		}
		// if all of the clusters does not pass, relax the PWM hit count criteria
		if (badClusters.size()==clusters.size()){
			for (int i=0;i<Math.min(5, clusters.size());i++)
				if (clusters.get(i).wm!=null )
					badClusters.remove(clusters.get(i));
		}
		clusters.removeAll(badClusters);
		
		// sort secondary clusters
		ArrayList<KmerCluster> secondaryClusters = new ArrayList<KmerCluster>();
		for (int i=1;i<clusters.size();i++)
			secondaryClusters.add(clusters.get(i));
		clusters.removeAll(secondaryClusters);
		Collections.sort(secondaryClusters);
		clusters.addAll(secondaryClusters);
		
		for (int i=0;i<clusters.size();i++){
			KmerCluster cluster = clusters.get(i);
			/** use all aligned sequences to find expected binding sites, set kmer offset */
	    	// average all the binding positions to decide the expected binding position
			StringBuilder sb = new StringBuilder();
			if (cluster.wm!=null){
				alignSequencesUsingPWM(seqList, cluster);
				if (config.refine_ksm)
					cluster.alignedKmers = getAlignedKmers (seqList, seed_range, new ArrayList<Kmer>());
				updateEngine(cluster.alignedKmers);
				cluster.ksmThreshold = optimizeKsmThreshold("", false);
			}
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
	    	int midPos=k_win/2;
			for (Sequence s : seqList){
				if (s.pos==UNALIGNED)
					continue;
				if (config.print_aligned_seqs)
					sb.append(String.format("%d\t%.1f\t%d\t%s\t%s%s\n", s.id, s.score, s.pos, s.isForward?"F":"R", CommonUtils.padding(-leftmost+s.pos, '.'), s.getSeq()));
				bs[count]=midPos+s.pos;
				count++;
			}
			
//			ArrayList<Integer> nums = new ArrayList<Integer>();
//			for (int j=0;j<count;j++)
//				nums.add((int)bs[j]);
//			Pair<int[], int[]> sorted = StatUtil.sortByOccurences (nums);
//			System.out.println("\nCluster: "+i);
//			System.out.println("Median: "+StatUtil.median(bs));
//			System.out.println("Mean:   "+StatUtil.mean(bs));
//			for (int j=sorted.car().length-1;j>=0;j--)
//				System.out.println(sorted.car()[j]+"\t"+sorted.cdr()[j]);
			
			// median BS position relative to seed k-mer start
			cluster.pos_BS_seed=(int)Math.ceil(StatUtil.median(bs));		
			if (config.print_aligned_seqs)
				CommonUtils.writeFile(outName+"_"+clusterID+"_seqs_aligned.txt", sb.toString());
			sb = null;
//			if (verbose>1 && cluster.wm!=null){
//				int pos_BS_midPWM = cluster.pos_BS_seed - cluster.pos_pwm_seed - cluster.wm.length()/2;
//				System.out.println(String.format("pos_BS_midPWM = %d, variance=%.2f", pos_BS_midPWM, StatUtil.std(bs)));
//			}
	    	
			ArrayList<Kmer> copy = new ArrayList<Kmer>();
			for (Kmer km: cluster.alignedKmers){			// set k-mer offset
				Kmer cp = km.clone();
				copy.add(cp);
				int shift = km.getShift();
				if (shift>RC/2){
					shift-=RC;
					cp.setKmerString(km.getKmerRC());
					cp.setShift(shift);
				}
				cp.setKmerStartOffset(shift-cluster.pos_BS_seed);
			}
			cluster.alignedKmers = copy;				// store all the aligned k-mers
			copy = null;
			
			if (i==0)
				primarySeed = clusters.get(0).seedKmer;
			clusters.get(i).clusterId = i;
			for (Kmer km: clusters.get(i).alignedKmers)
				km.setClusterId(i);
		}		
		
		ArrayList<Kmer> allAlignedKmers = new ArrayList<Kmer>();		// each kmer in diff cluster has been clone, would not overwrite
		for (KmerCluster c:clusters)
			allAlignedKmers.addAll(c.alignedKmers);
		
		Collections.sort(allAlignedKmers, new Comparator<Kmer>(){
		    public int compare(Kmer o1, Kmer o2) {
		    	return o1.compareByHGP(o2);
		    }
		});	
		
		// print PWM spatial distribtution
		printMotifDistanceDistribution(outName);
		
		outputClusters(allAlignedKmers, eventCounts, gemMsg);

		// print the clustered k-mers
		Collections.sort(kmers);
		KmerCluster pCluster = getPrimaryCluster();
		if (pCluster==null){			
			System.out.println("No motif found, exit here! "+CommonUtils.timeElapsed(tic));
			return null;
		}
		double score = pCluster.ksmThreshold==null?0:pCluster.ksmThreshold.score;
		Kmer.printKmers(allAlignedKmers, posSeqCount, negSeqCount, score, outName, false, true, false);
		
		System.out.println("\nFinish KMAC motif discovery, "+CommonUtils.timeElapsed(tic));
		System.out.println(StatUtil.cacheAccessCount);
		System.out.println(StatUtil.getCacheSize());
		return allAlignedKmers;
	}
	
	/** Index k-mers and sequences, remove un-enriched k-mers <br>
	 * 	keep track of the k-mer positions in the sequences so that we can use one to find the other
	 * */
	private void indexKmerSequences(ArrayList<Kmer> kmers, ArrayList<Sequence> seqList){

		/* Initialization, setup sequence list, and update kmers */
		// build the kmer search tree
		AhoCorasick oks = new AhoCorasick();
		for (Kmer km: kmers){
			oks.add(km.getKmerString().getBytes(), km);
	    }
		oks.prepare();
		
		// index kmer->seq, seq->kmer
		HashMap <Kmer, HashSet<Integer>> kmer2seq = new HashMap <Kmer, HashSet<Integer>>();
		for (Sequence s:seqList){
			String seq = s.seq;						// get the sequence from original strand
			s.reset();
			HashSet<Kmer> results = queryTree (seq, oks, config.strand_type==1);
			if (!results.isEmpty()){
				for (Kmer km: results){		
					if (!kmer2seq.containsKey(km)){
						kmer2seq.put(km, new HashSet<Integer>());
					}
					kmer2seq.get(km).add(s.id);

					// This is DIFFERENT from alignSequencesByCoOccurence(ArrayList<Kmer>)
					// farward strand
					ArrayList<Integer> pos = StringUtils.findAllOccurences(seq, km.getKmerString());
					for (int p:pos){
						if (!s.fPos.containsKey(km))
							s.fPos.put(km, new HashSet<Integer>());
						s.fPos.get(km).add(p);
					}
					pos = StringUtils.findAllOccurences(seq, km.getKmerRC());
					for (int p:pos){
						if (!s.fPos.containsKey(km))
							s.fPos.put(km, new HashSet<Integer>());
						s.fPos.get(km).add(p+RC);					// match kmer RC
					}
					if (config.strand_type != 1){
						// reverse strand
						pos = StringUtils.findAllOccurences(s.rc, km.getKmerString());
						for (int p:pos){
							if (!s.rPos.containsKey(km))
								s.rPos.put(km, new HashSet<Integer>());
							s.rPos.get(km).add(p);
						}
						pos = StringUtils.findAllOccurences(s.rc, km.getKmerRC());
						for (int p:pos){
							if (!s.rPos.containsKey(km))
								s.rPos.put(km, new HashSet<Integer>());
							s.rPos.get(km).add(p+RC);					// match kmer RC
						}
					}
				}
			}
		}
		oks = null;
		
		// update kmerCount, and hgp()
		ArrayList<Kmer> unenriched = new ArrayList<Kmer>();
		for (Kmer km:kmers){
			if (kmer2seq.containsKey(km)){
				km.setPosHits(kmer2seq.get(km));
				km.setHgp(computeHGP(km.getPosHitCount(), km.getNegHitCount()));	
				if (km.getHgp()>config.kmer_hgp || km.getPosHitCount()==0)
					unenriched.add(km);
			}
			else{
				unenriched.add(km);
			}
		}
		kmers.removeAll(unenriched);
		for (Sequence s:seqList){
			s.removeAllKmers(unenriched);
		}			
		unenriched = null;
		kmer2seq = null;
	}
	
//	private Kmer selectBestKmer0(ArrayList<Kmer> kmers){
//		ArrayList<Kmer> candidates = new ArrayList<Kmer>();
//		Collections.sort(kmers, new Comparator<Kmer>(){
//		    public int compare(Kmer o1, Kmer o2) {
//		    	return o1.compareByHGP(o2);
//		    }
//		});	
//		ArrayList<Kmer> newList = new ArrayList<Kmer>();
//		for (int i=0;i<kmers.size();i++)
//			newList.add(kmers.get(i));
//		for (Kmer km:kmers){
//			if (!newList.contains(km))		// if km has been selected as member of other kmer family
//				continue;
//			ArrayList<Kmer> family = new ArrayList<Kmer>();
//			family.add(km);
//			family.addAll(getMMKmers(newList, km.getKmerString(), 0));
//			newList.removeAll(family);
//			// compute KmerGroup hgp for the km family
//			KmerGroup kg = use_weighted_kmer ? new KmerGroup(family, 0, seq_weights) : new KmerGroup(family, 0);
//			km.familyHgp = computeHGP(kg.getGroupHitCount(), kg.getGroupNegHitCount());
//			candidates.add(km);
//			if (candidates.size()>3)
//				break;
//		}
//		Collections.sort(candidates, new Comparator<Kmer>(){
//		    public int compare(Kmer o1, Kmer o2) {
//		    	return o1.compareByFamilyHGP(o2);
//		    }
//		});	
//		return candidates.get(0);
//	}
	
	private Kmer selectBestKmer(ArrayList<Kmer> kmers){
		Kmer minHgpKmer = kmers.get(0);
		Kmer maxCountKmer = kmers.get(0);
		for (Kmer km:kmers){
			if (km.getHgp()<minHgpKmer.getHgp())
				minHgpKmer = km;
			if (km.getPosHitCount()>maxCountKmer.getPosHitCount())
				maxCountKmer = km;
		}
		if (minHgpKmer==maxCountKmer)
			return minHgpKmer;		
	
		ArrayList<Kmer> family1 = new ArrayList<Kmer>();
		family1.add(minHgpKmer);
		family1.addAll(getMMKmers(kmers, minHgpKmer.getKmerString(), 0));
		// compute KmerGroup hgp for the km family
		KmerGroup kg = config.use_weighted_kmer ? new KmerGroup(family1, 0, seq_weights, posSeqCount, negSeqCount) : new KmerGroup(family1, 0, posSeqCount, negSeqCount);
		minHgpKmer.familyHgp = computeHGP(kg.getGroupHitCount(), kg.getGroupNegHitCount());
		
		ArrayList<Kmer> family2 = new ArrayList<Kmer>();
		family2.add(maxCountKmer);
		family2.addAll(getMMKmers(kmers, maxCountKmer.getKmerString(), 0));
		// compute KmerGroup hgp for the km family
		kg = config.use_weighted_kmer ? new KmerGroup(family2, 0, seq_weights, posSeqCount, negSeqCount) : new KmerGroup(family2, 0, posSeqCount, negSeqCount);
		maxCountKmer.familyHgp = computeHGP(kg.getGroupHitCount(), kg.getGroupNegHitCount());
		
		if (minHgpKmer.familyHgp<=maxCountKmer.familyHgp)
			return minHgpKmer;
		else
			return maxCountKmer;
	}	
	
	private void alignSequencesUsingPWMmm(ArrayList<Sequence> seqList, KmerCluster cluster){
		String pwmStr = WeightMatrix.getMaxLetters(cluster.wm);
		AhoCorasick oks = new AhoCorasick();
		for (int p=0;p<pwmStr.length();p++){
			char[] pwmBases = pwmStr.toCharArray();
			for (char base:LETTERS){
				pwmBases[p]=base;
				String mm = new String(pwmBases);
				oks.add(mm.getBytes(), mm);
			}								
	    }
		oks.prepare();
		// align sequences using pwmStr with 1 mismatch */
		float[][] matrix = cluster.wm.matrix;
    	for (Sequence s : seqList){
    		if (s.pos!=UNALIGNED)
    			continue;
    		String seq = s.getSeq();
			HashSet<String> strFound = new HashSet<String>();	
			Iterator searcher = oks.search(seq.getBytes());
			while (searcher.hasNext()) {
				SearchResult result = (SearchResult) searcher.next();
				strFound.addAll(result.getOutputs());
			}
			String seq_rc = SequenceUtils.reverseComplement(seq);
			searcher = oks.search(seq_rc.getBytes());
			while (searcher.hasNext()) {
				SearchResult result = (SearchResult) searcher.next();
				strFound.addAll(result.getOutputs());
			}
			if (!strFound.isEmpty()){
				String bestStr = "";
				float bestScore = -100;
				for (String str: strFound){
					float score = 0;
					for (int p=0;p<str.length();p++)
						score += matrix[p][str.charAt(p)];
					if (bestScore<score)
						bestStr = str;
				}
				int pwm_seq = seq.indexOf(bestStr);
				if (pwm_seq<0){
					s.RC();
					pwm_seq = seq_rc.indexOf(bestStr);
					if (pwm_seq<0)
						continue;
				}
				s.pos = cluster.pos_pwm_seed - pwm_seq;		// i.e. seq_seed = pwm_seed-pwm_seq;
				continue;
			}
		}
	}
	
	private void outputClusters(ArrayList<Kmer> allAlignedKmers, int[] eventCounts, String gemMsg){
		// output cluster information, PFM, and PWM
		File f = new File(outName);
		String name = f.getName();
		f=null;
		
		// remove clusters with low hit count
		ArrayList<KmerCluster> toRemove = new ArrayList<KmerCluster>();
		for (int i=0;i<clusters.size();i++){
			KmerCluster c = clusters.get(i);
			double hitRatio = (double)c.pwmPosHitCount / posSeqCount;
			if (i>=10&&hitRatio<config.motif_hit_factor_report || hitRatio<config.motif_hit_factor)
					toRemove.add(c);
			if (config.evaluate_by_ksm && c.ksmThreshold.hgp>config.hgp)
				toRemove.add(c);
			if (!config.evaluate_by_ksm && c.pwmThresholdHGP>config.hgp)
					toRemove.add(c);
		}
		clusters.removeAll(toRemove);
		
		// output kmer alignments
		StringBuilder alignedKmer_sb = new StringBuilder();
		for (KmerCluster c:clusters){
			// print aligned k-mers of this cluster
			ArrayList<Kmer> alignedKmers = c.alignedKmers;
	    	if (!alignedKmers.isEmpty())
	    		alignedKmer_sb.append("Cluster #"+c.clusterId+", n="+alignedKmers.size()+"\n");
	    	// sort for output
			Collections.sort(alignedKmers, new Comparator<Kmer>(){
			    public int compare(Kmer o1, Kmer o2) {
			    	return o1.compareByHGP(o2);
			    }
			});	
	    	int leftmost_km = Integer.MAX_VALUE;
			for (Kmer km: alignedKmers){
				if (km.getKmerStartOffset()<leftmost_km)
					leftmost_km = km.getKmerStartOffset();
			}
			for (Kmer km: alignedKmers){
				alignedKmer_sb.append(km.getKmerStartOffset()+"\t"+CommonUtils.padding(-leftmost_km+km.getKmerStartOffset(), '.')
						+km.getKmerString()+"\t"+km.getPosHitCount()+"\t"+km.getNegHitCount()+"\t"+String.format("%.1f", km.getHgp())+"\t"+km.getAlignString()+"\n");
			}
		}
		CommonUtils.writeFile(outName+"_Alignement_k"+k+".txt", alignedKmer_sb.toString());
		alignedKmer_sb = null;
		
		// output PWM info
		System.out.println();
		StringBuilder pfm_sb = new StringBuilder();		// default TRANSFAC/STAMP format
		StringBuilder pfm_jasper_sb = new StringBuilder();		// JASPAR format
		StringBuilder pfm_meme_sb = new StringBuilder();		// MEME format
		StringBuilder pfm_homer_sb = new StringBuilder();		// HOMER format
		for (KmerCluster c:clusters){
     		WeightMatrix wm = c.wm;
    		System.out.println(String.format("--------------------------------------------------------------\n%s k-mer set #%d, aligned %d k-mers, %d sequences.", name, c.clusterId, c.alignedKmers.size(), c.total_aligned_seqs));
			int pos = c.pos_BS_seed-c.pos_pwm_seed;
    		if (pos>=0)
    			System.out.println(CommonUtils.padding(pos, ' ')+"|\n"+ WeightMatrix.printMatrixLetters(wm));
    		else
    			System.out.println(WeightMatrix.printMatrixLetters(wm));
    		System.out.println(String.format("PWM threshold: %.2f/%.2f, \thit=%d+/%d-, hgp=1e%.1f", c.pwmThreshold, c.wm.getMaxScore(), c.pwmPosHitCount, c.pwmNegHitCount, c.pwmThresholdHGP));
			pfm_sb.append(CommonUtils.makeTRANSFAC (c.pfm, c.pwmPosHitCount, 
					String.format("DE %s_%d_%d_c%d", name, c.clusterId, pos, c.pwmPosHitCount)));
			if (config.outputMEME)
				pfm_meme_sb.append(CommonUtils.makeMEME (c.pfm, c.pwmPosHitCount, 
						String.format("%s_%d_%d_c%d", name, c.clusterId, pos, c.pwmPosHitCount)));
			if (config.outputJASPAR)
				pfm_jasper_sb.append(CommonUtils.makeJASPAR (c.pfm, c.pwmPosHitCount, 
						String.format("%s_%d_%d_c%d", name, c.clusterId, pos, c.pwmPosHitCount)));
			if (config.outputHOMER)				
				pfm_homer_sb.append(CommonUtils.makeHOMER (c.pfm, c.pwmPosHitCount, 
						String.format("%s\t%s_%d_%d_c%d", WeightMatrix.getMaxLetters(c.wm), name, c.clusterId, pos, c.pwmPosHitCount)));
    		if (config.use_ksm && c.ksmThreshold!=null)
    			System.out.println(String.format("KSM threshold: %.2f, \thit=%d+/%d-, hgp=1e%.1f", c.ksmThreshold.score, c.ksmThreshold.posHit, c.ksmThreshold.negHit, c.ksmThreshold.hgp));
			
			// paint motif logo
			c.wm.setNameVerType(name, "#"+c.clusterId, "");
			CommonUtils.printMotifLogo(c.wm, new File(outName+"_"+c.clusterId+"_motif.png"), 75);
			
			WeightMatrix wm_rc = WeightMatrix.reverseComplement(wm);
			wm_rc.setNameVerType(name, "#"+c.clusterId, "rc");
			CommonUtils.printMotifLogo(wm_rc, new File(outName+"_"+c.clusterId+"_motif_rc.png"), 75);
		}
		CommonUtils.writeFile(outName+"_PFM.txt", pfm_sb.toString());
		if (config.outputMEME)
			CommonUtils.writeFile(outName+"_PFM_MEME.txt", pfm_meme_sb.toString());
		if (config.outputJASPAR)
			CommonUtils.writeFile(outName+"_PFM_JASPAR.txt", pfm_jasper_sb.toString());
		if (config.outputHOMER)
			CommonUtils.writeFile(outName+"_PFM_HOMER.txt", pfm_homer_sb.toString());

		// output HTML report
		StringBuffer html = new StringBuffer("<style type='text/css'>/* <![CDATA[ */ table, td{border-color: #600;border-style: solid;} table{border-width: 0 0 1px 1px; border-spacing: 0;border-collapse: collapse;} td{margin: 0;padding: 4px;border-width: 1px 1px 0 0;} /* ]]> */</style>");
		html.append("<script language='javascript' type='text/javascript'><!--\nfunction popitup(url) {	newwindow=window.open(url,'name','height=75,width=400');	if (window.focus) {newwindow.focus()}	return false;}// --></script>");
		html.append("<table><th bgcolor='#A8CFFF' colspan=2><font size='5'>");
		html.append(name).append("</font></th>");
		html.append("<tr><td valign='top'><br>");
		if (!this.standalone && eventCounts!=null){
			html.append("<b>Binding Event Predictions</b>:<p>");
			html.append("<a href='"+name+"_GEM_events.txt'>Significant Events</a>&nbsp;&nbsp;: "+eventCounts[0]);
			html.append("<br><a href='"+name+"_GEM_insignificant.txt'>Insignificant Events</a>: "+eventCounts[1]);
			html.append("<br><a href='"+name+"_GEM_filtered.txt'>Filtered Events</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;: "+eventCounts[2]);
			html.append("<p>Read distribution<br><img src='"+name.substring(0,name.length()-2)+"_All_Read_Distributions.png' width='350'><hr>");
		}
		html.append("<p><b>Motif Discovery Results</b>:<p>");
		html.append("<p>Total positive sequences: "+posSeqCount);
		html.append("<p><ul><li><a href='"+name+".KSM.txt'>Complete KSM (K-mer Set Motif) file.</a>");
		html.append("<li><a href='"+name+"_Alignement_k"+k+".txt'>K-mer alignment file.</a>");
		html.append("<li><a href='"+name+"_PFM.txt'>Motif PFMs</a></ul>");
		html.append("<p><table border=1><th>K-mer</th><th>Cluster</th><th>Offset</th><th>Pos Hit</th><th>Neg Hit</th><th>HGP</th>");
		
    	int leftmost_km = Integer.MAX_VALUE;
    	ArrayList<Kmer> outputs = new ArrayList<Kmer>();
    	for (int i=0;i<Math.min(20, allAlignedKmers.size());i++){
    		Kmer km = allAlignedKmers.get(i);
			if (km.getKmerStartOffset()<leftmost_km)
				leftmost_km = km.getKmerStartOffset();
			outputs.add(km);
		}
		Collections.sort(outputs, new Comparator<Kmer>(){
		    public int compare(Kmer o1, Kmer o2) {
		    		return o1.compareByClusterAndHGP(o2);
		    }
		});		
    	
		for (int i=0;i<outputs.size();i++){
			html.append("<tr><td>");
			html.append("<b><font size='4' face='Courier New'>");
			Kmer km = outputs.get(i);
			char[] kmStr = km.getKmerString().toCharArray();
			html.append(CommonUtils.padding(-leftmost_km+km.getKmerStartOffset(), '-'));
			for (char b:kmStr){
				switch(b){
				case 'A': html.append("<font color='green'>A</font>");break;
				case 'C': html.append("<font color='blue'>C</font>");break;
				case 'G': html.append("<font color='orange'>G</font>");break;
				case 'T': html.append("<font color='red'>T</font>");break;
				}
			}
			html.append("</font></b></td>");
			html.append(String.format("<td>%d</td><td>%d</td><td>%d</td><td>%d</td><td>%.1f</td></tr>", 
					km.getClusterId(), km.getKmerStartOffset(), km.getPosHitCount(), km.getNegHitCount(), km.getHgp()));
		}
		html.append("</table>");
		html.append("</td><td valign='top'><br>");
		html.append("<table border=0 align=center><th>Motif PWM</th><th>Motif spatial distribution (w.r.t. primary PWM)<br>Format: position,motif_occurences</th>");
		for (KmerCluster c:clusters){
    		html.append("<tr><td><img src='"+name+"_"+c.clusterId+"_motif.png"+"'><a href='#' onclick='return popitup(\""+name+"_"+c.clusterId+"_motif_rc.png\")'>rc</a><br>");
    		html.append(String.format("PWM: %.2f/%.2f, hit=%d+/%d-, hgp=1e%.1f<br>", 
    				c.pwmThreshold, c.wm.getMaxScore(), c.pwmPosHitCount, c.pwmNegHitCount, c.pwmThresholdHGP));
//    		html.append(String.format("KSM score: %.2f, \thit=%d+/%d-, hgp=1e%.1f<br><br>", 
//    				c.ksmThreshold.score, c.ksmThreshold.posHit, c.ksmThreshold.negHit, c.ksmThreshold.hgp));
    		String suffix = name+"_Spatial_dist_0_"+c.clusterId;
    		html.append("</td><td><a href='"+suffix+".txt'>"+"<img src='"+suffix+".png"+"' height='150'></a></td></tr>");
		}
		html.append("</table>");
		html.append("</td></tr></table>");
		html.append("<p><p>"+gemMsg);
		CommonUtils.writeFile(outName+"_result.htm", html.toString());
		
//		for (int i=0;i<Math.min(clusters.size(), 20);i++){
//			ArrayList<Kmer> clusterKmers = clusters.get(i).alignedKmers;
//			MotifThreshold t = this.estimateClusterKgsThreshold(clusterKmers);
//			if (t!=null)
//				System.out.println(String.format("%d\t%.2f\t%d\t%d\t%.1f", i, t.score, t.posHit, t.negHit, t.hgp ));
//		}
		
	}
	
	/** 
	 * Merge overlapped motif clusters<br>
	 * Assuming the clusters are sorted as its cluster id (0-based)
	 */
	@SuppressWarnings({ "unchecked", "rawtypes" })
	private void mergeOverlapClusters (String name, ArrayList<Sequence> seqList, int seed_range, 
			boolean use_KSM, boolean use_PWM_MM, boolean[][] checked, int depth){	
		if (depth>10)
			return;
		if (verbose>1)
    		System.out.println("\n"+CommonUtils.timeElapsed(tic)+": Merge overlapping motifs, iteration " + depth);
		boolean isChanged = false;
		int maxClusterId=0;
		for (KmerCluster c:clusters)
			if (c.clusterId > maxClusterId)
				maxClusterId = c.clusterId;
		
		ArrayList[][] hits = new ArrayList[seqs.length][maxClusterId+1];
		for (int j=0;j<clusters.size();j++){
			KmerCluster c = clusters.get(j);
			WeightMatrixScorer scorer = new WeightMatrixScorer(c.wm);
			for (int i=0;i<seqs.length;i++){
				hits[i][c.clusterId]=CommonUtils.getAllPWMHit(seqs[i], c.wm.length(), scorer, c.pwmThreshold);
			}
		}
				
		int seqLen = seqs[0].length();
		for (int m=0;m<clusters.size();m++){
			for (int j=m+1;j<clusters.size();j++){
				if (m>=clusters.size()||j>=clusters.size())			// removing a cluster may change the cluster index/size
					continue;
				KmerCluster cluster_main = clusters.get(m);
				KmerCluster cluster_junior = clusters.get(j);
				if (m!=0 && cluster_main.pwmPosHitCount<cluster_junior.pwmPosHitCount){
					cluster_main = clusters.get(j);
					cluster_junior = clusters.get(m);
				}
				if (checked[cluster_main.clusterId][cluster_junior.clusterId])
					continue;
				
				int range = seqLen - cluster_main.wm.length()/2 - cluster_junior.wm.length()/2 + 4;  // add 2 for rounding error
				int[] same = new int[range*2+1];
				int[] diff = new int[range*2+1];
				for (int i=0;i<seqs.length;i++){
					ArrayList<Integer> hitm = hits[i][cluster_main.clusterId];
					ArrayList<Integer> hitj = hits[i][cluster_junior.clusterId];
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
				
				int minOverlap = Math.min(cluster_main.wm.length(), cluster_junior.wm.length())/2;
				int minHitCount = Math.min(cluster_main.pwmPosHitCount, cluster_junior.pwmPosHitCount);
				int maxCount = 0;
				int maxDist = 0;
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
				
				if (maxCount>=minHitCount*0.3){				// if there is large enough overlap, try to merge 2 clusters
					if (verbose>1)
			    		System.out.println(String.format("\n%s: Trying to merge %s(#%d, %.1f) and %s(#%d, %.1f), dist=%d%s ... ", 
			    				CommonUtils.timeElapsed(tic), WeightMatrix.getMaxLetters(cluster_main.wm), cluster_main.clusterId, cluster_main.pwmThresholdHGP,
		    				WeightMatrix.getMaxLetters(cluster_junior.wm), cluster_junior.clusterId, cluster_junior.pwmThresholdHGP, maxDist, isRC?"rc":""));
					
					KmerCluster newCluster = cluster_main.clone();
					alignSequencesUsingPWM(seqList, newCluster);	// align with first PWM
					
					// align additional sequences matching second PWM
					WeightMatrix wm = isRC?WeightMatrix.reverseComplement(cluster_junior.wm):cluster_junior.wm;
			        WeightMatrixScorer scorer = new WeightMatrixScorer(wm);		
					int count_pwm_aligned=0;
					for (Sequence s:seqList){
			    	  String seq = s.getSeq();			// PWM to scan unaligned sequence, and align if pass threshold
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
			          if (maxSeqScore >= cluster_junior.pwmThreshold){
						if (maxScoringStrand =='-'){
							maxScoringShift = k_win-maxScoringShift-wm.length();
							s.RC();
							// i.e.  (seq.length()-1)-maxScoringShift-(wm.length()-1);
						}
						s.pos = cluster_main.pos_pwm_seed-(maxScoringShift+cluster_junior.wm.length()/2-maxDist-cluster_main.wm.length()/2);
						count_pwm_aligned ++;
			          }
			          else
			        	  s.pos = UNALIGNED;
			        }	// each sequence
					if (verbose>1)
			    		System.out.println(CommonUtils.timeElapsed(tic)+": PWM "+WeightMatrix.getMaxLetters(wm)+" align additional "+count_pwm_aligned+" sequences.");
					
			    	// build PWM
					buildPWM(seqList, newCluster, 0, tic, false);
					alignSequencesUsingPWM(seqList, newCluster);
					
					// Iteratively improve PWM and sequences alignment
					improvePWM (newCluster, seqList, seed_range, false, false);
					
					// merge if the new PWM is more enriched
					if ((newCluster.pwmThresholdHGP<cluster_main.pwmThresholdHGP && newCluster.pwmPosHitCount>=cluster_main.pwmPosHitCount)){		
						if (newCluster.pwmPosHitCount>newCluster.total_aligned_seqs)
							newCluster.total_aligned_seqs = newCluster.pwmPosHitCount;
						clusters.set(m, newCluster);
						isChanged = true;
						for (int d=0;d<checked.length;d++){
							checked[d][cluster_main.clusterId]=false;
							checked[cluster_main.clusterId][d]=false;
						}
						cluster_main = newCluster;
						if (verbose>1)
				    		System.out.println(CommonUtils.timeElapsed(tic)+": cluster #"+cluster_main.clusterId+" and cluster #"
				    				+cluster_junior.clusterId+" merge to new PWM "+WeightMatrix.getMaxLetters(newCluster.wm));	
					}
					else{	
						checked[cluster_main.clusterId][cluster_junior.clusterId]=true;
						if (verbose>1)
				    		System.out.println(String.format("%s: Merged PWM is not more enriched, do not merge", CommonUtils.timeElapsed(tic)));
					}
					// No matter successful merging or not, try the other PWM after removing the overlap hits
					if (verbose>1)
			    		System.out.println(String.format("%s: Testing the remaining cluster #%d.", 
			    			CommonUtils.timeElapsed(tic), cluster_junior.clusterId));						
					alignSequencesUsingPWM(seqList, cluster_main);
					ArrayList<Sequence> seqList_j = new ArrayList<Sequence>();
					for (Sequence s:seqList)
						if (s.pos == UNALIGNED)
							seqList_j.add(s);
					alignSequencesUsingPWM(seqList_j, cluster_junior);
					int aligned_seqs_count=0;
					for (Sequence s:seqList_j){
				    	  if (s.pos!=UNALIGNED)
				    		  aligned_seqs_count++;
					}
					// if aligned seq count is less than threshold, or if it contains less than half of total hit of the motif (i.e. majority of hits are still overlapped), remove it
					if (aligned_seqs_count<seqs.length*config.motif_hit_factor || aligned_seqs_count<cluster_junior.pwmPosHitCount/2){
						if (verbose>1)
				    		System.out.println(String.format("%s: Exclusive hits(%d) by cluster #%d is too few, remove it.", 
				    			CommonUtils.timeElapsed(tic), aligned_seqs_count, cluster_junior.clusterId));
						clusters.remove(cluster_junior);
					}
					else{
						WeightMatrix old_j = cluster_junior.wm;
						cluster_junior.wm = null;			// reset here, to get a new PWM
						int alignedSeqCount = 0;
						buildPWM(seqList_j, cluster_junior, 0, tic, false);
						if (cluster_junior.wm != null){		// got a new PWM
							alignSequencesUsingPWM(seqList_j, cluster_junior);
							improvePWM (cluster_junior, seqList_j, seed_range, false, false);
							alignedSeqCount = alignSequencesUsingPWM(seqList_j, cluster_junior);
						}
						if (alignedSeqCount<seqs.length*config.motif_hit_factor){
							if (verbose>1)
					    		System.out.println(String.format("%s: new PWM hit %d is too few, remove cluster #%d.", 
					    			CommonUtils.timeElapsed(tic), alignedSeqCount, cluster_junior.clusterId));
							clusters.remove(cluster_junior);
						}
						else{
							// update the # aligned_seqs using the number from all sequences
							if (cluster_junior.pwmPosHitCount>cluster_junior.total_aligned_seqs)
								cluster_junior.total_aligned_seqs = cluster_junior.pwmPosHitCount;
							if (old_j.isSame(cluster_junior.wm)){
								if (verbose>1)
						    		System.out.println(String.format("%s: Cluster #%d is unchanged.", 
						    			CommonUtils.timeElapsed(tic), cluster_junior.clusterId));
								if (!isChanged)		// both clusters are not changed
									checked[cluster_main.clusterId][cluster_junior.clusterId] = true;
								continue;
							}
							isChanged = true;
							if (verbose>1)
					    		System.out.println(String.format("%s: new PWM has sufficient hit %d, keep  cluster #%d.", 
					    			CommonUtils.timeElapsed(tic), alignedSeqCount, cluster_junior.clusterId));
							for (int d=0;d<checked.length;d++){
								checked[d][cluster_junior.clusterId]=false;
								checked[cluster_junior.clusterId][d]=false;
							}
						}
					}
				}		
				else{					// if overlap is not big enough
					checked[cluster_main.clusterId][cluster_junior.clusterId] = true;
				}
			}
		}
		
		if (isChanged)		// if merged, continue to merge, otherwise return
			mergeOverlapClusters (name, seqList, seed_range, use_KSM, use_PWM_MM, checked, depth+1);
	}
	
	/** Iteratively build PWM and align sequences <br>
	 *  The seqList should be aligned so that a new PWM can be built. */
	private void improvePWM (KmerCluster cluster, ArrayList<Sequence> seqList, int seed_range, boolean use_KSM, boolean use_PWM_MM){	
    	while(true){	
			int leftIdx = buildPWM(seqList, cluster, 0, tic, true);	
			if (leftIdx > -1){				    	
				alignSequencesUsingPWM(seqList, cluster);
				if (use_PWM_MM){
					alignSequencesUsingPWMmm(seqList, cluster);
				}
				if (use_KSM){
					cluster.alignedKmers = getAlignedKmers (seqList, seed_range, new ArrayList<Kmer>());
					if (cluster.alignedKmers.isEmpty())
						break;
					alignSequencesUsingKSM(seqList, cluster);
				}
			}
			else
				break;
    	}  // Iteratively improve PWM and sequences alignment
	}

	
	/**
	 * Compute the distance distributions between primary and secondary motifs
	 */
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public void printMotifDistanceDistribution (String name){		
		System.out.println("\nCompute motif distance distribution ...");
		
		ArrayList[][] hits = new ArrayList[seqs.length][clusters.size()];
		for (int j=0;j<clusters.size();j++){
			KmerCluster c = clusters.get(j);
			WeightMatrixScorer scorer = new WeightMatrixScorer(c.wm);
			for (int i=0;i<seqs.length;i++){
				hits[i][j]=CommonUtils.getAllForwardPWMHit(seqs[i], c.wm.length(), scorer, c.pwmThreshold);
			}
		}
		int seqLen = seqs[0].length();
		for (int m=0;m<1;m++){
			for (int j=0;j<clusters.size();j++){
				int range = seqLen - clusters.get(m).wm.length()/2 - clusters.get(j).wm.length()/2;
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
							boolean isForward_m = true;
							if (pm<0){
								pm=-pm;
								isForward_m = false;
							}
							for (int b=a;b<hitm.size();b++){
								int pj = hitm.get(b);
								boolean isForward_j = true;
								if (pj<0){
									pj=-pj;
									isForward_j = false;
								}
								int offset = pj-pm;
								if (Math.abs(offset)>range)
									continue;
								if(isForward_m){
									offset += range;		// shift to get array idx
									if (isForward_j)
										same[offset]++;
									else
										diff[offset]++;
								}
								else{
									offset = -offset;
									offset += range;		// shift to get array idx
									if (isForward_j)
										diff[offset]++;
									else
										same[offset]++;
								}
							}
						}
					}	// self
					else{
						for (int pm:hitm){
							boolean isForward_m = true;
							if (pm<0){
								pm=-pm;
								isForward_m = false;
							}
							for (int pj:hitj){
								boolean isForward_j = true;
								if (pj<0){
									pj=-pj;
									isForward_j = false;
								}
								int offset = pj-pm;
								if (Math.abs(offset)>range)
									continue;
								if(isForward_m){
									offset += range;		// shift to get array idx
									if (isForward_j)
										same[offset]++;
									else
										diff[offset]++;
								}
								else{
									offset = -offset;
									offset += range;		// shift to get array idx
									if (isForward_j)
										diff[offset]++;
									else
										same[offset]++;
								}
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
				String fileSuffix = name+"_Spatial_dist_"+clusters.get(m).clusterId+"_"+clusters.get(j).clusterId;
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
	
	private HashMap<Integer, PWMHit> findAllPWMHits(ArrayList<Sequence> seqList, KmerCluster cluster, double pwm_factor){
		
		WeightMatrix wm = cluster.wm;
	    WeightMatrixScorer scorer = new WeightMatrixScorer(wm);
	    HashMap<Integer, PWMHit> seq2hits = new HashMap<Integer, PWMHit>();
	    
		for (int i=0;i<seqList.size();i++){
		  Sequence s = seqList.get(i);
		  String seq = s.getSeqStrand(true);			// PWM to scan all sequence
		      	  
	      WeightMatrixScoreProfile profiler = scorer.execute(seq);
	      PwmMatch match = profiler.getBestMatch(config.strand_type);
	      
	      if (match.score >= wm.getMaxScore()*pwm_factor){
	    	  PWMHit hit = new PWMHit();
	    	  hit.clusterId = cluster.clusterId;
	    	  hit.wm = wm;
	    	  hit.seqId = i;
	    	  hit.start = match.shift;					// start is cordinate on + strand
	    	  hit.end = match.shift+wm.length()-1;		// end inclusive
			  hit.isForward = match.strand =='+';
			  String ss = seq.substring(match.shift, match.shift+wm.length());
			  if (match.strand =='-')
				  ss= SequenceUtils.reverseComplement(ss);
			  hit.str = ss;
			  
			  hit.weight = 1.0;
			  if (config.use_strength_weight)
				  hit.weight *= seq_weights[hit.seqId];
			  if (hit.clusterId==0 && config.use_pos_weight)
				  hit.weight *= profile[(hit.end+hit.start)/2];
			  
			  seq2hits.put(i, hit);
	      }
	    }	// each sequence
		return seq2hits;
	}
	/** refine PWMs by EM like iterations<br>
	 * given the un-masked positive sequences, and multiple PWMs<br>
	 * find the best PWMs to cluster the sequences
	 */
	public void refinePWMs(ArrayList<Sequence> seqList, String[] pos_seq_backup, String[] neg_seq_backup, 
			ArrayList<Kmer> kmers_in, int seed_range){
	
		ArrayList<KmerCluster> noPWM = new ArrayList<KmerCluster>();
		for (KmerCluster c:clusters){
			if (c.wm==null){
				noPWM.add(c);
				c.pi = c.pwmPosHitCount;
			}
		}
		clusters.removeAll(noPWM);
		if (clusters.size()<=1)
			return;
		
		File f = new File(outName);
		String name = f.getName();
		
		/** Resolve multiple-PWM-match conflict */
		for (int iter=0; iter<100; iter++){
			
			// record cluster PWM hgp for comparison
			double[] hgps = new double[clusters.size()];
			if (verbose>1)
				System.out.println("------------------------------------------------\n"+CommonUtils.timeElapsed(tic)+
						": Iteration #"+iter);
			
			for (int i=0; i<clusters.size(); i++){
				KmerCluster c = clusters.get(i);
				hgps[i] = c.pwmThresholdHGP;
				c.clusterId = i;
				c.seq2hits = findAllPWMHits(seqList, c, config.wm_factor);
				// paint motif logo
				c.wm.setNameVerType(name+"_i"+iter, "#"+c.clusterId, "");
				CommonUtils.printMotifLogo(c.wm, new File(outName+"_i"+iter+"_"+c.clusterId+"_motif.png"), 75);
			}	
			if (verbose>1){
				StringBuilder sb = new StringBuilder(CommonUtils.timeElapsed(tic)+": Cluster PWMs\t");
				for (KmerCluster c:clusters)
					sb.append(WeightMatrix.getMaxLetters(c.wm)).append(" ");
				sb.append("\n").append(CommonUtils.timeElapsed(tic)+": Cluster hgps\t"+CommonUtils.arrayToString(hgps));
				System.out.println(sb.toString());
			}
			int conflicts = resolveConflictHits(seqList);
			if (verbose>1)
				System.out.println(CommonUtils.timeElapsed(tic)+": "+conflicts+" sequences share hits by multiple PWMs");
			if (conflicts==0)
				break;
				
			buildAllClusterPWMfromHits(seqList);
			
			if (clusters.size()<=1)
				break;
			if (clusters.size()!=hgps.length)		// some PWM/cluster may be removed
				continue;
			
			boolean noChange = true;
			for (int i=0; i<clusters.size(); i++){
				if (hgps[i] != clusters.get(i).pwmThresholdHGP){
					noChange = false;
					break;
				}
			}
			if (noChange)
				break;
		}
	}
	
	public int resolveConflictHitsBG(ArrayList<Sequence> seqList){
		double alpha = seqs.length * config.motif_hit_factor_report;

		int conflict = 0;
		// set up mixing components: PWMs
		int wmLen[]=new int[clusters.size()];
		double totalHitSum=0;
		for (int j=0; j<clusters.size(); j++){
			KmerCluster cluster = clusters.get(j);
			wmLen[j]=cluster.wm.length();
			cluster.pi = cluster.pwmPosHitCount;
			totalHitSum += cluster.pwmPosHitCount;
			cluster.clusterId = j;
		}
		
		double bgHitCount = 0;
		// set up data points: PWM hit groups (overlapping hits are grouped as one data point)
		ArrayList<ArrayList<PWMHit>> data = new ArrayList<ArrayList<PWMHit>>();
		for (int i=0;i<seqList.size();i++){
			// get all hits and sort by starts and ends
			ArrayList<PWMHit> hits = new ArrayList<PWMHit>();
			for (int j=0; j<clusters.size(); j++){
				if (clusters.get(j).seq2hits.containsKey(i))
					hits.add(clusters.get(j).seq2hits.get(i));
			}
			if (hits.size()<1){
				bgHitCount += k_win;
				continue;
			}
			Collections.sort(hits, new Comparator<PWMHit>(){
			    public int compare(PWMHit o1, PWMHit o2) {
			    	return o1.compareByPosition(o2);
			    }
			});	
			
			// find independent hit groups
			HashSet<Integer> allOverlapHits = new HashSet<Integer>();
			int hitBases=0;
			for (int j=0;j<hits.size();j++){
				if (allOverlapHits.contains(j))	// if this hit has overlop with previous hit, it is used
					continue;
				PWMHit hit = hits.get(j);
				int start = hit.start;
				int end = hit.end;
				HashSet<Integer> overlapHits = new HashSet<Integer>();
				overlapHits.add(j);
				for (int k=j+1;k<hits.size();k++){
					PWMHit hit2 = hits.get(k);
					if (hit2.overlaps(start, end)){
						overlapHits.add(k);
						if (hit2.end>end)
							end = hit2.end;
					}
				}
				hitBases += end-start+1;
				
				allOverlapHits.addAll(overlapHits);

				ArrayList<PWMHit> hitGroup = new ArrayList<PWMHit>();	// group of overlapping hits
//				if (overlapHits.size()>1)
//					j=j;
				for (int h:overlapHits){
					PWMHit hh = hits.get(h);
					hitGroup.add(hh);
				}
				data.add(hitGroup);
				if (hitGroup.size()>1)
					conflict++;
			}
			bgHitCount += k_win-hitBases;
		}
		
		double bg_pi = bgHitCount/(seqList.size()*k_win);
		for (KmerCluster cluster : clusters){
			cluster.pi = cluster.pi/totalHitSum*(1-bg_pi);
		}
		bgHitCount = totalHitSum/(1-bg_pi)*bg_pi;			// scale BG hit count
		double[] bg_resp = new double[data.size()];

		if (verbose>1)
			System.out.println(String.format("%s: sparseness penalty = %.1f\t#HitGroup=%d", CommonUtils.timeElapsed(tic), alpha, data.size()));
		
		// init responsibilities
		for (int i=0;i<data.size();i++){
			ArrayList<PWMHit> hitGroup = data.get(i);
//			if (hitGroup.size()>1)
//				i=i;
			int maxPwmLen = 0;
			for (PWMHit h: hitGroup){
				int w = wmLen[h.clusterId];
				if (maxPwmLen<w)
					maxPwmLen=w;
			}
			for (PWMHit h: hitGroup){
				float[][]pfm = clusters.get(h.clusterId).pfm;
				double prob = 1;					// this need to be updated for each PWM update
				String seq = h.str;
				for (int p=0;p<seq.length();p++)
					prob *= pfm[p][seq.charAt(p)];
				int N = maxPwmLen - seq.length();
				if (N!=0)
					prob *= Math.pow(0.25,N);
				h.responsibility = prob*clusters.get(h.clusterId).pi;
			}
			bg_resp[i] = Math.pow(0.25,maxPwmLen) * bg_pi;
		}
		
		double logPosterior = Double.NEGATIVE_INFINITY;
		double currAlpha = 0;
		boolean eliminated = false;
		HashSet<Integer> zeroComp = new HashSet<Integer>();
		for(int iter=0;iter<200;iter++){
			// E-STEP for each hit group
			for (int i=0;i<data.size();i++){
				ArrayList<PWMHit> hitGroup = data.get(i);
				if (hitGroup.size()==1){					// only 1 hit, the cluster takes all responsibility
					for (PWMHit h: hitGroup){
						if (zeroComp.contains(h.clusterId))
							h.responsibility = 0;
						else{
							h.responsibility = h.responsibility/(h.responsibility+bg_resp[i]);
							bg_resp[i] = 1-h.responsibility;
						}
					}
					continue;
				}			
				// normalize data likelihood to get responsibility
				double totalResp=0;
				for (PWMHit h: hitGroup){
					totalResp += h.responsibility;
				}
				totalResp += bg_resp[i];
				for (PWMHit h: hitGroup){			
					h.responsibility = h.responsibility/totalResp;
				}
				bg_resp[i] = bg_resp[i] / totalResp;
			}
			
			// M Step
			StringBuilder sb = new StringBuilder("Cluster hit\t");
			for (int j=0; j<clusters.size(); j++){
				KmerCluster cluster = clusters.get(j);
				if(zeroComp.contains(j)){
					sb.append(String.format("%.1f\t", cluster.pi));
					continue;
				}
				float[][] pfm = new float[cluster.pfm.length][MAXLETTERVAL];
				for (int p=0;p<pfm.length;p++){
					for (char base:LETTERS)			// 0 count can cause log(0), set pseudo-count 0.375 to every pos, every base
						pfm[p][base]=0.375f; 		//http://www.ncbi.nlm.nih.gov.libproxy.mit.edu/pmc/articles/PMC2490743/
				} 

				double hitSum = 0;
				for (int i=0;i<seqList.size();i++){				// get all the hits from seqList, equivalent to hit groups list
					if (cluster.seq2hits.containsKey(i)){
						PWMHit hit = clusters.get(j).seq2hits.get(i);
						hitSum += hit.responsibility;
						for (int p=0;p<pfm.length;p++){
			    			char base = hit.str.charAt(p);
			    			pfm[p][base] += hit.responsibility*hit.weight;
			    		}					
					}
				}
				sb.append(String.format("%.1f\t", hitSum));
				cluster.pi = hitSum;
				
				for (int p=0;p<pfm.length;p++){
					float sum = 0;
					for (char base : LETTERS)
						sum += pfm[p][base];
					for (char base : LETTERS)
						pfm[p][base] /= sum;
				}
				cluster.pfm = pfm;
			}
			
			bg_pi=bgHitCount;
			sb.append(String.format("%.1f\t", bg_pi));
			
			// set alpha value, compute pi and BG pi
			if (iter<5)
				currAlpha = 0;
			else if (iter==5)
				currAlpha = 3;
			else if (iter>5 && !eliminated)	{			// if just eliminated, stay on same alpha
				currAlpha *= 3;
			}
			currAlpha = Math.min(currAlpha,alpha);
			
			totalHitSum = 0;
			int prePenaltyHitSum = 0;
			eliminated = false;
			for (KmerCluster cluster : clusters){
				if (zeroComp.contains(cluster.clusterId))
					continue;
				prePenaltyHitSum += cluster.pi;
				cluster.pi -= currAlpha;
				if (cluster.pi<=0){
					cluster.pi = 0;
					eliminated = true;
					zeroComp.add(cluster.clusterId);
				}
				totalHitSum += cluster.pi;
			}
			bg_pi = bg_pi/(bg_pi+prePenaltyHitSum);			
			for (KmerCluster cluster : clusters){
				cluster.pi = cluster.pi/totalHitSum*(1-bg_pi);
			}
			if (verbose>1)
				System.out.print(CommonUtils.timeElapsed(tic)+": i"+iter+" "+sb.toString());
			
			// Compute log likelihood
			double ll = 0;
			for (int i=0;i<data.size();i++){
				ArrayList<PWMHit> hitGroup = data.get(i);
				
				int maxPwmLen = 0;
				for (PWMHit h: hitGroup){
					if(zeroComp.contains(h.clusterId))				// skip the eliminated PWM components
						continue;
					int w = wmLen[h.clusterId];
					if (maxPwmLen<w)
						maxPwmLen=w;
				}
				if (maxPwmLen==0)
					continue;
				double resp_sum=0;
				for (PWMHit h: hitGroup){
					if(zeroComp.contains(h.clusterId))	{			// skip the eliminated PWM components
						h.responsibility = 0;
						continue;
					}
					float[][]pfm = clusters.get(h.clusterId).pfm;
					double prob = 1;								// this need to be updated for each PWM update
					String seq = h.str;
					for (int p=0;p<seq.length();p++)
						prob *= pfm[p][seq.charAt(p)];
					int N = maxPwmLen - seq.length();
					if (N!=0)
						prob *= Math.pow(0.25,N);		
					h.responsibility = prob*clusters.get(h.clusterId).pi;
					resp_sum += h.responsibility;
				}
				bg_resp[i] = Math.pow(0.25,maxPwmLen)*bg_pi;		// BG
				resp_sum += bg_resp[i];
				if (resp_sum!=0)
					ll += Math.log(resp_sum);
			}
			
            // log prior
            double LP=0;
            for(KmerCluster c: clusters)
                if (!zeroComp.contains(c.clusterId))
                    LP+=(-currAlpha)*Math.log(c.pi);		// sparse prior
            LP+=(-currAlpha)*Math.log(bg_pi);				// BG
            double lap = ll+LP;
            
			if (verbose>1)
				System.out.println(String.format("Alpha=%.1f\tLAP diff=%.5f", currAlpha, Math.abs(lap-logPosterior)));
			if (Math.abs(lap-logPosterior)>0.0001){
				logPosterior = lap;
			}
			else{		// normalize responsibility before exit the loop
				if (currAlpha<alpha)
					continue;
				for (int i=0;i<data.size();i++){
					ArrayList<PWMHit> hitGroup = data.get(i);
					if (hitGroup.size()==1){					// only 1 hit, the cluster takes all responsibility
						for (PWMHit h: hitGroup){
							if (zeroComp.contains(h.clusterId))
								h.responsibility = 0;
							else{
								h.responsibility = h.responsibility/(h.responsibility+bg_resp[i]);
								bg_resp[i] = 1-h.responsibility;
							}
						}
						continue;
					}			
					// normalize data likelihood to get responsibility
					double totalResp=0;
					for (PWMHit h: hitGroup){
						totalResp += h.responsibility;
					}
					totalResp += bg_resp[i];
					for (PWMHit h: hitGroup){			
						h.responsibility = h.responsibility/totalResp;
					}
					bg_resp[i] = bg_resp[i] / totalResp;
				}
				break;
			}
		}
		
		return conflict;
	}
	public int resolveConflictHits(ArrayList<Sequence> seqList){
		double alpha = seqs.length * config.motif_hit_factor_report;

		int conflict = 0;
		// set up mixing components: PWMs
		int wmLen[]=new int[clusters.size()];
		double totalHitSum=0;
		for (int j=0; j<clusters.size(); j++){
			KmerCluster cluster = clusters.get(j);
			wmLen[j]=cluster.wm.length();
			cluster.pi = cluster.pwmPosHitCount;
			totalHitSum += cluster.pwmPosHitCount;
			cluster.clusterId = j;
		}
		for (KmerCluster cluster : clusters){
			cluster.pi /= totalHitSum;
		}
		
		int bgBases = 0;
		// set up data points: PWM hit groups (overlapping hits are grouped as one data point)
		ArrayList<ArrayList<PWMHit>> data = new ArrayList<ArrayList<PWMHit>>();
		for (int i=0;i<seqList.size();i++){
			// get all hits and sort by starts and ends
			ArrayList<PWMHit> hits = new ArrayList<PWMHit>();
			for (int j=0; j<clusters.size(); j++){
				if (clusters.get(j).seq2hits.containsKey(i))
					hits.add(clusters.get(j).seq2hits.get(i));
			}
			if (hits.size()<1)
				continue;
			Collections.sort(hits, new Comparator<PWMHit>(){
			    public int compare(PWMHit o1, PWMHit o2) {
			    	return o1.compareByPosition(o2);
			    }
			});	
			
			// find independent hit groups
			HashSet<Integer> allOverlapHits = new HashSet<Integer>();
			for (int j=0;j<hits.size();j++){
				if (allOverlapHits.contains(j))	// if this hit has overlop with previous hit, it is used
					continue;
				PWMHit hit = hits.get(j);
				int start = hit.start;
				int end = hit.end;
				HashSet<Integer> overlapHits = new HashSet<Integer>();
				overlapHits.add(j);
				for (int k=j+1;k<hits.size();k++){
					PWMHit hit2 = hits.get(k);
					if (hit2.overlaps(start, end)){
						overlapHits.add(k);
						if (hit2.end>end)
							end = hit2.end;
					}
				}
				allOverlapHits.addAll(overlapHits);

				ArrayList<PWMHit> hitGroup = new ArrayList<PWMHit>();	// group of overlapping hits
//				if (overlapHits.size()>1)
//					j=j;
				for (int h:overlapHits){
					PWMHit hh = hits.get(h);
					hitGroup.add(hh);
				}
				data.add(hitGroup);
				if (hitGroup.size()>1)
					conflict++;
			}
		}
		if (verbose>1)
			System.out.println(String.format("%s: sparseness penalty = %.1f\t#HitGroup=%d", CommonUtils.timeElapsed(tic), alpha, data.size()));
		
		// init responsibilities
		for (int i=0;i<data.size();i++){
			ArrayList<PWMHit> hitGroup = data.get(i);
//			if (hitGroup.size()>1)
//				i=i;
			int maxPwmLen = 0;
			for (PWMHit h: hitGroup){
				int w = wmLen[h.clusterId];
				if (maxPwmLen<w)
					maxPwmLen=w;
			}
			for (PWMHit h: hitGroup){
				float[][]pfm = clusters.get(h.clusterId).pfm;
				double prob = 1;					// this need to be updated for each PWM update
				String seq = h.str;
				for (int p=0;p<seq.length();p++)
					prob *= pfm[p][seq.charAt(p)];
				int N = maxPwmLen - seq.length();
				if (N!=0)
					prob *= Math.pow(0.25,N);
				h.responsibility = prob*clusters.get(h.clusterId).pi;
			}
		}
		
		double logPosterior = Double.NEGATIVE_INFINITY;
		double currAlpha = 0;
		boolean eliminated = false;
		HashSet<Integer> zeroComp = new HashSet<Integer>();
		for(int iter=0;iter<200;iter++){
			// E-STEP for each hit group
			for (int i=0;i<data.size();i++){
				ArrayList<PWMHit> hitGroup = data.get(i);
				if (hitGroup.size()==1){					// only 1 hit, the cluster takes all responsibility
					for (PWMHit hit: hitGroup){
						if (zeroComp.contains(hit.clusterId))
							hit.responsibility = 0;
						else
							hit.responsibility = 1;
					}
					continue;
				}			
				// normalize data likelihood to get responsibility
				double totalResp=0;
				for (PWMHit c: hitGroup){
					totalResp += c.responsibility;
				}
				for (PWMHit c: hitGroup){			
					c.responsibility = c.responsibility/totalResp;
				}
			}
			
			// M Step
			StringBuilder sb = new StringBuilder("Cluster hit\t");
			totalHitSum = 0;
			double minHitSum = Double.MAX_VALUE;
			for (int j=0; j<clusters.size(); j++){
				KmerCluster cluster = clusters.get(j);
				if(zeroComp.contains(j)){
					sb.append(String.format("%.1f\t", cluster.pi));
					continue;
				}
				float[][] pfm = new float[cluster.pfm.length][MAXLETTERVAL];
				for (int p=0;p<pfm.length;p++){
					for (char base:LETTERS)			// 0 count can cause log(0), set pseudo-count 0.375 to every pos, every base
						pfm[p][base]=0.375f; 		//http://www.ncbi.nlm.nih.gov.libproxy.mit.edu/pmc/articles/PMC2490743/
				} 

				double hitSum = 0;
				for (int i=0;i<seqList.size();i++){				// get all the hits from seqList, equivalent to hit groups list
					if (cluster.seq2hits.containsKey(i)){
						PWMHit hit = clusters.get(j).seq2hits.get(i);
						hitSum += hit.responsibility;
						for (int p=0;p<pfm.length;p++){
			    			char base = hit.str.charAt(p);
			    			pfm[p][base] += hit.responsibility;
			    		}					
					}
				}
				if (minHitSum>hitSum && hitSum>10)				// some times the hitSum can be 0.0 or 1.0
					minHitSum = hitSum;
				sb.append(String.format("%.1f\t", hitSum));
				cluster.pi = hitSum;
				
				for (int p=0;p<pfm.length;p++){
					float sum = 0;
					for (char base : LETTERS)
						sum += pfm[p][base];
					for (char base : LETTERS)
						pfm[p][base] /= sum;
				}
				cluster.pfm = pfm;
			}

			if (iter<5)
				currAlpha = 0;
			else if (iter==5)
				currAlpha = 3;
			else if (iter>5 && !eliminated)	{			// if just eliminated, stay on same alpha
				currAlpha *= 3;
			}
			currAlpha = Math.min(currAlpha,alpha);
			
			eliminated = false;
			for (KmerCluster cluster : clusters){
				if (zeroComp.contains(cluster.clusterId))
					continue;
				cluster.pi -= currAlpha;
				if (cluster.pi<=0){
					cluster.pi = 0;
					eliminated = true;
					zeroComp.add(cluster.clusterId);
				}
				totalHitSum += cluster.pi;
			}
			for (KmerCluster cluster : clusters){
				cluster.pi /= totalHitSum;
			}
			if (verbose>1)
				System.out.print(CommonUtils.timeElapsed(tic)+": i"+iter+" "+sb.toString());
			
			// Compute log likelihood
			double ll = 0;
			for (int i=0;i<data.size();i++){
				ArrayList<PWMHit> hitGroup = data.get(i);
				
				int maxPwmLen = 0;
				for (PWMHit h: hitGroup){
					if(zeroComp.contains(h.clusterId))				// skip the eliminated PWM components
						continue;
					int w = wmLen[h.clusterId];
					if (maxPwmLen<w)
						maxPwmLen=w;
				}
				if (maxPwmLen==0)
					continue;
				double resp_sum=0;
				for (PWMHit h: hitGroup){
					if(zeroComp.contains(h.clusterId))	{			// skip the eliminated PWM components
						h.responsibility = 0;
						continue;
					}
					float[][]pfm = clusters.get(h.clusterId).pfm;
					double prob = 1;								// this need to be updated for each PWM update
					String seq = h.str;
					for (int p=0;p<seq.length();p++)
						prob *= pfm[p][seq.charAt(p)];
					int N = maxPwmLen - seq.length();
					if (N!=0)
						prob *= Math.pow(0.25,N);		
					h.responsibility = prob*clusters.get(h.clusterId).pi;
					resp_sum += h.responsibility;
				}
				if (resp_sum!=0)
					ll += Math.log(resp_sum);
			}
			
            // log prior
            double LP=0;
            for(KmerCluster c: clusters)
                if (!zeroComp.contains(c.clusterId))
                    LP+=(-currAlpha)*Math.log(c.pi);		// sparse prior
            double lap = ll+LP;
            
			if (verbose>1)
				System.out.println(String.format("Alpha=%.1f\tLAP diff=%.5f", currAlpha, Math.abs(lap-logPosterior)));
			if (Math.abs(lap-logPosterior)>0.0001){
				logPosterior = lap;
			}
			else{		// normalize responsibility before exit the loop
				if (currAlpha<alpha)
					continue;
				for (int i=0;i<data.size();i++){
					ArrayList<PWMHit> hitGroup = data.get(i);
					if (hitGroup.size()==1){					// only 1 hit, the cluster takes all responsibility
						for (PWMHit h: hitGroup){
							if(zeroComp.contains(h.clusterId))	{			// skip the eliminated PWM components
								h.responsibility = 0;
								continue;
							}
							h.responsibility = 1;
						}
						continue;
					}			
					// normalize data likelihood to get responsibility
					double totalResp=0;
					for (PWMHit h: hitGroup){
						totalResp += h.responsibility;
					}
					for (PWMHit h: hitGroup){			
						h.responsibility = h.responsibility/totalResp;
					}
				}
				break;
			}
		}
		
		return conflict;
	}	
	public int resolveConflictHits_old(ArrayList<Sequence> seqList){
		int conflict = 0;
		for (int i=0;i<seqList.size();i++){
			ArrayList<PWMHit> hits = new ArrayList<PWMHit>();
			for (int j=0; j<clusters.size(); j++){
				if (clusters.get(j).seq2hits.containsKey(i))
					hits.add(clusters.get(j).seq2hits.get(i));
			}
			if (hits.size()<=1)
				continue;
			Collections.sort(hits);
			ArrayList<Integer> toRemove= new ArrayList<Integer>();
			for (int j=0;j<hits.size()-1;j++){
				if (toRemove.contains(j))
					continue;
				PWMHit hit = hits.get(j);
				for (int k=j+1;k<hits.size();k++){
					if (toRemove.contains(k))
						continue;
					if (hit.overlaps(hits.get(k)))
						toRemove.add(k);
				}
			}
			if (toRemove.isEmpty())
				continue;
			for (int idx:toRemove){
				PWMHit hit = hits.get(idx);
				clusters.get(hit.clusterId).seq2hits.remove(hit.seqId);
				conflict++;
			}
		}
		return conflict;
	}
	public int resolveConflictHits_old2(ArrayList<Sequence> seqList){
		int conflict = 0;
		int wmLen[]=new int[clusters.size()];
		WeightMatrixScorer[] scorers=new WeightMatrixScorer[clusters.size()];
		for (int j=0; j<clusters.size(); j++){
			wmLen[j]=clusters.get(j).wm.length();
			scorers[j]=new WeightMatrixScorer(clusters.get(j).wm);
		}
		
		for (int i=0;i<seqList.size();i++){
			ArrayList<PWMHit> hits = new ArrayList<PWMHit>();
			for (int j=0; j<clusters.size(); j++){
				if (clusters.get(j).seq2hits.containsKey(i))
					hits.add(clusters.get(j).seq2hits.get(i));
			}
			if (hits.size()<=1)
				continue;
			Collections.sort(hits, new Comparator<PWMHit>(){
			    public int compare(PWMHit o1, PWMHit o2) {
			    	return o1.compareByPosition(o2);
			    }
			});	
			HashSet<Integer> toRemove= new HashSet<Integer>();
			for (int j=0;j<hits.size()-1;j++){
				PWMHit hit = hits.get(j);
				for (int k=j+1;k<hits.size();k++){
					PWMHit hit2 = hits.get(k);
					if (hit.overlaps(hit2)){
						int start = Math.max(hit.start, hit2.start);
						int end = Math.min(hit.end, hit.end);
						int width = Math.max(0,Math.max(wmLen[j], wmLen[k])-(end-start+1));
						String seq = CommonUtils.padding(width, 'N')+seqList.get(i).getSeqStrand(true).substring(start, end+1)+CommonUtils.padding(width, 'N');
						Pair<Integer, Double> score1 = CommonUtils.scanPWM(seq, wmLen[j], scorers[j]);
						Pair<Integer, Double> score2 = CommonUtils.scanPWM(seq, wmLen[k], scorers[k]);
						if (score1.cdr()>score2.cdr())
							toRemove.add(k);
						else if (score1.cdr()<score2.cdr())
							toRemove.add(j);
						else{						// if equal, remove all
							toRemove.add(k);
							toRemove.add(j);
						}
					}
				}
			}
			if (toRemove.isEmpty())
				continue;
			for (int idx:toRemove){
				PWMHit hit = hits.get(idx);
				clusters.get(hit.clusterId).seq2hits.remove(hit.seqId);				//TODO: fix out of range error
			}
			conflict+=toRemove.size();
		}
		return conflict;
	}
	public int removeAllConflictHits(ArrayList<Sequence> seqList){
		int conflict = 0;
		for (int i=0;i<seqList.size();i++){
			ArrayList<PWMHit> hits = new ArrayList<PWMHit>();
			for (int j=0; j<clusters.size(); j++){
				if (clusters.get(j).seq2hits.containsKey(i))
					hits.add(clusters.get(j).seq2hits.get(i));
			}
			if (hits.size()<=1)
				continue;
			Collections.sort(hits, new Comparator<PWMHit>(){
			    public int compare(PWMHit o1, PWMHit o2) {
			    	return o1.compareByPosition(o2);
			    }
			});	
			HashSet<Integer> toRemove= new HashSet<Integer>();
			for (int j=0;j<hits.size()-1;j++){
				PWMHit hit = hits.get(j);
				for (int k=j+1;k<hits.size();k++){
					if (hit.overlaps(hits.get(k))){
						toRemove.add(j);
						toRemove.add(k);
					}
				}
			}
			if (toRemove.isEmpty())
				continue;
			for (int idx:toRemove){
				PWMHit hit = hits.get(idx);
				clusters.get(hit.clusterId).seq2hits.remove(hit.seqId);				//TODO: fix out of range error
			}
			conflict+=toRemove.size();
		}
		return conflict;
	}

	private int buildPWM(ArrayList<Sequence> seqList, KmerCluster cluster, double noiseRatio, long tic, boolean onlyBetter){	    	
		ArrayList<String> alignedSeqs = new ArrayList<String>();
		ArrayList<Double> weights = new ArrayList<Double>();
		for (Sequence seq:seqList){
			if (seq.pos==UNALIGNED)
				continue;
//			if (scorer!=null){
//				Pair<Integer, Double> hit = CommonUtils.scanPWM(seq.getSeq(), cluster.wm.length(), scorer);
//				if (hit.cdr()<cluster.pwmThreshold*0.5)
//					continue;
//				if (hit.car()+seq.pos-cluster.pos_pwm_seed!=0)
//					continue;
//			}
			// align window = k_win
//			int start = k/2-k_win/2-seq.pos;
//			int end = start+k_win;
//			String startPadding = "";
//			String endPadding = "";
//			if (start<0){
//				startPadding = CommonUtils.padding(-start, "N");
//				start=0;
//			}
//			if (end>k_win){
//				endPadding = CommonUtils.padding(end-k_win, "N");
//				end=k_win;
//			}
			// align window = k*2+1		
			int start = k/2-k-seq.pos;
			int end = start+2*k+1;
			String startPadding = "";
			String endPadding = "";
			if (start<0){
				startPadding = CommonUtils.padding(-start, "N");
				start=0;
			}
			int seqLen = seq.getSeq().length();
			if (end>seqLen){
				endPadding = CommonUtils.padding(end-seqLen, "N");
				end=seqLen;
			}
			if (start>=end)
				continue;
 			String s = startPadding+seq.getSeq().substring(start, end)+endPadding;
 			alignedSeqs.add(s);
 			
 			// determine the weight for the sequence
 			double weight = 1.0;
			if (config.use_strength_weight)
				weight *= seq_weights[seq.id];
			if (cluster.clusterId==0 && config.use_pos_weight){
				int prof_pos = k/2-seq.pos;
				if (prof_pos<0)
					prof_pos = 0;
				else if (prof_pos>seqLen-1)
					prof_pos = seqLen-1;
				weight *=profile[prof_pos];
			}
			weights.add(weight);
    	}
		if (alignedSeqs.size()<seqs.length*config.motif_hit_factor){
			if (verbose>1)
	    		System.out.println(String.format("%s: Number of sequences %d is too few, stop building PWM.", CommonUtils.timeElapsed(tic), alignedSeqs.size()));
			return -1;
		}
		
		int bestLeft = buildPWMfromAlignedSequences(alignedSeqs, weights, cluster, noiseRatio, onlyBetter);
		
    	return bestLeft;
    }
	
	public int buildPWMfromHits(ArrayList<Sequence> seqList, KmerCluster cluster, Iterator<PWMHit> hits){
		ArrayList<String> alignedSeqs = new ArrayList<String>();
		ArrayList<Double> weights = new ArrayList<Double>();
		while(hits.hasNext()){
			PWMHit hit = hits.next();
			String seq = seqList.get(hit.seqId).getSeqStrand(true);
			String s = seq.substring(hit.start, hit.end+1);
			if (!hit.isForward)
				s = SequenceUtils.reverseComplement(s);
			alignedSeqs.add(s);
			weights.add(hit.weight*hit.responsibility);
		}
		// make PWM from aligned sequence segament
		return buildPWMfromAlignedSequences(alignedSeqs, weights, cluster, 0, true);	
	}
	
	public void buildAllClusterPWMfromHits(ArrayList<Sequence> seqList){
		int extend = 0;
		ArrayList<KmerCluster> badClusters = new ArrayList<KmerCluster>();
		for (int j=0; j<clusters.size(); j++){
			KmerCluster cluster = clusters.get(j);
			if(cluster.pi==0)	{						// remove the eliminated PWM components
				badClusters.add(cluster);
				continue;
			}
			Iterator<PWMHit> hits = cluster.seq2hits.values().iterator();
			ArrayList<String> alignedSeqs = new ArrayList<String>();
			ArrayList<Double> weights = new ArrayList<Double>();
			while(hits.hasNext()){
				PWMHit hit = hits.next();
				String seq = seqList.get(hit.seqId).getSeqStrand(true);
				int beginIndex = hit.start-extend;
				int left = 0;
				if (beginIndex<0){
					left=-beginIndex;
					beginIndex = 0;
				}
				int endIndex = hit.end+extend+1;		// hit.end is inclusive, endIndex is exclusive
				int right=0;
				if (endIndex>k_win){
					right = endIndex-k_win;
					endIndex = k_win;
				}
				String s = CommonUtils.padding(left, 'N')+
							seq.substring(beginIndex, endIndex)+
							CommonUtils.padding(right, 'N');
				if (!hit.isForward)
					s = SequenceUtils.reverseComplement(s);
				alignedSeqs.add(s);
				weights.add(hit.weight*hit.responsibility);
			}
			// make PWM from aligned sequence segament
			int result = buildPWMfromAlignedSequences(alignedSeqs, weights, cluster, 0, false);
			if (result==-1){
				badClusters.add(cluster);
				if (verbose>1)
					System.out.println("Cluster "+j+" failed to form a PWM, remvoe it.");
			}
		}
		clusters.removeAll(badClusters);
	}
	
	private int buildPWMfromAlignedSequences(ArrayList<String> alignedSeqs, ArrayList<Double> weights, 
		KmerCluster cluster, double noiseRatio, boolean onlyBetter){
		if (alignedSeqs.isEmpty())
			return -1;
		double[][] pfm = new double[alignedSeqs.get(0).length()][MAXLETTERVAL];
    	if (verbose>1)
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
    		if (ic[p]>=ic_trim){
    			leftIdx = p;
    			break;
    		}
    	}
    	int rightIdx=0;
    	for (int p=ic.length-1;p>=0;p--){
    		if (ic[p]>=ic_trim){    			
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
			if (verbose>1)
				System.out.println("PWM is too short, W="+(rightIdx-leftIdx+1));
			return -1;
		}
		
		/* try all pwm length with the most IC-rich columns, find the best PWM */
//TODO:		if (progressive_PWM_trim){
		int[] left, right;
		int tooShort = k-2;
		if (k<=7)
			tooShort = k-1;
		if (rightIdx-leftIdx+1>tooShort){		// length > k-1
			left=new int[rightIdx-leftIdx+1-tooShort];
			right=new int[rightIdx-leftIdx+1-tooShort];
			for (int i=0;i<left.length;i++){
				int bestLeft = -1;
				double bestSumIC = 0;
				for(int p=leftIdx;p<=rightIdx-tooShort-i;p++){
					int end = tooShort+i+p;
					if (ic[p]<ic_trim || ic[end]<ic_trim)			// if the ends have low ic, skip
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
				right[i]=bestLeft+tooShort+i;
			}
		}
		else{				// if it is not very long
			left=new int[1];
			right=new int[1];
			left[0]=leftIdx;
			right[0]=rightIdx;
		}
		
		MotifThreshold bestEstimate = null;
    	double bestHGP = -0.1;
    	WeightMatrix bestWM = null;
    	int bestLeft=0;
    	int bestRight=0;    
		for (int i=0;i<left.length;i++){
			if (left[i]==-1)
				continue;
	    	float[][] matrix = new float[right[i]-left[i]+1][MAXLETTERVAL];   
	    	for(int p=left[i];p<=right[i];p++){
	    		for (char base:LETTERS){							// ignore 'N' count
	    			matrix[p-left[i]][base]=(float) pwm[p][base];
	    		}
	    	}

	    	WeightMatrix wm = new WeightMatrix(matrix);

//	    	if (bmverbose>1)
//	    		System.out.println(String.format("%s: got PWM.", CommonUtils.timeElapsed(tic)));
	    	// Check the quality of new PWM: hyper-geometric p-value test using the positive and negative sequences
	    	MotifThreshold estimate = estimatePwmThreshold(wm, wm.getMaxScore()*config.wm_factor);
//		    	MotifThreshold estimate = optimizePwmThreshold(wm, "", false, wm.getMaxScore()*wm_factor);
	    	double pwmThreshold = estimate.score;
	    	double pwmThresholdHGP = estimate.hgp;
    		if (verbose>1)
    			if (pwmThresholdHGP==0)
    				System.out.println(String.format("%s: PWM %s is not enriched", CommonUtils.timeElapsed(tic), WeightMatrix.getMaxLetters(wm)));
        		else
        			System.out.println(String.format("%s: PWM score %.2f/%.2f\tmatch %d+/%d- seqs\thgp=1e%.1f\t%s", CommonUtils.timeElapsed(tic), 
        					pwmThreshold, wm.getMaxScore(), estimate.posHit, estimate.negHit, pwmThresholdHGP, WeightMatrix.getMaxLetters(wm)));
    		if (pwmThresholdHGP<=bestHGP){
    			bestWM = wm;
    			bestHGP = pwmThresholdHGP;
    			bestEstimate = estimate;
    			bestLeft=left[i];
    			bestRight=right[i];
    		}
		}
		if (bestEstimate==null) {
			if (verbose>1)
				System.out.println(CommonUtils.timeElapsed(tic)+": None of PWM is enriched.");
			return -1;
		}
		
		// test if we want to accept the new PWM
		if (cluster.wm!=null){
			if( cluster.pwmThresholdHGP<=bestHGP && onlyBetter){		// previous pwm is more enriched, stop here
				return -2;
			}
		}

    	cluster.wm = bestWM;
    	cluster.pwmGoodQuality = (bestEstimate.posHit>seqs.length*config.motif_hit_factor);
    	cluster.pwmThreshold = bestEstimate.score;
    	cluster.pwmThresholdHGP = bestHGP;
    	cluster.pwmPosHitCount = bestEstimate.posHit;
    	cluster.pwmNegHitCount = bestEstimate.negHit;
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
    	cluster.pfm = pfm_trim;
    	cluster.pos_pwm_seed = bestLeft-(alignedSeqs.get(0).length()/2-k/2);		// pwm_seed = pwm_seqNew-seed_seqNew
    	
    	return bestLeft;
	}
	
	private void alignSequencesUsingKSM(ArrayList<Sequence> seqList, KmerCluster cluster){
		ArrayList<Kmer> alignedKmers = cluster.alignedKmers;
		updateEngine(alignedKmers);
		
		MotifThreshold threshold = new MotifThreshold();
		if (config.estimate_ksm_threshold){
		threshold = optimizeKsmThreshold("", false);
			if (threshold==null)
				return;
			if (verbose>1)
				System.out.println(String.format("%s: KSM score %.2f\tmatch %d+/%d- seqs\thgp=1e%.1f", 
						CommonUtils.timeElapsed(tic), threshold.score, threshold.posHit, threshold.negHit, threshold.hgp));
		}
//		if (threshold.hgp<cluster.ksmThreshold.hgp){
//			cluster.ksmThreshold = threshold;
			
			// scan all sequences, align them using kmer hit results
			for (Sequence s : seqList){
				s.reset();			
				KmerGroup[] kgs = queryS(s.getSeq());
				//Check if there are matches on the same position from both strand??
				for (int i=0;i<kgs.length-1;i++){
					KmerGroup kg=kgs[i];
					for (int j=i;i<kgs.length;i++){
						KmerGroup other=kgs[j];
						if (Math.abs(other.bs-kg.bs)==RC){	// if so, add all kmers to the stronger kmer
							System.out.println("Kmer match on same position:"+kg.toString()+" "+other.toString()); //TODO: Comment out on release
							if (kg.hgp<other.hgp){
								kg.kmers.addAll(other.kmers);
								kg.setHgp(computeHGP(kg.getGroupHitCount(), kg.getGroupNegHitCount()));
							}
							else{
								other.kmers.addAll(kg.kmers);
								other.setHgp(computeHGP(other.getGroupHitCount(), other.getGroupNegHitCount()));										
							}
						}
					}
				}
				if (kgs.length!=0){
					Arrays.sort(kgs);
					KmerGroup kg = kgs[0];
					if (kg.hgp<=-threshold.score){
						s.score = -kg.hgp;		// score = -log10 hgp, becomes positive value
						if (kg.bs>RC/2){		// match on reverse strand
							s.RC();
							s.pos = -(kg.bs-RC);
						}
						else
							s.pos = -kg.bs;
					}
				}
			}
//		}
	}
	/** align all sequences using a PWM<br>
	 * sequences have a PWM hit is aligned according hit position, if no PWM hit, set it to unaligned
	 */
	private int alignSequencesUsingPWM(ArrayList<Sequence> seqList, KmerCluster cluster){
			    	
    		WeightMatrix wm = cluster.wm;
	        WeightMatrixScorer scorer = new WeightMatrixScorer(wm);		
			int count_pwm_aligned=0;
			for (Sequence s:seqList){
	    	  String seq = s.getSeq();			// PWM to scan all sequence, and align if pass threshold
//	    	  if (s.pos!=UNALIGNED || seq.length()<wm.length())
//	    		  continue;
	    	      	  
	          WeightMatrixScoreProfile profiler = scorer.execute(seq);
	          PwmMatch match = profiler.getBestMatch(config.strand_type);	          
	          // if a sequence pass the motif score, align with PWM hit
	          if (match.score >= cluster.pwmThreshold){
				if (match.strand =='-'){
					match.shift = k_win-match.shift-wm.length();
					s.RC();
					// i.e.  (seq.length()-1)-maxScoringShift-(wm.length()-1);
				}
				s.pos = cluster.pos_pwm_seed-match.shift;
				count_pwm_aligned ++;
	          }
	          else
	        	  s.pos = UNALIGNED;
	        }	// each sequence

	    	if (verbose>1)
	    		System.out.println(CommonUtils.timeElapsed(tic)+": PWM "+WeightMatrix.getMaxLetters(wm)+" align "+count_pwm_aligned+" sequences.");
	    	
	    	return count_pwm_aligned;
	}
	
	/**
	 * Get all the kmers that has k/4 (k>=8) or 1 (k<8) mismatch to the kmerStr<br>
	 * and set shift and alignString for the kmers 
	 * @param kmers
	 * @param kmerStr the length should be the same as kmers
	 * @param shift
	 * @return
	 */
	private ArrayList<Kmer> getMMKmers(ArrayList<Kmer> kmers, String kmerStr, int shift) {
		ArrayList<Kmer> family = new ArrayList<Kmer>();
		if (kmers.isEmpty())
			return family;
		ArrayList<Kmer> copy = (ArrayList<Kmer>) kmers.clone();
		ArrayList<Kmer> toRemove = new ArrayList<Kmer>();
		String kmerRC = SequenceUtils.reverseComplement(kmerStr);
		int mm = 1;
		if (k>=8 && this.use_smart_mm)
			mm = k/4;
		// progressively allow more mismatch
		for (int i=1;i<=mm;i++){
			for (Kmer kmer: copy){
		    	if (CommonUtils.mismatch(kmerStr, kmer.getKmerString())==i)
		    		kmer.setAlignString(kmer.getKmerString());
		    	else if (config.strand_type!=1 && CommonUtils.mismatch(kmerRC, kmer.getKmerString())==i)	// if not single strand
		    		kmer.setAlignString(kmer.getKmerRC());		// do not RC kmer, set the string for alignment
		    	else 
		    		continue;
	    		kmer.setShift(0);
	    		family.add(kmer);
		    	toRemove.add(kmer);
			}
			copy.removeAll(toRemove);
			toRemove.clear();
		}
		return family;
	}
	/**
	 * Get aligned kmers within seed_range using the aligned sequences
	 * @param seqList
	 * @param seed_range
	 * @param kmer_inRange_fraction
	 * @return
	 */
	private ArrayList<Kmer> getAlignedKmers (ArrayList<Sequence> seqList, int seed_range, ArrayList<Kmer> excludes){

    	/** set kmer consensus position */
		HashMap<Kmer, ArrayList<Integer>> kmer2pos = new HashMap<Kmer, ArrayList<Integer>>();
		for (Sequence s:seqList){
			if (s.pos != UNALIGNED){		// aligned seqs
				HashMap<Kmer, HashSet<Integer>> kmerPos = s.getKmerPos();
				for (Kmer km:kmerPos.keySet()){
					if (excludes.contains(km))
						continue;
					if (!kmer2pos.containsKey(km))
						kmer2pos.put(km, new ArrayList<Integer>());
					HashSet<Integer> hits = kmerPos.get(km);
					if (hits==null)
						continue;
					for (int pRef:hits){
						if (pRef>RC/2){
							int pos = (pRef-RC)+s.pos;
							if (pos>=-seed_range && pos<=seed_range)
								kmer2pos.get(km).add(pos+RC);		// use kmerRC 
						}
						else{
							int pos = pRef+s.pos;
							if (pos>=-seed_range && pos<=seed_range)
								kmer2pos.get(km).add(pos);		
						}
					}
				}
			}
		}

		ArrayList<Kmer> alignedKmers = new ArrayList<Kmer>();
		for (Kmer km:kmer2pos.keySet()){
			ArrayList<Integer> posKmer = kmer2pos.get(km);		// all in_sequence positions of this kmer
			if (posKmer==null || posKmer.size() < km.getPosHitCount()*config.kmer_inRange_fraction){			// The kmer hit in the 2*k region should be at least 1/2 of total hit
				km.setAlignString("Too few hit "+posKmer.size());
				continue;
			}	
			// find the most frequent kmerPos
			Pair<int[], int[]> sorted = StatUtil.sortByOccurences(posKmer);
			int counts[] = sorted.cdr();
			int posSorted[] = sorted.car();
			int maxCount = counts[counts.length-1];
			if (maxCount < Math.min(posKmer.size(),km.getPosHitCount()) * config.kmer_inRange_fraction)	// for palindromic kmer, posKmer>hitCount
				continue;
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
						pos-=RC;
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
//					
			km.setShift(shift);
			
			km.setAlignString(maxCount+"/"+posKmer.size());
			km.setKmerStartOffset(km.getShift());
			alignedKmers.add(km);
		}	
		kmer2pos = null;
		alignedKmers.trimToSize();
		return alignedKmers;
	}


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
    
	public class Sequence implements Comparable<Sequence>{
		int id;				// original input id
		String seq;
		String rc;
		int pos=UNALIGNED;
		private boolean isForward = true;
		
		HashMap<Kmer, HashSet<Integer>> fPos = new HashMap<Kmer, HashSet<Integer>>();		// forward
		HashMap<Kmer, HashSet<Integer>> rPos = new HashMap<Kmer, HashSet<Integer>>();		// reverse
		int totalCount;
		int maxCount;
		int kmerPosCount;
		double score;
		
		public Sequence(String seq, int id){
			this.id = id;
			this.seq = seq;
			this.rc = SequenceUtils.reverseComplement(seq);
		}
		
		public void RC(){
			isForward = !isForward;
		}
		public String getSeq(){
			return isForward?seq:rc;
		}
		public void setSeq(String seq){
			if (isForward){
				this.seq = seq;
				this.rc = SequenceUtils.reverseComplement(seq);
			}
			else{
				this.rc = seq;
				this.seq = SequenceUtils.reverseComplement(seq);
			}
		}

		public void reset(){
			pos = UNALIGNED;
		}
		public void removeAllKmers(Collection<Kmer> kmers){
			for (Kmer km:kmers){
				fPos.remove(km);
				rPos.remove(km);
			}
		}
		public String getSeqStrand(boolean isForward){
			return isForward?seq:rc;
		}
		public HashMap<Kmer, HashSet<Integer>> getKmerPosStrand(boolean isForward){
			return isForward?fPos:rPos;
		}
		public HashMap<Kmer, HashSet<Integer>> getKmerPos(){
			return isForward?fPos:rPos;
		}
		public HashMap<Kmer, HashSet<Integer>> getKmerRCPos(){
			return (!isForward)?fPos:rPos;
		}
		public int getKmerPosCount(){
			return kmerPosCount;
		}		
		public void setCount(){
			maxCount=0;
			totalCount=0;
			kmerPosCount=0;
			for (Kmer km:fPos.keySet()){
				int count = km.getPosHitCount();
				if (maxCount<count)
					maxCount=count;
				totalCount+=count;
				kmerPosCount+=fPos.get(km).size();
			}
		}
	
		public int compareTo(Sequence s) {					// descending count
			int max = compareByMaxCount(s);
			if (max!=0)
				return max;
			if(totalCount<s.totalCount){return(1);}
			else if(totalCount>s.totalCount){return(-1);}
			else return(0);
		}
		
		public int compareByMaxCount(Sequence s) {					// descending count
			if(maxCount<s.maxCount){return(1);}
			else if(maxCount>s.maxCount){return(-1);}
			else return(0);
		}
		
		public String toString(){
			return String.format("%d\t%d\t%s\t%d\t%s", id, maxCount, isForward?"F":"R", pos, isForward?seq:rc);
		}
	}

	private int getRcCoor(int original){
		return numPos-1-original;
	}
	
	private Pair<Double, Integer> findBestPosition(Sequence top, ArrayList<Sequence> seqAligned){
		//estimate the potential sequence positions
		HashSet<Integer> allPos = new HashSet<Integer>();
		// based on each kmer's most frequent position in the alignment
		HashMap<Kmer, HashSet<Integer>> kmer2pos = top.getKmerPos();
		if (kmer2pos.isEmpty())
			return new Pair<Double, Integer>(Double.NEGATIVE_INFINITY, 0);
		
		for (Kmer km:kmer2pos.keySet()){
			HashSet<Integer> kmerPos = kmer2pos.get(km);				
			ArrayList<Integer> possiblePos = new ArrayList<Integer>(); //kmer's all positions in the alignment
			for (Sequence sa:seqAligned){
				if (!sa.getKmerPos().containsKey(km))
					continue;
				for (int pRef:sa.getKmerPos().get(km)){
					if (pRef>RC/2)
						possiblePos.add((pRef-RC)+sa.pos+RC);		// kmer on rc of sa, use kmerRC coor
					else
						possiblePos.add(pRef+sa.pos);
				}
			}
			if (possiblePos.isEmpty())
				continue;
			
			Pair<int[],int[]> sorted = StatUtil.sortByOccurences(possiblePos);
			int maxCount = sorted.cdr()[sorted.cdr().length-1];		// get the most frequent position
			for (int i=sorted.cdr().length-1;i>=0;i--){
				if (sorted.cdr()[i]==maxCount){
					int goodPos = sorted.car()[i];
					if (goodPos>RC/2){				// kmer on rc of sa
						for (int p: kmerPos){
							if (p>RC/2)				// kmer on rc of top, same strand, use kmerRC coor
								allPos.add((goodPos-RC)-(p-RC));
							else					// kmer on seq of top
								allPos.add((goodPos-RC)-p+RC);
						}	
					}		
					else{							// kmer on seq of sa
						for (int p: kmerPos){
							if (p>RC/2)				// kmer on rc of top, diff strand, use kmerRC coor
								allPos.add(goodPos-(p-RC)+RC);
							else
								allPos.add(goodPos-p);
						}
					}
				}
			}
		}
		if (allPos.isEmpty())
			return new Pair<Double, Integer>(Double.NEGATIVE_INFINITY, 0);
		
		Integer[] allPosArray = new Integer[allPos.size()];
		allPos.toArray(allPosArray);
		double[][] scores = new double[allPosArray.length][seqAligned.size()];
		for (int i=0;i<allPosArray.length;i++){
			int p=allPosArray[i];
			for (int j=0;j<seqAligned.size();j++){
				Sequence ref = seqAligned.get(j);
				scores[i][j] = scoreSequenceAlignment(ref, top, p);
			}
		}
		// top 10 scores
		int topSeqNum = Math.min(seqAligned.size(), 10);
		double[] topScores = new double[allPosArray.length];
		for (int i=0;i<allPosArray.length;i++){
			double scorePos[] = scores[i];
			Arrays.sort(scorePos);
			for (int j=0;j<topSeqNum;j++){
				topScores[i]+=scorePos[scorePos.length-1-j];
			}
		}
		Pair<Double, TreeSet<Integer>> max = StatUtil.findMax(topScores);
		return new Pair<Double, Integer>(max.car(), allPosArray[max.cdr().first()]);
	}
	
	public int scoreSequenceAlignment(Sequence ref, Sequence top, int sPos){
		int score=0;
		HashMap<Kmer, HashSet<Integer>> sKmerPos = (sPos>RC/2)?top.getKmerRCPos():top.getKmerPos();
		if (sPos>RC/2)
			sPos-=RC;
		for (Kmer km:ref.getKmerPos().keySet()){
			if (!top.getKmerPos().containsKey(km))
				continue;
			boolean match = false;
			findMatch: for (int pRef:ref.getKmerPos().get(km))
				for (int pS:sKmerPos.get(km)){
					if ((pRef+ref.pos)-(pS+sPos)==0){
						match = true;
						break findMatch;
					}
				}
			if (match)
				score+=km.getPosHitCount();
			else
				score-=km.getPosHitCount();
		}
		return score;
	}
	
	public class KmerCluster implements Comparable<KmerCluster>{
		public int pos_pwm_seed;
		public int pos_BS_seed;
		public MotifThreshold ksmThreshold = new MotifThreshold();
		public WeightMatrix wm;
		public ArrayList<Kmer> alignedKmers;			// The K-mer set motif, a set of aligned k-mers
		
		int clusterId;
		Kmer seedKmer;
		float[][] pfm;
		boolean pwmGoodQuality = false;
		public double pwmThreshold;
		public double pwmThresholdHGP;
		public int pwmNegHitCount;
		public int pwmPosHitCount;
		int total_aligned_seqs;
		HashMap<Integer, PWMHit> seq2hits = null;
		double pi;
		
		KmerCluster(){}	// empty constructor;

		protected KmerCluster clone(){
			KmerCluster cluster = new KmerCluster();
			cluster.clusterId = this.clusterId;
			cluster.seedKmer = this.seedKmer;
			cluster.pfm = this.pfm.clone();
			cluster.wm = this.wm;
			cluster.pwmGoodQuality = this.pwmGoodQuality;
			cluster.pwmThreshold = this.pwmThreshold;
			cluster.pwmThresholdHGP = this.pwmThresholdHGP;
			cluster.pwmNegHitCount = this.pwmNegHitCount;
			cluster.pos_pwm_seed = this.pos_pwm_seed;
			cluster.pos_BS_seed = this.pos_BS_seed;
			cluster.pwmPosHitCount = this.pwmPosHitCount;
			cluster.alignedKmers = this.alignedKmers;
			return cluster;
		}
		
		public int compareTo(KmerCluster c) {					// ascending pwmThresholdHGP
			if(pwmThresholdHGP<c.pwmThresholdHGP){return(-1);}
			else if(pwmThresholdHGP>c.pwmThresholdHGP){return(+1);}
			else return(0);
		}
		/**
		 * This method is only used for selecting K when the hgp of the clusters are very close (<5%)<br>
		 * It first select for k that gives pwm width most close to k, then select for pwm that have a highest pwm score per position
		 * @param c
		 * @return
		 */
		public int compareForSelectingK(KmerCluster c) {					
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
			return String.format("PWM %d: %s, pi=%.2f", clusterId, wm!=null?WeightMatrix.getMaxLetters(wm):"----", pi);
		}
	}
	
	/**
	 * Compute hgp (log10) using the positive/negative sequences
	 */	
	public double computeHGP(int posHitCount, int negHitCount){
		if (posHitCount==0)
			return 0;
		return computeHGP(posSeqCount, negSeqCount, posHitCount, negHitCount);
	}
	/**
	 * Compute hgp (log10) using the positive/negative sequences
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
	private static double computeHGP_TINY(int posSeq, int negSeq, int posHit, int negHit){
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
	 * Compute hyper-geometric p-values (log10) for the kmer strings<br>
	 */
	public double[] computeHGPs(ArrayList<String> kmerStrings){
		AhoCorasick tree = new AhoCorasick();
		for (int i=0;i<kmerStrings.size();i++){
			tree.add(kmerStrings.get(i).getBytes(), i);
			tree.add(SequenceUtils.reverseComplement(kmerStrings.get(i)).getBytes(), i);
	    }
	    tree.prepare();
	    int[] posHitCount = new int[kmerStrings.size()];
	    int[] negHitCount = new int[kmerStrings.size()];
	    for (String seq: seqs){
			Iterator searcher = tree.search(seq.getBytes());
			while (searcher.hasNext()) {
				SearchResult result = (SearchResult) searcher.next();
				Set<Integer> idxs = result.getOutputs();
				for (int idx:idxs)
					posHitCount[idx]++;
			}
	    }
	    for (String seq: seqsNegList){
			Iterator searcher = tree.search(seq.getBytes());
			while (searcher.hasNext()) {
				SearchResult result = (SearchResult) searcher.next();
				Set<Integer> idxs = result.getOutputs();
				for (int idx:idxs)
					negHitCount[idx]++;
			}
	    }
	    double hgps[] = new double[kmerStrings.size()];
	    for (int i=0;i<kmerStrings.size();i++){
			hgps[i] = computeHGP(posSeqCount, negSeqCount, posHitCount[i], negHitCount[i]);
//	    	hgps[i] = 1-StatUtil.hyperGeometricCDF_cache(posHitCount[i], totalSeqCount, posHitCount[i]+negHitCount[i], posSeqCount);
	    }
	    return hgps;
	}
	
	/**
	 * Update k-mer hit counts and compute hyper-geometric p-values<br>
	 */
	public void updateKmerCounts(ArrayList<Kmer> kmers, ArrayList<ComponentFeature> events){
		AhoCorasick tree = new AhoCorasick();
		for (int i=0;i<kmers.size();i++){
			String kmStr = kmers.get(i).getKmerString();
			tree.add(kmStr.getBytes(), i);
			tree.add(SequenceUtils.reverseComplement(kmStr).getBytes(), i);
	    }
	    tree.prepare();
	    ArrayList<HashSet<Integer>> posHits = new ArrayList<HashSet<Integer>>(kmers.size());
	    for (int i=0;i<kmers.size();i++)
	    	posHits.add(new HashSet<Integer>());
	    ArrayList<HashSet<Integer>> negHits = new ArrayList<HashSet<Integer>>(kmers.size());
	    for (int i=0;i<kmers.size();i++)
	    	negHits.add(new HashSet<Integer>());
//	    double[] kmerStrength = new double[kmers.size()];		// TODO: ignore kmer Strength for now
	    for (int i=0;i<posSeqCount;i++){
	    	String seq = seqs[i];
			Iterator searcher = tree.search(seq.getBytes());
			HashSet<Integer> idxs = new HashSet<Integer>();
			while (searcher.hasNext()) {
				SearchResult result = (SearchResult) searcher.next();
				idxs.addAll(result.getOutputs());
			}
			double kmerSum = 0;
//			for (int idx:idxs){
//				kmerSum += kmers.get(idx).getPosHitCount();
//			}
			for (int idx:idxs){
				posHits.get(idx).add(i);
//				kmerStrength[idx] += kmerSum==0?0:kmers.get(idx).getPosHitCount()/kmerSum*events.get(i).getTotalEventStrength();
			}
	    }
	    for (int i=0;i<negSeqCount;i++){
	    	String seq = seqsNegList.get(i);
			Iterator searcher = tree.search(seq.getBytes());
			HashSet<Integer> idxs = new HashSet<Integer>();
			while (searcher.hasNext()) {
				SearchResult result = (SearchResult) searcher.next();
				idxs.addAll(result.getOutputs());
			}
			for (int idx:idxs)
				negHits.get(idx).add(i);
	    }
	    for (int i=0;i<kmers.size();i++){
	    	Kmer km = kmers.get(i);
	    	km.setPosHits(posHits.get(i));
	    	km.setNegHits(negHits.get(i));
			km.setHgp( computeHGP(km.getPosHitCount(), km.getNegHitCount()));
			km.setStrength(seq_weights[i]);
	    }
	}
	
	/**
	 * Optimize the threshold of a PWM (larger than startingScore) using the positive/negative sequences<br>
	 * Multi-thread to compute HGP<br>
	 * Approximate grid search to find best HGP, to reduce run time
	 */
	public MotifThreshold optimizePwmThreshold(WeightMatrix wm, String outName, boolean printPwmHgp, double startingScore){
		double[] posSeqScores = new double[posSeqCount];
		double[] negSeqScores = new double[negSeqCount];
		for (int i=0;i<posSeqCount;i++){
			posSeqScores[i]=WeightMatrixScorer.getMaxSeqScore(wm, seqs[i], config.strand_type==1);
		}
//		Arrays.sort(posSeqScores);
		int[] posIdx = StatUtil.findSort(posSeqScores);		// index of sequence after sorting the scores
		
		int startIdx = Arrays.binarySearch(posSeqScores, startingScore);
		if( startIdx < 0 ) { startIdx = -startIdx - 1; }
		
		for (int i=0;i<negSeqCount;i++){
			negSeqScores[i]=WeightMatrixScorer.getMaxSeqScore(wm, seqsNegList.get(i), config.strand_type==1);
		}
		Arrays.sort(negSeqScores);
		
		TreeSet<Double> posScoreUnique = new TreeSet<Double>();
		for (int i=startIdx;i<posSeqScores.length;i++)
			posScoreUnique.add(posSeqScores[i]);
		if (posScoreUnique.isEmpty()){						// this could happen if all pwm scores are less than 0
			MotifThreshold score = new MotifThreshold();
			score.score = 0;
			score.hgp = 0;
			score.posHit = 0;
			score.negHit = 0;
			return score;
		}
		
		// count hits at each score, compute hgp
		Double[] posScores_u = new Double[posScoreUnique.size()];
		posScoreUnique.toArray(posScores_u);
		int[] poshits = new int[posScoreUnique.size()];
		int[] neghits = new int[posScoreUnique.size()];
		double[] hgps = new double[posScoreUnique.size()];
		for (int i=0;i<posScores_u.length;i++){
			double key = posScores_u[i];
			int index = CommonUtils.findKey(posSeqScores, key);
			if (config.use_weighted_kmer){
				double weightedHit = 0;
				for (int s=index; s<posSeqScores.length; s++)
					weightedHit += seq_weights[ posIdx[s] ];
				poshits[i] = (int) weightedHit;
			}
			else
				poshits[i] = posSeqScores.length-index;
			index = CommonUtils.findKey(negSeqScores, key);
			neghits[i] = negSeqScores.length-index;
		}

		ArrayList<Integer> idxs = new ArrayList<Integer>();			// the score ids to compute HGP
		for (int i=posScores_u.length-1;i>=0;i--)
			if (poshits[i]>neghits[i]*2.0*posSeqCount/negSeqCount)	// posHit should be at least 2 fold
				idxs.add(i);
		
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
			
			best = findBestScore(idxCoarse, poshits, neghits, hgps);
			
			// finer resolution search
			int bestIdx = idxs.indexOf(best.cdr());
			int start = Math.max(bestIdx-gridStep+1, 0);
			int end = Math.min(bestIdx+gridStep-1,idxs.size()-1) ;
			ArrayList<Integer> idxFine = new ArrayList<Integer>();	
			for (int i=start;i<=end;i++){
				idxFine.add(idxs.get(i));
			}
			
			best = findBestScore(idxFine, poshits, neghits, hgps);
		}
		else
			best = findBestScore(idxs, poshits, neghits, hgps);
		
		MotifThreshold score = new MotifThreshold();
		score.score = posScores_u[best.cdr()];
		score.hgp = best.car();
		score.posHit = poshits[best.cdr()];
		score.negHit = neghits[best.cdr()];
//		if (printPwmHgp)
//			CommonUtils.writeFile(outName+"_"+WeightMatrix.getMaxLetters(wm)+"_PwmHgp.txt", sb.toString());
		return score;
	}
	/**
	 * Estimate threshold of a PWM (larger than and closest to wmScore) using the positive/negative sequences<br>
	 * This is different from optimizePwmThreshold in that it only compute 1 HGP
	 */
	public MotifThreshold estimatePwmThreshold(WeightMatrix wm, double wmScore){
		double[] posSeqScores = new double[posSeqCount];
		double[] negSeqScores = new double[negSeqCount];
		for (int i=0;i<posSeqCount;i++){
			posSeqScores[i]=WeightMatrixScorer.getMaxSeqScore(wm, seqs[i], config.strand_type==1);
		}
		Arrays.sort(posSeqScores);
		int startIdx = Arrays.binarySearch(posSeqScores, wmScore);
		if( startIdx < 0 ) { startIdx = -startIdx - 1; }
		
		for (int i=0;i<negSeqCount;i++){
			negSeqScores[i]=WeightMatrixScorer.getMaxSeqScore(wm, seqsNegList.get(i), config.strand_type==1);
		}
		Arrays.sort(negSeqScores);
		
		TreeSet<Double> posScoreUnique = new TreeSet<Double>();
		for (int i=startIdx;i<posSeqScores.length;i++)
			posScoreUnique.add(posSeqScores[i]);
		if (posScoreUnique.isEmpty()){						// this could happen if all pwm scores are less than 0
			MotifThreshold score = new MotifThreshold();
			score.score = 0;
			score.hgp = 0;
			score.posHit = 0;
			score.negHit = 0;
			return score;
		}
		// count hits at each score, compute hgp
		Double[] posScores_u = new Double[posScoreUnique.size()];
		posScoreUnique.toArray(posScores_u);
		
		MotifThreshold score = new MotifThreshold();
		for (int i=0;i<posScores_u.length;i++){
			score.score = posScores_u[i];
			int index = CommonUtils.findKey(posSeqScores, score.score);
			score.posHit = posSeqScores.length-index;
			index = CommonUtils.findKey(negSeqScores, score.score);
			score.negHit = negSeqScores.length-index;
			if (score.posHit>=score.negHit*2.0*posSeqCount/negSeqCount){	// posHit should be at least 2 fold
				score.hgp = computeHGP(posSeqCount, negSeqCount, score.posHit, score.negHit);
				return score;
			}			
		}
		// if we can not find a threshold that is larger than wmScore and have 2 fold enrichment, return with hgp=0
		score.hgp = 0;
		return score;
	}
	
	private Pair<Double, Integer> findBestScore(ArrayList<Integer> idxs, int[] poshits, int[] neghits, double[] hgps){

		int numThread = Math.min(config.maxThreads, java.lang.Runtime.getRuntime().availableProcessors());
		Thread[] threads = new Thread[numThread];
		for (int i=0;i<numThread;i++){
            Thread t = new Thread(new HGPThread(idxs, posSeqCount, negSeqCount, poshits, neghits, hgps));
            t.start();
            threads[i] = t;
		}
		boolean anyrunning = true;
        while (anyrunning) {
            anyrunning = false;
            try {
                Thread.sleep(100);
            } catch (InterruptedException e) { }
            for (int i = 0; i < threads.length; i++) {
                if (threads[i].isAlive()) {
                    anyrunning = true;
                    break;
                }
            }   
        }
		
		Pair<Double, TreeSet<Integer>> minHgp = StatUtil.findMin(hgps);
		int minIdx = minHgp.cdr().last();
		
		return new Pair<Double, Integer>(minHgp.car(),minIdx);
	}
	
	public class MotifThreshold{
		/** k-mer group score,  score = -log10 hgp, becomes positive value */
		public double score;
		public int posHit;
		public int negHit;
		/* the significance of the k-mer group score */
		public double hgp=0;		
	}

	/**
	 * Grid search threshold of a Kmer Group Score using the positive/negative sequences<br>
	 * Compute the hyper-geometric p-value from number of pos/neg sequences that have the scores higher than the considered score.<br>
	 * Return the score gives the most significant p-value.
	 */
	public MotifThreshold optimizeKsmThreshold(String outName, boolean printKgcHgp){
		if (! engineInitialized)
			return null;
		
		double[] posSeqScores = new double[posSeqCount];
		double[] negSeqScores = new double[negSeqCount];
		for (int i=0;i<posSeqCount;i++){
			KmerGroup[] kgs = query(seqs[i]);
			if (kgs.length==0)
				posSeqScores[i]=0;
			else{
				Arrays.sort(kgs);
				posSeqScores[i]=-kgs[0].getHgp();	// use first kg, the best match	// score = -log10 hgp, becomes positive value
			}
		}
//		Arrays.sort(posSeqScores);
		int[] posIdx = StatUtil.findSort(posSeqScores);		// index of sequence after sorting the scores
		for (int i=0;i<negSeqCount;i++){
			KmerGroup[] kgs = query(seqsNegList.get(i));
			if (kgs.length==0)
				negSeqScores[i]=0;
			else{
				Arrays.sort(kgs);
				negSeqScores[i]=-kgs[0].getHgp();
			}
		}
		Arrays.sort(negSeqScores);
		
		// find the threshold motif score
		TreeSet<Double> posScoreUnique = new TreeSet<Double>();
		for (double s:posSeqScores)
			posScoreUnique.add(s);
		Double[] posScores_u = new Double[posScoreUnique.size()];
		posScoreUnique.toArray(posScores_u);
		int[] poshits = new int[posScoreUnique.size()];
		int[] neghits = new int[posScoreUnique.size()];
		double[] hgps = new double[posScoreUnique.size()];
		for (int i=0;i<posScores_u.length;i++){
			double key = posScores_u[i];
			int index = CommonUtils.findKey(posSeqScores, key);
			if (config.use_weighted_kmer){
				double weightedHit = 0;
				for (int s=index; s<posSeqScores.length; s++)
					weightedHit += seq_weights[ posIdx[s] ];
				poshits[i] = (int) weightedHit;
			}
			else
				poshits[i] = posSeqScores.length-index;
			index = CommonUtils.findKey(negSeqScores, key);
			neghits[i] = negSeqScores.length-index;
		}

		ArrayList<Integer> idxs = new ArrayList<Integer>();			// the score ids to compute HGP
		for (int i=posScores_u.length-1;i>=0;i--)
			if (poshits[i]>neghits[i]*2.0*posSeqCount/negSeqCount)	// posHit should be at least 2 fold
				idxs.add(i);
		
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
			
			best = findBestScore(idxCoarse, poshits, neghits, hgps);
			
			// finer resolution search
			int bestIdx = idxs.indexOf(best.cdr());
			int start = Math.max(bestIdx-gridStep+1, 0);
			int end = Math.min(bestIdx+gridStep-1,idxs.size()-1) ;
			ArrayList<Integer> idxFine = new ArrayList<Integer>();	
			for (int i=start;i<=end;i++){
				idxFine.add(idxs.get(i));
			}
			
			best = findBestScore(idxFine, poshits, neghits, hgps);
		}
		else
			best = findBestScore(idxs, poshits, neghits, hgps);
		
		MotifThreshold score = new MotifThreshold();
		score.score = posScores_u[best.cdr()];
		score.hgp = best.car();
		score.posHit = poshits[best.cdr()];
		score.negHit = neghits[best.cdr()];
//		if (printPwmHgp)
//			CommonUtils.writeFile(outName+"_"+WeightMatrix.getMaxLetters(wm)+"_PwmHgp.txt", sb.toString());
		return score;
	}

	/**
	 * Optimize threshold of a Kmer Group Score using the positive/negative sequences<br>
	 * Compute the hyper-geometric p-value from number of pos/neg sequences that have the scores higher than the considered score.<br>
	 * Return the score gives the most significant p-value.
	 */
	public MotifThreshold estimateKsmThreshold(String outName, boolean printKgcHgp){
		if (! engineInitialized)
			return null;
		
		double[] posSeqScores = new double[posSeqCount];
		double[] negSeqScores = new double[negSeqCount];
		for (int i=0;i<posSeqCount;i++){
			KmerGroup[] kgs = query(seqs[i]);
			if (kgs.length==0)
				posSeqScores[i]=0;
			else{
				Arrays.sort(kgs);
				posSeqScores[i]=-kgs[0].getHgp();				// score = -log10 hgp, becomes positive value
			}
		}
		Arrays.sort(posSeqScores);
		for (int i=0;i<negSeqCount;i++){
			KmerGroup[] kgs = query(seqsNegList.get(i));
			if (kgs.length==0)
				negSeqScores[i]=0;
			else{
				Arrays.sort(kgs);
				negSeqScores[i]=-kgs[0].getHgp();
			}
		}
		Arrays.sort(negSeqScores);
		
		// find the threshold motif score
		TreeSet<Double> posScoreUnique = new TreeSet<Double>();
		for (double s:posSeqScores)
			posScoreUnique.add(s);
		Double[] posScores_u = new Double[posScoreUnique.size()];
		posScoreUnique.toArray(posScores_u);
		int[] poshits = new int[posScoreUnique.size()];
		int[] neghits = new int[posScoreUnique.size()];
		double[] hgps = new double[posScoreUnique.size()];
		StringBuilder sb = new StringBuilder();
		for (int i=0;i<posScores_u.length;i++){
			double key = posScores_u[i];
			int index = CommonUtils.findKey(posSeqScores, key);
			poshits[i] = posSeqScores.length-index;
			index = CommonUtils.findKey(negSeqScores, key);
			neghits[i] = negSeqScores.length-index;
		}

		int numThread = Math.min(config.maxThreads, java.lang.Runtime.getRuntime().availableProcessors());
		Thread[] threads = new Thread[numThread];
		ArrayList<Integer> idxs = new ArrayList<Integer>();
		for (int i=posScores_u.length-1;i>=0;i--)
			if (poshits[i]>neghits[i]*2.0*posSeqCount/negSeqCount)	// posHit should be at least 2 fold
				idxs.add(i);
		for (int i=0;i<numThread;i++){
            Thread t = new Thread(new HGPThread(idxs, posSeqCount, negSeqCount, poshits, neghits, hgps));
            t.start();
            threads[i] = t;
		}
		boolean anyrunning = true;
        while (anyrunning) {
            anyrunning = false;
            try {
                Thread.sleep(100);
            } catch (InterruptedException e) { }
            for (int i = 0; i < threads.length; i++) {
                if (threads[i].isAlive()) {
                    anyrunning = true;
                    break;
                }
            }   
        }
		hgps[0]=0;		// the lowest threshold will match all positive sequences, lead to hgp=0
		
		if (printKgcHgp){
			for (int i=0;i<posScores_u.length;i++)
				sb.append(String.format("%d\t%.3f\t%d\t%d\t%.1f\n", i, posScores_u[i], poshits[i], neghits[i], hgps[i] ));
		}
		Pair<Double, TreeSet<Integer>> minHgp = StatUtil.findMin(hgps);
		int minIdx = minHgp.cdr().last();
		MotifThreshold score = new MotifThreshold();
		score.score = posScores_u[minIdx];
		score.hgp = minHgp.car();
		score.posHit = poshits[minIdx];
		score.negHit = neghits[minIdx];
		if (printKgcHgp)
			CommonUtils.writeFile(outName+"_KgsHgp.txt", sb.toString());
		return score;
	}	
	/**
	 * Check if a k-mer is a k-mer that was not enriched in positive sets
	 */
	public boolean isNegativeKmer(String kmerStr){
		// is kmer in the negative k-mer set
		Iterator found = tree_negatives.search(kmerStr.getBytes());
		if (found.hasNext())
			return true;
		// try reverse compliment
		found = tree_negatives.search(SequenceUtils.reverseComplement(kmerStr).getBytes());
		if (found.hasNext())
			return true;
		return false;
	}
	
	/** load Kmers and prepare the search Engine, print k-mer list<br>
	 *  assuming the kmers are unique
	 * 
	 * @param kmers List of kmers (with kmerString, sequence hit count)
	 */
	public void updateEngine(ArrayList<Kmer> kmers, String outPrefix){
		if (kmers.isEmpty()){
			engineInitialized = false;
			return;
		}
		Collections.sort(kmers);
//		Kmer.printKmers(kmers, posSeqCount, negSeqCount, outPrefix, false, true);
		
		//Aho-Corasick for searching Kmers in sequences
		//ahocorasick_java-1.1.tar.gz is an implementation of Aho-Corasick automata for Java. BSD license.
		//from <http://hkn.eecs.berkeley.edu/~dyoo/java/index.html> 
		tree = new AhoCorasick();
		str2kmer.clear();
		for (Kmer km: kmers){
			str2kmer.put(km.getKmerString(), km);
			tree.add(km.getKmerString().getBytes(), km.getKmerString());
	    }
	    tree.prepare();
	    engineInitialized = true;
	}	
	/** load Kmers and prepare the search Engine, do not print k-mer list<br>
	 *  assuming the kmers are unique
	 * 
	 * @param kmers List of kmers (with kmerString, sequence hit count)
	 */
	public void updateEngine(ArrayList<Kmer> kmers){
		if (kmers.isEmpty()){
			engineInitialized = false;
			return;
		}		
		
		//Init Aho-Corasick alg. for searching multiple Kmers in sequences
		//ahocorasick_java-1.1.tar.gz is an implementation of Aho-Corasick automata for Java. BSD license.
		//from <http://hkn.eecs.berkeley.edu/~dyoo/java/index.html> 
		tree = new AhoCorasick();
		str2kmer.clear();
		for (Kmer km: kmers){
			int offset = km.getKmerStartOffset();
			String kmerStr;
			if (offset>RC/2){
				offset-=RC;
				kmerStr=km.getKmerRC();
				km.setKmerStartOffset(offset);
			}
			else
				kmerStr = km.getKmerString();
			str2kmer.put(kmerStr, km);
			tree.add(kmerStr.getBytes(), kmerStr);
	    }
	    tree.prepare();
	    engineInitialized = true;
	}		
	
	/** 
	 * Search all k-mers in the sequence
	 * @param seq sequence string to search k-mers, config.strand_type determines whether to search on RC strand
	 * @return an array of KmerGroups:<br>
	 * Each k-mer group maps to a binding position in the sequence
	 * Note: matches on negative strand are combined with matches on positive strand
	 */
	public KmerGroup[] query (String seq){
		seq = seq.toUpperCase();
		HashSet<Object> kmerFound = new HashSet<Object>();	// each kmer is only used 
		//Search for all kmers in the sequences using Aho-Corasick algorithms (initialized)
		//ahocorasick_java-1.1.tar.gz is an implementation of Aho-Corasick automata for Java. BSD license.
		//from <http://hkn.eecs.berkeley.edu/~dyoo/java/index.html> 

		Iterator searcher = tree.search(seq.getBytes());
		while (searcher.hasNext()) {
			SearchResult result = (SearchResult) searcher.next();
			kmerFound.addAll(result.getOutputs());
		}
		if (config.strand_type!=1){
			// the reverse compliment
			String seq_rc = SequenceUtils.reverseComplement(seq);
			searcher = tree.search(seq_rc.getBytes());
			while (searcher.hasNext()) {
				SearchResult result = (SearchResult) searcher.next();
				kmerFound.addAll(result.getOutputs());
			}
		}
		
		// Aho-Corasick only gives the patterns (kmers) matched, need to search for positions
		// matches on negative strand are combined with matches on positive strand
		HashMap<Integer, ArrayList<Kmer>> result = new HashMap<Integer, ArrayList<Kmer>> ();
		String seqRC = SequenceUtils.reverseComplement(seq);
		for (Object o: kmerFound){
			String kmerStr = (String) o;
			Kmer kmer = str2kmer.get(kmerStr);
			ArrayList<Integer> pos = StringUtils.findAllOccurences(seq, kmerStr);
			for (int p: pos){
				int x = p-kmer.getKmerStartOffset();	// minus kmerShift to get the motif position
				if (!result.containsKey(x))
					result.put(x, new ArrayList<Kmer>());
				result.get(x).add(kmer);	
			}
			if (config.strand_type!=1){
				ArrayList<Integer> pos_rc = StringUtils.findAllOccurences(seqRC, kmerStr);
				for (int p: pos_rc){
					int x = p-kmer.getKmerStartOffset();	// motif position in seqRC
					x = seq.length()-1-x;		// convert to position in Seq
					if (!result.containsKey(x))
						result.put(x, new ArrayList<Kmer>());
					result.get(x).add(kmer);	
				}
			}
		}
		KmerGroup[] matches = new KmerGroup[result.keySet().size()];
		int idx = 0;
		for (int p:result.keySet()){
			KmerGroup kg = config.use_weighted_kmer ? new KmerGroup(result.get(p), p, seq_weights, posSeqCount, negSeqCount) : new KmerGroup(result.get(p), p, posSeqCount, negSeqCount);
			matches[idx]=kg;
			kg.setHgp(computeHGP(kg.getGroupHitCount(), kg.getGroupNegHitCount()));
			idx++;
		}
		Arrays.sort(matches);
		return matches;
	}
	
	/** 
	 * Search all k-mers (loaded in the AhoCorasick tree) in the sequence, both strand
	 * @param seq sequence string to search k-mers
	 * @return a set of kmers found
	 */
	public static HashSet<Kmer> queryTree (String seq, AhoCorasick tree, boolean isForwardOnly){
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
		if (!isForwardOnly){
			String seq_rc = SequenceUtils.reverseComplement(seq);
			searcher = tree.search(seq_rc.getBytes());
			while (searcher.hasNext()) {
				SearchResult result = (SearchResult) searcher.next();
				kmerFound.addAll(result.getOutputs());
			}
		}
		
		HashSet<Kmer> result = new HashSet<Kmer>();
		for (Object km: kmerFound)
			result.add((Kmer)km);
		return result;
	}

	/** 
	 * Search all k-mers in the sequence
	 * @param seq sequence string to search k-mers
	 * @return an array of KmerGroups:<br>
	 * Each k-mer group maps to a binding position in the sequence<br>
	 * Note: the return value is different from query() in that match on RC strand is labeled (pos+RC)
	 */
	public KmerGroup[] queryS (String seq){
		seq = seq.toUpperCase();
		HashSet<Object> kmerFound = new HashSet<Object>();	// each kmer is only used 
		//Search for all kmers in the sequences using Aho-Corasick algorithms (initialized)
		//ahocorasick_java-1.1.tar.gz is an implementation of Aho-Corasick automata for Java. BSD license.
		//from <http://hkn.eecs.berkeley.edu/~dyoo/java/index.html> 

		Iterator searcher = tree.search(seq.getBytes());
		while (searcher.hasNext()) {
			SearchResult result = (SearchResult) searcher.next();
			kmerFound.addAll(result.getOutputs());
		}
		if (config.strand_type!=1){
			// the reverse compliment
			String seq_rc = SequenceUtils.reverseComplement(seq);
			searcher = tree.search(seq_rc.getBytes());
			while (searcher.hasNext()) {
				SearchResult result = (SearchResult) searcher.next();
				kmerFound.addAll(result.getOutputs());
			}
		}
		// Aho-Corasick only gives the patterns (kmers) matched, need to search for positions
		// matches on negative strand are combined with matches on positive strand
		TreeMap<Integer, ArrayList<Kmer>> result = new TreeMap<Integer, ArrayList<Kmer>> ();
		String seqRC = SequenceUtils.reverseComplement(seq);
		for (Object o: kmerFound){
			String kmerStr = (String) o;
			Kmer kmer = str2kmer.get(kmerStr);
			ArrayList<Integer> pos = StringUtils.findAllOccurences(seq, kmerStr);
			for (int p: pos){
				int x = p-kmer.getKmerStartOffset();	// minus kmerShift to get the motif position
				if (!result.containsKey(x))
					result.put(x, new ArrayList<Kmer>());
				result.get(x).add(kmer);	
			}
			if (config.strand_type!=1){
				ArrayList<Integer> pos_rc = StringUtils.findAllOccurences(seqRC, kmerStr);
				for (int p: pos_rc){
					int x = p-kmer.getKmerStartOffset();	// motif position in seqRC
					x += RC;									// label it as "found on RC"
					if (!result.containsKey(x))
						result.put(x, new ArrayList<Kmer>());
					result.get(x).add(kmer);	
				}
			}
		}
		KmerGroup[] matches = new KmerGroup[result.keySet().size()];
		int idx = 0;
		for (int p:result.keySet()){
			KmerGroup kg = config.use_weighted_kmer ? new KmerGroup(result.get(p), p, seq_weights, posSeqCount, negSeqCount) : new KmerGroup(result.get(p), p, posSeqCount, negSeqCount);
			matches[idx]=kg;
			kg.setHgp(computeHGP(kg.getGroupHitCount(), kg.getGroupNegHitCount()));
			idx++;
		}
		return matches;
	}
		
	public String getSequenceUppercase(Region r){
		return seqgen.execute(r).toUpperCase();
	}
	
	public void indexKmers(List<File> files){
		long tic = System.currentTimeMillis();
		ArrayList<Kmer> kmers = Kmer.loadKmers(files);
		if (kmers.isEmpty())
			return;
		int step = 100000000;
		this.k = kmers.get(0).getK();
		HashMap<String, Integer> map = new HashMap<String, Integer>();
		for (Kmer kmer:kmers){
			map.put(kmer.getKmerString(), 0);
		}
		for (String chr : genome.getChromList()){
			System.out.println(chr);
			int chrLen = genome.getChromLengthMap().get(chr);
			for (int l=0;l<chrLen;l+=step-k+2){		// the step size is set so that the overlap is k-1
				int end = Math.min(l+step-1, chrLen-1);
				String seq = seqgen.execute(new Region(genome, chr, l, end)).toUpperCase();
				for (int i=0;i<seq.length()-k;i++){
					String s = seq.substring(i, i+k);
					if (map.containsKey(s)){			// only count known kmers, save memory space
						 map.put(s, (map.get(s)+1));
					}
					else {		// try the other strand
						String rc = SequenceUtils.reverseComplement(s);
						if (map.containsKey(rc)){			
							 map.put(rc, (map.get(rc)+1));
						}
					}						
				}
			}
		}
		Collections.sort(kmers);
		StringBuilder sb = new StringBuilder();
		for (Kmer km:kmers){
			sb.append(km.getKmerString()).append("\t").append(map.get(km.getKmerString())).append("\n");
		}
		CommonUtils.writeFile(genome.getVersion()+"_kmers_"+k+".txt", sb.toString());
		System.out.println(CommonUtils.timeElapsed(tic));
	}

	private ArrayList<Kmer> getMismatchKmers(ArrayList<Kmer> aligned, ArrayList<Kmer> candidates){
			ArrayList<Kmer> mmaligned = new ArrayList<Kmer>();			// kmers aligned by using mismatch
			Collections.sort(aligned);
			for (Kmer km: candidates){
				String seq = km.getKmerString();
				for (Kmer akm: aligned){
	//				if (Math.abs(akm.getShift())>config.k/2)		// the mismatch must be proximal to seed kmer
					if (km.getPosHitCount() > akm.getPosHitCount())		// the mismatch kmer should have less count than the reference kmer, avoiding a bad weak kmer misguiding a strong kmer
						continue;
					if (CommonUtils.mismatch(akm.getKmerString(), seq)==1){
						km.setShift(akm.getShift());
						km.setAlignString("MM:"+akm.getKmerString());
						mmaligned.add(km);
						break;
					}
					if (CommonUtils.mismatch(akm.getKmerRC(), seq)==1){
						int shift = akm.getShift();
						km.setShift(shift>RC/2?shift:shift-RC);
						km.setAlignString("MM:"+akm.getKmerString());
						mmaligned.add(km);
						break;
					}
				}
			}
			return mmaligned;
		}
	private boolean mergeClusters(double setOverlapRatio){
		int originalCount = clusters.size();
		HashSet<Integer> merged = new HashSet<Integer>();
		for (int i=0; i<originalCount-1; i++){
			for (int j=i+1; j<originalCount; j++){
				if (clusters.get(j).alignedKmers.size()<2){		// if a cluster contains less than 2 kmers, it is removed
					merged.add(j);	
					continue;
				}
				if (merged.contains(j))
					continue;
				SetTools<Kmer> setTools = new SetTools<Kmer>();
				Set<Kmer> first = new HashSet<Kmer>();
				first.addAll(clusters.get(i).alignedKmers);
				Set<Kmer> second = new HashSet<Kmer>();
				second.addAll(clusters.get(j).alignedKmers);
				Set<Kmer> common = setTools.intersection(first, second);
				if (common.size()>Math.min(first.size(), second.size())*setOverlapRatio){
					Set<Kmer>extra = setTools.subtract(second, common);
					clusters.get(i).alignedKmers.addAll(extra);
					merged.add(j);
				}
			}
		}
		Integer[] ids = new Integer[merged.size()];
		merged.toArray(ids);
		Arrays.sort(ids);
		for (int i=ids.length-1;i>=0;i--){
			clusters.remove(ids[i].intValue());
		}
		clusters.trimToSize();
		
		if (verbose>1)
			System.out.println(CommonUtils.timeElapsed(tic)+": Merge "+originalCount+" clusters to "+clusters.size()+" clusters.");
	
		return !merged.isEmpty();
	}


	class HGPThread implements Runnable {
		ArrayList<Integer> idxs;
		int posTotal;
		int negTotal;
		int[] posHits;
		int[] negHits;
		double[] hgps;
		HGPThread(ArrayList<Integer> idxs, int posTotal, int negTotal, int[] posHits, int[] negHits, double[] hgps){
			this.idxs = idxs;
			this.posTotal = posTotal;
			this.negTotal = negTotal;
			this.posHits = posHits;
			this.negHits = negHits;
			this.hgps = hgps;
		}
		public void run() {
			int i;
			while (!idxs.isEmpty()) {
				synchronized (idxs){
	            	if (!idxs.isEmpty()){
	            		i = idxs.get(0);
	            		idxs.remove(0);
	        		}
	            	else
	            		break;
	        	}		
				if (i==0 && posHits[0]==posSeqCount)
	    			hgps[0]=0;
				hgps[i]=computeHGP(posTotal, negTotal, posHits[i], negHits[i]);
			}
		}
	}
	public static void main1(String[] args){
		ArrayList<Integer> x_list = new ArrayList<Integer>();
		ArrayList<Integer> same_list = new ArrayList<Integer>();
		ArrayList<Integer> diff_list = new ArrayList<Integer>();
		try {	
			BufferedReader bin = new BufferedReader(new InputStreamReader(new FileInputStream(new File(args[0]))));
	        String line;
	        String[]f = null;
	        while((line = bin.readLine()) != null) { 
	            f = line.trim().split("\t");
	            if (f.length==3){
	            	x_list.add(Integer.parseInt(f[0]));
	            	same_list.add(Integer.parseInt(f[1]));
	            	diff_list.add(Integer.parseInt(f[2]));
	            }
	        }			
	        if (bin != null) {
	            bin.close();
	        }
	    } catch (IOException e) {
	    	System.err.println("Error when processing "+args[0]);
	        e.printStackTrace(System.err);
	    }
	    int[] x = new int[x_list.size()];
	    int[] same = new int[x_list.size()];
	    int[] diff = new int[x_list.size()];
	    for (int i=0;i<x.length;i++)
	    	x[i] = x_list.get(i);
	    for (int i=0;i<x.length;i++)
	    	same[i] = same_list.get(i);
	    for (int i=0;i<x.length;i++)
	    	diff[i] = diff_list.get(i);
//	    KMAC kmf = new KMAC();
//	    kmf.plotMotifDistanceDistribution(x, same, diff, args[0]+".png");
	}
	
	// options cMyc_cMyc cMyc_HeLa_61bp_GEM.fasta cMyc_HeLa_61bp_GEM_neg.fasta 5 8 CCACGTG
	public static void main(String[] args){
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
		
		String name = Args.parseString(args, "out_name", null);
		File outFolder = new File(name+"_outputs");
		outFolder.mkdir();
		name = new File(outFolder, name).getAbsolutePath();
		
		String pos_file = Args.parseString(args, "pos_seq", null);
		String neg_file = Args.parseString(args, "neg_seq", null);
		if (pos_file==null || (config.k==-1&&config.k_min==-1)){
			System.err.println("Example: KMAC --pos_seq c-Myc_Crawford_HeLa-S3_61bp_GEM.fasta [--neg_seq c-Myc_Crawford_HeLa-S3_61bp_GEM_neg.fasta] --k_min 5 --k_max 8 --out_name cMyc_cMyc --seed CACGTG");
			System.exit(-1);
		}
		if (config.seed!=null){
			config.k=config.seed.length();
			System.out.println("Starting seed k-mer is "+config.seed+".\n");
		}
		String format = Args.parseString(args, "format", "fasta");
		ArrayList<String> strs = CommonUtils.readTextFile(pos_file);
        String[]f = null;
		for (String line: strs){
			if (format.equals("fasta")){
	            if (line.startsWith(">")){
	        		f = line.split(" ");
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
	        		int left = line.length()/2-config.k_win/2;
	        		if (left<0)
	        			continue;
	        		pos_seqs.add(line.toUpperCase().substring(left, left+config.k_win));
	        	}
			}
			else{
				f = line.split("\t");
        		if (f.length>1){
	            	seq_w.add(Double.parseDouble(f[1]));
	            	pos_seqs.add(f[0].toUpperCase());
	            }
			}
		}
		
		ArrayList<String> neg_seqs = new ArrayList<String>();
		if (neg_file!=null){
			strs = CommonUtils.readTextFile(neg_file);
			for (String line: strs){
				if (format.equals("fasta")){
					if (!line.startsWith(">")){
		        		int left = line.length()/2-config.k_win/2;
		        		if (left<0)
		        			continue;
		        		neg_seqs.add(line.toUpperCase().substring(left, left+config.k_win));
	        		}
				}
				else{
		            	neg_seqs.add(line.substring(0,config.k_win).toUpperCase());
				}
			}
		}
        
        KMAC kmf = new KMAC(config, name);
        kmf.setSequences(pos_seqs, neg_seqs, seq_w);
        kmf.setStandalone();
		
		if (config.strand_type==1)
			System.out.println("Running single-strand KMAC motif discovery ...\n");
		
        System.out.println(String.format("%d input positive sequences, use top %d center sequences (%dbp) to find motif ...", pos_seqs.size(), kmf.seqs.length, config.k_win));
        if (config.k==-1){
        	if (config.selectK_byTopKmer)
        		config.k = kmf.selectK_byTopKmer(config.k_min, config.k_max, null);
        	else
        		config.k = kmf.selectK(config.k_min, config.k_max, null);
        }
        ArrayList<Kmer>kmers = kmf.selectEnrichedKmers(config.k);
        StringBuilder sb=new StringBuilder();
		for (String arg:args){
			if (arg.trim().indexOf(" ")!=-1)
				sb.append("\"").append(arg).append("\" ");
			else
				sb.append(arg).append(" ");
		}
        kmf.KmerMotifAlignmentClustering(kmers, -1, false, null, sb.toString());
	}
}


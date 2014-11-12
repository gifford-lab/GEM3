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
import java.util.*;

import javax.imageio.ImageIO;

import com.ctc.wstx.util.StringUtil;

import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.SetTools;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.utils.strings.StringUtils;
import edu.mit.csail.cgs.utils.strings.multipattern.*;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.features.ComponentFeature;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.utils.stats.StatUtil;
import edu.mit.csail.cgs.utils.stats.StatUtil.DensityClusteringPoint;
import edu.mit.csail.cgs.deepseq.discovery.Config;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KMAC2.KmerCluster;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KMAC2.KmerGroup;

public class KMAC2WK {
	private static final int RC=100000;		// extra bp add to indicate negative strand match of kmer
	private static final int UNALIGNED=9999;	// the special shift for unaligned kmer
	public static final char[] LETTERS = {'A','C','G','T'};
	public static final int MAXLETTERVAL = Math.max(Math.max(Math.max('A','C'),Math.max('T','G')),
            Math.max(Math.max('a','c'),Math.max('t','g'))) + 1;
	
	private boolean standalone = false;
	public void setStandalone(){
		standalone = true;
	}
	Config config = new Config();	

	private Genome genome;
	private boolean engineInitialized =false;
	private int k;
	private int minHitCount = 2;
	private int numPos;
	private double[] bg= new double[4];	// background frequency based on GC content
	private double ic_trim = 0.4;
	private String outName;
	private boolean use_smart_mm = false;	

	
	private int seqLen;
	private double[] profile;
	private ArrayList<Sequence> seqList;
	private ArrayList<Sequence> seqListNeg;
	private String[] seqs;			// DNA sequences around binding sites
	private double[] seq_weights;	// sequence hit count weighted by binding strength, corresponding to the seqs[]
	private double totalWeight;
	private String[] seqsNeg;		// DNA sequences in negative sets
	private ArrayList<String> seqsNegList=new ArrayList<String>(); // Effective negative sets, excluding overlaps in positive sets
	public int getNegSeqCount(){return negSeqCount;}
    private int posSeqCount;
    private int negSeqCount;
    public void setTotalSeqCount(int pos, int neg){
    	posSeqCount = pos;
    	negSeqCount = neg;
    }
	private int negRegionDistance;
	/** region-->index for negative sequences */
	private TreeMap<Region, Integer> neg_region_map;
	public String[] getPositiveSeqs(){return seqs;};
	public double get_NP_ratio(){return (double)negSeqCount/posSeqCount;}
	
	private ArrayList<ArrayList<Kmer>> forward = new ArrayList<ArrayList<Kmer>>();
	private ArrayList<ArrayList<Kmer>> reverse = new ArrayList<ArrayList<Kmer>>();
	private int negPositionPadding = 0;

	private TreeMap<Integer, ArrayList<Kmer>> k2kmers = new TreeMap<Integer, ArrayList<Kmer>>();

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
	public KMAC2WK(){
	}

	public void setTotalSeqCounts(int posSeqCount, int negSeqCount){
		this.posSeqCount = posSeqCount;
		this.negSeqCount = negSeqCount;
	}
	
	public void setConfig(Config config, String outName){
		this.config = config;
	    this.outName = outName;
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
		
		seq_weights = new double[seqNum];
		totalWeight=0;
		for (int i=0;i<seqNum;i++){
			if (config.weight_by_sqrt_strength)
				seq_weights[i]=Math.sqrt(pos_w.get(i));
			else
				seq_weights[i]=pos_w.get(i);
			totalWeight += seq_weights[i];
		}
		for (int i=0;i<seq_weights.length;i++){
			seq_weights[i] = seq_weights[i]*seqs.length/totalWeight;	// scale weights with total sequence count, and total weight
		}
		if (config.use_weighted_kmer)
			Kmer.set_seq_weights(seq_weights);
		
		if (config.k_neg_shuffle){
			System.out.println("Use shuffled sequences as negative sequences.\n");
			Random randObj = new Random(config.rand_seed);
			for (int i=0;i<seqNum;i++)
				seqsNegList.add(SequenceUtils.shuffle(seqs[i], randObj));
		}
		else if (config.k_neg_dinu_shuffle){
			System.out.println("Use di-nucleotide shuffled sequences as negative sequences.\n");
			Random randObj = new Random(config.rand_seed);
			for (int i=0;i<seqNum;i++)
				seqsNegList.add(SequenceUtils.dinu_shuffle(seqs[i], randObj));
		}
		else{
			if (neg_seqs.size()<seqNum)
				seqNum = neg_seqs.size();
			for (int i=0;i<seqNum;i++){
				String seq = neg_seqs.get(i);
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
				seqsNegList.add(seq.toUpperCase());
			}

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
		seqLen = seqs[0].length();
		
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
	   	
	    // count cg-content
		int gcCount = 0;
		for (String seq:seqsNegList){
			for (char c:seq.toCharArray())
				if (c=='C'||c=='G')
					gcCount ++;
		}
		double gcRatio = (double)gcCount/negSeqCount/seqLen;
		bg[0]=0.5-gcRatio/2; 
    	bg[1]=gcRatio/2; 
    	bg[2]=bg[1]; 
    	bg[3]=bg[0];
	}
	
	public KMAC2WK(Genome g, boolean useCache, boolean use_db_genome, String genomePath){
//		setUseKmerWeight();

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
	public KMAC2WK(ArrayList<Kmer> kmers, String outPrefix){
		if (!kmers.isEmpty()){
			if (outPrefix!=null)
				updateEngine(kmers, outPrefix);
			else
				updateEngine(kmers);
			k=kmers.get(0).getK();
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
				posSeqs.add(seq.toUpperCase());		// if repeat_fraction>=1, allow all repeats, convert to upper case
				if (config.weight_by_sqrt_strength)
					posSeqWeights.add(Math.sqrt(events.get(i).getTotalEventStrength()));	// use sq root of event strength as weight
				else
					posSeqWeights.add(events.get(i).getTotalEventStrength());		// use event strength as weight
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
			if (config.k_neg_shuffle){
				System.out.println("Use shuffled sequences as negative sequences.\n");
				Random randObj = new Random(config.rand_seed);
				for (int i=0;i<seqs.length;i++)
					seqsNegList.add(SequenceUtils.shuffle(seqs[i], randObj));
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
		
		ArrayList<KmerCluster> allClusters = new ArrayList<KmerCluster>();
		
		/** Initialization of the sequences */
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
		
		// init list to keep track of kmer matches in a sequence
		negPositionPadding = k_max*2;
		for (int i=0;i<config.k_win+negPositionPadding*2;i++){
			forward.add(new ArrayList<Kmer>());
			reverse.add(new ArrayList<Kmer>());
		}
		forward.trimToSize();
		reverse.trimToSize();
		
		for (int i=0;i<k_max-k_min+1;i++){
			int k = i+k_min;
			StringBuilder sb = new StringBuilder("\n------------------------- k = "+ k +" ----------------------------\n");
			System.out.println("\n----------------------------------------------------------\nTrying k="+k+" ...\n");
			ArrayList<Kmer> kmers = generateEnrichedKmers(k);
			if (config.use_gapped)
				kmers.addAll(generateEnrichedGappedKmers(k,1));
			
			/** Index sequences and kmers, for each k, need to index this only once */
			indexKmerSequences(kmers, seqList, seqListNeg, config.kmer_hgp);
			Collections.sort(kmers);
			k2kmers.put(k, kmers);

	        ArrayList<KmerCluster> tmp = new ArrayList<KmerCluster>();
			if (!kmers.isEmpty()){
				double[][]distanceMatrix = computeDistanceMatrix(kmers);
				ArrayList<Kmer> centerKmers = selectCenterKmersByDensityClustering(kmers, distanceMatrix);
		        for (int j=0;j<centerKmers.size();j++){	
		        	Kmer seedKmer = centerKmers.get(j);
		    		if (config.verbose>1)
		    			System.out.println("------------------------------------------------\n"+
		    					"Aligning k-mers with "+seedKmer.getKmerString()+",   #"+j);

		        	KmerCluster c=KmerMotifAlignmentClustering(kmers, seedKmer);
		        	if (c!=null){
			        	c.k = k;
		        		tmp.add(c);
		        	}
		        	if (tmp.size()==config.k_top)
		        		break;
		        }
			}
			if (tmp.size()<=1){
				clusters = tmp;
			}
			else{
				// sort clusters, set clusterid
				if (config.evaluate_by_ksm){
					Collections.sort(tmp, new Comparator<KmerCluster>() {
		                public int compare(KmerCluster o1, KmerCluster o2) {
		                    return o1.compareToByKsmHGP(o2);
		                }
		            });
				}
				else{
					Collections.sort(tmp);
				}
				for (int j=0;j<tmp.size();j++){
					tmp.get(j).clusterId = j;
				}
				// print all the motifs for k, before merging
				for (KmerCluster c:tmp){
					if (config.evaluate_by_ksm){
						sb.append(String.format("k=%d\thit=%d\thgp=1e%.1f\tTopKmer=%s\n", k, c.ksmThreshold.posHit, 
								c.ksmThreshold.hgp, c.seedKmer.getKmerString()));
					}
					else if (c.wm!=null){
						sb.append(String.format("k=%d\thit=%d\thgp=1e%.1f\tW=%d\tPWM=%s.\n", k, c.pwmPosHitCount, 
								c.pwmThresholdHGP, c.wm.length(), WeightMatrix.getMaxLetters(c.wm)));
					}
				}
				sb.append("\n");
				clusters = tmp;
	
				if (config.verbose>1)
					System.out.println("\nMerge overlapping motif clusters ...\n");
				boolean[][] checked = new boolean[clusters.size()][clusters.size()];	// a table indicating whether a pair has been checked
				
				tic = System.currentTimeMillis();
				mergeOverlapClusters (seqList, checked, config.use_ksm, 0);
			}
			// print motifs for k, after merging
			for (KmerCluster c:clusters){
				if (config.evaluate_by_ksm){
					sb.append(String.format("k=%d\thit=%d+/%d-\thgp=1e%.1f\tTopKmer=%s\n", k, c.ksmThreshold.posHit, c.ksmThreshold.negHit,
							c.ksmThreshold.hgp, c.seedKmer.getKmerString()));
				}
				else if (c.wm!=null){
					sb.append(String.format("k=%d\thit=%d+/%d-\thgp=1e%.1f\tW=%d\tPWM=%s.\n", k, c.pwmPosHitCount, c.pwmNegHitCount,
							c.pwmThresholdHGP, c.wm.length(), WeightMatrix.getMaxLetters(c.wm)));
				}
			}
//			sb.append("\n-----------------------------------------------------------\n");
			System.out.println(sb.toString());
			
			allClusters.addAll(clusters);
		} // for each k
		
		// merge motifs from all K values
		// sort clusters, set clusterid
		Collections.sort(allClusters);
		for (int j=0;j<allClusters.size();j++){
			allClusters.get(j).clusterId = j;
		}
		clusters = allClusters;
		
		if (config.verbose>1)
			System.out.println("\n---------------------------------------------\nMerge all overlapping motif clusters ...\n");

		StringBuilder sb_all = new StringBuilder("\n------------------------- "+ new File(outName).getName() +" ----------------------------\n");
		for (KmerCluster c:clusters){
			if (config.evaluate_by_ksm){
				sb_all.append(String.format("k=%d\thit=%d+/%d-\thgp=1e%.1f\tTopKmer=%s\n", c.k, c.ksmThreshold.posHit, c.ksmThreshold.negHit, 
						c.ksmThreshold.hgp, c.seedKmer.getKmerString()));
			}
			else if (c.wm!=null){
				sb_all.append(String.format("k=%d\thit=%d+/%d-\thgp=1e%.1f\tW=%d\tPWM=%s.\n", c.k, c.pwmPosHitCount,  c.pwmNegHitCount,
						c.pwmThresholdHGP, c.wm.length(), WeightMatrix.getMaxLetters(c.wm)));
			}
//			if (c.wm!=null){
//				pfm_sb.append(CommonUtils.makeTRANSFAC (c.pfm, c.pwmPosHitCount, 
//						String.format("DE %s_%d_k%d_c%d", config.out_name, c.clusterId, c.k, c.pwmPosHitCount)));
//			}
		}
		sb_all.append("\n");
		System.out.print(sb_all.toString());
		
		tic = System.currentTimeMillis();
		boolean[][] checked = new boolean[clusters.size()][clusters.size()];	// a table indicating whether a pair has been checked
		mergeOverlapClusters (seqList, checked, config.use_ksm, 0);
		if (config.verbose>1)
			System.out.println("\nFinished merging motifs.\n");

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
		
		if (config.evaluate_by_ksm){
			Collections.sort(clusters, new Comparator<KmerCluster>() {
                public int compare(KmerCluster o1, KmerCluster o2) {
                    return o1.compareToByKsmHGP(o2);
                }
            });
		}
		else{
			Collections.sort(clusters);
		}
		for (int j=0;j<clusters.size();j++){
			clusters.get(j).clusterId = j;
		}
		
		// Extract final KSM from PWM alignment and output KSM (need to output it here, o.w. the k-mers may be overwritten by other motifs
		// TODO: if no PWM, it is the seed_family alignment
		for (KmerCluster cluster: clusters){
			if (cluster.wm!=null){
				indexKmerSequences(k2kmers.get(cluster.k), seqList, seqListNeg, config.kmer_hgp);  // need this to get KSM
	    		alignSequencesUsingPWM(seqList, cluster);
	    		if (cluster.pwmPosHitCount>cluster.total_aligned_seqs)
	    			cluster.total_aligned_seqs = cluster.pwmPosHitCount;
		    	NewKSM newKSM = extractKSM (seqList, cluster.seedKmer.getK(), new ArrayList<Kmer>());
		    	if (newKSM==null)
		    		continue;
				if (newKSM.threshold!=null && newKSM.threshold.hgp <= cluster.ksmThreshold.hgp){
					cluster.alignedKmers = newKSM.kmers;
			    	cluster.ksmThreshold = newKSM.threshold;
				}
				GappedKmer.printGappedKmers(cluster.alignedKmers, cluster.k, posSeqCount, negSeqCount, cluster.ksmThreshold.score, outName+".m"+cluster.clusterId, false, true, false);
	    	}
		}
		
		// show all the motifs for k, after merging
		sb_all = new StringBuilder("\n------------------------- "+ new File(outName).getName() +" ----------------------------\n");
//		StringBuilder pfm_sb = new StringBuilder();
		for (KmerCluster c:clusters){
			if (config.evaluate_by_ksm){
				sb_all.append(String.format("k=%d\thit=%d+/%d-\thgp=1e%.1f\tTopKmer=%s\n", c.k, c.ksmThreshold.posHit, c.ksmThreshold.negHit, 
						c.ksmThreshold.hgp, c.seedKmer.getKmerString()));
			}
			else if (c.wm!=null){
				sb_all.append(String.format("k=%d\thit=%d+/%d-\thgp=1e%.1f\tW=%d\tPWM=%s.\n", c.k, c.pwmPosHitCount,  c.pwmNegHitCount,
						c.pwmThresholdHGP, c.wm.length(), WeightMatrix.getMaxLetters(c.wm)));
			}
//			if (c.wm!=null){
//				pfm_sb.append(CommonUtils.makeTRANSFAC (c.pfm, c.pwmPosHitCount, 
//						String.format("DE %s_%d_k%d_c%d", config.out_name, c.clusterId, c.k, c.pwmPosHitCount)));
//			}
		}
		sb_all.append("\n");
		System.out.print(sb_all.toString());
		
//		CommonUtils.writeFile(outName+".all.PFM.txt", pfm_sb.toString());

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
				html.append("<a href='"+name+"_GEM_events.txt'>Significant Events</a>&nbsp;&nbsp;: "+eventCounts[0]);
				html.append("<br><a href='"+name+"_GEM_insignificant.txt'>Insignificant Events</a>: "+eventCounts[1]);
				html.append("<br><a href='"+name+"_GEM_filtered.txt'>Filtered Events</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;: "+eventCounts[2]);
			}
			html.append("<p><p>Motif can not be found!<p>");
			html.append("</td></tr></table>");
			CommonUtils.writeFile(outName+"_result.htm", html.toString());
			return;
		}
		else{
			// print PWM spatial distribtution
			computeMotifDistanceDistribution(outName);
			outputClusters(eventCounts);
		}
	}
	
	/** 
	 * Index k-mers from the positive sequences, select enriched k-mers
	 * */
	public ArrayList<Kmer> generateEnrichedKmers(int k){
		this.k = k;
		// expected count of kmer = total possible unique occurences of kmer in sequence / total possible kmer sequence permutation
		tic = System.currentTimeMillis();
		numPos = seqLen-k+1;

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
		
		// Merge kmer and its reverse compliment (RC)	
		ArrayList<Kmer> kms = new ArrayList<Kmer>();
		ArrayList<String> kmerStrings = new ArrayList<String>();
		kmerStrings.addAll(kmerstr2seqs.keySet());
		
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
			String seq_rc = SequenceUtils.reverseComplement(seq);
			searcher = tmp.search(seq_rc.getBytes());
			while (searcher.hasNext()) {
				SearchResult result = (SearchResult) searcher.next();
				kmerHits.addAll(result.getOutputs());
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
	 * Select gapped k-mers
	 * */
	public ArrayList<Kmer> generateEnrichedGappedKmers(int k, int numgapped){
		k += numgapped;
		// expected count of kmer = total possible unique occurences of kmer in sequence / total possible kmer sequence permutation
		tic = System.currentTimeMillis();
		numPos = seqLen-k+1;

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
		
		// Merge kmer and its reverse compliment (RC)	
		ArrayList<String> kmerStrings = new ArrayList<String>();
		kmerStrings.addAll(kmerstr2seqs.keySet());
		
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
					kmerstr2seqs.get(winner).addAll(kmerstr2seqs.get(loser));	// winner takes all
					kmerstr2seqs.remove(loser);					// remove the loser kmer because it is represented by its RC
				}
			}
		}
		System.out.println(String.format("k=%d+%d, mapped %d k-mers, %s", k-numgapped, numgapped, 
				kmerstr2seqs.keySet().size(), CommonUtils.timeElapsed(tic)));

		// create the kmer object
		HashMap<String, Kmer> kmMap = new HashMap<String, Kmer>();
		for (String s:kmerstr2seqs.keySet()){	
			Kmer kmer = new Kmer(s, kmerstr2seqs.get(s));
			kmMap.put(kmer.getKmerString(),kmer);
		}
		kmerstr2seqs=null;	// clean up
		System.gc();
		
		// form gapped kmers by mutating each non-edge base of the k-mers
		HashMap<String, GappedKmer> gkMap = new HashMap<String, GappedKmer>();
		for (String s:kmMap.keySet()){
			char[] cs = s.toCharArray();
			if (numgapped==1){
				for (int i=1; i<k-1; i++){	// only mutate the non-edge bases
					char ci = s.charAt(i);
					// check whether this gapped has been made
					cs[i]='N';
					String wkStr = String.valueOf(cs);
					if (gkMap.containsKey(wkStr)||gkMap.containsKey(SequenceUtils.reverseComplement(wkStr)))
						continue;	// this gapped has been made, skip to next position
					GappedKmer gk = new GappedKmer(wkStr);
					for (char c: LETTERS){
						cs[i]=c;
						String m = String.valueOf(cs);
						if (kmMap.containsKey(m))
							gk.addSubKmer(kmMap.get(m), true);
						else{
							String mrc = SequenceUtils.reverseComplement(m);
							if (kmMap.containsKey(mrc)){
								gk.addSubKmer(kmMap.get(mrc), false); 
							}
						}
					}
					cs[i]=ci;		// restore position i					
					if (gk.getSubKmers().size()>1){			// ignore if a singleton longer kmer, b/c it will be covered by a basic k-mer
						gk.mergePosHits();
						gkMap.put(gk.getKmerString(), gk);
					}
				}
			}
		}
		
		// compute the smallest PosCount needed to be significant, with negCount=0
		int smallestPosCount;
		for (smallestPosCount=minHitCount;smallestPosCount<posSeqCount;smallestPosCount++){
			double hgp = computeHGP(posSeqCount, negSeqCount, smallestPosCount, 0);
			if (hgp<config.kmer_hgp){
				break;
			}
		}
		int expectedCount = (int) Math.max( smallestPosCount, Math.round(seqs.length*2*(seqs[0].length()-k+1) / Math.pow(4, k)));
		// the purpose of expectedCount is to limit kmers to run HGP test, if total kmer number is low, it can be relaxed a little
//		if (gkMap.keySet().size()<10000){
//			expectedCount = Math.min(smallestPosCount, expectedCount);
//		}
		ArrayList<GappedKmer> gks = new ArrayList<GappedKmer>();
		HashSet<Kmer> subKmers = new HashSet<Kmer>();
		for (String key:gkMap.keySet()){	
			if (gkMap.get(key).getPosHitCount()< expectedCount)
				continue;	// skip low count kmers 
			GappedKmer gk = gkMap.get(key);
			gks.add(gk);
			subKmers.addAll(gk.getSubKmers());
		}		
		System.out.println("Expected kmer hit count="+expectedCount + ", gapped kmer count="+gks.size()+
				", sub-kmers count="+subKmers.size());
		
		/**
		 * Select significantly over-representative kmers 
		 * Count kmer hits in the negative sequences, then compute hgp using positive/negative hit counts
		 */
		tic = System.currentTimeMillis();
		//Aho-Corasick for searching Kmers in negative sequences
		//ahocorasick_java-1.1.tar.gz is an implementation of Aho-Corasick automata for Java. BSD license.
		//from <http://hkn.eecs.berkeley.edu/~dyoo/java/index.html> 
		AhoCorasick tmp = new AhoCorasick();
		for (Kmer km: subKmers){
			String s = km.getKmerString();
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
		for (Kmer kmer:subKmers){
			if (kmerstr2negSeqs.containsKey(kmer.getKmerString())){
				kmer.setNegHits(kmerstr2negSeqs.get(kmer.getKmerString()));				
			}
		}
		
		// score the kmers, hypergeometric p-value, select significant k-mers
		ArrayList<Kmer> results = new ArrayList<Kmer>();
		for (GappedKmer gk:gks){
			gk.mergeNegHits();		// aggregate sub-kmer negative hits to the WC kmer
			if (gk.getPosHitCount() >= gk.getNegHitCount()/get_NP_ratio() * config.k_fold ){
				// optimize the subkmers to get best hgp
				ArrayList<Kmer> subKmerList = new ArrayList<Kmer>();
				subKmerList.addAll(gk.getSubKmers());
				optimizeKSM(subKmerList);
				HashSet<Kmer> toremove = new HashSet<Kmer>();
				for (Kmer subkm:gk.getSubKmers())
					if (!subKmerList.contains(subkm))
						toremove.add(subkm);
				for (Kmer km: toremove)
					gk.removeSubkmers(km);

				if (gk.getSubKmers().size()<=1)		// skip if only has 1 subkmer
					continue;
				gk.mergePosHits();
				gk.mergeNegHits();
				gk.setHgp(computeHGP(posSeqCount, negSeqCount, gk.getPosHitCount(), gk.getNegHitCount()));
				if (gk.getHgp() <= config.kmer_hgp){		// these GK are significant
					gk.linkSubKmers();
					results.add(gk);
				}
			}
		}
		results.trimToSize();
		Collections.sort(results);
		System.out.println(String.format("k=%d, selected %d gapped k-mers from %d+/%d- sequences, %s", k-numgapped, results.size(), posSeqCount, negSeqCount, CommonUtils.timeElapsed(tic)));
		
		ArrayList<Kmer> kms = new ArrayList<Kmer>();
		for (GappedKmer gk:gks)
			kms.add(gk);
		if (config.print_all_kmers){
			Collections.sort(kms);
			GappedKmer.printGappedKmers(kms, k-numgapped, posSeqCount, negSeqCount, 0, outName+"_g_all_w"+seqs[0].length(), true, false, true);
		}
		
		return results;
	}

	private double[][] computeDistanceMatrix(ArrayList<Kmer> kmers) {
		boolean print_dist_matrix = false;
		int kmerCount = kmers.size();
		double[][]distanceMatrix = new double[kmerCount][kmerCount];
		for (int i=0;i<kmerCount;i++){
			String ks_i = kmers.get(i).getKmerString();
			int cutoff = ks_i.length()/2+1;
//			String krc_i = kmers.get(i).getKmerRC();
			for (int j=0;j<=i;j++){
//				String ks_j = kmers.get(j).getKmerString();
				distanceMatrix[i][j]=CommonUtils.strMinDistanceWithCutoff(ks_i, kmers.get(j).getKmerString(), cutoff);
//				distanceMatrix[i][j]=Math.min(CommonUtils.strMinDistance(ks_i, ks_j), CommonUtils.strMinDistance(krc_i, ks_j));  // much slower
			}
		}
		for (int i=0;i<kmerCount;i++){
			for (int j=i+1;j<kmerCount;j++){
				distanceMatrix[i][j]=distanceMatrix[j][i];
			}
		}
		System.out.println("computeDistanceMatrix, "+CommonUtils.timeElapsed(tic));

		if (print_dist_matrix){
	        StringBuilder output = new StringBuilder();
	        for (int j=0;j<kmerCount;j++){
	        	output.append(String.format("%s\t",kmers.get(j).getKmerString()));
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
		}
		return distanceMatrix;
	}
	
	private ArrayList<Kmer> selectCenterKmersByDensityClustering(ArrayList<Kmer> kmers, double[][]distanceMatrix) {
		long tic = System.currentTimeMillis();
		int k = kmers.get(0).getK();
		
//		double weights[] = new double[kmerCount];
		ArrayList<BitSet> posHitList = new ArrayList<BitSet>();
		ArrayList<BitSet> negHitList = new ArrayList<BitSet>();
		for (Kmer km:kmers){
			posHitList.add(km.posBits);
			negHitList.add(km.negBits);
        }
		
//		ArrayList<DensityClusteringPoint> centers = StatUtil.weightedDensityClustering(distanceMatrix, weights, config.kd, config.delta);
		int dc = k>=8?config.dc+1:config.dc;
		int delta = dc + (k>=11?2:1);
		ArrayList<DensityClusteringPoint> centers = StatUtil.hitWeightedDensityClustering(distanceMatrix, posHitList, negHitList, dc, delta);
		ArrayList<Kmer> results = new ArrayList<Kmer>();
		//TODO: for each cluster, select the strongest kmer that is far away from other stronger cluster center kmers
		for (DensityClusteringPoint p:centers){
			// find the best representative kmer
			Kmer km = kmers.get(p.id);
//			System.out.println(String.format("\n%s    \t%d\t%.1f\t%.1f\t%.1f\t%d", km.toShortString(), p.id, p.density, p.delta, p.gamma, p.members.size()));
			for (int id : p.members){
				if (kmers.get(id).getHgp()<km.getHgp()){ 
//					System.out.println(km.toShortString()+"\t"+kmers.get(id).toShortString());
					if (results.isEmpty()){
						km = kmers.get(id);
					}
					else{
						boolean isNotSimilar = true;
						for (Kmer kmer:results){
//							System.out.println(distanceMatrix[id][p.id]+"\t"+CommonUtils.strMinDistanceWithCutoff(kmers.get(id).getKmerString(), kmer.getKmerString(), km.getKmerString().length()));
							if (CommonUtils.strMinDistance(kmers.get(id).getKmerString(), kmer.getKmerString())<dc)
								isNotSimilar = false;
						}
						if (isNotSimilar)
							km = kmers.get(id);
					}
				}
			}
			results.add(km);

		}

		System.out.println(String.format("cluster_num=%d, kmer_distance<=%d, delta>=%d", centers.size(), dc, delta));
		System.out.println(Kmer.toShortHeader(k)+"\tId\tDensity\tDelta\tGamma\tCluster_size");
		int displayCount = Math.min(config.k_top, centers.size());
		for (int i=0;i<displayCount;i++){
			DensityClusteringPoint p = centers.get(i);
			System.out.println(String.format("%s    \t%d\t%.1f\t%.1f\t%.1f\t%d",
					results.get(i).toShortString(), p.id, p.density, p.delta, p.gamma, p.members.size()));
		}
		
		
		System.out.println(CommonUtils.timeElapsed(tic));
		
		return results;
	}
	
	/** Index k-mers for a set of sequences <br>
	 * 	record hit sequence id in the k-mers<br>
	 * 	keep track of the k-mer positions in the sequences so that we can use one to find the other<br>
	 * 	Sequences are modified to index the kmer match positions, and are set as unaligned.<br>
	 *  @return returns the mapping of kmer--[sequence ids]
	 * */
	private static HashMap <Kmer, HashSet<Integer>> mapBasicKmerSequences(HashSet<Kmer> kmers, ArrayList<Sequence> seqList){
		/** Next, k-mer search and index is done with basicKmers */
		// build the kmer search tree
		AhoCorasick oks = new AhoCorasick();
		for (Kmer km: kmers){
			oks.add(km.getKmerString().getBytes(), km);
	    }
		oks.prepare();
		
		// index kmer->seq id, seq->(kmer and position)
		// it start with the basicKmers, then the sub-kmers will be replaced by the WC-kmers
		HashMap <Kmer, HashSet<Integer>> kmer2seq = new HashMap <Kmer, HashSet<Integer>>();
		for (Sequence s:seqList){
			String seq = s.seq;						// get the sequence from original strand
			s.clearKmerIndex();
			s.resetAlignment();
			HashSet<Kmer> results = queryTree (seq, oks);
			if (!results.isEmpty()){
				for (Kmer km: results){		
					if (!kmer2seq.containsKey(km)){
						kmer2seq.put(km, new HashSet<Integer>());
					}
					kmer2seq.get(km).add(s.id);

					// This is DIFFERENT from alignSequencesByCoOccurence(ArrayList<Kmer>)
					// forward strand
					ArrayList<Integer> pos = StringUtils.findAllOccurences(seq, km.getKmerString());
					if (!s.fPos.containsKey(km))
						s.fPos.put(km, new HashSet<Integer>());
					for (int p:pos){
						s.fPos.get(km).add(p);
					}
					pos = StringUtils.findAllOccurences(seq, km.getKmerRC());
					for (int p:pos){
						s.fPos.get(km).add(p+RC);					// match kmer RC
					}
					// reverse strand
					pos = StringUtils.findAllOccurences(s.rc, km.getKmerString());
					if (!s.rPos.containsKey(km))
						s.rPos.put(km, new HashSet<Integer>());
					for (int p:pos){
						s.rPos.get(km).add(p);
					}
					pos = StringUtils.findAllOccurences(s.rc, km.getKmerRC());
					for (int p:pos){
						s.rPos.get(km).add(p+RC);					// match kmer RC
					}
				}
			}
		}
		
		// for each sequence, replace the sub-kmer hits with its WC kmers 
		// sub-kmer to WC-kmer mapping is many to many, need to collect all the positions for each WK, from all sub-kmers, 
		// then add WK-mers and remove sub-kmers.
		HashMap<Kmer, HashSet<Integer>> wkPosMap = new HashMap<Kmer, HashSet<Integer>>();
		HashSet<Kmer> replacedSubKmers = new HashSet<Kmer>();
		for (Sequence s:seqList){
			for (Kmer km: s.fPos.keySet()){		
				if (km.getGappedKmers()!=null){	// only process sub-kmers, ignore exact kmers
					for (GappedKmer wk: km.getGappedKmers()){
						if (! wkPosMap.containsKey(wk)){
							wkPosMap.put(wk, new HashSet<Integer>());
						}
						if (wk.getSubKmerOrientation(km)){		// same orientation, directly copy
							wkPosMap.get(wk).addAll(s.fPos.get(km));
						}
						else{		// diff orientation, reverse strand for the match kmer position
							for (int pos: s.fPos.get(km)){
								if (pos>RC/2)
									wkPosMap.get(wk).add(pos-RC);
								else
									wkPosMap.get(wk).add(pos+RC);
							}
						}
					}
					replacedSubKmers.add(km);
				}
			}
//			k=k+0;
			for (Kmer rk:replacedSubKmers)
				s.fPos.remove(rk);
			replacedSubKmers.clear(); 	// clear for each sequence
			for (Kmer wk: wkPosMap.keySet())
				s.fPos.put(wk, wkPosMap.get(wk));
			wkPosMap.clear(); 			// clear for each sequence
			
			// same for the reverse sequence
			for (Kmer km: s.rPos.keySet()){		
				if (km.getGappedKmers()!=null){	// only process sub-kmers, ignore exact kmers
					for (GappedKmer wk: km.getGappedKmers()){
						if (! wkPosMap.containsKey(wk)){
							wkPosMap.put(wk, new HashSet<Integer>());
						}
						if (wk.getSubKmerOrientation(km)){		// same orientation, directly copy
							wkPosMap.get(wk).addAll(s.rPos.get(km));
						}
						else{		// diff orientation, reverse strand for the match kmer position
							for (int pos: s.rPos.get(km)){
								if (pos>RC/2)
									wkPosMap.get(wk).add(pos-RC);
								else
									wkPosMap.get(wk).add(pos+RC);
							}
						}
					}
					replacedSubKmers.add(km);
				}
			}
			for (Kmer rk:replacedSubKmers)
				s.rPos.remove(rk);
			replacedSubKmers.clear(); 	// clear for each sequence
			for (Kmer wk: wkPosMap.keySet())
				s.rPos.put(wk, wkPosMap.get(wk));
			wkPosMap.clear(); 			// clear for each sequence
		}
		return kmer2seq;
	}
	
	/** Index k-mers and pos/neg sequences, remove un-enriched k-mers <br>
	 * 	record hit sequence id in the k-mers
	 * 	keep track of the k-mer positions in the sequences so that we can use one to find the other<br>
	 * 	This will set the sequences as unaligned.
	 * */
	private static void indexKmerSequences(ArrayList<Kmer> kmers, ArrayList<Sequence> seqList, 
			ArrayList<Sequence> seqListNeg, double kmer_hgp){

		/* Initialization, setup sequence list, and update kmers */
		/** basicKmers include the exact kmers and the sub-kmers of the gapped kmers */
		HashSet<Kmer> basicKmers = new HashSet<Kmer>();
		for (Kmer km: kmers){
			km.addBasicKmersToSet(basicKmers);
	    }
		int totalPosCount = seqList.size();
		int totalNegCount = seqListNeg.size();

		HashMap <Kmer, HashSet<Integer>> kmer2seq = mapBasicKmerSequences(basicKmers, seqList);
		HashMap <Kmer, HashSet<Integer>> kmer2seqNeg = mapBasicKmerSequences(basicKmers, seqListNeg);
		
		// update kmerCount, and hgp()
		HashSet<Kmer> unenriched = new HashSet<Kmer>();
		for (Kmer km:basicKmers){
			if (kmer2seq.containsKey(km))
				km.setPosHits(kmer2seq.get(km));
			else
				km.setPosHits(new HashSet<Integer>());
			if (kmer2seqNeg.containsKey(km))
				km.setNegHits(kmer2seqNeg.get(km));
			else
				km.setNegHits(new HashSet<Integer>());
		}
		for (Kmer km:kmers){
			if (km instanceof GappedKmer){	
				((GappedKmer)km).update();	// the pos hits of the sub-kmers has been updated
				if (km.getPosHitCount()==0)
					continue;
				km.setHgp(computeHGP(totalPosCount, totalNegCount, km.getPosHitCount(), km.getNegHitCount()));	// neg hit count is not change
				if (km.getHgp()>kmer_hgp)
					unenriched.add(km);
			}
			else{		// for exact kmers
				if (kmer2seq.containsKey(km)){
					if (km.getPosHitCount()==0)
						continue;
					km.setHgp(computeHGP(totalPosCount, totalNegCount, km.getPosHitCount(), km.getNegHitCount()));	
					if (km.getHgp()>kmer_hgp)
						unenriched.add(km);
				}
				else{
					unenriched.add(km);
				}
			}
		}
		kmers.removeAll(unenriched);
		
		// update seqList with GK hits
		HashSet<Kmer> basicKmersUsed = new HashSet<Kmer>();
		for (Kmer km: kmers){
			km.addBasicKmersToSet(basicKmersUsed);
	    }
		basicKmers.removeAll(basicKmersUsed); // now basicKmers contains all the un-used k-mers
		for (Sequence s:seqList){
			s.removeAllKmers(basicKmers);
		}
	}
	/**
	 * This is the main method for KMAC2 (density clustering) motif discovery
	 */
	public KmerCluster KmerMotifAlignmentClustering (ArrayList<Kmer> kmers, Kmer seed){
//		if (seed.getKmerString().equals("AGATTA"))
//			k+=0;
		tic = System.currentTimeMillis();
		int seed_range = k;
		KmerCluster cluster = new KmerCluster();
		
		if (kmers.size()==0)
			return null;
			
		cluster.seedKmer = seed;
		
		/** init kmerSet with seed family of seed kmer, by adding mismatch k-mers, order by #mm */
		ArrayList<Kmer> seedFamily = new ArrayList<Kmer>();
		seedFamily.add(seed);
		kmers.remove(seed);
		seed.setShift(0);
		seed.setSeedOrientation(true);
		seed.setAlignString(seed.getKmerString());
		if (config.use_seed_family){
			seedFamily.addAll(getFamilyKmers(kmers, cluster.seedKmer.getKmerString()));
			KmerGroup kg = config.use_weighted_kmer ? new KmerGroup(seedFamily, 0, seq_weights) : new KmerGroup(seedFamily, 0);
			cluster.seedKmer.familyHgp = computeHGP(kg.getGroupHitCount(), kg.getGroupNegHitCount());
			if (config.verbose>1)
				System.out.println(CommonUtils.timeElapsed(tic)+": Seed family hgp = "+cluster.seedKmer.familyHgp);
		}
		
		/** get the first sequence aligment using seedFamily kmer positions */
    	NewKSM newksm = new NewKSM(seedFamily);
    	cluster.alignedKmers = newksm.kmers;
    	cluster.ksmThreshold = newksm.threshold;
		alignSequencesUsingKmers (seedFamily, cluster);
    	seedFamily = null;
    	
    	// build first PWM
    	if (config.verbose>1 && config.noise!=0)
			System.out.println(CommonUtils.timeElapsed(tic)+ ": PWM noise = " + config.noise);
    	NewPWM newPWM = buildPWM(seqList, cluster, config.noise, tic, true);	
		
    	if (newPWM!=null){
			newPWM.updateClusterPwmInfo(cluster);
			alignSequencesUsingPWM(seqList, cluster);
			
			// extract first KSM
	    	if (config.use_ksm){
	    		NewKSM newKSM = extractKSM (seqList, seed_range, null);
				if (newKSM==null){
					return null;
				}		    	
				cluster.alignedKmers = newKSM.kmers;
		    	cluster.ksmThreshold = newKSM.threshold;

		    	alignSequencesUsingKmers(cluster.alignedKmers, cluster);
	    	}
    	
			/** Iteratively build PWM and KSM */
			iteratePWMandKSM (cluster, seqList, seed_range, config.use_ksm);
			
			if (config.evaluate_by_ksm && cluster.ksmThreshold!=null){
				if (cluster.ksmThreshold.posHit>cluster.total_aligned_seqs)
	    			cluster.total_aligned_seqs = cluster.ksmThreshold.posHit;
			}
			else{
				if (cluster.pwmPosHitCount>cluster.total_aligned_seqs)
	    			cluster.total_aligned_seqs = cluster.pwmPosHitCount;
			}
		}

		return cluster;
	}
	
	private void alignSequencesUsingKmers (ArrayList<Kmer> kmers, KmerCluster cluster){

		HashSet<Kmer> kmerSet = new HashSet<Kmer>();
		kmerSet.addAll(kmers);
		
		HashSet<Integer> kmerLinkedSeqIds = new HashSet<Integer>();
		for (Kmer km:kmers)
			kmerLinkedSeqIds.addAll(km.getPosHits());
		
		for (Sequence s : seqList){
    		s.resetAlignment();	
			if (!kmerLinkedSeqIds.contains(s.id))
				continue;
    		KmerGroup[] matches = scanKsmMatches(s, kmerSet);
    		if (matches==null)
    			continue;
			Arrays.sort(matches);
			KmerGroup kg = matches[0];	// take the top kg as match
			if (kg.hgp<=-cluster.ksmThreshold.score){
				s.score = -kg.hgp;		// score = -log10 hgp, becomes positive value
				if (kg.bs>RC/2){		// match on reverse strand
					s.RC();
					s.pos = -(kg.bs-RC); // seq_seed = - seed_seq
				}
				else
					s.pos = -kg.bs;
			}
		}
//		if (config.verbose>1)
//			System.out.println(CommonUtils.timeElapsed(tic)+ ": alignSequencesUsingKmers, done scan pos sequences.");
	}
	
	/** 
	 * Return KmerGroup matches in the sequence<br>
	 * Using the input kmers list as the KSM (passed in KSM with a threshold??)
	 * @param s
	 * @param kmers
	 * @return return KmerGroup object, return null if not match is found
	 */
	private KmerGroup[] scanKsmMatches(Sequence s, HashSet<Kmer> kmers){
		// need to add negPositionPadding to seed_seq for indexing the arrayList because seed_seq may be slightly negative.
		long t=System.currentTimeMillis();
		for (Kmer km: s.fPos.keySet()){
			if(!kmers.contains(km))
				continue;
			// the string contains this k-mer
			if (km.isSeedOrientation()){
				for (int km_seq: s.fPos.get(km)){
					if (km_seq<RC/2){		// if string match kmer
						int seed_seq = km_seq-km.getShift();		// seed_seq = km_seq - km_seed 
						forward.get(seed_seq+negPositionPadding).add(km);
					}
				}
				for (int km_seq: s.rPos.get(km)){
					if (km_seq<RC/2){		// if RC string match kmer, use reverse list
						int seed_seq = km_seq-km.getShift();		// seed_seq = km_seq - km_seed 
						reverse.get(seed_seq+negPositionPadding).add(km);
					}
				}
			}
			else{
				for (int km_seq: s.fPos.get(km)){
					if (km_seq>RC/2){		// if string match kmer RC (seed orientation)
						km_seq -= RC;		// remove RC_kmer
						int seed_seq = km_seq-km.getShift();		// seed_seq = km_seq - km_seed 
						forward.get(seed_seq+negPositionPadding).add(km);
					}
				}
				for (int km_seq: s.rPos.get(km)){
					if (km_seq>RC/2){		// if RC string match kmer RC (seed orientation), use reverse list
											// no change, use RC_kmer as RC_string
						int seed_seq = km_seq-km.getShift();		// seed_seq = km_seq - km_seed 
						reverse.get(seed_seq-RC+negPositionPadding).add(km);
					}
				}
			}
		}
//		System.out.println("Index k-mer Pos "+CommonUtils.timeElapsed(t));
		// count position hits in both strand
		int forwardCount=0,reverseCount=0;
		for (ArrayList<Kmer>kms:forward)
			if (!kms.isEmpty())
				forwardCount++;
		for (ArrayList<Kmer>kms:reverse)
			if (!kms.isEmpty())
				reverseCount++;
		
		if (forwardCount+reverseCount==0)
			return null;
		
		KmerGroup[] matches = new KmerGroup[forwardCount+reverseCount];
		int idx = 0;
		for (int p=0;p<forward.size();p++){
			if (!forward.get(p).isEmpty()){
				KmerGroup kg = config.use_weighted_kmer ? new KmerGroup(forward.get(p), p-negPositionPadding, seq_weights) : new KmerGroup(forward.get(p), p-negPositionPadding);
				matches[idx]=kg;
//				System.out.println("Made KG "+CommonUtils.timeElapsed(t));
				kg.setHgp(computeHGP(kg.getGroupHitCount(), kg.getGroupNegHitCount()));
				idx++;
				forward.get(p).clear();
//				System.out.println("HGP "+CommonUtils.timeElapsed(t));
			}
		}	
//		System.out.println(CommonUtils.timeElapsed(t));
		for (int p=0;p<reverse.size();p++){
			if (!reverse.get(p).isEmpty()){
				KmerGroup kg = config.use_weighted_kmer ? new KmerGroup(reverse.get(p), p-negPositionPadding+RC, seq_weights) : new KmerGroup(reverse.get(p), p-negPositionPadding+RC);
				matches[idx]=kg;
//				System.out.println("Made KG "+CommonUtils.timeElapsed(t));
				kg.setHgp(computeHGP(kg.getGroupHitCount(), kg.getGroupNegHitCount()));
				idx++;
				reverse.get(p).clear();
//				System.out.println("HGP "+CommonUtils.timeElapsed(t));
			}
		}
//		System.out.println("Done "+CommonUtils.timeElapsed(t));

		return matches;
	}

	private void outputClusters(int[] eventCounts){
		// output cluster information, PFM, and PWM
		File f = new File(outName);
		String name = f.getName();
		f=null;
			
		// output kmer alignments
		StringBuilder alignedKmer_sb = new StringBuilder();
		for (KmerCluster c:clusters){
			// print aligned k-mers of this cluster
			ArrayList<Kmer> alignedKmers = c.alignedKmers;
	    	if (!alignedKmers.isEmpty())
	    		alignedKmer_sb.append("Motif m"+c.clusterId+", n="+alignedKmers.size()+"\n");
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
		CommonUtils.writeFile(outName+".Alignement_k"+k+".txt", alignedKmer_sb.toString());
		alignedKmer_sb = null;

		// output PWM info
		System.out.println();
		StringBuilder pfm_sb = new StringBuilder();		// default TRANSFAC/STAMP format
		StringBuilder pfm_jasper_sb = new StringBuilder();		// JASPAR format
		StringBuilder pfm_meme_sb = new StringBuilder();		// MEME format
		StringBuilder pfm_homer_sb = new StringBuilder();		// HOMER format
		for (KmerCluster c:clusters){
     		WeightMatrix wm = c.wm;
    		System.out.println(String.format("--------------------------------------------------------------\n%s k-mer cluster #%d, aligned %d k-mers, %d sequences.", name, c.clusterId, c.alignedKmers.size(), c.total_aligned_seqs));
			int pos = c.pos_BS_seed-c.pos_pwm_seed;
    		if (pos>=0)
    			System.out.println(CommonUtils.padding(pos, ' ')+"|\n"+ WeightMatrix.printMatrixLetters(wm));
    		else
    			System.out.println(WeightMatrix.printMatrixLetters(wm));
    		System.out.println(String.format("PWM threshold: %.2f/%.2f, \thit=%d+/%d-, hgp=1e%.1f", c.pwmThreshold, c.wm.getMaxScore(), c.pwmPosHitCount, c.pwmNegHitCount, c.pwmThresholdHGP));
    		pfm_sb.append(CommonUtils.makeTRANSFAC (c.pfm, c.pwmPosHitCount, 
    				String.format("DE %s_m%d_k%d_p%d_c%d", name, c.clusterId, c.k, pos, c.pwmPosHitCount)));
			if (config.outputMEME)
				pfm_meme_sb.append(CommonUtils.makeMEME (c.pfm, c.pwmPosHitCount, 
						String.format("DE %s_m%d_k%d_p%d_c%d", name, c.clusterId, c.k, pos, c.pwmPosHitCount)));
			if (config.outputJASPAR)
				pfm_jasper_sb.append(CommonUtils.makeJASPAR (c.pfm, c.pwmPosHitCount, 
						String.format("DE %s_m%d_k%d_p%d_c%d", name, c.clusterId, c.k, pos, c.pwmPosHitCount)));
			if (config.outputHOMER)				
				pfm_homer_sb.append(CommonUtils.makeHOMER (c.pfm, c.pwmPosHitCount, 
						String.format("DE %s_m%d_k%d_p%d_c%d", name, c.clusterId, c.k, pos, c.pwmPosHitCount)));
    		if (config.use_ksm && c.ksmThreshold!=null)
    			System.out.println(String.format("KSM threshold: %.2f, \thit=%d+/%d-, hgp=1e%.1f", c.ksmThreshold.score, c.ksmThreshold.posHit, c.ksmThreshold.negHit, c.ksmThreshold.hgp));
			
			// paint motif logo
			c.wm.setNameVerType(name, "m"+c.clusterId, null);
			CommonUtils.printMotifLogo(c.wm, new File(outName+".m"+c.clusterId+".PWM.png"), 75);
			
			WeightMatrix wm_rc = WeightMatrix.reverseComplement(wm);
			wm_rc.setNameVerType(name, "m"+c.clusterId, "rc");
			CommonUtils.printMotifLogo(wm_rc, new File(outName+".m"+c.clusterId+".PWM_rc.png"), 75);
		}
		CommonUtils.writeFile(outName+".all.PFM.txt", pfm_sb.toString());
		if (config.outputMEME)
			CommonUtils.writeFile(outName+".all.PFM_MEME.txt", pfm_meme_sb.toString());
		if (config.outputJASPAR)
			CommonUtils.writeFile(outName+".all.PFM_JASPAR.txt", pfm_jasper_sb.toString());
		if (config.outputHOMER)
			CommonUtils.writeFile(outName+".all.PFM_HOMER.txt", pfm_homer_sb.toString());

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
		html.append("<p><ul>");
		for (KmerCluster c:clusters)
			html.append("<li><a href='"+name+".m"+c.clusterId+".KSM.txt'>Motif m"+c.clusterId+" KSM (K-mer Set Motif) file.</a>");
		html.append("<li><a href='"+name+".Alignement_k"+k+".txt'>K-mer alignment file.</a>");
		html.append("<li><a href='"+name+".all.PFM.txt'>All motif PFMs</a></ul>");
		html.append("<p>K-mers of the motifs<p>");
		html.append("<table border=1><th>K-mer</th><th>Motif</th><th>Offset</th><th>Pos Hit</th><th>Neg Hit</th><th>HGP</th>");
		
    	int leftmost_km = Integer.MAX_VALUE;
    	ArrayList<Kmer> outputs = new ArrayList<Kmer>();		
    	for (int j=0;j<clusters.size();j++){
    		ArrayList<Kmer> alignedKmers = Kmer.cloneKmerList(clusters.get(j).alignedKmers);		// each kmer in diff cluster has been clone, would not overwrite
	    	for (int i=0;i<Math.min(20, alignedKmers.size());i++){
	    		Kmer km = alignedKmers.get(i);
	    		km.setClusterId(i);
				if (km.getKmerStartOffset()<leftmost_km)
					leftmost_km = km.getKmerStartOffset();
				outputs.add(km);
			}
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
			String kmString = km.isSeedOrientation()?km.getKmerString():km.getKmerRC();
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
			html.append(String.format("<td>%d</td><td>%d</td><td>%d</td><td>%d</td><td>%.1f</td></tr>", 
					km.getClusterId(), km.getKmerStartOffset(), km.getPosHitCount(), km.getNegHitCount(), km.getHgp()));
		}
		html.append("</table>");
		html.append("</td><td valign='top'><br>");
		html.append("<table border=0 align=center><th>Motif PWM</th><th>Motif spatial distribution (w.r.t. primary PWM)<br>Format: position,motif_occurences</th>");
		for (KmerCluster c:clusters){
    		html.append("<tr><td><img src='"+name+".m"+c.clusterId+".PWM.png"+"'><a href='#' onclick='return popitup(\""+name+".m"+c.clusterId+".PWM_rc.png\")'>rc</a><br>");
    		html.append(String.format("PWM: %.2f/%.2f, hit=%d+/%d-, hgp=1e%.1f<br>", 
    				c.pwmThreshold, c.wm.getMaxScore(), c.pwmPosHitCount, c.pwmNegHitCount, c.pwmThresholdHGP));
//    		html.append(String.format("KSM score: %.2f, \thit=%d+/%d-, hgp=1e%.1f<br><br>", 
//    				c.ksmThreshold.score, c.ksmThreshold.posHit, c.ksmThreshold.negHit, c.ksmThreshold.hgp));
    		String suffix = name+".Spatial_dist.m0_m"+c.clusterId;
    		html.append("</td><td><a href='"+suffix+".txt'>"+"<img src='"+suffix+".png"+"' height='150'></a></td></tr>");
		}
		html.append("</table>");
		html.append("</td></tr></table>");
		CommonUtils.writeFile(outName+"_result.htm", html.toString());
		
//		for (int i=0;i<Math.min(clusters.size(), 20);i++){
//			ArrayList<Kmer> clusterKmers = clusters.get(i).alignedKmers;
//			MotifThreshold t = this.estimateClusterKgsThreshold(clusterKmers);
//			if (t!=null)
//				System.out.println(String.format("%d\t%.2f\t%d\t%d\t%.1f", i, t.score, t.posHit, t.negHit, t.hgp ));
//		}
		
	}
	
	/** 
	 * Recursively merge overlapped motif clusters using PWMs<br>
	 * Assuming the clusters are sorted by cluster id (0-based)<br>
	 * This method consider the motif hit distance, thus avoiding merging dis-similar motifs that co-bind many sequences<br>
	 * checked matrix in indexed by the clusterid, thus is ok to sort the clusters without changing clusterid 
	 */
	@SuppressWarnings({ "unchecked", "rawtypes" })
	private void mergeOverlapClusters (ArrayList<Sequence> seqList, boolean[][] checked, 
			boolean useKSM, int recursion){	
		
		if (recursion>10)
			return;
		if (config.verbose>1)
    		System.out.println("\n"+CommonUtils.timeElapsed(tic)+": Merge overlapping motifs, iteration " + recursion);
		boolean isChanged = false;
		int maxClusterId=0;
		for (KmerCluster c:clusters)
			if (c.clusterId > maxClusterId)
				maxClusterId = c.clusterId;
		
		// hits matrix contains the sequence-motif hit, indexed by clusterid value
		ArrayList[][] hits = new ArrayList[seqs.length][maxClusterId+1];
		for (int j=0;j<clusters.size();j++){
			KmerCluster c = clusters.get(j);
			if (c.wm==null)
				continue;
			WeightMatrixScorer scorer = new WeightMatrixScorer(c.wm);
			for (int i=0;i<seqs.length;i++){
				hits[i][c.clusterId]=CommonUtils.getAllPWMHit(seqs[i], c.wm.length(), scorer, c.pwmThreshold);
			}
		}
				
		int seqLen = seqs[0].length();
		ArrayList<KmerCluster> toRemove = new ArrayList<KmerCluster>();
		for (int m=0;m<clusters.size();m++){
			for (int j=m+1;j<clusters.size();j++){
//				if (m>=clusters.size()||j>=clusters.size())			// removing a cluster may change the cluster index/size
//					continue;
				KmerCluster cluster1 = clusters.get(m);
				KmerCluster cluster2 = clusters.get(j);
				if (cluster1.pwmThresholdHGP > cluster2.pwmThresholdHGP){
					cluster1 = clusters.get(j);
					cluster2 = clusters.get(m);
				}
				if (cluster1.wm==null||cluster2.wm==null||checked[cluster1.clusterId][cluster2.clusterId])
					continue;
				
				int range = seqLen - cluster1.wm.length()/2 - cluster2.wm.length()/2 + 4;  // add 2 for rounding error
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
				int minHitCount = Math.min(cluster1.pwmPosHitCount, cluster2.pwmPosHitCount);
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
				
				if (maxCount>=minHitCount*config.motif_hit_overlap_ratio){				// if there is large enough overlap, try to merge 2 clusters
					if (config.verbose>1)
			    		System.out.println(String.format("\n%s: Motifs to be merged, %d overlap hits, motif hit distance = %d%s:\n#%d %s\t(hit=%d, hgp=1e%.1f)\n#%d %s\t(hit=%d, hgp=1e%.1f)\n", 
			    				CommonUtils.timeElapsed(tic), maxCount, maxDist, isRC?"rc":"", 
			    				cluster1.clusterId, WeightMatrix.getMaxLetters(cluster1.wm), cluster1.pwmPosHitCount, cluster1.pwmThresholdHGP, 
			    				cluster2.clusterId, WeightMatrix.getMaxLetters(cluster2.wm), cluster2.pwmPosHitCount, cluster2.pwmThresholdHGP));
					
					KmerCluster newCluster = cluster1.clone();
					alignSequencesUsingPWM(seqList, newCluster);	// align with first PWM
					
					// align additional sequences matching second PWM
					WeightMatrix wm = isRC?WeightMatrix.reverseComplement(cluster2.wm):cluster2.wm;
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
			        	  double score = profiler.getMaxScore(p);
			        	  if (maxSeqScore<score || (maxSeqScore==score && maxScoringStrand=='-')){	// equal score, prefer on '+' strand
			        		  maxSeqScore = score;
			        		  maxScoringShift = p;
			        		  maxScoringStrand = profiler.getMaxStrand(p);
			        	  }
			          }
			          // if a sequence pass the motif score, align with PWM hit
			          if (maxSeqScore >= cluster2.pwmThreshold){
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
						checked[cluster1.clusterId][cluster2.clusterId]=true;
						checked[cluster2.clusterId][cluster1.clusterId]=true;
			    		System.out.println(CommonUtils.timeElapsed(tic)+": PWM "+WeightMatrix.getMaxLetters(wm)+" is removed.");
						continue;
					}
						
			    	// build PWM
					NewPWM newPWM = buildPWM(seqList, newCluster, 0, tic, false);
					indexKmerSequences(k2kmers.get(newCluster.k), seqList, seqListNeg, config.kmer_hgp);  // need this to get KSM
					if (newPWM!=null){
						newPWM.updateClusterPwmInfo(newCluster);
						alignSequencesUsingPWM(seqList, newCluster);
					}
					
					// Iteratively improve PWM and sequences alignment
					iteratePWMandKSM (newCluster, seqList, newCluster.seedKmer.getK(), useKSM);
					
					// merge if the new PWM is more enriched
					if ((newCluster.pwmThresholdHGP<cluster1.pwmThresholdHGP)){		
						if (newCluster.pwmPosHitCount>newCluster.total_aligned_seqs)
							newCluster.total_aligned_seqs = newCluster.pwmPosHitCount;
						clusters.set(m, newCluster);
						isChanged = true;
						for (int d=0;d<checked.length;d++){
							checked[d][cluster1.clusterId]=false;
							checked[cluster1.clusterId][d]=false;
						}
						cluster1 = newCluster;
						if (config.verbose>1)
				    		System.out.println(String.format("\n%s: motif #%d and motif #%d merge to new PWM %s, hgp=1e%.1f", 
				    				CommonUtils.timeElapsed(tic), cluster1.clusterId, cluster2.clusterId,
				    				WeightMatrix.getMaxLetters(newCluster.wm), newCluster.pwmThresholdHGP));	
					}
					else{	
						checked[cluster1.clusterId][cluster2.clusterId]=true;
						checked[cluster2.clusterId][cluster1.clusterId]=true;
						if (config.verbose>1)
				    		System.out.println(String.format("%s: Merged PWM is not more enriched, do not merge", CommonUtils.timeElapsed(tic)));
					}
					
					// No matter successful merging or not, try the other PWM after removing the overlap hits
					if (config.verbose>1)
			    		System.out.println(String.format("%s: Testing the remaining motif #%d.", 
			    			CommonUtils.timeElapsed(tic), cluster2.clusterId));						
					alignSequencesUsingPWM(seqList, cluster1);
					ArrayList<Sequence> seqList_j = new ArrayList<Sequence>();
					for (Sequence s:seqList)
						if (s.pos == UNALIGNED)
							seqList_j.add(s);
					int aligned_seqs_count = alignSequencesUsingPWM(seqList_j, cluster2);
					
					// if aligned seq count is less than threshold, or if it contains less than half of total hit of the motif (i.e. majority of hits are still overlapped), remove it
					if (aligned_seqs_count<seqs.length*config.motif_hit_factor || aligned_seqs_count<cluster2.pwmPosHitCount/3){
						if (config.verbose>1)
				    		System.out.println(String.format("%s: Motif #%d has too few (%d) non-shared hits, remove it.", 
				    			CommonUtils.timeElapsed(tic), cluster2.clusterId, aligned_seqs_count));
						cluster2.wm = null;
						toRemove.add(cluster2);
						checked[cluster1.clusterId][cluster2.clusterId]=true;
						checked[cluster2.clusterId][cluster1.clusterId]=true;
					}
					else{
						// for the second motif, build PWM using the non-cluster1-hit sequences
						// do not interate to improve motif, just learn and stop
						WeightMatrix old_j = cluster2.wm;
						cluster2.wm = null;			// reset here, to get a new PWM
						int alignedSeqCount = 0;
						newPWM = buildPWM(seqList_j, cluster2, 0, tic, false);
						if (newPWM!=null){
							newPWM.updateClusterPwmInfo(cluster2);
							// TODO: only clone the list, not the k-mers
							indexKmerSequences( (ArrayList<Kmer>)(k2kmers.get(cluster2.k).clone()), seqList_j, seqListNeg, config.kmer_hgp);  // need this to get KSM
							alignedSeqCount = alignSequencesUsingPWM(seqList_j, cluster2);
							if (config.evaluate_by_ksm){
								NewKSM newKSM = extractKSM (seqList, cluster2.k, null);
								if (newKSM!=null&&newKSM.threshold!=null){
									cluster2.alignedKmers = newKSM.kmers;
									cluster2.ksmThreshold = newKSM.threshold;
								}
								else{	// if no good KSM
									cluster2.wm = null;
									toRemove.add(cluster2);
									checked[cluster1.clusterId][cluster2.clusterId]=true;
									checked[cluster2.clusterId][cluster1.clusterId]=true;
						    		System.out.println(CommonUtils.timeElapsed(tic)+": PWM "+WeightMatrix.getMaxLetters(wm)+" is removed.");
									continue;
								}
							}
							// check if new PWM has sufficient hit
							if (alignedSeqCount<seqs.length*config.motif_hit_factor){
								if (config.verbose>1)
						    		System.out.println(String.format("%s: too few new PWM hits (%d) , remove motif #%d.", 
						    			CommonUtils.timeElapsed(tic), alignedSeqCount, cluster2.clusterId));
								cluster2.wm = null;
								toRemove.add(cluster2);
								checked[cluster1.clusterId][cluster2.clusterId]=true;
								checked[cluster2.clusterId][cluster1.clusterId]=true;
							}
							else{
								if (config.verbose>1)
						    		System.out.println(String.format("%s: new PWM has sufficient hit %d, keep motif #%d.", 
						    			CommonUtils.timeElapsed(tic), alignedSeqCount, cluster2.clusterId));
								// update the # aligned_seqs using the number from all sequences
								cluster2.total_aligned_seqs = cluster2.pwmPosHitCount;
								if (!isChanged){		// if motif 1 is not changed, set this pair as checked.
									checked[cluster1.clusterId][cluster2.clusterId] = true;
									checked[cluster2.clusterId][cluster1.clusterId] = true;
								}
							}
						}
					}
				}		
				else{					// if overlap is not big enough
					checked[cluster1.clusterId][cluster2.clusterId] = true;
					checked[cluster2.clusterId][cluster1.clusterId] = true;
				}
			}
		}
		clusters.removeAll(toRemove);
		
		// if merged, recursively merge, otherwise return
		if (isChanged){		
			// checked matrix is indexed by the clusterid, thus is ok to sort the clusters without changing clusterid 
			Collections.sort(clusters);			
			mergeOverlapClusters (seqList, checked, useKSM, recursion+1);
		}
	}
	
	/** Iteratively build PWM and KSM <br>
	 *  The seqList should have been aligned so that a new PWM can be built. <br>
	 *  The seqList should have been indexed with kmers. <br>
	 *  Need to have use_KSM flag here because of PWM only refinement */
	private void iteratePWMandKSM (KmerCluster cluster, ArrayList<Sequence> seqList, int seed_range, boolean use_KSM){	
    	while(true){	
			NewPWM newPWM = buildPWM(seqList, cluster, 0, tic, true);	
			if (newPWM==null)
				return;
			// test if we want to accept the new PWM
			if( (!config.evaluate_by_ksm || config.evaluate_by_ksm&&!use_KSM) && cluster.pwmThresholdHGP<=newPWM.hgp){
				return;		// previous pwm is more enriched, stop here
			}
			else
				newPWM.updateClusterPwmInfo(cluster);
			alignSequencesUsingPWM(seqList, cluster);
			
			if (use_KSM){
				NewKSM newKSM = extractKSM (seqList, seed_range, null);
				if (newKSM==null ||newKSM.threshold==null || config.evaluate_by_ksm && (newKSM.threshold.hgp >= cluster.ksmThreshold.hgp))
					return;
				cluster.alignedKmers = newKSM.kmers;
				cluster.ksmThreshold = newKSM.threshold;

				alignSequencesUsingKmers(cluster.alignedKmers, cluster);
			}
    	}  // Iteratively improving the HGP of the motif
	}

	
	/**
	 * Compute the distance distributions between primary and secondary motifs
	 */
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public void computeMotifDistanceDistribution (String name){		
		System.out.println("\nCompute motif distance distribution ...");
		
		ArrayList[][] hits = new ArrayList[seqs.length][clusters.size()];
		for (int j=0;j<clusters.size();j++){
			KmerCluster c = clusters.get(j);
//			if (c.wm!=null)
//				System.err.println("Cluster "+j+" PWM length="+c.wm.length());
			if (c.wm!=null){
				WeightMatrixScorer scorer = new WeightMatrixScorer(c.wm);
				for (int i=0;i<seqs.length;i++){
					hits[i][j]=CommonUtils.getAllPWMHit(seqs[i], c.wm.length(), scorer, c.pwmThreshold);
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
			for (int j=0;j<clusters.size();j++){
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
	private HashMap<Integer, PWMHit> findAllPWMHits(ArrayList<Sequence> seqList, KmerCluster cluster, double pwm_factor){
		
		WeightMatrix wm = cluster.wm;
	    WeightMatrixScorer scorer = new WeightMatrixScorer(wm);
	    HashMap<Integer, PWMHit> seq2hits = new HashMap<Integer, PWMHit>();
	    
		for (int i=0;i<seqList.size();i++){
		  Sequence s = seqList.get(i);
		  String seq = s.getSeqStrand(true);			// PWM to scan all sequence
		      	  
	      WeightMatrixScoreProfile profiler = scorer.execute(seq);
	      double maxSeqScore = Double.NEGATIVE_INFINITY;
	      int maxScoringShift = 0;
	      char maxScoringStrand = '+';
	      for (int j=0;j<profiler.length();j++){
	    	  double score = profiler.getMaxScore(j);
	    	  if (maxSeqScore<score || (maxSeqScore==score && maxScoringStrand=='-')){	// equal score, prefer on '+' strand
	    		  maxSeqScore = score;
	    		  maxScoringShift = j;
	    		  maxScoringStrand = profiler.getMaxStrand(j);
	    	  }
	      }
	      // if a sequence pass the motif score, store it
//	      if (maxSeqScore > 0){
	      if (maxSeqScore >= wm.getMaxScore()*pwm_factor){
	    	  PWMHit hit = new PWMHit();
	    	  hit.clusterId = cluster.clusterId;
	    	  hit.wm = wm;
	    	  hit.seqId = i;
	    	  hit.start = maxScoringShift;					// start is cordinate on + strand
	    	  hit.end = maxScoringShift+wm.length()-1;		// end inclusive
			  hit.isForward = maxScoringStrand =='+';
			  String ss = seq.substring(maxScoringShift, maxScoringShift+wm.length());
			  if (maxScoringStrand =='-')
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
//	
//	/** refine PWMs by EM like iterations<br>
//	 * given the un-masked positive sequences, and multiple PWMs<br>
//	 * find the best PWMs to cluster the sequences
//	 */
//	private void refinePWMs(ArrayList<Sequence> seqList, String[] pos_seq_backup, String[] neg_seq_backup, 
//			ArrayList<Kmer> kmers_in, int seed_range){
//	
//		ArrayList<KmerCluster> noPWM = new ArrayList<KmerCluster>();
//		for (KmerCluster c:clusters){
//			if (c.wm==null){
//				noPWM.add(c);
//				c.pi = c.pwmPosHitCount;
//			}
//		}
//		clusters.removeAll(noPWM);
//		if (clusters.size()<=1)
//			return;
//		
//		File f = new File(outName);
//		String name = f.getName();
//		
//		/** Resolve multiple-PWM-match conflict */
//		for (int iter=0; iter<100; iter++){
//			
//			// record cluster PWM hgp for comparison
//			double[] hgps = new double[clusters.size()];
//			if (config.verbose>1)
//				System.out.println("------------------------------------------------\n"+CommonUtils.timeElapsed(tic)+
//						": Iteration #"+iter);
//			
//			for (int i=0; i<clusters.size(); i++){
//				KmerCluster c = clusters.get(i);
//				hgps[i] = c.pwmThresholdHGP;
//				c.clusterId = i;
//				c.seq2hits = findAllPWMHits(seqList, c, config.wm_factor);
//				// paint motif logo
//				c.wm.setNameVerType(name+"_i"+iter, "#"+c.clusterId, "");
//				CommonUtils.printMotifLogo(c.wm, new File(outName+"_i"+iter+"_"+c.clusterId+"_motif.png"), 75);
//			}	
//			if (config.verbose>1){
//				StringBuilder sb = new StringBuilder(CommonUtils.timeElapsed(tic)+": Cluster PWMs\t");
//				for (KmerCluster c:clusters)
//					sb.append(WeightMatrix.getMaxLetters(c.wm)).append(" ");
//				sb.append("\n").append(CommonUtils.timeElapsed(tic)+": Cluster hgps\t"+CommonUtils.arrayToString(hgps));
//				System.out.println(sb.toString());
//			}
//			int conflicts = resolveConflictHits(seqList);
//			if (config.verbose>1)
//				System.out.println(CommonUtils.timeElapsed(tic)+": "+conflicts+" sequences share hits by multiple PWMs");
//			if (conflicts==0)
//				break;
//				
//			buildAllClusterPWMfromHits(seqList);
//			
//			if (clusters.size()<=1)
//				break;
//			if (clusters.size()!=hgps.length)		// some PWM/cluster may be removed
//				continue;
//			
//			boolean noChange = true;
//			for (int i=0; i<clusters.size(); i++){
//				if (hgps[i] != clusters.get(i).pwmThresholdHGP){
//					noChange = false;
//					break;
//				}
//			}
//			if (noChange)
//				break;
//		}
//	}
	

	private NewPWM buildPWM(ArrayList<Sequence> seqList, KmerCluster cluster, double noiseRatio, long tic, boolean onlyBetter){	    	
		ArrayList<String> alignedSeqs = new ArrayList<String>();
		ArrayList<Double> weights = new ArrayList<Double>();
		for (Sequence seq:seqList){
			if (seq.pos==UNALIGNED)
				continue;
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
		alignedSeqs.trimToSize();
		
		if (alignedSeqs.size()<seqs.length*config.motif_hit_factor){
			if (config.verbose>1)
	    		System.out.println(String.format("%s: Seq_Count %d is too few, stop building PWM.", CommonUtils.timeElapsed(tic), alignedSeqs.size()));
			return null;
		}
		
		return buildPWMfromAlignedSequences(alignedSeqs, weights, cluster, noiseRatio, onlyBetter);
    }
	
	public NewPWM buildPWMfromHits(ArrayList<Sequence> seqList, KmerCluster cluster, Iterator<PWMHit> hits){
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
	

	private NewPWM buildPWMfromAlignedSequences(ArrayList<String> alignedSeqs, ArrayList<Double> weights, 
		KmerCluster cluster, double noiseRatio, boolean onlyBetter){
		if (alignedSeqs.isEmpty())
			return null;
		double[][] pfm = new double[alignedSeqs.get(0).length()][MAXLETTERVAL];
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
			if (config.verbose>1)
				System.out.println("PWM is too short, W="+(rightIdx-leftIdx+1));
			return null;
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

//	    	if (bmconfig.verbose>1)
//	    		System.out.println(String.format("%s: got PWM.", CommonUtils.timeElapsed(tic)));
	    	
	    	// Check the quality of new PWM: hyper-geometric p-value test using the positive and negative hit counts
	    	MotifThreshold estimate = null;
	    	if (config.optimize_pwm_threshold)
	    		estimate = optimizePwmThreshold(wm, "", wm.getMaxScore()*config.wm_factor);
	    	else
	    		estimate = estimatePwmThreshold(wm, wm.getMaxScore()*config.wm_factor);
	    	double pwmThreshold = estimate.score;
	    	double pwmThresholdHGP = estimate.hgp;
    		if (config.verbose>1)
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
			if (config.verbose>1)
				System.out.println(CommonUtils.timeElapsed(tic)+": None of PWM is enriched.");
			return null;
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
		newPWM.pwmGoodQuality = (bestEstimate.posHit>seqs.length*config.motif_hit_factor);
		newPWM.threshold = bestEstimate.score;
		newPWM.hgp = bestHGP;
		newPWM.pwmPosHitCount = bestEstimate.posHit;
		newPWM.pwmNegHitCount = bestEstimate.negHit;		
		newPWM.pfm = pfm_trim;
		newPWM.pos_pwm_seed = bestLeft-(alignedSeqs.get(0).length()/2-k/2);		// pwm_seed = pwm_seqNew-seed_seqNew

    	return newPWM;
	}
	
	private class NewPWM{
		double threshold = 0;
    	double hgp = -0.1;
    	WeightMatrix wm = null;
    	boolean pwmGoodQuality = false;
    	int pwmPosHitCount = 0;
    	int pwmNegHitCount = 0;
    	float[][] pfm = null;
    	int pos_pwm_seed = 0;   
    	
    	private void updateClusterPwmInfo(KmerCluster cluster){
    		if (wm==null)
    			return;
    		cluster.wm = wm;
        	cluster.pwmGoodQuality = pwmGoodQuality;
        	cluster.pwmThreshold = threshold;
        	cluster.pwmThresholdHGP = hgp;
        	cluster.pwmPosHitCount = pwmPosHitCount;
        	cluster.pwmNegHitCount = pwmNegHitCount;
        	cluster.pfm = pfm;
        	cluster.pos_pwm_seed = pos_pwm_seed;
    	}
	}
	
	private void alignSequencesUsingKSM_old(ArrayList<Sequence> seqList, KmerCluster cluster){
		if (cluster==null)
			return;
		
		// scan all sequences, align them using kmer hit results
		for (Sequence s : seqList){
			s.resetAlignment();			
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
				if (kg.hgp<=-cluster.ksmThreshold.score){
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
	}
	
	private void alignSequencesUsingKSM(ArrayList<Sequence> seqList, KmerCluster cluster){
		if (cluster==null)
			return;
		
		// scan all sequences, align them using kmer hit results
		for (Sequence s : seqList){
			s.resetAlignment();			
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
				if (kg.hgp<=-cluster.ksmThreshold.score){
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
	}
	/** align seqList using a PWM<br>
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
	          double maxSeqScore = Double.NEGATIVE_INFINITY;
	          int maxScoringShift = 0;
	          char maxScoringStrand = '+';
	          for (int j=0;j<profiler.length();j++){
	        	  double score = profiler.getMaxScore(j);
	        	  if (maxSeqScore<score || (maxSeqScore==score && maxScoringStrand=='-')){	// equal score, prefer on '+' strand
	        		  maxSeqScore = score;
	        		  maxScoringShift = j;
	        		  maxScoringStrand = profiler.getMaxStrand(j);
	        	  }
	          }
	          // if a sequence pass the motif score, align with PWM hit
	          if (maxSeqScore >= cluster.pwmThreshold){
				if (maxScoringStrand =='-'){
					maxScoringShift = seqLen-maxScoringShift-wm.length();
					s.RC();
					// i.e.  (seq.length()-1)-maxScoringShift-(wm.length()-1);
				}
				s.pos = cluster.pos_pwm_seed-maxScoringShift;
				count_pwm_aligned ++;
	          }
	          else
	        	  s.pos = UNALIGNED;
	        }	// each sequence

	    	if (config.verbose>1)
	    		System.out.println(CommonUtils.timeElapsed(tic)+": PWM "+WeightMatrix.getMaxLetters(wm)+" align "+count_pwm_aligned+" sequences.");
	    	
	    	return count_pwm_aligned;
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
	 * Get all the kmers that has k/4 (k>=8) or 1 (k<8) mismatch to the kmerStr<br>
	 * and set alignString and its shift for the kmers, the shift is wrt the input kmerStr orientation
	 * @param kmers
	 * @param kmerStr the length should be the same as kmers
	 * @param shift
	 * @return
	 */
	private ArrayList<Kmer> getFamilyKmers(ArrayList<Kmer> kmers, String kmerStr) {
		ArrayList<Kmer> family = new ArrayList<Kmer>();
		if (kmers.isEmpty())
			return family;
		ArrayList<Kmer> kmerListCopy = Kmer.copyKmerList(kmers);
		ArrayList<Kmer> toRemove = new ArrayList<Kmer>();
		int mm = 1;
		if (k>=8 && this.use_smart_mm)
			mm = k/4;
		// progressively allow more mismatch, this will give the more similar kmer priorty for sequence matching
		for (int i=1;i<=mm;i++){
			for (Kmer kmer: kmerListCopy){
		    	if (CommonUtils.strMinDistance(kmerStr, kmer.getKmerString())==i){
		    		Pair<Integer,Integer> p = CommonUtils.strMinDistanceAndShift(kmerStr, kmer.getKmerString());
		    		kmer.setAlignString(kmer.getKmerString());
		    		kmer.setShift(p.cdr());
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
	
	/** Create a NewKSM will update the ksm Engine */
	private class NewKSM{
		private ArrayList<Kmer> kmers=null;
		private MotifThreshold threshold = null;
		public ArrayList<Kmer> getKmers(){
			return kmers;
		}
		public MotifThreshold getThreshold(){
			return threshold;
		}
		/** Create a NewKSM will update the ksm Engine */
		NewKSM(ArrayList<Kmer> kmers){
			this.kmers = kmers;
			this.setThreshold();
		}
		void setThreshold(){	
			threshold = optimizeKsmThreshold(seqList, seqListNeg, kmers);
			if (threshold==null)
				return;
			if (config.verbose>1)
				System.out.println(String.format("%s: KSM score %.2f\tmatch %d+/%d- seqs\thgp=1e%.1f", 
						CommonUtils.timeElapsed(tic), threshold.score, threshold.posHit, threshold.negHit, threshold.hgp));
		}
	}
	/**
	 * Get aligned kmers within seed_range using the aligned sequences
	 * @param seqList
	 * @param seed_range
	 * @param excludes
	 * @return
	 */
	private NewKSM extractKSM (ArrayList<Sequence> seqList, int seed_range, ArrayList<Kmer> excludes){

		/** kmer2pos: record all the hit positions (in reference to the seed position) of all k-mers in the alignment */
		HashMap<Kmer, ArrayList<Integer>> kmer2pos_seed = new HashMap<Kmer, ArrayList<Integer>>();
		for (Sequence s:seqList){
			if (s.pos != UNALIGNED){		// aligned seqs
				HashMap<Kmer, HashSet<Integer>> kmerPos_seq = s.getKmerPos();	// get from the aligned strand orientation
				for (Kmer km:kmerPos_seq.keySet()){
					if (excludes!=null && excludes.contains(km))
						continue;
					if (!kmer2pos_seed.containsKey(km))
						kmer2pos_seed.put(km, new ArrayList<Integer>());
					HashSet<Integer> hits_seq = kmerPos_seq.get(km);
					if (hits_seq==null)
						continue;
					for (int km_seq:hits_seq){
						if (km_seq>RC/2){		// if it is kmerRC match on the aligned strand
							int km_seed = (km_seq-RC)+s.pos;	// km_seed = km_seq + seq_seed
							if (km_seed>=-seed_range && km_seed<=seed_range)
								kmer2pos_seed.get(km).add(km_seed+RC);		// RC to label that kmer is in opposite orientation of seed kmer
						}
						else{
							int pos = km_seq+s.pos;
							if (pos>=-seed_range && pos<=seed_range)
								kmer2pos_seed.get(km).add(pos);		
						}
					}
				}
			}
		}

    	/** find k-mers that are consistently aligned, set kmer consensus position */
		ArrayList<Kmer> alignedKmers = new ArrayList<Kmer>();
		for (Kmer km:kmer2pos_seed.keySet()){
			ArrayList<Integer> posKmer = kmer2pos_seed.get(km);		// all km_seed positions of this kmer
			// The kmer hit in the 2*k region should be at least 1/2 of total hit
			if (posKmer==null || posKmer.size() < km.getPosHitCount()*config.kmer_aligned_fraction){			
				km.setAlignString("Too few hit "+posKmer.size());
				continue;
			}	
			// find the most frequent kmerPos
			Pair<int[], int[]> sorted = StatUtil.sortByOccurences(posKmer);
			int counts[] = sorted.cdr();
			int posSorted[] = sorted.car();
			int maxCount = counts[counts.length-1];
			// posKmer.size() only count in seqs in the alignment, getPosHitCount() count all seqs
			if (maxCount < Math.min(posKmer.size(),km.getPosHitCount()) * config.kmer_aligned_fraction){
				km.setAlignString("Low consistent hit count "+maxCount);
				continue;
			}
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
			if (shift>RC/2){		// if kmer is in opposite orientation of seed kmer
				km.setSeedOrientation(false);
				shift -= RC;
			}
			else
				km.setSeedOrientation(true);
			km.setShift(shift);			
			km.setAlignString(maxCount+"/"+posKmer.size());
			km.setKmerStartOffset(km.getShift());
			alignedKmers.add(km);
		}	
		
		kmer2pos_seed = null;
//		System.out.println(String.format("%s: Extracted %d k-mers.", CommonUtils.timeElapsed(tic), alignedKmers.size()));

		if (alignedKmers.isEmpty())
			return null;
		
		NewKSM newKSM = null;
		if (config.optimize_kmer_set){
			int tmp = alignedKmers.size();
			optimizeKSM(alignedKmers);
			if (config.verbose>1)
				System.out.println(String.format("%s: Optimized %d to %d k-mers.", CommonUtils.timeElapsed(tic), tmp, alignedKmers.size()));

			newKSM = new NewKSM(alignedKmers);
		}
		else{
			alignedKmers.trimToSize();
			newKSM = new NewKSM(alignedKmers);
		}		
		return newKSM;
	}

	/**
	 * Optimize the kmer set by removing non-essential k-mers to improve the overall HGP
	 * @param alignedKmers
	 */
	private void optimizeKSM(ArrayList<Kmer> kmers){
		if (kmers.size()<=1)
			return;

		Collections.sort(kmers, new Comparator<Kmer>(){
		    public int compare(Kmer o1, Kmer o2) {
		    	return o1.compareByHGP(o2);
		    }
		});	
		Collections.reverse(kmers);		// reverse so that weaker hit are remove before stronger hit
		boolean changed=true;
		
		// mapping from sequence id to kmers
		HashMap<Integer, HashSet<Kmer>> seq2kmers = new HashMap<Integer, HashSet<Kmer>>();
		for (Kmer km: kmers){
			HashSet<Integer> hits = km.getPosHits();
			for (int h:hits){
				if (!seq2kmers.containsKey(h))
					seq2kmers.put(h, new HashSet<Kmer>());
				seq2kmers.get(h).add(km);
			}
		}
		HashMap<Integer, HashSet<Kmer>> seq2kmers_neg = new HashMap<Integer, HashSet<Kmer>>();
		for (Kmer km: kmers){
			HashSet<Integer> hits = km.getNegHits();
			for (int h:hits){
				if (!seq2kmers_neg.containsKey(h))
					seq2kmers_neg.put(h, new HashSet<Kmer>());
				seq2kmers_neg.get(h).add(km);
			}
		}
		
		int posHitCount = seq2kmers.size();
		int negHitCount = seq2kmers_neg.size();
		double hgp_all = computeHGP(posHitCount, negHitCount);
		
		while(changed){
			changed = false;
			ArrayList<Kmer> kmers_toRemove = new ArrayList<Kmer>();

			for (int i=0;i<kmers.size();i++){
				
				Kmer km = kmers.get(i);
				int count_with_single_kmer = 0;
				int count_with_single_kmer_neg = 0;
				HashSet<Integer> hits_to_remove = new HashSet<Integer>();
				HashSet<Integer> hits_to_remove_neg = new HashSet<Integer>();
				HashSet<Integer> hits = km.getPosHits();
				HashSet<Integer> hits_neg = km.getNegHits();
				for (int h:hits){
					if (seq2kmers.get(h).size()==1){
						count_with_single_kmer++;
						hits_to_remove.add(h);
					}
				}
				for (int h:hits_neg){
					if (seq2kmers_neg.get(h).size()==1){
						count_with_single_kmer_neg++;
						hits_to_remove_neg.add(h);
					}
				}
				if (count_with_single_kmer_neg>0){
					double hgp_remove_this = computeHGP(posHitCount-count_with_single_kmer, negHitCount-count_with_single_kmer_neg);
					// test whether removing this k-mer will improve enrichment significance hgp
					if (hgp_remove_this<hgp_all){
//						System.err.println(String.format("%s: p%d n%d p-%d n-%d hgp: %.1f-->%.1f", 
//								km.toString(), posHitCount, negHitCount, count_with_single_kmer, count_with_single_kmer_neg,
//								hgp_all, hgp_remove_this));
						hgp_all = hgp_remove_this;
						// remove this k-mer 
						changed = true;
						kmers_toRemove.add(km);
						// remove the sequence hit containing only this kmer
						for (int h: hits_to_remove)
							seq2kmers.remove(h);
						for (int h: hits_to_remove_neg)
							seq2kmers_neg.remove(h);						
					}
				}
			}
			kmers.removeAll(kmers_toRemove);
		}
		kmers.trimToSize();
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
    
    private class Sequence implements Comparable<Sequence>{
		int id;				// original input id
		String seq;
		String rc;
		/** seq_seed: positin of sequence relative to the seed kmer position, the squence should be RC() to match the seed orientation */
		int pos=UNALIGNED;
		private boolean isForward = true;
		/** forward strand (as the original orientation) matches */
		HashMap<Kmer, HashSet<Integer>> fPos = new HashMap<Kmer, HashSet<Integer>>();		// forward
		/** reverse strand (as the original orientation) matches */
		HashMap<Kmer, HashSet<Integer>> rPos = new HashMap<Kmer, HashSet<Integer>>();		// reverse
		int totalCount;
		int maxCount;
		double score;
		
		private Sequence(String seq, int id){
			this.id = id;
			this.seq = seq;
			this.rc = SequenceUtils.reverseComplement(seq);
		}
		
		private void RC(){
			isForward = false;
		}
		private String getSeq(){
			return isForward?seq:rc;
		}
		private void resetAlignment(){
			pos = UNALIGNED;
			isForward = true;
		}
		private void removeAllKmers(Collection<Kmer> kmers){
			for (Kmer km:kmers){
				fPos.remove(km);
				rPos.remove(km);
			}
		}
		private void clearKmerIndex(){
			fPos.clear();
			rPos.clear();
		}
		private String getSeqStrand(boolean isForward){
			return isForward?seq:rc;
		}

		private HashMap<Kmer, HashSet<Integer>> getKmerPos(){
			return isForward?fPos:rPos;
		}
	
		public int compareTo(Sequence s) {					// descending count
			int max = compareByMaxCount(s);
			if (max!=0)
				return max;
			if(totalCount<s.totalCount){return(1);}
			else if(totalCount>s.totalCount){return(-1);}
			else return(0);
		}
		
		private int compareByMaxCount(Sequence s) {					// descending count
			if(maxCount<s.maxCount){return(1);}
			else if(maxCount>s.maxCount){return(-1);}
			else return(0);
		}
		
		public String toString(){
			return String.format("%d\t%d\t%s\t%d\t%s", id, maxCount, isForward?"F":"R", pos, isForward?seq:rc);
		}
	}

	
	public class KmerCluster implements Comparable<KmerCluster>{
		public int pos_pwm_seed;
		public int pos_BS_seed;
		public MotifThreshold ksmThreshold = new MotifThreshold();
		public WeightMatrix wm;
		
		int clusterId;
		int k;
		Kmer seedKmer;
		float[][] pfm;
		boolean pwmGoodQuality = false;
		double pwmThreshold;
		/** log10(hg p-value)*/
		double pwmThresholdHGP;
		int pwmNegHitCount;
		int pwmPosHitCount;
		ArrayList<Kmer> alignedKmers;			// The K-mer set motif, a set of aligned k-mers
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
			cluster.k = this.k;
			return cluster;
		}
		
		/**
		 * Sort by PWM HGP, then by KSM HGP
		 */
		public int compareTo(KmerCluster c) {					// ascending pwmThresholdHGP
			if(pwmThresholdHGP<c.pwmThresholdHGP){return(-1);}
			else if(pwmThresholdHGP>c.pwmThresholdHGP){return(+1);}
			else // if same pwmHGP
				if(ksmThreshold.hgp<c.ksmThreshold.hgp){return(-1);}
				else if(ksmThreshold.hgp>c.ksmThreshold.hgp){return(+1);}
				else return(0);
		}
		
		/**
		 * Sort by KSM HGP, then by PWM HGP
		 */
		public int compareToByKsmHGP(KmerCluster c) {					// ascending ksmThresholdHGP
			if(ksmThreshold.hgp<c.ksmThreshold.hgp){return(-1);}
			else if(ksmThreshold.hgp>c.ksmThreshold.hgp){return(+1);}
			else 	// if same ksmHGP
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
			return String.format("PWM %d: %s, hgp=1e%.1f", clusterId, wm!=null?WeightMatrix.getMaxLetters(wm):"----", pwmThresholdHGP);
		}
	}
	
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
	 * Option to use multi-thread to compute HGP<br>
	 * Approximate grid search to find best HGP, to reduce run time
	 */
	private MotifThreshold optimizePwmThreshold(WeightMatrix wm, String outName, double startingScore){
		double[] posSeqScores = new double[posSeqCount];
		double[] negSeqScores = new double[negSeqCount];
		for (int i=0;i<posSeqCount;i++){
			posSeqScores[i]=WeightMatrixScorer.getMaxSeqScore(wm, seqs[i]);
		}
		for (int i=0;i<negSeqCount;i++){
			negSeqScores[i]=WeightMatrixScorer.getMaxSeqScore(wm, seqsNegList.get(i));
		}

		
		int[] posIdx = StatUtil.findSort(posSeqScores);		// index of sequence after sorting the scores
		
		int startIdx = Arrays.binarySearch(posSeqScores, startingScore);
		if( startIdx < 0 ) { startIdx = -startIdx - 1; }
		
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
	private MotifThreshold estimatePwmThreshold(WeightMatrix wm, double wmScore){
		double[] posSeqScores = new double[posSeqCount];
		double[] negSeqScores = new double[negSeqCount];
		for (int i=0;i<posSeqCount;i++){
			posSeqScores[i]=WeightMatrixScorer.getMaxSeqScore(wm, seqs[i]);
		}
		Arrays.sort(posSeqScores);
		int startIdx = Arrays.binarySearch(posSeqScores, wmScore);
		if( startIdx < 0 ) { startIdx = -startIdx - 1; }
		
		for (int i=0;i<negSeqCount;i++){
			negSeqScores[i]=WeightMatrixScorer.getMaxSeqScore(wm, seqsNegList.get(i));
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
	
	/**
	 * Multi-thread implementaton to compute HGP scores in pararell<br>
	 * @param idxs	The indices of elements to compute
	 * @param poshits
	 * @param neghits
	 * @param hgps
	 * @return
	 */
	private Pair<Double, Integer> findBestScore(ArrayList<Integer> idxs, int[] poshits, int[] neghits, double[] hgps){
		if (config.multl_thread_hgp){
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
		}
		else{
			for (int i:idxs){
				if (i==0 && poshits[0]==posSeqCount)	//why?
	    			hgps[0]=0;
				hgps[i]=computeHGP(posSeqCount, negSeqCount, poshits[i], neghits[i]);
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
		/** the significance of the k-mer group score */
		public double hgp=0;		
	}
	
	/**
	 * Find the best score cutoff given the positive scores and negative scores<br>
	 * This should be applicable to both PWM and KSM
	 * @param posSeqScores	motif scanning score of positive sequences, should be in the same order
	 * @param negSeqScores motif scanning score of negative sequences
	 * @return
	 */
	private MotifThreshold optimizeThreshold(double[] posSeqScores, double[] negSeqScores){
		int[] posIdx = StatUtil.findSort(posSeqScores);		
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
		for (int i=posScores_u.length-1;i>=0;i--){
			if (poshits[i]*1.0 >= neghits[i]*1.5*posSeqScores.length/negSeqScores.length)	// posHit should be at least 2 fold
				idxs.add(i);		
		}
		if (idxs.isEmpty())
			idxs.add(posScores_u.length-1);
		
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
		
		return score;		
	}
	/**
	 * Grid search threshold of a Kmer Group Score using the positive/negative sequences<br>
	 * Compute the hyper-geometric p-value from number of pos/neg sequences that have the scores higher than the considered score.<br>
	 * The KSM k-mers are assumed to have been loaded into the Engine
	 * @returns the KSM score gives the most significant p-value.
	 */
	private MotifThreshold optimizeKsmThreshold(ArrayList<Sequence> seqList, ArrayList<Sequence> seqListNeg, ArrayList<Kmer> kmers){
//		if (config.verbose>1)
//			System.out.println(CommonUtils.timeElapsed(tic)+ ": Ksm threshold, start.");

		HashSet<Kmer> kmerSet = new HashSet<Kmer>();
		kmerSet.addAll(kmers);
		HashSet<Integer> kmerLinkedPosSeqIds = new HashSet<Integer>();
		HashSet<Integer> kmerLinkedNegSeqIds = new HashSet<Integer>();
		for (Kmer km:kmers){
			kmerLinkedPosSeqIds.addAll(km.getPosHits());
			kmerLinkedNegSeqIds.addAll(km.getNegHits());
		}
		
		double[] posSeqScores = new double[seqList.size()];
		double[] negSeqScores = new double[seqListNeg.size()];
		for (int i=0;i<seqList.size();i++){
			Sequence s = seqList.get(i);
			if (!kmerLinkedPosSeqIds.contains(s.id))
				continue;
			KmerGroup[] kgs = scanKsmMatches(s, kmerSet);
			if (kgs==null)
				posSeqScores[i]=0;
			else{
				Arrays.sort(kgs);
				posSeqScores[i]=-kgs[0].getHgp();	// use first kg, the best match	// score = -log10 hgp, becomes positive value
			}
		}
//		if (config.verbose>1)
//			System.out.println(CommonUtils.timeElapsed(tic)+ ": Ksm threshold, scanned Pos seqs.");

		for (int i=0;i<seqListNeg.size();i++){
			Sequence s = seqListNeg.get(i);
			if (!kmerLinkedNegSeqIds.contains(s.id))
				continue;
			KmerGroup[] kgs = scanKsmMatches(s, kmerSet);
			if (kgs==null)
				negSeqScores[i]=0;
			else{
				Arrays.sort(kgs);
				negSeqScores[i]=-kgs[0].getHgp();
			}
		}
//		if (config.verbose>1)
//			System.out.println(CommonUtils.timeElapsed(tic)+ ": Ksm threshold, done scan neg sequences.");
		MotifThreshold score = optimizeThreshold(posSeqScores, negSeqScores);
		return score;
	}

	/**
	 * Optimize threshold of a KSM Score using the positive/negative sequences<br>
	 * Compute the hyper-geometric p-value from number of pos/neg sequences that have the scores higher than the considered score.<br>
	 * Return the score gives the most significant p-value.<br>
	 * Multi-thread version
	 */
	// currently NOT USED!!!
	private MotifThreshold estimateKsmThreshold(String outName, boolean printKgcHgp){
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
		if (config.multl_thread_hgp){
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
		}
		else{
			for (int i:idxs){
				if (i==0 && poshits[0]==posSeqCount)	//why?
	    			hgps[0]=0;
				hgps[i]=computeHGP(posSeqCount, negSeqCount, poshits[i], neghits[i]);
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
			if (km instanceof GappedKmer){
				GappedKmer gk = (GappedKmer)km;
				for (Kmer sk: gk.getSubKmers()){
					String kmerStr = gk.isSeedOrientation()&&gk.getSubKmerOrientation(sk)?sk.kmerString:sk.kmerRC;
					sk.kmerStartOffset=gk.kmerStartOffset;
					str2kmer.put(kmerStr, sk);			// use individual sub-kmers
					tree.add(kmerStr.getBytes(), kmerStr);
				}
			}
			else{
				String kmerStr = km.isSeedOrientation()?km.kmerString:km.kmerRC;				
				str2kmer.put(kmerStr, km);
				tree.add(kmerStr.getBytes(), kmerStr);
			}
	    }
	    tree.prepare();
	    engineInitialized = true;
	}		
	
	/** 
	 * Search all k-mers in the sequence
	 * @param seq sequence string to search k-mers
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
		// the reverse compliment
		String seq_rc = SequenceUtils.reverseComplement(seq);
		searcher = tree.search(seq_rc.getBytes());
		while (searcher.hasNext()) {
			SearchResult result = (SearchResult) searcher.next();
			kmerFound.addAll(result.getOutputs());
		}
		
		// Aho-Corasick only gives the patterns (kmers) matched, need to search for positions
		// matches on negative strand are combined with matches on positive strand
		HashMap<Integer, ArrayList<Kmer>> result = new HashMap<Integer, ArrayList<Kmer>> ();
		for (Object o: kmerFound){
			String kmerStr = (String) o;
			Kmer kmer = str2kmer.get(kmerStr);
			ArrayList<Integer> pos = StringUtils.findAllOccurences(seq, kmerStr);
			for (int p: pos){
				int x = p-kmer.getKmerStartOffset();	// minus kmerShift to get the motif position
				if (!result.containsKey(x))
					result.put(x, new ArrayList<Kmer>());
				result.get(x).add(kmer);	
//				System.out.println(String.format("%s\toff=%d\tp=%d\tx=%d", kmerStr, kmer.getKmerStartOffset(), p,x));
			}
			ArrayList<Integer> pos_rc = StringUtils.findAllOccurences(seq_rc, kmerStr);
			for (int p: pos_rc){
				int x = p-kmer.getKmerStartOffset();	// motif position in seqRC
				x = seq.length()-1-x;		// convert to position in Seq
				if (!result.containsKey(x))
					result.put(x, new ArrayList<Kmer>());
				result.get(x).add(kmer);	
//				System.out.println(String.format("%s\toff=%d\tp=%d\tx=%d", kmerStr, kmer.getKmerStartOffset(), p,x));
			}
			k=k+0;
		}
		KmerGroup[] matches = new KmerGroup[result.keySet().size()];
		int idx = 0;
		for (int p:result.keySet()){
			KmerGroup kg = config.use_weighted_kmer ? new KmerGroup(result.get(p), p, seq_weights) : new KmerGroup(result.get(p), p);
			matches[idx]=kg;
			kg.setHgp(computeHGP(kg.getGroupHitCount(), kg.getGroupNegHitCount()));
			idx++;
		}
		return matches;
	}
	
	/** 
	 * Search all k-mers (loaded in the AhoCorasick tree) in the sequence, both strand
	 * @param seq sequence string to search k-mers
	 * @return a set of kmers found
	 */
	public static HashSet<Kmer> queryTree (String seq, AhoCorasick tree){
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
	 * Search all k-mers in the sequence, strand-specific<br>
	 * Assuming the object variable "kmer tree" has been constructed.
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
		// the reverse compliment
		String seq_rc = SequenceUtils.reverseComplement(seq);
		searcher = tree.search(seq_rc.getBytes());
		while (searcher.hasNext()) {
			SearchResult result = (SearchResult) searcher.next();
			kmerFound.addAll(result.getOutputs());
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
			ArrayList<Integer> pos_rc = StringUtils.findAllOccurences(seqRC, kmerStr);
			for (int p: pos_rc){
				int x = p-kmer.getKmerStartOffset();	// motif position in seqRC
				x += RC;									// label it as "found on RC"
				if (!result.containsKey(x))
					result.put(x, new ArrayList<Kmer>());
				result.get(x).add(kmer);	
			}
		}
		KmerGroup[] matches = new KmerGroup[result.keySet().size()];
		int idx = 0;
		for (int p:result.keySet()){
			KmerGroup kg = config.use_weighted_kmer ? new KmerGroup(result.get(p), p, seq_weights) : new KmerGroup(result.get(p), p);
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
		
		if (config.verbose>1)
			System.out.println(CommonUtils.timeElapsed(tic)+": Merge "+originalCount+" clusters to "+clusters.size()+" clusters.");
	
		return !merged.isEmpty();
	}

	/**
	 * This KmerGroup class is used for recording the overlapping kmer instances mapped to the same binding position in a sequence
	 * @author yuchun
	 */
	public class KmerGroup implements Comparable<KmerGroup>{
		ArrayList<Kmer> kmers;
		int bs = 999;
		int clusterId = -1;
		int posHitGroupCount;
		int negHitGroupCount;
		boolean isKmersSorted=false;
		/** hgp (log10) using the positive/negative sequences */
		double hgp;
//		public KmerGroup(ArrayList<Kmer> kmers, int bs, int old){
//			this.bs = bs;
//			this.kmers = kmers;
//    		HashSet<Integer> allPosHits = new HashSet<Integer>();
//     		HashSet<Integer> allNegHits = new HashSet<Integer>();
//     		for (Kmer km:kmers){
//        		allPosHits.addAll(km.getPosHits());
//        		allNegHits.addAll(km.getNegHits());
//    		}
//    		posHitGroupCount = allPosHits.size();
//    		negHitGroupCount = allNegHits.size();
//		}
		public KmerGroup(ArrayList<Kmer> kmers, int bs){
			this.bs = bs;
			this.kmers = kmers;
			BitSet b_pos = new BitSet(posSeqCount);
			BitSet b_neg = new BitSet(negSeqCount);
     		for (Kmer km:kmers){
     			b_pos.or(km.posBits);
     			b_neg.or(km.negBits);
    		}
     		posHitGroupCount = b_pos.cardinality();
    		negHitGroupCount = b_pos.cardinality();
		}		
		public KmerGroup(ArrayList<Kmer> kmers, int bs, double[]weights){
			this.bs = bs;
			this.kmers = kmers;			
			BitSet b_pos = new BitSet(posSeqCount);
			BitSet b_neg = new BitSet(negSeqCount);
     		for (Kmer km:kmers){
     			b_pos.or(km.posBits);
     			b_neg.or(km.negBits);
    		}

    		if (weights==null){
        		posHitGroupCount = b_pos.cardinality();
    		}
    		else{
	    		double weight=0;
	    		for (int i = b_pos.nextSetBit(0); i >= 0; i = b_pos.nextSetBit(i+1))
	    			weight+=weights[i];
	    		posHitGroupCount = (int)(weight);
    		}
    		negHitGroupCount = b_neg.cardinality();
		}
//		public KmerGroup(ArrayList<Kmer> kmers, int bs, double[]weights){
//			this.bs = bs;
//			this.kmers = kmers;
//			
//    		HashSet<Integer> allPosHits = new HashSet<Integer>();
//     		HashSet<Integer> allNegHits = new HashSet<Integer>();
//     		for (Kmer km:kmers){
//        		allPosHits.addAll(km.getPosHits());
//        		allNegHits.addAll(km.getNegHits());
//    		}
//    		
//    		if (weights==null){
//        		posHitGroupCount = allPosHits.size();
//    		}
//    		else{
//	    		double weight=0;
//	    		for (int i: allPosHits)
//	    			weight+=weights[i];
//	    		posHitGroupCount = (int)(weight);
//    		}
//    		negHitGroupCount = allNegHits.size();
//		}
		
		/** hgp (log10) using the positive/negative sequences */
		public double getHgp() {return hgp;	}
		/** hgp (log10) using the positive/negative sequences */
		public void setHgp(double hgp) {this.hgp = hgp;	}
		public int getClusterId(){return clusterId;}
		public void setClusterId(int id){clusterId = id;}
		
		public ArrayList<Kmer> getKmers(){
			return kmers;
		}
		public Kmer getBestKmer(){
			if (!isKmersSorted){
				Collections.sort(this.kmers);
				isKmersSorted = true;
			}
			return kmers.get(0);
		}
//		public int getTotalKmerCount(){
//    		int kmerCountSum = 0;
//    		for (Kmer kmer:kmers){
//        		kmerCountSum+=kmer.getPosHitCount();	
//    		}
//    		return kmerCountSum;
//		}
		/** Get the number of sequences hit by any kmer in the group */
		public int getGroupHitCount(){
			return posHitGroupCount;
		}
		public int getGroupNegHitCount(){
			return negHitGroupCount;
		}
		public double getTotalKmerStrength(){
    		double total = 0;
    		for (Kmer kmer:kmers){
    			// on first kpp round, kmers do not have strength value, use count here
    			total+=kmer.getStrength()>1?kmer.getStrength():kmer.getPosHitCount();	
    		}
    		return total;
		}	
		/** Get the weighted kmer strength<cr>
		 *  The weight is 1 for top kmer, 1/k for other kmer 
		 *  Note: this is only approximate */
		public double getWeightedKmerStrength(){
    		double total = kmers.get(0).getWeightedHitCount();
    		double k = kmers.get(0).getK();
    		for (int i=1;i<kmers.size();i++){
    			total+=kmers.get(i).getWeightedHitCount()/k;	
    		}
    		return total;
		}			
		public int getPosBS(){
			return bs;
		}
		public int compareToByPosHitCount(KmerGroup kg) {		// descending pos hit count
			if(posHitGroupCount>kg.getGroupHitCount()){return(-1);}
			else if(posHitGroupCount<kg.getGroupHitCount()){return(1);}
			else return(0);
		}
		public int compareTo(KmerGroup kg) {					// ascending hgp
			if(hgp<kg.getHgp()){return(-1);}
			else if(hgp>kg.getHgp()){return(1);}
			else return(0);
		}
		public String toString(){
			return String.format("%s|%b %d: %d+/%d-, hpg=%.2f", getBestKmer().getKmerStrRC(), getBestKmer().isSeedOrientation(), bs, posHitGroupCount, negHitGroupCount, hgp);
		}
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
	    KMAC2WK kmf = new KMAC2WK();
	    kmf.plotMotifDistanceDistribution(x, same, diff, args[0]+".png");
	}
	
	// options cMyc_cMyc cMyc_HeLa_61bp_GEM.fasta cMyc_HeLa_61bp_GEM_neg.fasta 5 8 CCACGTG
	public static void main(String[] args){
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
			System.err.println("Example: KMAC2 --pos_seq c-Myc_Crawford_HeLa-S3_61bp_GEM.fasta [--neg_seq c-Myc_Crawford_HeLa-S3_61bp_GEM_neg.fasta] --k_min 5 --k_max 8 --out_name cMyc_cMyc --seed CACGTG");
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
	        				seq_w.add(1.0);
	        			}
	        		}
		            else
		            	seq_w.add(1.0);
	        	}
	        	else{
	        		if (config.k_win==-1){
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
						if (config.k_win==-1){
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
		KMAC2WK kmf = new KMAC2WK();

        kmf.setConfig(config, out_prefix);
        kmf.setSequences(pos_seqs, neg_seqs, seq_w);
        kmf.setStandalone();

        System.out.println(String.format("%d input positive sequences, use top %d center sequences (%dbp) to find motif ...", pos_seqs.size(), kmf.seqs.length, config.k_win));

        kmf.discoverMotifs(config.k_min, config.k_max, null);
        
		System.out.println(StatUtil.cacheAccessCount);
		System.out.println(StatUtil.getCacheSize());

        System.out.println("Done: "+CommonUtils.timeElapsed(tic));
	}
	

	
} 


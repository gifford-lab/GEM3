package edu.mit.csail.cgs.ewok.verbs.motifs;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import net.sf.samtools.util.SequenceUtil;

import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.ArrayUtils;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.utils.strings.StringUtils;
import edu.mit.csail.cgs.utils.strings.multipattern.*;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.features.ComponentFeature;
import edu.mit.csail.cgs.deepseq.features.Feature;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public class KmerEngine {
	private Genome genome;
	private boolean engineInitialized =false;
	private int k=10;
	private int minHitCount = 3;
	private int numPos;
	
	private String[] seqs;		// DNA sequences around binding sites
	private String[] seqsNeg;	// DNA sequences in negative sets
	private ArrayList<String> seqsNegList=new ArrayList<String>(); // Effective negative sets, excluding overlaps in positive sets
	public int getNegSeqCount(){return seqsNegList.size();}
	private int negRegionDistance;
	/** region-->index for negative sequences */
	private TreeMap<Region, Integer> neg_region_map;
	public String[] getPositiveSeqs(){return seqs;};
	public double get_NP_ratio(){return (double)seqsNegList.size()/seqs.length;}
	
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
		
	public KmerEngine(Genome g, boolean useCache){
		genome = g;
		seqgen = new SequenceGenerator<Region>();
		if (useCache)
			seqgen.useCache(true);		
	}
	/* 
	 * Contruct a Kmer Engine from a list of Kmers
	 */
	public KmerEngine(ArrayList<Kmer> kmers, String outPrefix){
		if (!kmers.isEmpty()){
			if (outPrefix!=null)
				updateEngine(kmers, outPrefix);
			else
				updateEngine(kmers);
			k=kmers.get(0).k;
		}
	}
	/**
	 * Set up the light weight genome cache. Only load the sequences for the specified regions.<br>
	 * At the same time, retrieve negative sequences (for only once)
	 * @param regions
	 */
	public double setupRegionCache(ArrayList<Region> cacheRegions, ArrayList<Region> negativeRegions, int negRegionDistance){
		this.negRegionDistance = negRegionDistance;
		double gcRatio=0;
		if (!seqgen.isRegionCached()){
			seqsNeg = seqgen.setupRegionCache(cacheRegions, negativeRegions);
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
		}
		return gcRatio;
	}
	/*
	 * Find significant Kmers that have high HyperGeometric p-value
	 * in sequence around the binding events 
	 * and build the kmer AhoCorasick engine
	 */
	public void buildEngine(int k, ArrayList<Point> events, int winSize, double hgp, double k_fold, String outPrefix){
		ArrayList<Kmer> kms = selectEnrichedKmers(k, events, winSize, hgp, k_fold, outPrefix);
		updateEngine(kms, outPrefix);		
	}
	
	public ArrayList<Kmer> selectEnrichedKmers(int k, ArrayList<Point> events, int winSize, double hgp, double k_fold, String outPrefix){
		cern.jet.random.engine.RandomEngine randomEngine = new cern.jet.random.engine.MersenneTwister();
		this.k = k;
		numPos = (winSize+1)-k+1;
		tic = System.currentTimeMillis();
		int eventCount = events.size();
		Collections.sort(events);		// sort by location
		// expected count of kmer = total possible unique occurence of kmer in sequence / total possible kmer sequence permutation
		int expectedCount = (int) Math.round(eventCount / Math.pow(4, k));

		// collect pos/neg test sequences based on event positions
		loadTestSequences(events, winSize);
		
		HashMap<String, HashSet<Integer>> kmerstr2seqs = new HashMap<String, HashSet<Integer>>();
		for (int seqId=0;seqId<seqs.length;seqId++){
			String seq = seqs[seqId];
			HashSet<String> uniqueKmers = new HashSet<String>();			// only count repeated kmer once in a sequence
			for (int i=0;i<numPos;i++){
				if ((i+k)>seq.length()) // endIndex of substring is exclusive
					break;
				uniqueKmers.add(seq.substring(i, i+k));
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

		// create the kmer object
		for (String key:kmerstr2seqs.keySet()){	
			if (kmerstr2seqs.get(key).size()< Math.max(expectedCount, minHitCount))
				continue;	// skip low count kmers 
			Kmer kmer = new Kmer(key, kmerstr2seqs.get(key));
			kms.add(kmer);
		}
		kmerstr2seqs=null;
		System.gc();
		System.out.println("k="+k+", mapped "+kms.size()+" k-mers, "+CommonUtils.timeElapsed(tic));
		
		/**
		 * Select significantly over-representative kmers 
		 * Search the kmer counts in the negative sequences, then compare to positive counts
		 */
		tic = System.currentTimeMillis();
		//Aho-Corasick for searching Kmers in negative sequences
		//ahocorasick_java-1.1.tar.gz is an implementation of Aho-Corasick automata for Java. BSD license.
		//from <http://hkn.eecs.berkeley.edu/~dyoo/java/index.html> 
		AhoCorasick tmp = new AhoCorasick();
		for (Kmer km: kms){
			tmp.add(km.kmerString.getBytes(), km.kmerString);
	    }
		tmp.prepare();
		
		// count hits in the negative sequences
		HashMap<String, HashSet<Integer>> kmerstr2negSeqs = new HashMap<String, HashSet<Integer>>();
		for (int negSeqId=0; negSeqId<seqsNegList.size();negSeqId++){
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
		int posSeq = seqs.length;
		int negSeq = seqsNegList.size();
		
		ArrayList<Kmer> toRemove = new ArrayList<Kmer>();
		ArrayList<Kmer> highHgpKmers = new ArrayList<Kmer>();
		for (Kmer kmer:kms){
			if (kmer.getPosHitCount()<=1){
				toRemove.add(kmer);	
				continue;
			}
			if (kmerstr2negSeqs.containsKey(kmer.kmerString)){
				kmer.setNegHits(kmerstr2negSeqs.get(kmer.kmerString));
			}
			if (kmer.getPosHitCount() < kmer.getNegHitCount()/get_NP_ratio() * k_fold ){
				highHgpKmers.add(kmer);	
				continue;
			}
			kmer.hgp_lg10 = computeHGP(posSeq, negSeq, kmer.getPosHitCount(), kmer.getNegHitCount());
			if (kmer.hgp_lg10>hgp)
				highHgpKmers.add(kmer);		
		}
		// remove un-enriched kmers		
		kms.removeAll(toRemove);
		Collections.sort(kms);
		Kmer.printKmers(kms, outPrefix+"_all_w"+winSize, true);
		
		kms.removeAll(highHgpKmers);
		System.out.println(String.format("k=%d, selected %d k-mers from %d+/%d- sequences, %s", k, kms.size(), posSeq, negSeq, CommonUtils.timeElapsed(tic)));
		
		return kms;
	}

	/**
	 * Load pos/neg test sequences based on event positions
	 * @param events
	 * @param winSize
	 * @param winShift
	 */
	public void loadTestSequences(ArrayList<Point> events, int winSize){
		int eventCount = events.size();
		seqs = new String[eventCount];	// DNA sequences around binding sites
		ArrayList<Region> posImpactRegion = new ArrayList<Region>();			// to make sure negative region is not within negRegionDistance of positive regions.

		for(int i=0;i<eventCount;i++){
			Region posRegion = events.get(i).expand(winSize/2);
			seqs[i] = seqgen.execute(posRegion).toUpperCase();
			posImpactRegion.add(events.get(i).expand(negRegionDistance));
		}
		/** Negative sequences has been retrieved when setting up region caches */
		// need to filter out those overlap with positive sequences
		ArrayList<Region> negRegions = new ArrayList<Region>();
		negRegions.addAll(neg_region_map.keySet());
		Region.filterOverlapRegions(negRegions, posImpactRegion);				// to make sure negative region is not within negRegionDistance of positive regions.
		seqsNegList.clear();
		for (Region r:negRegions){
			seqsNegList.add(seqsNeg[neg_region_map.get(r)].substring(0, winSize+1));
		}
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
	
	public double computeHGP(int posHitCount, int negHitCount){
	    int posSeqCount = seqs.length;
	    int negSeqCount = seqsNegList.size();
		return computeHGP(posSeqCount, negSeqCount, posHitCount, negHitCount);
	}
	/**
	 * Compute hgp using the positive/negative sequences
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
	    int posSeqCount = seqs.length;
	    int negSeqCount = seqsNegList.size();
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
	    for (int i=0;i<seqs.length;i++){
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
	    for (int i=0;i<seqsNegList.size();i++){
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
//			km.setStrength(kmerStrength[i]);
	    }
	}
	
	/**
	 * Estimate threshold of a PWM using the positive/negative sequences<br>
	 * Multi-thread to compute HGP<br>
	 * Do not consider negative score. (to reduce run time; negative pwm score means the PWM is of bad quality anyway)
	 */
	public MotifThreshold estimatePwmThreshold(WeightMatrix wm, String outName, boolean printPwmHgp){
		WeightMatrixScorer scorer = new WeightMatrixScorer(wm);
		double[] posSeqScores = new double[seqs.length];
		double[] negSeqScores = new double[seqsNegList.size()];
		for (int i=0;i<seqs.length;i++){
			posSeqScores[i]=scorer.getMaxSeqScore(wm, seqs[i]);
		}
		Arrays.sort(posSeqScores);
		int zeroIdx = Arrays.binarySearch(posSeqScores, 0);
		if( zeroIdx < 0 ) { zeroIdx = -zeroIdx - 1; }
		
		for (int i=0;i<seqsNegList.size();i++){
			negSeqScores[i]=scorer.getMaxSeqScore(wm, seqsNegList.get(i));
		}
		Arrays.sort(negSeqScores);
		
		TreeSet<Double> posScoreUnique = new TreeSet<Double>();
		for (int i=zeroIdx;i<posSeqScores.length;i++)
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
			poshits[i] = posSeqScores.length-index;
			index = CommonUtils.findKey(negSeqScores, key);
			neghits[i] = negSeqScores.length-index;
		}

		int numThread = java.lang.Runtime.getRuntime().availableProcessors();
		Thread[] threads = new Thread[numThread];
		ArrayList<Integer> idxs = new ArrayList<Integer>();
		for (int i=posScores_u.length-1;i>=0;i--)
			idxs.add(i);
		for (int i=0;i<numThread;i++){
            Thread t = new Thread(new HGPThread(idxs, seqs.length, seqsNegList.size(), poshits, neghits, hgps));
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
		
		StringBuilder sb = new StringBuilder();
		if (printPwmHgp){
			for (int i=0;i<posScores_u.length;i++)
				sb.append(String.format("%d\t%.2f\t%d\t%d\t%.1f\n", i, posScores_u[i], poshits[i], neghits[i], hgps[i] ));
		}
		
		Pair<Double, TreeSet<Integer>> minHgp = StatUtil.findMin(hgps);
		int minIdx = minHgp.cdr().last();
		MotifThreshold score = new MotifThreshold();
		score.score = posScores_u[minIdx];
		score.hgp = minHgp.car();
		score.posHit = poshits[minIdx];
		score.negHit = neghits[minIdx];
		if (printPwmHgp)
			CommonUtils.writeFile(outName+"_"+WeightMatrix.getMaxLetters(wm)+"_PwmHgp.txt", sb.toString());
		return score;
	}

	public class MotifThreshold{
		public double score;
		public int posHit;
		public int negHit;
		public double hgp;		
	}

	/**
	 * Estimate threshold of a Kmer Group Score using the positive/negative sequences<br>
	 * Compute the hyper-geometric p-value from number of pos/neg sequences that have the scores higher than the considered score.<br>
	 * Return the score gives the most significant p-value.
	 */
	public MotifThreshold estimateKgsThreshold(String outName, boolean printKgcHgp){
		double[] posSeqScores = new double[seqs.length];
		double[] negSeqScores = new double[seqsNegList.size()];
		for (int i=0;i<seqs.length;i++){
			KmerGroup[] kgs = query(seqs[i]);
			if (kgs.length==0)
				posSeqScores[i]=0;
			else{
				Arrays.sort(kgs);
				posSeqScores[i]=-kgs[0].getHgp();				// score = -log10 hgp, becomes positive value
			}
		}
		Arrays.sort(posSeqScores);
		for (int i=0;i<seqsNegList.size();i++){
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

		int numThread = java.lang.Runtime.getRuntime().availableProcessors();
		Thread[] threads = new Thread[numThread];
		ArrayList<Integer> idxs = new ArrayList<Integer>();
		for (int i=posScores_u.length-1;i>=0;i--)
			idxs.add(i);
		for (int i=0;i<numThread;i++){
            Thread t = new Thread(new HGPThread(idxs, seqs.length, seqsNegList.size(), poshits, neghits, hgps));
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
		Kmer.printKmers(kmers, outPrefix, false);
		
		//Aho-Corasick for searching Kmers in sequences
		//ahocorasick_java-1.1.tar.gz is an implementation of Aho-Corasick automata for Java. BSD license.
		//from <http://hkn.eecs.berkeley.edu/~dyoo/java/index.html> 
		tree = new AhoCorasick();
		for (Kmer km: kmers){
			str2kmer.put(km.kmerString, km);
			tree.add(km.kmerString.getBytes(), km.kmerString);
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
		for (Kmer km: kmers){
			str2kmer.put(km.kmerString, km);
			tree.add(km.kmerString.getBytes(), km.kmerString);
	    }
	    tree.prepare();
	    engineInitialized = true;
	}		
	
	/** 
	 * Search all k-mers in the sequence
	 * @param seq sequence string to search k-mers
	 * @return an array of KmerGroups:<br>
	 * Each k-mer group maps to a binding position in the sequence
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
				x = seq.length()-1-x;		// convert to position in Seq
				if (!result.containsKey(x))
					result.put(x, new ArrayList<Kmer>());
				result.get(x).add(kmer);	
			}
		}
		KmerGroup[] matches = new KmerGroup[result.keySet().size()];
		int idx = 0;
		for (int p:result.keySet()){
			KmerGroup kg = new KmerGroup(result.get(p), p);
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

	public String getSequence(Region r){
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
			map.put(kmer.kmerString, 0);
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
						String rc = SequenceUtil.reverseComplement(s);
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
			sb.append(km.kmerString).append("\t").append(map.get(km.kmerString)).append("\n");
		}
		CommonUtils.writeFile(genome.getVersion()+"_kmers_"+k+".txt", sb.toString());
		System.out.println(CommonUtils.timeElapsed(tic));
	}

	/**
	 * This KmerGroup class is used for recording the overlapping kmer instances mapped to the same binding position in a sequence
	 * @author yuchun
	 */
	public class KmerGroup implements Comparable<KmerGroup>{
		ArrayList<Kmer> kmers;
		int bs = 999;
		int posHitGroupCount;
		int negHitGroupCount;
		double hgp;
		public KmerGroup(ArrayList<Kmer> kmers, int bs){
			this.bs = bs;
			this.kmers = kmers;
			Collections.sort(this.kmers);
			
    		HashSet<Integer> allPosHits = new HashSet<Integer>();
    		for (int i=0;i<kmers.size();i++){
        		allPosHits.addAll(kmers.get(i).getPosHits());
    		}
    		posHitGroupCount = allPosHits.size();
    		
    		HashSet<Integer> allNegHits = new HashSet<Integer>();
    		for (int i=0;i<kmers.size();i++){
        		allNegHits.addAll(kmers.get(i).getNegHits());
    		}
    		negHitGroupCount = allNegHits.size();
		}
		public double getHgp() {return hgp;	}
		public void setHgp(double hgp) {this.hgp = hgp;	}
		
		public ArrayList<Kmer> getKmers(){
			return kmers;
		}
		public Kmer getBestKmer(){
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
		 *  The weight is 1 for top kmer, 1/k for other kmer */
		public double getWeightedKmerStrength(){
    		double total = kmers.get(0).getStrength();
    		double k = kmers.get(0).k;
    		for (int i=1;i<kmers.size();i++){
    			total+=kmers.get(i).getStrength()/k;	
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
				hgps[i]=computeHGP(posTotal, negTotal, posHits[i], negHits[i]);
			}
		}
	}
	public static void main0(String[] args){
		Genome g = null;
		ArgParser ap = new ArgParser(args);
	    try {
		      Pair<Organism, Genome> pair = Args.parseGenome(args);
		      if(pair==null){
		        if(ap.hasKey("geninfo")){
		          g = new Genome("Genome", new File(ap.getKeyValue("geninfo")));
		        }else{
	              System.err.println("No genome provided; provide a Gifford lab DB genome name or a file containing chromosome name/length pairs.");;System.exit(1);
	            }
		      }else{
		        g = pair.cdr();
		      }
	    } catch (NotFoundException e) {
	      e.printStackTrace();
	    }
	    KmerEngine kEngine = new KmerEngine(g, false);
	    List<File> files = Args.parseFileHandles(args, "kmers_file");
	    kEngine.indexKmers(files);
	}
	public static void main1(String[] args){
		KmerEngine ke = new KmerEngine(new ArrayList<Kmer>(),"");
		System.err.println(ke.computeHGP(652,665,652,665));
		System.err.println(ke.computeHGP(50,1112,40,357));
	}
}


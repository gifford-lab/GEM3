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
	private HashMap<String, Integer> negKmerHitCounts;
	
	String[] seqs;		// DNA sequences around binding sites
	String[] seqsNeg;	// DNA sequences in negative sets
	ArrayList<String> seqsNegList=new ArrayList<String>(); // Effective negative sets, excluding overlaps in positive sets
	/** region-->index for negative sequences */
	TreeMap<Region, Integer> neg_region_map;
	public String[] getPositiveSeqs(){return seqs;};
	public double get_NP_ratio(){return (double)seqsNegList.size()/seqs.length;}
	
	private ArrayList<Kmer> allKmers = new ArrayList<Kmer>();	// all the kmers in the sequences
	public ArrayList<Kmer> getAllKmers() {
		return allKmers;
	}
	private HashMap<String, Kmer> str2kmer = new HashMap<String, Kmer>();
	
	// AhoCorasick algorithm for multi-pattern search
	// Pre-processing is to build the tree with all the patters (kmers)
	// Then each individual search can be done in scan()
	private AhoCorasick tree;
	private AhoCorasick tree_negatives;
	
	private int maxShift;		// max kmer flanking Shift of whole dataset
	public int getMaxShift() {return maxShift;}
	public void setMaxShift(int maxShift) {this.maxShift = maxShift;}
	private int minShift;		// min kmer flanking Shift of whole dataset
	public int getMinShift() {return minShift;}
	public void setMinShift(int minShift) {this.minShift = minShift;}
	
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
				updateEngine(kmers, outPrefix, false);
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
	public double setupRegionCache(ArrayList<Region> cacheRegions, ArrayList<Region> negativeRegions){
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
		updateEngine(kms, outPrefix, true);		
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
		
		HashMap<String, Integer> map = new HashMap<String, Integer>();
		for (int seqId=0;seqId<seqs.length;seqId++){
			String seq = seqs[seqId];
			HashSet<String> uniqueKmers = new HashSet<String>();			// only count repeated kmer once in a sequence
			for (int i=0;i<numPos;i++){
				if ((i+k)>seq.length()) // endIndex of substring is exclusive
					break;
				uniqueKmers.add(seq.substring(i, i+k));
			}
			for (String s: uniqueKmers){
				if (map.containsKey(s)){
					 map.put(s, (map.get(s)+1));
				}
				else{
					 map.put(s, 1);
				}
			}
		}
		
		// sort the kmer strings, to make a kmer to also represent its reverse compliment (RC)	
		ArrayList<Kmer> kms = new ArrayList<Kmer>();
		ArrayList<String> sortedKeys = new ArrayList<String>();
		sortedKeys.addAll(map.keySet());
		Collections.sort(sortedKeys);	
		
		// create kmers from its and RC's counts
		for (String key:sortedKeys){
			if (!map.containsKey(key))		// this kmer has been removed, represented by RC
				continue;

			// consolidate kmer and its reverseComplment kmer
			String key_rc = SequenceUtils.reverseComplement(key);				
			if (!key_rc.equals(key)){	// if it is not reverse compliment itself
				if (map.containsKey(key_rc)){
					map.put(key, (map.get(key)+map.get(key_rc)));
					map.remove(key_rc);		// remove this kmer because it is represented by its RC
				}
			}
			
			if (map.get(key)< Math.max(expectedCount, minHitCount))
				continue;	// skip low count kmers, 
			
			// create the kmer object
			Kmer kmer = new Kmer(key, map.get(key));
			kms.add(kmer);
		}
		allKmers = new ArrayList<Kmer>(kms);		//TODO: make sure it does not take too much memory
		map=null;
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
		HashMap<String, Integer> negHitCounts = new HashMap<String, Integer>();
		for (String seq: seqsNegList){
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
				if (negHitCounts.containsKey(kmer))
					negHitCounts.put(kmer, negHitCounts.get(kmer)+1);
				else
					negHitCounts.put(kmer, 1);
			}
		}
		
		// score the kmers, hypergeometric p-value
		int posSeq = seqs.length;
		int negSeq = seqsNegList.size();
		
		ArrayList<Kmer> toRemove = new ArrayList<Kmer>();
		ArrayList<Kmer> highHgpKmers = new ArrayList<Kmer>();
		for (Kmer kmer:kms){
			if (kmer.seqHitCount<=1){
				toRemove.add(kmer);	
				continue;
			}
			if (negHitCounts.containsKey(kmer.kmerString)){
				kmer.negCount = negHitCounts.get(kmer.kmerString);
			}
			if (kmer.seqHitCount < kmer.negCount/get_NP_ratio() * k_fold ){
				highHgpKmers.add(kmer);	
				continue;
			}
			kmer.hgp = computeHGP(posSeq, negSeq, kmer.seqHitCount, kmer.negCount);
			if (kmer.hgp>hgp)
				highHgpKmers.add(kmer);		
		}
		// remove un-enriched kmers		
		kms.removeAll(toRemove);
		Collections.sort(kms);
		Kmer.printKmers(kms, outPrefix+"_all_w"+winSize, true);
		
		kms.removeAll(highHgpKmers);
		System.out.println(String.format("k=%d, selected %d k-mers from %d+/%d- sequences, %s", k, kms.size(), posSeq, negSeq, CommonUtils.timeElapsed(tic)));
		
		negKmerHitCounts = new HashMap<String, Integer>();
		for (Kmer km: highHgpKmers){
			negKmerHitCounts.put(km.kmerString, km.negCount);
		}
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
		ArrayList<Region> seqCoors = new ArrayList<Region>();

		for(int i=0;i<eventCount;i++){
			Region posRegion = events.get(i).expand(winSize/2);
			seqCoors.add(posRegion);
			seqs[i] = seqgen.execute(posRegion).toUpperCase();
		}
		/** Negative sequences has been retrieved when setting up region caches */
		// need to filter out those overlap with positive sequences
		ArrayList<Region> negRegions = new ArrayList<Region>();
		negRegions.addAll(neg_region_map.keySet());
		Region.filterOverlapRegions(negRegions, seqCoors);
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
	
	/**
	 * Compute hgp using the positive/negative sequences
	 */
	private double computeHGP(int posSeq, int negSeq, int posHit, int negHit){
		int allHit = posHit + negHit;
		int allSeq = posSeq + negSeq;
		// add one pseudo-count for negative set (zero count in negative set leads to tiny p-value)
		if (posHit<negHit){		// select smaller x for hyperGeometricCDF_cache(), to reduce # of x sum operations
			double hgcdf = StatUtil.hyperGeometricCDF_cache(posHit, allSeq, allHit, posSeq);
			return Math.log(1-hgcdf);
		}
		else{				// flip the problem, compute cdf of negative count
			double hgcdf=0;
			if (negHit==0)
				hgcdf = StatUtil.hyperGeometricCDF_cache(0, allSeq, allHit+2+1, negSeq);
			else
				hgcdf = StatUtil.hyperGeometricCDF_cache(negHit-1, allSeq, allHit, negSeq);
			if (hgcdf<Double.MIN_VALUE)
				return computeHGP_BIG(posSeq, negSeq, posHit, negHit);
			else
				return Math.log10(hgcdf);
		}
	}
	/**
	 * Compute hgp using the positive/negative sequences, high precision version
	 * Runs much slower, but more accurate for very small p-value (<MIN_VALUE, 2^-1074)
	 */
	private double computeHGP_BIG(int posSeq, int negSeq, int posHit, int negHit){
		int allHit = posHit + negHit;
		int allSeq = posSeq + negSeq;
		// add one pseudo-count for negative set (zero count in negative set leads to tiny p-value)
		if (posHit<negHit){		// select smaller x for hyperGeometricCDF_cache(), to reduce # of x sum operations
			double hgcdf = StatUtil.hyperGeometricCDF_cache(posHit, allSeq, allHit, posSeq);
			return Math.log(1-hgcdf);
		}
		else{				// flip the problem, compute cdf of negative count
			double hgcdf_log10=0;
			if (negHit==0)
				hgcdf_log10 = StatUtil.log10_hyperGeometricCDF_cache(0, allSeq, allHit+2+1, negSeq);
			else
				hgcdf_log10 = StatUtil.log10_hyperGeometricCDF_cache(negHit-1, allSeq, allHit, negSeq);
			return hgcdf_log10;
		}
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
	    int[] posHitCount = new int[kmers.size()];
	    int[] negHitCount = new int[kmers.size()];
	    double[] kmerStrength = new double[kmers.size()];
	    for (int i=0;i<seqs.length;i++){
	    	String seq = seqs[i];
			Iterator searcher = tree.search(seq.getBytes());
			while (searcher.hasNext()) {
				SearchResult result = (SearchResult) searcher.next();
				Set<Integer> idxs = result.getOutputs();
				int kmerSum = 0;
				for (int idx:idxs){
					kmerSum += kmers.get(idx).getSeqHitCount();
				}
				for (int idx:idxs){
					posHitCount[idx]++;
					kmerStrength[idx] += kmerSum==0?0:kmers.get(idx).getSeqHitCount()/kmerSum*events.get(i).getTotalEventStrength();
				}
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
	    for (int i=0;i<kmers.size();i++){
	    	Kmer km = kmers.get(i);
	    	km.setSeqHitCount(posHitCount[i]);
	    	km.setNegCount(negHitCount[i]);
			km.setHgp( updateHGP(posHitCount[i], negHitCount[i]));
			km.setStrength(kmerStrength[i]);
	    }
	}
	
	public double updateHGP(int posHitCount, int negHitCount){
	    int posSeqCount = seqs.length;
	    int negSeqCount = seqsNegList.size();
		return computeHGP(posSeqCount, negSeqCount, posHitCount, negHitCount);
	}
	
	/**
	 * Estimate threshold of a PWM using the positive/negative sequences<br>
	 * If a PWM is good PWM, the curve of the difference between the number of positive and negative sequences match vs score
	 * should have a peak value, then we set PWM threshold = the largest score corresponding to 0.9*peak_value
	 */
	public Pair<Double, Double> estimatePwmThreshold(WeightMatrix wm, String outName, boolean printFDR){
		WeightMatrixScorer scorer = new WeightMatrixScorer(wm);
		double[] posSeqScores = new double[seqs.length];
		double[] negSeqScores = new double[seqsNegList.size()];
		for (int i=0;i<seqs.length;i++){
			posSeqScores[i]=scorer.getMaxSeqScore(wm, seqs[i].toCharArray());
		}
		Arrays.sort(posSeqScores);
		int zeroIdx = Arrays.binarySearch(posSeqScores, 0);
		if( zeroIdx < 0 ) { zeroIdx = -zeroIdx - 1; }
		
		for (int i=0;i<seqsNegList.size();i++){
			negSeqScores[i]=scorer.getMaxSeqScore(wm, seqsNegList.get(i).toCharArray());
		}
		Arrays.sort(negSeqScores);
		
		// find the threshold motif score
		double[] hgps = new double[seqs.length];
		double diffs[] = new double[seqs.length];
		double fdrs[] = new double[seqs.length];
		StringBuilder sb = new StringBuilder();
		for (int i=zeroIdx;i<seqs.length;i++){
			int index = Arrays.binarySearch(negSeqScores, posSeqScores[i]);
			if( index < 0 ) { index = -index - 1; }
			int positiveCount = posSeqScores.length-i;
			long negativeCount = Math.round((double)(negSeqScores.length-index)*seqs.length/seqsNegList.size());		// scale the count
			hgps[i]=computeHGP(seqs.length, seqsNegList.size(), posSeqScores.length-i, negSeqScores.length-index);
			fdrs[i] = (double)negativeCount/(positiveCount+negativeCount);
			diffs[i] = positiveCount-negativeCount;
			if (printFDR)
				sb.append(String.format("%d\t%.2f\t%d\t%d\t%.0f\t%.4f\t%.1f\n", i, posSeqScores[i], positiveCount, negativeCount, diffs[i], fdrs[i], hgps[i] ));
		}	
		hgps[0]=0;		// the lowest threshold will match all positive sequences, lead to hgp=0
		Pair<Double, TreeSet<Integer>> maxDiff = StatUtil.findMin(hgps);
		double minHGP = maxDiff.car();
		int minIdx = maxDiff.cdr().last();
		double threshold = posSeqScores[minIdx];
//		for (int i=maxIdx; i<diffs.length;i++){
//			if (diffs[i]<max*0.9){
//				threshold = posSeqScores[i];
//				break;
//			}
//		}
		System.out.println(String.format("%.2f\t%.0f\t%.4f\t%.1f\n", threshold, diffs[minIdx], fdrs[minIdx], hgps[minIdx] ));
		if (printFDR)
			CommonUtils.writeFile(outName+"_"+WeightMatrix.getMaxLetters(wm)+"_fdr.txt", sb.toString());
		return new Pair<Double, Double>(threshold, minHGP);
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
	public void updateEngine(ArrayList<Kmer> kmers, String outPrefix, boolean isOverlappedKmer){
		if (kmers.isEmpty()){
			engineInitialized = false;
			return;
		}
		Collections.sort(kmers);
		Kmer.printKmers(kmers, outPrefix, isOverlappedKmer);
		
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
	 * @return an array of KmerMatches:<br>
	 * pos is the binding site position in the sequence<br>
	 * kmers are the kmers that map to this position on both positive and negative strands<br>
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
			matches[idx]=new KmerGroup(result.get(p), p);
			idx++;
		}
		return matches;
	}
	
	/** 
	 * Search all k-mers in the sequence
	 * @param seq sequence string to search k-mers
	 * @return a kmers list<br>
	 * pos is the binding site position in the sequence<br>
	 * kmers are the kmers that map to this position<br>
	 * if pos is negative, then the kmer match is on the reverse compliment seq string
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
	 * This KmerMatch class is used for recording kmer instances in sequences in the conventional motif finding setting
	 * @author yuchun
	 */
	class KmerMatch {
		int seqId;			// sequence id in the dataset
		int pos;			// position in the se
		int isMinus;		// the kmer is on positive strand (0), or negative strand (1). use int (not boolean) for iteration.
		KmerMatch(int seqId, int pos, int isMinus){
			this.seqId = seqId;
			this.pos = pos;
			this.isMinus = isMinus;
		}
	}
	
	/**
	 * This KmerMatches class is used for recording the overlapping kmer instances mapped to the same binding position in a sequence
	 * @author yuchun
	 */
	public class KmerGroup {
		ArrayList<Kmer> kmers;
		int bs = 999;
		public KmerGroup(ArrayList<Kmer> kmers, int bs){
			this.bs = bs;
			this.kmers = kmers;
			Collections.sort(this.kmers);
		}
		public ArrayList<Kmer> getKmers(){
			return kmers;
		}
		public Kmer getBestKmer(){
			return kmers.get(0);
		}
		public int getTotalKmerCount(){
    		int kmerCountSum = 0;
    		for (Kmer kmer:kmers){
        		kmerCountSum+=kmer.getSeqHitCount();	
    		}
    		return kmerCountSum;
		}
		/** Get the weighted kmer count<cr>
		 *  The weight is 1 for top kmer, 1/k for other kmer */
		public int getWeightedKmerCount(){
    		float kmerCountSum = kmers.get(0).getSeqHitCount();
    		float k = kmers.get(0).k;
    		for (int i=1;i<kmers.size();i++){
        		kmerCountSum+=kmers.get(i).getSeqHitCount()/k;	
    		}
    		return Math.round(kmerCountSum);
		}
		public double getTotalKmerStrength(){
    		double total = 0;
    		for (Kmer kmer:kmers){
    			// on first kpp round, kmers do not have strength value, use count here
    			total+=kmer.getStrength()>1?kmer.getStrength():kmer.getSeqHitCount();	
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
}


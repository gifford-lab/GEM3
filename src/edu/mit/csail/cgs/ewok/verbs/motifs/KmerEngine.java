package edu.mit.csail.cgs.ewok.verbs.motifs;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.TreeSet;
import java.util.Vector;

import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.utils.strings.StringUtils;
import edu.mit.csail.cgs.utils.strings.multipattern.*;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.features.ComponentFeature;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public class KmerEngine {
	private Genome genome;
	private boolean engineInitialized =false;
	private int k=10;
	private int minHitCount = 3;
	private int numPos;
	
	private ArrayList<Kmer> kmers = new ArrayList<Kmer>();		// current set of kmers
	private ArrayList<Kmer> allKmers = new ArrayList<Kmer>();	// all the kmers in the sequences
	private HashMap<String, Kmer> str2kmer = new HashMap<String, Kmer>();
	
	// AhoCorasick algorithm for multi-pattern search
	// Pre-processing is to build the tree with all the patters (kmers)
	// Then each individual search can be done in scan()
	private AhoCorasick tree;
	
	private int maxCount;		// max kmer hit count of whole dataset
	public int getMaxCount() {return maxCount;}
	private int minCount;		// min kmer hit count of whole dataset
	public int getMinCount() {return minCount;}
	
	private int maxShift;		// max kmer flanking Shift of whole dataset
	public int getMaxShift() {return maxShift;}
	public void setMaxShift(int maxShift) {this.maxShift = maxShift;}
	private int minShift;		// min kmer flanking Shift of whole dataset
	public int getMinShift() {return minShift;}
	public void setMinShift(int minShift) {this.minShift = minShift;}
	
	public boolean isInitialized(){ return engineInitialized;}
	
	// The average profile/density of kmers along the sequence positions
	private double[] positionProbs;
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
			updateEngine(kmers, outPrefix);
			k=kmers.get(0).k;
		}
	}
	public void compactRegionCache(ArrayList<Region> regions){
		seqgen.compactRegionCache(regions);
	}
	/*
	 * Find significant Kmers that have high HyperGeometric p-value
	 * in sequence around the binding events 
	 * and build the kmer AhoCorasick engine
	 */
	public void buildEngine(int k, ArrayList<Point> events, int winSize, int winShift, double hgp, double k_fold, String outPrefix){
		this.k = k;
		numPos = (winSize+1)-k+1;
		tic = System.currentTimeMillis();
		int eventCount = events.size();
		Collections.sort(events);		// sort by location
		// input data
		String[]  seqs = new String[eventCount];	// DNA sequences around binding sites
		String[]  seqsNeg = new String[eventCount];	// DNA sequences in negative sets
//		double[] seqProbs = new double[eventCount];	// Relative strength of that binding site, i.e. ChIP-Seq read count
		Region[] seqCoors = new Region[eventCount];

		// expected count of kmer = total possible unique occurence of kmer in sequence / total possible kmer sequence permutation
		int expectedCount = (int) (eventCount / Math.pow(4, k) + 0.5);
		ArrayList<Region> negRegions = new ArrayList<Region>();
		// prepare for progress reporting
		int displayStep = (int) Math.pow(10, (int) (Math.log10(eventCount)));
		TreeSet<Integer> reportTriggers = new TreeSet<Integer>();
		for (int i=1;i<=eventCount/displayStep; i++){
			reportTriggers.add(i*displayStep);
		}
		reportTriggers.add(100);
		reportTriggers.add(1000);
		reportTriggers.add(10000);
		System.out.println("Retrieving sequences from "+eventCount+" binding event regions ... ");
		for(int i=0;i<eventCount;i++){
			Region posRegion = events.get(i).expand(winSize/2);
			seqCoors[i] = posRegion;
			seqs[i] = seqgen.execute(seqCoors[i]).toUpperCase();
			// getting negative sequences
			// exclude negative regions that overlap with positive regions, or exceed start of chrom
			// it is OK if we lose a few sequences here, so some entries of the seqsNeg will be null
			int start = posRegion.getStart()-winShift;
			int end = posRegion.getEnd()-winShift;
			if (start < 0)
				continue;
			if (i>0){
				if (seqCoors[i-1].overlaps( new Region(genome, posRegion.getChrom(), start, end)))
					continue;
			}
			Region negRegion = new Region(genome, posRegion.getChrom(), start, end);			
			negRegions.add(negRegion);
			int trigger = eventCount;
            if (!reportTriggers.isEmpty())
            	trigger = reportTriggers.first();
            if (i>trigger){
//				System.out.println(trigger+"\t/"+eventCount+"\t"+CommonUtils.timeElapsed(tic));
				reportTriggers.remove(reportTriggers.first());
            }
            seqsNeg[i] = seqgen.execute(negRegion).toUpperCase();
		}
//		System.out.println(eventCount+"\t/"+eventCount+"\t"+CommonUtils.timeElapsed(tic));
		
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
//		System.out.println("\nKmers indexed "+CommonUtils.timeElapsed(tic));
//		tic = System.currentTimeMillis();
		
		// sort the kmer strings, to make a kmer to also represent its reverse compliment (RC)	
		ArrayList<Kmer> kmers = new ArrayList<Kmer>();
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
			kmers.add(kmer);
		}
		allKmers = new ArrayList<Kmer>(kmers);		//TODO: make sure it does not take too much memory
		map=null;
		System.gc();
		System.out.println("Kmers("+kmers.size()+") mapped, "+CommonUtils.timeElapsed(tic));
		
		/**
		 * Select significantly over-representative kmers 
		 */
		tic = System.currentTimeMillis();
		/*
		Aho-Corasick for searching Kmers in negative sequences
		ahocorasick_java-1.1.tar.gz is an implementation of Aho-Corasick automata for Java. BSD license.
		from <http://hkn.eecs.berkeley.edu/~dyoo/java/index.html> 
		 */		
		AhoCorasick tmp = new AhoCorasick();
		for (Kmer km: kmers){
			tmp.add(km.kmerString.getBytes(), km.kmerString);
	    }
		tmp.prepare();
		
		// count hits in the negative sequences
		HashMap<String, Integer> negHitCounts = new HashMap<String, Integer>();
		int negSeqCount = 0;
		for (String seq: seqsNeg){
			if (seq==null)			// some neg seq may be null, if overlap with positive sequences
				continue;
			negSeqCount++;
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
		int n = seqs.length;
		int N = n + negSeqCount;
		
		ArrayList<Kmer> toRemove = new ArrayList<Kmer>();
		for (Kmer kmer:kmers){
			if (kmer.seqHitCount<=1){
				toRemove.add(kmer);	
				continue;
			}
			// add one pseudo-count for negative set (zero count in negative set leads to tiny p-value)
			int kmerAllHitCount = kmer.seqHitCount+1;
			if (negHitCounts.containsKey(kmer.kmerString)){
				kmer.negCount = negHitCounts.get(kmer.kmerString);
				kmerAllHitCount += kmer.negCount;
			}
			if (kmer.seqHitCount < kmer.negCount * k_fold){
				toRemove.add(kmer);	
				continue;
			}
			kmer.hg = 1-StatUtil.hyperGeometricCDF(kmer.seqHitCount, N, kmerAllHitCount, n);
//			System.out.println(String.format("%s\t%d\t%.4f", kmer.kmerString, kmer.seqHitCount, kmer.hg));
			if (kmer.hg>hgp)
				toRemove.add(kmer);		
		}
		// remove un-enriched kmers		
		kmers.removeAll(toRemove);
		System.out.println(String.format("Kmers(%d) selected from %d positive vs %d negative sequences, %s", kmers.size(), n, negSeqCount, CommonUtils.timeElapsed(tic)));

		// set Kmers and prepare the search Engine
		updateEngine(kmers, outPrefix);
	}
	
	/** load Kmers and prepare the search Engine
	 *  assuming the kmers are unique
	 * 
	 * @param kmers List of kmers (with kmerString, sequence hit count)
	 */
	public void updateEngine(ArrayList<Kmer> kmers, String outPrefix){
		if (kmers.isEmpty()){
			engineInitialized = false;
			return;
		}
		tic = System.currentTimeMillis();
		this.kmers = kmers;
		Collections.sort(kmers);
		this.maxCount = kmers.get(0).seqHitCount;
		this.minCount = kmers.get(kmers.size()-1).seqHitCount;
		Kmer.printKmers(kmers, outPrefix);
		
		/*
		Init Aho-Corasick alg. for searching multiple Kmers in sequences
		ahocorasick_java-1.1.tar.gz is an implementation of Aho-Corasick automata for Java. BSD license.
		from <http://hkn.eecs.berkeley.edu/~dyoo/java/index.html> 
		 */		
		tree = new AhoCorasick();
		for (Kmer km: kmers){
			str2kmer.put(km.kmerString, km);
			tree.add(km.kmerString.getBytes(), km.kmerString);
	    }
	    tree.prepare();
	    engineInitialized = true;
//	    System.out.println("Kmers("+kmers.size()+") loaded to the Kmer Engine, "+CommonUtils.timeElapsed(tic));
	}	
	
	/** 
	 * Search all kmers in the sequence
	 * @param seq
	 * @return a map (pos-->kmers)
	 * pos is the binding site position in the sequence
	 * kmers are the kmers that map to this position
	 * if pos is negative, then the kmer match is on the reverse compliment seq string
	 */
	public HashMap<Integer, ArrayList<Kmer>> query (String seq){
		seq = seq.toUpperCase();
		HashSet<Object> kmerFound = new HashSet<Object>();	// each kmer is only used 
		/*
		Search for all kmers in the sequences using Aho-Corasick algorithms (initialized)
		ahocorasick_java-1.1.tar.gz is an implementation of Aho-Corasick automata for Java. BSD license.
		from <http://hkn.eecs.berkeley.edu/~dyoo/java/index.html> 
		 */			
		Iterator searcher = tree.search(seq.getBytes());
//		System.out.println(seq);
		while (searcher.hasNext()) {
			SearchResult result = (SearchResult) searcher.next();
			kmerFound.addAll(result.getOutputs());
//			System.out.println(result.getOutputs()+"\tat index: " + result.getLastIndex());
		}
		// the reverse compliment
		String seq_rc = SequenceUtils.reverseComplement(seq);
		searcher = tree.search(seq_rc.getBytes());
		while (searcher.hasNext()) {
			SearchResult result = (SearchResult) searcher.next();
			kmerFound.addAll(result.getOutputs());
		}
		
		// Aho-Corasick only gives the patterns (kmers) matched, need to search for positions
		// negative postion --> the start position matches the rc of kmer
		HashMap<Integer, ArrayList<Kmer>> result = new HashMap<Integer, ArrayList<Kmer>> ();
		String seqRC = SequenceUtils.reverseComplement(seq);
		for (Object o: kmerFound){
			String kmerStr = (String) o;
			Kmer kmer = str2kmer.get(kmerStr);
			ArrayList<Integer> pos = StringUtils.findAllOccurences(seq, kmerStr);
			for (int p: pos){
				int x = p-kmer.getKmerShift();	// minus kmerShift to get the motif position
				if (!result.containsKey(x))
					result.put(x, new ArrayList<Kmer>());
				result.get(x).add(kmer);	
			}
			ArrayList<Integer> pos_rc = StringUtils.findAllOccurences(seqRC, kmerStr);
			for (int p: pos_rc){
				int x = p-kmer.getKmerShift();	// motif position in seqRC
				x = -(seq.length()-1-x);		// convert to position in Seq, "-" for reverse strand
				if (!result.containsKey(x))
					result.put(x, new ArrayList<Kmer>());
				result.get(x).add(kmer);	
			}
		}

		return result;
	}

	public String getSequence(Region r){
		return seqgen.execute(r).toUpperCase();
	}

	
	private float scorePWM(String seq){
		float score = 0;
//		char[] chars = seq.toCharArray();
//        for (int j = 0; j < motif.length(); j++) {
//            score += motif.matrix[j][chars[j]];
//        }
        return score;
	}
	

	/**
	 * This KmerMatch class is used for recording kmer instances in sequences in the conventional motif finding setting
	 * @author yuchun
	 *
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
}


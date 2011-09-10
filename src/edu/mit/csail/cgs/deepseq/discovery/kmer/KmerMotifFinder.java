package edu.mit.csail.cgs.deepseq.discovery.kmer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
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
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.SetTools;
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
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public class KmerMotifFinder {
	private final int RC=100000;		// extra bp add to indicate negative strand match of kmer
	private final int UNALIGNED=9999;	// the special shift for unaligned kmer
	private final char[] LETTERS = {'A','C','G','T'};
	private final int MAXLETTERVAL = Math.max(Math.max(Math.max('A','C'),Math.max('T','G')),
            Math.max(Math.max('a','c'),Math.max('t','g'))) + 1;
	
	private int verbose = 2;
	private Genome genome;
	private boolean engineInitialized =false;
	private int k;
	private int minHitCount = 3;
	private int numPos;
	private double[] bg= new double[4];	// background frequency based on GC content
	private double ic_trim = 0.4;
	private double motif_hit_factor=0.01;			// KSM or PWM
	private double motif_hit_factor_report=0.05;			// KSM or PWM
	private String outName;
	private boolean use_grid_search=true;
	private double seed_search_fraction = 0.2;
	private double kmer_set_overlap_ratio = 0.5;
	
	private double k_fold;
	private double hgp = -3; 	// p-value threshold of hyper-geometric test for enriched kmer 
	private Kmer bestSeed = null;
	private boolean select_seed=false;
	private int max_cluster = 30;
	
	private int k_win;
	private String[] seqs;		// DNA sequences around binding sites
	private String[] seqsNeg;	// DNA sequences in negative sets
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
		
	public KmerMotifFinder(){ }
	
	public void setParameters(double hgp, double k_fold, double motif_hit_factor, double motif_hit_factor_report, String outName, 
			boolean select_seed, boolean use_grid_search, int verbose, double seed_search_fraction, double kmer_set_overlap_ratio){
	    this.hgp = hgp;
	    this.k_fold = k_fold;	
	    this.motif_hit_factor = motif_hit_factor;
	    this.outName = outName;
	    this.select_seed = select_seed;
	    this.verbose = verbose;
	    this.seed_search_fraction = seed_search_fraction;
	    this.motif_hit_factor_report = motif_hit_factor_report;
	    this.kmer_set_overlap_ratio = kmer_set_overlap_ratio;
	    this.use_grid_search = use_grid_search;
	}
	
	public void setSequences(ArrayList<String> pos_seqs, ArrayList<String> neg_seqs){
		seqs = new String[pos_seqs.size()];	
		pos_seqs.toArray(seqs);
		seqsNegList = neg_seqs;
		posSeqCount = seqs.length;
	    negSeqCount = seqsNegList.size();
	    updateSequenceInfo();
	}
	private void updateSequenceInfo(){
		k_win = seqs[0].length();
	    
	    // count cg-content
		int gcCount = 0;
		for (String s:seqsNegList){
			for (char c:s.toCharArray())
				if (c=='C'||c=='G')
					gcCount ++;
		}
		double gcRatio = (double)gcCount/negSeqCount/k_win;
		bg[0]=0.5-gcRatio/2; 
    	bg[1]=gcRatio/2; 
    	bg[2]=bg[1]; 
    	bg[3]=bg[0];
	}
	
	public KmerMotifFinder(Genome g, boolean useCache){
		genome = g;
		seqgen = new SequenceGenerator<Region>();
		if (useCache)
			seqgen.useCache(true);		
	}
	/* 
	 * Contruct a Kmer Engine from a list of Kmers
	 */
	public KmerMotifFinder(ArrayList<Kmer> kmers, String outPrefix){
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
			bg[0]=0.5-gcRatio/2; 
        	bg[1]=gcRatio/2; 
        	bg[2]=bg[1]; 
        	bg[3]=bg[0];
		}
		return gcRatio;
	}
	/*
	 * Find significant Kmers that have high HyperGeometric p-value
	 * in sequence around the binding events 
	 * and build the kmer AhoCorasick engine
	 */
	public void buildEngine(int k, ArrayList<Point> events, int winSize, double hgp, double k_fold, String outPrefix, boolean print_kmer_hits){
		ArrayList<Kmer> kms = selectEnrichedKmers(k, events, winSize, hgp, k_fold, outPrefix);
		updateEngine(kms, outPrefix);		
	}
	
	public ArrayList<Kmer> selectEnrichedKmers(int k, ArrayList<Point> events, int winSize, double hgp, double k_fold, String outName){
		this.hgp = hgp;
		this.outName = outName;
		this.k_fold = k_fold;

		// collect pos/neg test sequences based on event positions
		loadTestSequences(events, winSize);
		
		return selectEnrichedKmers(k);
	}
	
	/** 
	 * Select the value of k by the coverage of enriched k-mers<br>
	 * Count the exact base of coverage from every sequence
	 * @param k_min
	 * @param k_max
	 * @return
	 */
	public int selectK(int k_min, int k_max){
		// compare different values of k to select most enriched k value			
		ArrayList<Integer> widths = new ArrayList<Integer>(); 
		// compare different values of k to select most enriched k value
		int bestK = 0;
		int bestCoverage = 0;
		for (int i=0;i<k_max-k_min+1;i++){
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
				HashSet<Kmer> results = queryTree (seq, oks);
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
				HashSet<Kmer> results = queryTree (seq, oks);
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
	public ArrayList<Kmer> selectEnrichedKmers(int k){
		this.k = k;
		// expected count of kmer = total possible unique occurence of kmer in sequence / total possible kmer sequence permutation
		tic = System.currentTimeMillis();
		numPos = k_win-k+1;
		int expectedCount = (int) Math.round(seqs.length / Math.pow(4, k));
		HashMap<String, HashSet<Integer>> kmerstr2seqs = new HashMap<String, HashSet<Integer>>();
		for (int seqId=0;seqId<posSeqCount;seqId++){
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
			tmp.add(km.getKmerString().getBytes(), km.getKmerString());
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
		ArrayList<Kmer> toRemove = new ArrayList<Kmer>();
		ArrayList<Kmer> highHgpKmers = new ArrayList<Kmer>();
		for (Kmer kmer:kms){
			if (kmer.getPosHitCount()<=1){
				toRemove.add(kmer);	
				continue;
			}
			if (kmerstr2negSeqs.containsKey(kmer.getKmerString())){
				kmer.setNegHits(kmerstr2negSeqs.get(kmer.getKmerString()));
			}
			if (kmer.getPosHitCount() < kmer.getNegHitCount()/get_NP_ratio() * k_fold ){
				highHgpKmers.add(kmer);	
				continue;
			}
			kmer.setHgp(computeHGP(posSeqCount, negSeqCount, kmer.getPosHitCount(), kmer.getNegHitCount()));
			if (kmer.getHgp()>hgp)
				highHgpKmers.add(kmer);		
		}
		// remove un-enriched kmers		
		kms.removeAll(toRemove);
		Collections.sort(kms);
		Kmer.printKmers(kms, posSeqCount, negSeqCount, outName+"_all_w"+seqs[0].length(), true, false);
		
		kms.removeAll(highHgpKmers);
		System.out.println(String.format("k=%d, selected %d k-mers from %d+/%d- sequences, %s", k, kms.size(), posSeqCount, negSeqCount, CommonUtils.timeElapsed(tic)));
		
		return kms;
	}
//	public ArrayList<Kmer> alignByKmerScan (ArrayList<Kmer> kmers_in, int seed_rage, double kmer_aligned_fraction, boolean print_aligned_seqs){
//		long tic = System.currentTimeMillis();
//		if (kmers_in.size()==0)
//			return kmers_in;
//		System.out.println("\nAlign and cluster k-mers ...");
//		boolean bestSeed_is_reset = false;					// clone to modify locally
//		ArrayList<Kmer> kmers = new ArrayList<Kmer>();
//		for (Kmer km:kmers_in)
//			kmers.add(km.clone());
//		kmers.trimToSize();
//		
//		ArrayList<Sequence> seqList = new ArrayList<Sequence>();
//		for (int i=0;i<seqs.length;i++){
//			Sequence s = new Sequence(seqs[i], i);
//			seqList.add(s);
//		}
//		seqList.trimToSize();
//		
//    	clusters.clear();
//    	StringBuilder alignedKmer_sb = new StringBuilder();
//    	
//		for (int clusterID = 0; clusterID<max_cluster; clusterID++){
//			
//			/** Initialization of new cluster and the remaining kmers */
//			KmerCluster cluster = new KmerCluster();
//			cluster.clusterId = clusterID;
//			clusters.add(cluster);
//			
//			/* Initialization, setup sequence list, and update kmers */
//			// build the kmer search tree
//			AhoCorasick oks = new AhoCorasick();
//			for (Kmer km: kmers){
//				oks.add(km.getKmerString().getBytes(), km);
//		    }
//			oks.prepare();
//			
//			// index kmer->seq, seq->kmer
//			HashMap <Kmer, HashSet<Integer>> kmer2seq = new HashMap <Kmer, HashSet<Integer>>();
//			for (Sequence s:seqList){
//				String seq = s.seq;						// get the sequence from original strand
//				HashSet<Kmer> results = queryTree (seq, oks);
//				if (!results.isEmpty()){
//					for (Kmer km: results){		
//						if (!kmer2seq.containsKey(km)){
//							kmer2seq.put(km, new HashSet<Integer>());
//						}
//						kmer2seq.get(km).add(s.id);
//
//						// This is DIFFERENT from alignSequencesByCoOccurence(ArrayList<Kmer>)
//						// farward strand
//						ArrayList<Integer> pos = StringUtils.findAllOccurences(seq, km.getKmerString());
//						for (int p:pos){
//							if (!s.fPos.containsKey(km))
//								s.fPos.put(km, new HashSet<Integer>());
//							s.fPos.get(km).add(p);
//						}
//						pos = StringUtils.findAllOccurences(seq, km.getKmerRC());
//						for (int p:pos){
//							if (!s.fPos.containsKey(km))
//								s.fPos.put(km, new HashSet<Integer>());
//							s.fPos.get(km).add(p+RC);					// match kmer RC
//						}
//						// reverse strand
//						pos = StringUtils.findAllOccurences(s.rc, km.getKmerString());
//						for (int p:pos){
//							if (!s.rPos.containsKey(km))
//								s.rPos.put(km, new HashSet<Integer>());
//							s.rPos.get(km).add(p);
//						}
//						pos = StringUtils.findAllOccurences(s.rc, km.getKmerRC());
//						for (int p:pos){
//							if (!s.rPos.containsKey(km))
//								s.rPos.put(km, new HashSet<Integer>());
//							s.rPos.get(km).add(p+RC);					// match kmer RC
//						}
////						s.setCount();
//					}
//				}
//			}
//			oks = null;
//			
//			// update kmerCount, and hgp()
//			ArrayList<Kmer> unenriched = new ArrayList<Kmer>();
//			for (Kmer km:kmers){
//				if (kmer2seq.containsKey(km)){
//					km.setPosHits(kmer2seq.get(km));
//					km.setHgp(computeHGP(km.getPosHitCount(), km.getNegHitCount()));	
//					if (km.getHgp()>hgp)
//						unenriched.add(km);
//				}
//				else{
//					unenriched.add(km);
//				}
//			}
//			kmers.removeAll(unenriched);
//			for (Sequence s:seqList){
//				s.removeAllKmers(unenriched);
//			}			
//			unenriched = null;
//			kmer2seq = null;
//			
//			if (kmers.isEmpty()){
//				clusters.remove(clusterID);
//				break;
//			}
//			
//			if (verbose>1)
//				System.out.println("------------------------------------------------\n"+CommonUtils.timeElapsed(tic)+
//						": Aligning cluster #"+clusterID+",   n="+kmers.size());
//
//			Collections.sort(kmers, new Comparator<Kmer>(){
//			    public int compare(Kmer o1, Kmer o2) {
//			    	return o1.compareByHGP(o2);
//			    }
//			});			
//	
//			/** get seed kmer */
//			Kmer seed=null;
//			if (bestSeed!=null && clusterID==0){
//				String bestStr = bestSeed.getKmerString();
//				String bestRc = bestSeed.getKmerRC();
//				for (Kmer km:kmers){
//					if (km.getKmerString().equals(bestStr)||
//							km.getKmerString().equals(bestRc)){
//						seed = km;
//						break;
//					}
//				}
//				if (seed==null){				// not found
//					for (Kmer km:kmers){
//						if (km.getHgp()<hgp){
//							seed = km;
//							break;
//						}
//					}
//					if (seed==null){
//						if (kmers.get(0).getHgp()>hgp){		// if the top kmer do not pass HGP test, stop
//							clusters.remove(clusterID);
//							break;
//						}
//						seed = kmers.get(0);
//					}
//				}
//			}	
//			else{				
//				for (Kmer km:kmers){
//					if (km.getHgp()<hgp){
//						seed = km;
//						break;
//					}
//				}
//				if (seed==null){
//					if (kmers.get(0).getHgp()>hgp){
//						clusters.remove(clusterID);
//						break;
//					}
//					seed = kmers.get(0);
//				}
//			}
//			
//			// print out top kmer information
//			if (clusterID==0){
//				System.out.println("\nTop 5 k-mers");
//				System.out.println(Kmer.toShortHeader(k));
//				for (int i=0;i<Math.min(5,kmers.size());i++){
//					System.out.println(kmers.get(i).toShortString());
//				}
//				System.out.println("Seed k-mer:\n"+seed.toShortString()+"\n");
//			}
//			else
//				if (verbose>1)
//					System.out.println(CommonUtils.timeElapsed(tic)+": Seed k-mer: "+seed.toShortString());
//			
//			cluster.seedKmer = seed;
//			
//			ArrayList<Kmer> seedFamily = getSeedKmerFamily(kmers, seed);
//			
//			/** align sequences using kmer positions */
//	    	for (Sequence s : seqList){
//	    		for (Kmer km:seedFamily){
//					int seed_seq = s.getSeq().indexOf(km.getAlignString());
//					if (seed_seq<0){
//						s.RC();
//						seed_seq = s.getSeq().indexOf(km.getAlignString());
//						if (seed_seq<0)
//							continue;
//					}
////					if (km.getKmerString().equals("CCACGCG")||km.getKmerRC().equals("CCACGCG"))
////						km.getK();
//					s.pos = -seed_seq;
//					break;				// seq is aligned, do not try with weaker k-mers
//	    		}
//			}
//	    	seedFamily =null;
//	    	
//	    	/** Iteratively align kmers and align sequences */
//	    	ArrayList<Kmer> alignedKmers=null;
//	    	while(true){
//		    	/** set kmer consensus position */
//	    		HashMap<Kmer, ArrayList<Integer>> kmer2pos = new HashMap<Kmer, ArrayList<Integer>>();
//				for (Sequence s:seqList){
//					if (s.pos != UNALIGNED){		// aligned seqs
//						HashMap<Kmer, HashSet<Integer>> kmerPos = s.getKmerPos();
//						for (Kmer km:kmerPos.keySet()){
////							if (km.getKmerString().equals("CCACGCG")||km.getKmerRC().equals("CCACGCG"))
////								km.getK();;
//							if (!kmer2pos.containsKey(km))
//								kmer2pos.put(km, new ArrayList<Integer>());
//							HashSet<Integer> hits = kmerPos.get(km);
//							if (hits==null)
//								continue;
//							for (int pRef:hits){
//								if (pRef>RC/2){
//									int pos = (pRef-RC)+s.pos;
//									if (pos>=-seed_rage && pos<=seed_rage)
//										kmer2pos.get(km).add(pos+RC);		// use kmerRC 
//								}
//								else{
//									int pos = pRef+s.pos;
//									if (pos>=-seed_rage && pos<=seed_rage)
//										kmer2pos.get(km).add(pos);		
//								}
//							}
//						}
//					}
//				}
//	
//				alignedKmers = new ArrayList<Kmer>();
//				for (Kmer km:kmer2pos.keySet()){
////					if (km.getKmerString().equals("CCACGCG")||km.getKmerRC().equals("CCACGCG"))
////						km.getK();
//					ArrayList<Integer> posKmer = kmer2pos.get(km);
//					if (posKmer==null || posKmer.isEmpty()){
//						km.setAlignString("NOT in 2k region");
//						continue;
//					}
//					if (posKmer.size()<km.getPosHitCount()*kmer_aligned_fraction){			// The kmer hit in seed_range region should be at least 1/2 of total hit
//						km.setAlignString("Too few hit "+posKmer.size());
//						continue;
//					}
//						
//					// find the most frequent kmerPos
//					Pair<int[], int[]> sorted = StatUtil.sortByOccurences(posKmer);
//					int counts[] = sorted.cdr();
//					int posSorted[] = sorted.car();
//					int maxCount = counts[counts.length-1];
//					ArrayList<Integer> maxPos = new ArrayList<Integer>();
//					for (int i=counts.length-1;i>=0;i--){
//						if (counts[i]==maxCount){
//							int p = posSorted[i];
//							maxPos.add(p);
//						}	
//						else		// do not need to count for non-max elements
//							break;
//					}
//					int shift=0;
//					if (maxPos.size()>1){		// if tie with 1+ positions, get the one closest to seed kmer
//						int min = Integer.MAX_VALUE;
//						int minIdx = 0;
//						for (int i=0;i<maxPos.size();i++){
//							int pos = maxPos.get(i);
//							if (pos>RC/2)
//								pos-=RC;
//							int distance = Math.abs(pos);
//							if (distance<min){
//								minIdx = i;
//								min = distance;
//							}
//						}
//						shift=maxPos.get(minIdx);
//					}
//					else{
//						shift=maxPos.get(0);
//					}
////						
//					km.setShift(shift);
//					
//					km.setAlignString(maxCount+"/"+posKmer.size());
//					km.setKmerStartOffset(km.getShift());
//					alignedKmers.add(km);
//				}	
//				kmer2pos = null;
//				alignedKmers.trimToSize();
//				if (verbose>1)
//					System.out.println(String.format("%s: %d kmers aligned.", 
//							CommonUtils.timeElapsed(tic), alignedKmers.size()));
//				if (alignedKmers.size()<=1)
//					break;
//
//				// print out aligned k-mers
//				/*	
//				Collections.sort(alignedKmers, new Comparator<Kmer>(){
//				    public int compare(Kmer o1, Kmer o2) {
//				    	return o1.compareByHGP(o2);
//				    }
//				});	
//				
//				int leftmost_km = Integer.MAX_VALUE;
//				for (Kmer km: alignedKmers){
//					int offset = km.getKmerStartOffset();
//					if (offset>RC/2)
//						offset-=RC;
//					if (offset<leftmost_km)
//						leftmost_km = offset;
//				}	   
//
//				for (Kmer km: alignedKmers){
//					int offset = km.getKmerStartOffset();
//					String kmerStr;
//					if (offset>RC/2){
//						offset-=RC;
//						kmerStr=km.getKmerRC();
//					}
//					else
//						kmerStr = km.getKmerString();
//					System.out.print(offset+"\t"+CommonUtils.padding(-leftmost_km+offset, '.')
//							+kmerStr+"\t"+km.getPosHitCount()+"\t"+km.getNegHitCount()+"\t"+String.format("%.1f", km.getHgp())+"\t"+km.getAlignString()+"\n");
//				}
//				System.out.println("******************************************************");
//					*/
//				
//				// set kmers to search engine, estimate threshold
//				updateEngine(alignedKmers);
//				MotifThreshold threshold = estimateKsmThreshold("", false);
//				if (verbose>1)
//					System.out.println(String.format("%s: KSM score %.2f\tmatch %d+/%d- seqs\thgp=1e%.1f", 
//							CommonUtils.timeElapsed(tic), threshold.score, threshold.posHit, threshold.negHit, threshold.hgp));
//				
//				if (threshold.hgp<cluster.ksmThreshold.hgp){
//					cluster.ksmThreshold = threshold;
//					
//					// scan all sequences, align them using kmer hit results
//					for (Sequence s : seqList){
//						s.reset();			
//						KmerGroup[] kgs = queryS(s.getSeq());
//						//Check if there are matches on the same position from both strand??
//						for (int i=0;i<kgs.length-1;i++){
//							KmerGroup kg=kgs[i];
//							for (int j=i;i<kgs.length;i++){
//								KmerGroup other=kgs[j];
//								if (Math.abs(other.bs-kg.bs)==RC){	// if so, add all kmers to the stronger kmer
//									System.out.println("Kmer match on same position:"+kg.toString()+" "+other.toString()); //TODO: Comment out on release
//									if (kg.hgp<other.hgp){
//										kg.kmers.addAll(other.kmers);
//										kg.setHgp(computeHGP(kg.getGroupHitCount(), kg.getGroupNegHitCount()));
//									}
//									else{
//										other.kmers.addAll(kg.kmers);
//										other.setHgp(computeHGP(other.getGroupHitCount(), other.getGroupNegHitCount()));										
//									}
//								}
//							}
//						}
//						if (kgs.length!=0){
//							Arrays.sort(kgs);
//							KmerGroup kg = kgs[0];
//							if (kg.hgp<=-threshold.score){
//								s.score = -kg.hgp;		// score = -log10 hgp, becomes positive value
//								if (kg.bs>RC/2){		// match on reverse strand
//									s.RC();
//									s.pos = -(kg.bs-RC);
//								}
//								else
//									s.pos = -kg.bs;
//							}
//						}
//					}
//					/*
//			    	int leftmost = Integer.MAX_VALUE;
//			    	for (Sequence s : seqList){
//						if (s.pos==UNALIGNED)
//							continue;
//						if (s.pos < leftmost )
//							leftmost = s.pos;		
//					}
//			
//					StringBuilder sb = new StringBuilder();
//					for (Sequence s : seqList){
//						if (s.pos==UNALIGNED)
//							continue;
//						sb.append(String.format("%d\t%.1f\t%d\t%s\t%s%s\n", s.id, s.score, s.pos, s.isForward?"F":"R", CommonUtils.padding(-leftmost+s.pos, '.'), s.getSeq()));
//					}
//					CommonUtils.writeFile("seqs_aligned.txt", sb.toString());
//					sb=null;
//					*/
//				}
//				else{					// if kmer-scan already converge
//					int aligned_seqs_count=0;
//					for (Sequence s:seqList){
//				    	  if (s.pos!=UNALIGNED)
//				    		  aligned_seqs_count++;
//					}
//					if (aligned_seqs_count<seed.getPosHitCount()*2)	{		// if most of sequences are aligned by seed kmer, stop
//						if (verbose>1)
//							System.out.println(CommonUtils.timeElapsed(tic)+": Sequence: "+aligned_seqs_count+", seed k-mer hit: "+seed.getPosHitCount());
//						break;
//					}
//					
//			    	// build PWM
//					int leftIdx = buildPWM(seqList, cluster, tic);	
//					if (leftIdx > -1){				    	
//						alignSequencesUsingPWM(seqList, cluster);
//					}
//					else
//						break;
////				}
//	    	}  // Iteratively improve k-mer and sequences alignment
//	    	
//			// compare pwm Hgp to previous cluster Hgp, so that the primary cluster will have the best Hgp
//	    	if (select_seed){
//				if (bestSeed==null)						// first cluster, save the seed
//					bestSeed = seed;
//				else if (cluster.ksmThreshold.hgp<clusters.get(0).ksmThreshold.hgp*1.1 && !bestSeed_is_reset){		// this cluster is at least 10% better
//					// reset kmers and posSeqs, to start over with this new bestSeed
//					seqList = new ArrayList<Sequence>();
//					for (int i=0;i<seqs.length;i++){
//						Sequence s = new Sequence(seqs[i], i);
//						seqList.add(s);
//					}
//					seqList.trimToSize();
//					
//					kmers.clear();
//					for (Kmer km:kmers_in)
//						kmers.add(km.clone());	
//					kmers.trimToSize();
//					
//			    	clusters.clear();
//					clusterID = 0;
//					alignedKmer_sb = new StringBuilder();
//					bestSeed = seed;
//					System.out.println("Current motif is more enriched, set primary seed="+seed.getKmerString());
//					bestSeed_is_reset=true;							// marked as "reset", only once, avoid potential infinite loop
//					continue;										// start over with new seed
//				}
//			}
//	    	/** use all aligned sequences to find expected binding sites, set kmer offset */
//	    	// average all the binding positions to decide the expected binding position
//			StringBuilder sb = new StringBuilder();
//			double sum_offset = 0;
//	    	double count = 0;
//	    	int leftmost = Integer.MAX_VALUE;
//	    	int total_aligned_seqs = 0;
//	    	for (Sequence s : seqList){
//				if (s.pos==UNALIGNED)
//					continue;
//				if (s.pos < leftmost )
//					leftmost = s.pos;		
//				total_aligned_seqs++;
//			}
//	    	cluster.total_aligned_seqs = total_aligned_seqs;
//	    	int midPos=k_win/2;
//			for (Sequence s : seqList){
//				if (s.pos==UNALIGNED)
//					continue;
//				if (print_aligned_seqs)
//					sb.append(String.format("%d\t%.1f\t%d\t%s\t%s%s\n", s.id, s.score, s.pos, s.isForward?"F":"R", CommonUtils.padding(-leftmost+s.pos, '.'), s.getSeq()));
//				sum_offset+=midPos+s.pos;
//				count++;
//			}
//			cluster.pos_BS_seed=StatUtil.round(sum_offset/count);		// mean BS position relative to seed k-mer start
//			if (print_aligned_seqs)
//				CommonUtils.writeFile(outName+"_"+clusterID+"_seqs_aligned.txt", sb.toString());
//			sb = null;
//	    	
//	    	// set k-mer offset, print aligned k-mers of this cluster
//	    	if (!alignedKmers.isEmpty())
//	    		alignedKmer_sb.append("Cluster #"+clusterID+", n="+alignedKmers.size()+"\n");
//	    	// sort for output
//			Collections.sort(alignedKmers, new Comparator<Kmer>(){
//			    public int compare(Kmer o1, Kmer o2) {
//			    	return o1.compareByHGP(o2);
//			    }
//			});	
//	    	int leftmost_km = Integer.MAX_VALUE;
//			for (Kmer km: alignedKmers){
//				int shift = km.getShift();
//				if (shift>RC/2){
//					shift-=RC;
//					km.RC();
//				}
//				km.setKmerStartOffset(shift-cluster.pos_BS_seed);
//				km.setClusterId(clusterID);
//				if (km.getKmerStartOffset()<leftmost_km)
//					leftmost_km = km.getKmerStartOffset();
//			}	    	
//			for (Kmer km: alignedKmers){
//				alignedKmer_sb.append(km.getKmerStartOffset()+"\t"+CommonUtils.padding(-leftmost_km+km.getKmerStartOffset(), '.')
//						+km.getKmerString()+"\t"+km.getPosHitCount()+"\t"+km.getNegHitCount()+"\t"+String.format("%.1f", km.getHgp())+"\t"+km.getAlignString()+"\n");
//			}
//	    	
//			/** store aligned kmers in the cluster */
//			ArrayList<Kmer> copy = new ArrayList<Kmer>();
//			ArrayList<Kmer> zero_shift = new ArrayList<Kmer>();
//			for (Kmer km:alignedKmers){
//				Kmer cp = km.clone();
//				copy.add(cp);
//				if (km.getShift()==0)
//					zero_shift.add(km);
//			}
//			cluster.alignedKmers = copy;				// store all the aligned k-mers
//			copy = null;
//			alignedKmers = null;
//			
////			kmers.removeAll(zero_shift);
////			for (Sequence s:seqList){					
////				s.removeAllKmers(zero_shift);
////			}	
////			zero_shift = null;
//			
//			/** mask aligned sequences at seed kmer region */
//			for (Sequence s : seqList){
//				if (s.pos==UNALIGNED)
//					continue;
//				String seq = s.getSeq();
//				if (-s.pos<0 || (-s.pos+k)>seq.length())
//					continue;
//				seq=seq.substring(0, -s.pos)
//				.concat(CommonUtils.padding(k, 'N'))
//				.concat(seq.substring(-s.pos+k, seq.length()));
//				s.setSeq(seq);
//			}
//		} // Loop for each cluster
//		
//		CommonUtils.writeFile(outName+"_Alignement_k"+k+".txt", alignedKmer_sb.toString());
//		alignedKmer_sb = null;
//		
//		return processClusters();
//	}
	
	public ArrayList<Kmer> alignBySimplePWM (ArrayList<Kmer> kmers_in, int seed_range, double kmer_aligned_fraction, 
			boolean print_aligned_seqs, boolean re_train){
		long tic = System.currentTimeMillis();
		if (kmers_in.size()==0)
			return kmers_in;
		System.out.println("\nAlign and cluster k-mers ...");
//		boolean bestSeed_is_reset = false;					
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
		while (!kmers.isEmpty()){
			
			/** Initialization of new cluster and the remaining kmers */
			KmerCluster cluster = new KmerCluster();
			clusters.add(cluster);
			
			indexKmerSequences(kmers, seqList);
			
			if (kmers.isEmpty()){
				clusters.remove(clusterID);
				break;
			}
			
			if (verbose>1)
				System.out.println("------------------------------------------------\n"+CommonUtils.timeElapsed(tic)+
						": Aligning cluster #"+clusterID+",   n="+kmers.size());

			Collections.sort(kmers, new Comparator<Kmer>(){
			    public int compare(Kmer o1, Kmer o2) {
			    	return o1.compareByHGP(o2);
			    }
			});			
	
			/** get seed kmer */
			Kmer seed=null;
			if (bestSeed!=null && clusterID==0){
				String bestStr = bestSeed.getKmerString();
				String bestRc = bestSeed.getKmerRC();
				for (Kmer km:kmers){
					if (km.getKmerString().equals(bestStr)||
							km.getKmerString().equals(bestRc)){
						seed = km;
						break;
					}
				}
				if (seed==null){				// not found
					for (Kmer km:kmers){
						if (km.getHgp()<hgp){
							seed = km;
							break;
						}
					}
					if (seed==null){
						if (kmers.get(0).getHgp()>hgp){		// if the top kmer do not pass HGP test, stop
							clusters.remove(clusterID);
							break;
						}
						seed = kmers.get(0);
					}
				}
			}	
			else{				
				for (Kmer km:kmers){
					if (km.getHgp()<hgp){
						seed = km;
						break;
					}
				}
				if (seed==null){
					if (kmers.get(0).getHgp()>hgp){
						clusters.remove(clusterID);
						break;
					}
					seed = kmers.get(0);
				}
			}
			
			// print out top kmer information
			if (clusterID==0){
				System.out.println("\nTop 5 k-mers");
				System.out.println(Kmer.toShortHeader(k));
				for (int i=0;i<Math.min(5,kmers.size());i++){
					System.out.println(kmers.get(i).toShortString());
				}
				System.out.println("Seed k-mer:\n"+seed.toShortString()+"\n");
			}
			else
				if (verbose>1)
					System.out.println(CommonUtils.timeElapsed(tic)+": Seed k-mer: "+seed.toShortString());
			
			cluster.seedKmer = seed;
			
			ArrayList<Kmer> seedFamily = getSeedKmerFamily(kmers, seed);
			
			/** align sequences using seedFamily kmer positions */
	    	for (Sequence s : seqList){
	    		for (Kmer km:seedFamily){
					int seed_seq = s.getSeq().indexOf(km.getAlignString());
					if (seed_seq<0){
						s.RC();
						seed_seq = s.getSeq().indexOf(km.getAlignString());
						if (seed_seq<0)
							continue;
					}
//					if (km.getKmerString().equals("CCACGCG")||km.getKmerRC().equals("CCACGCG"))
//						km.getK();
					s.pos = -seed_seq;
					break;				// seq is aligned, do not try with weaker k-mers
	    		}
			}
	    	
	    	cluster.alignedKmers = getAlignedKmers (seqList, seed_range, kmer_aligned_fraction, new ArrayList<Kmer>());
	    	
	    	alignSequencesUsingKSM(seqList, cluster);
	    	
	    	if (re_train)	{
		    	/** build PWM, do not iterate */
				int aligned_seqs_count=0;
				for (Sequence s:seqList){
			    	  if (s.pos!=UNALIGNED)
			    		  aligned_seqs_count++;
				}
				if (aligned_seqs_count<seed.getPosHitCount()*1.5)	{		// if most of sequences are aligned by seed kmer, stop
					if (verbose>1)
						System.out.println(CommonUtils.timeElapsed(tic)+": Sequence: "+aligned_seqs_count+", seed k-mer hit: "+seed.getPosHitCount());
					kmers.remove(seed);
					continue;
				}
				
		    	// build PWM
				int leftIdx = buildPWM(seqList, cluster, tic);	
				if (leftIdx <= -1){		
					kmers.remove(seed);
					continue;
				}
	    	}
	    	else{
		    	/** Iteratively build PWM and align sequences */
		    	while(true){			
					int aligned_seqs_count=0;
					for (Sequence s:seqList){
				    	  if (s.pos!=UNALIGNED)
				    		  aligned_seqs_count++;
					}
					if (aligned_seqs_count<seed.getPosHitCount()*1.5)	{		// if most of sequences are aligned by seed kmer, stop
						if (verbose>1)
							System.out.println(CommonUtils.timeElapsed(tic)+": Sequence: "+aligned_seqs_count+", seed k-mer hit: "+seed.getPosHitCount());
						break;
					}
					
			    	// build PWM
					int leftIdx = buildPWM(seqList, cluster, tic);	
					if (leftIdx > -1){				    	
						alignSequencesUsingPWM(seqList, cluster);
					}
					else
						break;
		    	}  // Iteratively improve PWM and sequences alignment
		    	
		    	/** use all aligned sequences to find expected binding sites, set kmer offset */
		    	// average all the binding positions to decide the expected binding position
				StringBuilder sb = new StringBuilder();
				double sum_offset = 0;
		    	double count = 0;
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
		    	int midPos=k_win/2;
				for (Sequence s : seqList){
					if (s.pos==UNALIGNED)
						continue;
					if (print_aligned_seqs)
						sb.append(String.format("%d\t%.1f\t%d\t%s\t%s%s\n", s.id, s.score, s.pos, s.isForward?"F":"R", CommonUtils.padding(-leftmost+s.pos, '.'), s.getSeq()));
					sum_offset+=midPos+s.pos;
					count++;
				}
				cluster.pos_BS_seed=StatUtil.round(sum_offset/count);		// mean BS position relative to seed k-mer start
				if (print_aligned_seqs)
					CommonUtils.writeFile(outName+"_"+clusterID+"_seqs_aligned.txt", sb.toString());
				sb = null;
				if (verbose>1 && cluster.wm!=null){
					int pos_BS_midPWM = cluster.pos_BS_seed - cluster.pos_pwm_seed - cluster.wm.length()/2;
					System.out.println("pos_BS_midPWM = "+pos_BS_midPWM);
				}
	    	}
	    	
			// get Aligned Kmers either from PWM or seed family alignement
	    	ArrayList<Kmer> alignedKmers = getAlignedKmers (seqList, seed_range, kmer_aligned_fraction, new ArrayList<Kmer>());

	    	// set k-mer offset
			for (Kmer km: alignedKmers){
				int shift = km.getShift();
				if (shift>RC/2){
					shift-=RC;
					km.RC();
					km.setShift(shift);
				}
				km.setKmerStartOffset(shift-cluster.pos_BS_seed);
			}	    	
			
			if (!re_train)	{	
				/** store aligned kmers in the cluster */
				ArrayList<Kmer> copy = new ArrayList<Kmer>();
				for (Kmer km:alignedKmers){
					Kmer cp = km.clone();
					copy.add(cp);
				}
				cluster.alignedKmers = copy;				// store all the aligned k-mers
				copy = null;
			}
			
			/** mask aligned sequences */
			float k_mask_f=1;
			if (cluster.wm!=null){				
		        if (cluster.pwmGoodQuality){			// if PWM quality is not too bad, mask
			        WeightMatrixScorer scorer = new WeightMatrixScorer(cluster.wm);
			        int left = Math.round(cluster.wm.length()*(1-k_mask_f)/2);
			        int right = Math.round(cluster.wm.length()*(1-(1-k_mask_f)/2));
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
							// replace with N
							seq=seq.substring(0, start+left)
									.concat(CommonUtils.padding(right-left, 'N'))
									.concat(seq.substring(start+right, seq.length()));
						}
						if (found){
							seqsNegList.set(i, seq);
						}
					}
		        }
			}
			for (Kmer km: alignedKmers)
				if (km.getShift()<=k/2)
					kmers.remove(km);
			kmers.remove(seed);
//				for (Sequence s : seqList){
//					if (s.pos==UNALIGNED)
//						continue;
//					String seq = s.getSeq();
//					if (-s.pos<0 || (-s.pos+k)>seq.length())
//						continue;
//					seq=seq.substring(0, -s.pos)
//					.concat(CommonUtils.padding(k, 'N'))
//					.concat(seq.substring(-s.pos+k, seq.length()));
//					s.setSeq(seq);
//				}
			clusterID++;
		} // Loop for each cluster
		
		if (re_train)	{
			ArrayList<KmerCluster> noPWM = new ArrayList<KmerCluster>();
			for (KmerCluster c:clusters)
				if (c.wm==null)
					noPWM.add(c);
			clusters.removeAll(noPWM);
			
			// reload sequences and kmers
			for (int i=0;i<seqs.length;i++)
				seqList.set(i, new Sequence(seqs[i], i));
			
			kmers.clear();
			for (Kmer km:kmers_in)
				kmers.add(km.clone());
			kmers.trimToSize();
			
			indexKmerSequences(kmers, seqList);
				
			/** Resolve multiple-PWM-match conflict */
			if (clusters.size()>1){
				for (int iter=0; iter<100; iter++){
					
					// record cluster PWM hgp for comparison
					double[] hgps = new double[clusters.size()];
					if (verbose>1)
						System.out.println("------------------------------------------------\n"+CommonUtils.timeElapsed(tic)+
								": Iteration #"+iter);
	
					Collections.sort(clusters);
					for (int i=0; i<clusters.size(); i++){
						hgps[i] = clusters.get(i).pwmThresholdHGP;
						clusters.get(i).clusterId = i;
						clusters.get(i).seq2hits = findAllPWMHits(seqList, clusters.get(i));
					}				
					System.out.println(CommonUtils.arrayToString(hgps));
					int conflicts = resolveConflictHits(seqList);
					if (verbose>1)
						System.out.println(conflicts+" sequence hits are matched by multiple PWMs");
					if (conflicts==0)
						break;
						
					buildPWMfromHits(seqList);
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
			
			for (int i=0; i<clusters.size(); i++){
				/** use all aligned sequences to find expected binding sites, set kmer offset */
		    	// average all the binding positions to decide the expected binding position
				KmerCluster cluster = clusters.get(i);
				alignSequencesUsingPWM(seqList, cluster);
				
				StringBuilder sb = new StringBuilder();
				double sum_offset = 0;
		    	double count = 0;
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
		    	int midPos=k_win/2;
				for (Sequence s : seqList){
					if (s.pos==UNALIGNED)
						continue;
					if (print_aligned_seqs)
						sb.append(String.format("%d\t%.1f\t%d\t%s\t%s%s\n", s.id, s.score, s.pos, s.isForward?"F":"R", CommonUtils.padding(-leftmost+s.pos, '.'), s.getSeq()));
					sum_offset+=midPos+s.pos;
					count++;
				}
				cluster.pos_BS_seed=StatUtil.round(sum_offset/count);		// mean BS position relative to seed k-mer start
				if (print_aligned_seqs)
					CommonUtils.writeFile(outName+"_"+i+"_seqs_aligned.txt", sb.toString());
				sb = null;
		    	
				// get Aligned Kmers either from PWM or seed family alignement
				ArrayList<Kmer> alignedKmers = getAlignedKmers (seqList, seed_range, kmer_aligned_fraction, new ArrayList<Kmer>());
			    	
				/** store aligned kmers in the cluster */
				ArrayList<Kmer> copy = new ArrayList<Kmer>();
				for (Kmer km:alignedKmers){
					Kmer cp = km.clone();
					int shift = cp.getShift();
					if (shift>RC/2){
						shift-=RC;
						cp.RC();
					}
					cp.setKmerStartOffset(shift-cluster.pos_BS_seed);
					copy.add(cp);
				}
				cluster.alignedKmers = copy;				// store all the aligned k-mers
				copy = null;			
			}
		}
		return processClusters();
	}
	
    /**
	 * Compute the distance distributions between primary and secondary motifs
	 */
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public void printMotifDistanceDistribution (String name){		
		System.out.println("Compute motif distance distribution ...");

		ArrayList[][] hits = new ArrayList[seqs.length][clusters.size()];
		for (int j=0;j<clusters.size();j++){
			KmerCluster c = clusters.get(j);
			WeightMatrixScorer scorer = new WeightMatrixScorer(c.wm);
			for (int i=0;i<seqs.length;i++){
				hits[i][j]=CommonUtils.getAllPWMHit(seqs[i], c.wm.length(), scorer, c.pwmThreshold);
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
				for (int i=0;i<same.length;i++){
					sb.append(String.format("%d\t%d\t%d\n", i-range, same[i], diff[i]));
				}
				CommonUtils.writeFile(name+"_Spatial_dist_"+clusters.get(m).clusterId+"_"+clusters.get(j).clusterId+".txt", sb.toString());
			}
		}
	}
	
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
			HashSet<Kmer> results = queryTree (seq, oks);
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
//					s.setCount();
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
				if (km.getHgp()>hgp)
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
	private ArrayList<Kmer> processClusters(){
		// output cluster information, PFM, and PWM
		File f = new File(outName);
		String name = f.getName();
		f=null;
		StringBuilder pfm_sb = new StringBuilder();	
		
		ArrayList<KmerCluster> badClusters = new ArrayList<KmerCluster>();
		for (KmerCluster c:clusters){
    		if (c.wm==null || (!c.pwmGoodQuality) || c.total_aligned_seqs<seqs.length*motif_hit_factor_report)
    			badClusters.add(c);
		}
		clusters.removeAll(badClusters);
		Collections.sort(clusters);
		for (int i=0;i<clusters.size();i++){
			clusters.get(i).clusterId = i;
			for (Kmer km: clusters.get(i).alignedKmers)
				km.setClusterId(i);
		}
		
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
		for (KmerCluster c:clusters){
     		WeightMatrix wm = c.wm;
    		System.out.println(String.format("------------------------------------------------\n%s k-mer cluster #%d, aligned %d k-mers, %d sequences.", outName, c.clusterId, c.alignedKmers.size(), c.total_aligned_seqs));
//    		System.out.println(String.format("KSM threshold: %.2f, \thit=%d+/%d-, hgp=1e%.1f", c.ksmThreshold.score, c.ksmThreshold.posHit, c.ksmThreshold.negHit, c.ksmThreshold.hgp));
			int pos = c.pos_BS_seed-c.pos_pwm_seed;
    		if (pos>0)
    			System.out.println(CommonUtils.padding(pos, ' ')+"|\n"+ WeightMatrix.printMatrixLetters(wm));
    		else
    			System.out.println(WeightMatrix.printMatrixLetters(wm));
    		System.out.println(String.format("PWM threshold: %.2f/%.2f, \thit=%d+/%d-, hgp=1e%.1f", c.pwmThreshold, c.wm.getMaxScore(), c.pwmPosHitCount, c.pwmNegHitCount, c.pwmThresholdHGP));
			pfm_sb.append(c.pfmString);
			
			// paint motif logo
			c.wm.setNameVerType(name, "#"+c.clusterId, "");
			CommonUtils.printMotifLogo(c.wm, new File(outName+"_"+c.clusterId+"_motif.png"), 75);
		}
		CommonUtils.writeFile(outName+"_PFM_k"+k+".txt", pfm_sb.toString());
		
//		// get all the kmers, if duplicates, the kmer from earlier cluster has higher priority
//		HashMap<String, Kmer> allKmers = new HashMap<String, Kmer>();
//		for (KmerCluster c:clusters){
//			for (Kmer km: c.alignedKmers){
//				if ( !(allKmers.containsKey(km.getKmerString())||allKmers.containsKey(km.getKmerRC())))
//					allKmers.put(km.getKmerString(), km);
//			}
//		}
//		ArrayList<Kmer> allAlignedKmers = new ArrayList<Kmer>();
//		allAlignedKmers.addAll(allKmers.values());
//		allKmers = null;

		ArrayList<Kmer> allAlignedKmers = new ArrayList<Kmer>();		// each kmer in diff cluster has been clone, would not overwrite
		for (KmerCluster c:clusters)
			allAlignedKmers.addAll(c.alignedKmers);
		
		// output HTML report
		Collections.sort(allAlignedKmers, new Comparator<Kmer>(){
		    public int compare(Kmer o1, Kmer o2) {
		    	return o1.compareByHGP(o2);
		    }
		});	
		StringBuffer html = new StringBuffer("<style type='text/css'>/* <![CDATA[ */ table, td{border-color: #600;border-style: solid;} table{border-width: 0 0 1px 1px; border-spacing: 0;border-collapse: collapse;} td{margin: 0;padding: 4px;border-width: 1px 1px 0 0;} /* ]]> */</style>");
		html.append("<table><th bgcolor='#A8CFFF' colspan=2><font size='5'>");
		html.append(name).append("</font></th>");
		html.append("<tr><td>");
		html.append("<a href='"+name+"_GPS_significant.txt'>Significant Events</a>&nbsp;&nbsp;: "+seqs.length);
//		html.append("<br><a href='"+name+"_GPS_insignificant.txt'>Insignificant Events</a>: "+insignificantFeatures.size());
//		html.append("<br><a href='"+name+"_GPS_filtered.txt'>Filtered Events</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;: "+filteredFeatures.size());
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
		html.append("</table><br><a href='"+name+"_kmer_k"+k+".txt'>The complete K-mer list.</a>");
		html.append("<br><a href='"+name+"_Alignement_k"+k+".txt'>The k-mer alignment file.</a>");
		html.append("<br><a href='"+name+"_PFM_k"+k+".txt'>The PFM of motifs</a>");
		html.append("</td><td><br>");
		for (KmerCluster c:clusters){
    		WeightMatrix wm = c.wm;
    		if (wm==null || (!c.pwmGoodQuality)|| c.total_aligned_seqs<seqs.length*motif_hit_factor_report)
    			continue;
    		html.append("<img src='"+name+"_"+c.clusterId+"_motif.png"+"'><br>");
    		html.append(String.format("PWM score: %.2f/%.2f, hit=%d+/%d-, hgp=1e%.1f<br>", 
    				c.pwmThreshold, c.wm.getMaxScore(), c.pwmPosHitCount, c.pwmNegHitCount, c.pwmThresholdHGP));
//    		html.append(String.format("KSM score: %.2f, \thit=%d+/%d-, hgp=1e%.1f<br><br>", 
//    				c.ksmThreshold.score, c.ksmThreshold.posHit, c.ksmThreshold.negHit, c.ksmThreshold.hgp));
		}
		html.append("</td></tr></table>");
		CommonUtils.writeFile(outName+"_result.htm", html.toString());
		
//		for (int i=0;i<Math.min(clusters.size(), 20);i++){
//			ArrayList<Kmer> clusterKmers = clusters.get(i).alignedKmers;
//			MotifThreshold t = this.estimateClusterKgsThreshold(clusterKmers);
//			if (t!=null)
//				System.out.println(String.format("%d\t%.2f\t%d\t%d\t%.1f", i, t.score, t.posHit, t.negHit, t.hgp ));
//		}
		System.out.print("\nAlignment done, "+CommonUtils.timeElapsed(tic)+"\n");
		
		return allAlignedKmers;
	}
	
	public ArrayList<Kmer> clusterKmers (ArrayList<Kmer> kmers_in, int seed_range, double kmer_aligned_fraction, boolean print_aligned_seqs){

		tic = System.currentTimeMillis();
		if (kmers_in.size()==0)
			return kmers_in;
		System.out.println("\nCluster k-mers ...");
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
		
		/* Initialization, setup sequence list, and update kmers */
		// build the kmer search tree
		AhoCorasick oks = new AhoCorasick();
		for (Kmer km: kmers){
			oks.add(km.getKmerString().getBytes(), km);
	    }
		oks.prepare();
		
		// index seq->kmer
//		HashMap <Kmer, HashSet<Integer>> kmer2seq = new HashMap <Kmer, HashSet<Integer>>();
		for (Sequence s:seqList){
			String seq = s.seq;						// get the sequence from original strand
			HashSet<Kmer> results = queryTree (seq, oks);
			if (!results.isEmpty()){
				for (Kmer km: results){		
//					if (!kmer2seq.containsKey(km)){
//						kmer2seq.put(km, new HashSet<Integer>());
//					}
//					kmer2seq.get(km).add(s.id);

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
//					s.setCount();
				}
			}
		}
		oks = null;
		
		// Make the seedCandidates list
		Collections.sort(kmers, new Comparator<Kmer>(){
		    public int compare(Kmer o1, Kmer o2) {
		    	return o1.compareByHGP(o2);
		    }
		});		
		ArrayList<Kmer> seedCandidates = new ArrayList<Kmer>();
		for (int i=0;i<kmers.size()*seed_search_fraction;i++){
			seedCandidates.add(kmers.get(i));
		}
		
    	clusters.clear();
    	
    	/** Construct inital clusters and merge similar clusters */
    	int clusterID = 0;
    	ArrayList<Kmer> alignedCandidates = new ArrayList<Kmer>(); 
		while ( !seedCandidates.isEmpty()){
			
			/** Initialization of new cluster and the remaining kmers */
			KmerCluster cluster = new KmerCluster();
			cluster.clusterId = clusterID;
			clusters.add(cluster);

			Collections.sort(kmers, new Comparator<Kmer>(){
			    public int compare(Kmer o1, Kmer o2) {
			    	return o1.compareByHGP(o2);
			    }
			});			
	
			/** get seed kmer */
			Kmer seed = seedCandidates.get(0);
			cluster.seedKmer = seed;
			
			alignSequencesUsingSeedFamily(seqList, seedCandidates, seed);
	    	
			cluster.alignedKmers = getAlignedKmers (seqList, seed_range, kmer_aligned_fraction, alignedCandidates);
			
			alignedCandidates.addAll(cluster.alignedKmers);
			
			seedCandidates.removeAll(cluster.alignedKmers);

			clusterID++;
		}
		clusters.trimToSize();
		
		// Consolidate the clusters
//		mergeClusters(kmer_set_overlap_ratio);
		
		/** Use seed-derived kmer sets to align sequences, and grow PWM by progressively improve its enrichment score */
		for (int i=0; i<clusters.size(); i++){
			// re-align the k-mers in each cluster
			KmerCluster cluster = clusters.get(i);
			cluster.clusterId = i;

			if (verbose>1)
				System.out.println("------------------------------------------------\n"+CommonUtils.timeElapsed(tic)+
						": Aligning cluster #"+cluster.clusterId+",   seed="+cluster.seedKmer.toShortString());

//			alignSequencesUsingSeedFamily(seqList, kmers, seed);
//			ArrayList<Kmer> alignedKmers = getAlignedKmers (seqList, 10000, kmer_aligned_fraction);				// no range limit
//			
//			SetTools<Kmer> setTools = new SetTools<Kmer>();
//			Set<Kmer> first = new HashSet<Kmer>();
//			first.addAll(cluster.alignedKmers);
//			Set<Kmer> second = new HashSet<Kmer>();
//			second.addAll(alignedKmers);
//			Set<Kmer> unaligned = setTools.subtract(first, second);			// ignore those can not be aligned			
//			cluster.alignedKmers.removeAll(unaligned);

			// aligned sequences with aliged kmer set (KSM)
			alignSequencesUsingKSM(seqList, cluster);

			// grow PWM
			while(true){
		    	// build PWM
				int leftIdx = buildPWM(seqList, cluster, tic);	
				if (leftIdx <= -1)
					break;
				
				// align sequences with PWM
				alignSequencesUsingPWM(seqList, cluster);
			}
			
			cluster.alignedKmers = getAlignedKmers (seqList, seed_range, kmer_aligned_fraction, new ArrayList<Kmer>());
		}
		
		ArrayList<KmerCluster> noPWM = new ArrayList<KmerCluster>();
		for (KmerCluster c:clusters)
			if (c.wm==null)
				noPWM.add(c);
		clusters.removeAll(noPWM);

		// Consolidate the clusters
		mergeClusters(kmer_set_overlap_ratio);


		/** Resolve multiple-PWM-match conflict */
		if (clusters.size()>1){
			for (int iter=0; iter<10; iter++){
				// record cluster PWM hgp for comparison
				double[] hgps = new double[clusters.size()];
				if (verbose>1)
					System.out.println("------------------------------------------------\n"+CommonUtils.timeElapsed(tic)+
							": Iteration #"+iter);

				Collections.sort(clusters);
				for (int i=0; i<clusters.size(); i++){
					hgps[i] = clusters.get(i).pwmThresholdHGP;
					clusters.get(i).clusterId = i;
					clusters.get(i).seq2hits = findAllPWMHits(seqList, clusters.get(i));
				}				
				System.out.println(CommonUtils.arrayToString(hgps));
				int conflicts = removeAllConflictHits(seqList);
				if (verbose>1)
					System.out.println(conflicts+" sequence hits are matched by multiple PWMs");
				if (conflicts==0)
					break;
					
				buildPWMfromHits(seqList);
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
		
		if (verbose>1)
			System.out.println("------------------------------------------------\n"+CommonUtils.timeElapsed(tic)+
					": Collect aligned k-mers based on PWM hit aligned sequences.");
		
		Collections.sort(clusters);
		
		StringBuilder alignedKmer_sb = new StringBuilder();
		for (int i=0; i<clusters.size(); i++){
			KmerCluster cluster = clusters.get(i);
			cluster.clusterId = i;
			
			// align sequences with PWM
			alignSequencesUsingPWM(seqList, clusters.get(i));
			
			// get aligned k-mers from PWM aligned sequences
			ArrayList<Kmer> alignedKmers = getAlignedKmers (seqList, seed_range, kmer_aligned_fraction, new ArrayList<Kmer>());
			updateEngine(alignedKmers);
			cluster.ksmThreshold = estimateKsmThreshold("", false);		// only estimate KSM threshold, do not align sequences
			
			/** use all aligned sequences to find expected binding sites, set kmer offset */
	    	// average all the binding positions to decide the expected binding position
			StringBuilder sb = new StringBuilder();
			double sum_offset = 0;
	    	double count = 0;
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
	    	int midPos=k_win/2;
			for (Sequence s : seqList){
				if (s.pos==UNALIGNED)
					continue;
				if (print_aligned_seqs)
					sb.append(String.format("%d\t%.1f\t%d\t%s\t%s%s\n", s.id, s.score, s.pos, s.isForward?"F":"R", CommonUtils.padding(-leftmost+s.pos, '.'), s.getSeq()));
				sum_offset+=midPos+s.pos;
				count++;
			}
			cluster.pos_BS_seed=StatUtil.round(sum_offset/count);		// mean BS position relative to seed k-mer start
			if (print_aligned_seqs)
				CommonUtils.writeFile(outName+"_"+i+"_seqs_aligned.txt", sb.toString());
			sb = null;
	    	
	    	// set k-mer offset, print aligned k-mers of this cluster
	    	if (!alignedKmers.isEmpty())
	    		alignedKmer_sb.append("Cluster #"+i+", n="+alignedKmers.size()+"\n");
	    	// sort for output
			Collections.sort(alignedKmers, new Comparator<Kmer>(){
			    public int compare(Kmer o1, Kmer o2) {
			    	return o1.compareByHGP(o2);
			    }
			});	
	    	int leftmost_km = Integer.MAX_VALUE;
			for (Kmer km: alignedKmers){
				km.setClusterId(i);
				int shift = km.getShift();
				if (shift>RC/2){
					shift-=RC;
					km.RC();
				}
				km.setKmerStartOffset(shift-cluster.pos_BS_seed);
				km.setClusterId(i);
				if (km.getKmerStartOffset()<leftmost_km)
					leftmost_km = km.getKmerStartOffset();
			}	    	
			for (Kmer km: alignedKmers){
				alignedKmer_sb.append(km.getKmerStartOffset()+"\t"+CommonUtils.padding(-leftmost_km+km.getKmerStartOffset(), '.')
						+km.getKmerString()+"\t"+km.getPosHitCount()+"\t"+km.getNegHitCount()+"\t"+String.format("%.1f", km.getHgp())+"\t"+km.getAlignString()+"\n");
			}
			
			ArrayList<Kmer> copy = new ArrayList<Kmer>();
			for (Kmer km:alignedKmers){
				Kmer cp = km.clone();
				copy.add(cp);
			}
			clusters.get(i).alignedKmers = copy;
		}
		
		CommonUtils.writeFile(outName+"_Alignement_k"+k+".txt", alignedKmer_sb.toString());
		alignedKmer_sb = null;
		
		return processClusters();
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
	
	public int resolveConflictHits(ArrayList<Sequence> seqList){
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

	public void buildPWMfromHits(ArrayList<Sequence> seqList){
		int extend = 2;
		ArrayList<KmerCluster> badClusters = new ArrayList<KmerCluster>();
		for (int j=0; j<clusters.size(); j++){
			KmerCluster cluster = clusters.get(j);
			Iterator<PWMHit> hits = cluster.seq2hits.values().iterator();
			ArrayList<String> alignedSeqs = new ArrayList<String>();
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
			}
			// make PWM from aligned sequence segament
			int result = buildPWMfromAlignedSequences(alignedSeqs, cluster, false);
			if (result==-1){
				badClusters.add(cluster);
				if (verbose>1)
					System.out.println("Cluster "+j+" failed to form a PWM, remvoe it.");
			}
		}
		clusters.removeAll(badClusters);
	}
	
	private int buildPWMfromAlignedSequences(ArrayList<String> alignedSeqs, KmerCluster cluster, boolean onlyBetter){
			if (alignedSeqs.isEmpty())
				return -1;
			double[][] pfm = new double[alignedSeqs.get(0).length()][MAXLETTERVAL];
	    	if (verbose>1)
	    		System.out.println(String.format("%s: %d seqs to build PWM.", CommonUtils.timeElapsed(tic), alignedSeqs.size()));
			// count base frequencies
			for (int p=0;p<pfm.length;p++){
				for (char base:LETTERS)			// 0 count can cause log(0), set pseudo-count 0.375 to every pos, every base
					pfm[p][base]=0.375; 		//http://www.ncbi.nlm.nih.gov.libproxy.mit.edu/pmc/articles/PMC2490743/
			} 
	    	for (String s:alignedSeqs){
				for (int p=0;p<pfm.length;p++){
	    			char base = s.charAt(p);
	    			pfm[p][base] ++;
	    		}
	    	}
	    	
			double[][] pwm = pfm.clone();
			for (int i=0;i<pfm.length;i++)
				pwm[i]=pfm[i].clone();
			
	    	// make the PWM
	    	// normalize, compare to background, and log2
	    	double[] ic = new double[pwm.length];						// information content
	    	for (int p=0;p<pwm.length;p++){						// for each position
	    		double countN = pwm[p]['N'];
	    		if (countN!=0){
	    			for (int b=0;b<LETTERS.length;b++){
	        			char base = LETTERS[b];						// add the fraction of 'N' according to bg dist
		    			pwm[p][base] += countN*bg[b];
		    		}   
	    		}
	    		int sum=0;
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
			int[] left, right;
			if (rightIdx-leftIdx+1>(k-2)){
				left=new int[rightIdx-leftIdx+1-(k-2)];
				right=new int[rightIdx-leftIdx+1-(k-2)];
				for (int i=0;i<left.length;i++){
					int bestLeft = -1;
					double bestSumIC = 0;
					for(int p=leftIdx;p<=rightIdx-(k-2)-i;p++){
						int end = (k-2)+i+p;
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
					right[i]=bestLeft+(k-2)+i;
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
		    	MotifThreshold estimate = estimatePwmThreshold(wm, "", false, wm.getMaxScore()/2);
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
	    	cluster.pwmGoodQuality = (bestEstimate.posHit>seqs.length*motif_hit_factor);
	    	cluster.pwmThreshold = bestEstimate.score;
	    	cluster.pwmThresholdHGP = bestHGP;
	    	cluster.pwmPosHitCount = bestEstimate.posHit;
	    	cluster.pwmNegHitCount = bestEstimate.negHit;
	    	// record pfm
	    	float[][] pfm_trim = new float[bestRight-bestLeft+1][MAXLETTERVAL];   
	    	for(int p=bestLeft;p<=bestRight;p++){
	    		for (char base:LETTERS){
	    			pfm_trim[p-bestLeft][base]=(float) pfm[p][base];
	    		}
	    	}
	    	cluster.pfmString = makeTRANSFAC (pfm_trim, String.format("DE %s_%d_c%d\n", outName, cluster.clusterId, cluster.pwmPosHitCount));
	    	cluster.pos_pwm_seed = bestLeft-(alignedSeqs.get(0).length()/2-k/2);		// pwm_seed = pwm_seqNew-seed_seqNew
	    	
	    	return bestLeft;
		}
	private void alignSequencesUsingKSM(ArrayList<Sequence> seqList, KmerCluster cluster){
		ArrayList<Kmer> alignedKmers = cluster.alignedKmers;
		updateEngine(alignedKmers);
		MotifThreshold threshold = estimateKsmThreshold("", false);
		if (verbose>1)
			System.out.println(String.format("%s: KSM score %.2f\tmatch %d+/%d- seqs\thgp=1e%.1f", 
					CommonUtils.timeElapsed(tic), threshold.score, threshold.posHit, threshold.negHit, threshold.hgp));
		
		if (threshold.hgp<cluster.ksmThreshold.hgp){
			cluster.ksmThreshold = threshold;
			
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
		}
	}
	/** align all sequences using a PWM<br>
	 * sequences have a PWM hit is aligned according hit position, if no PWM hit, set it to unaligned
	 */
	private void alignSequencesUsingPWM(ArrayList<Sequence> seqList, KmerCluster cluster){
			    	
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
					maxScoringShift = k_win-maxScoringShift-wm.length();
					s.RC();
					// i.e.  (seq.length()-1)-maxScoringShift-(wm.length()-1);
				}
				s.pos = cluster.pos_pwm_seed-maxScoringShift;
				count_pwm_aligned ++;
	          }
	          else
	        	  s.pos = UNALIGNED;
	        }	// each sequence

	    	if (verbose>1)
	    		System.out.println(CommonUtils.timeElapsed(tic)+": PWM "+WeightMatrix.getMaxLetters(wm)+" align "+count_pwm_aligned+" sequences.");
	}
	
	private HashMap<Integer, PWMHit> findAllPWMHits(ArrayList<Sequence> seqList, KmerCluster cluster){
    	
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
          if (maxSeqScore >= cluster.pwmThreshold){
        	  PWMHit hit = new PWMHit();
        	  hit.clusterId = cluster.clusterId;
        	  hit.wm = wm;
        	  hit.seqId = i;
        	  hit.start = maxScoringShift;					// start is cordinate on + strand
        	  hit.end = maxScoringShift+wm.length()-1;		// end inclusive
			  hit.isForward = maxScoringStrand =='+';
			  hit.score = maxSeqScore/wm.getMaxScore();			// TODO: score normalized by PWM width
			  seq2hits.put(i, hit);
          }
        }	// each sequence
		return seq2hits;
	}
	
	private void alignSequencesUsingSeedFamily(ArrayList<Sequence> seqList, ArrayList<Kmer> kmers, Kmer seed){
		ArrayList<Kmer> seedFamily = getSeedKmerFamily(kmers, seed);	// from all kmers
		
		/** align sequences using kmer positions */
    	for (Sequence s : seqList){						// use all sequences
    		s.reset();
    		for (Kmer km:seedFamily){
				int seed_seq = s.getSeq().indexOf(km.getAlignString());
				if (seed_seq<0){
					s.RC();
					seed_seq = s.getSeq().indexOf(km.getAlignString());
					if (seed_seq<0)
						continue;
				}
//				if (km.getKmerString().equals("CCACGCG")||km.getKmerRC().equals("CCACGCG"))
//					km.getK();
				s.pos = -seed_seq;
				break;				// seq is aligned, do not try with weaker k-mers
    		}
		}
	}

	// 
	private ArrayList<Kmer> getSeedKmerFamily(ArrayList<Kmer> kmers, Kmer seed) {
		ArrayList<Kmer> family = new ArrayList<Kmer>();
		family.add(seed);
		String seedKmerStr = seed.getKmerString();
		String seedKmerRC = seed.getKmerRC();
		seed.setShift(0);
		seed.setAlignString(seedKmerStr);
		for (Kmer kmer: kmers){
	    	if (CommonUtils.mismatch(seedKmerStr, kmer.getKmerString())==1){
	    		kmer.setShift(0);
	    		family.add(kmer);
	    		kmer.setAlignString(kmer.getKmerString());
	    	}
	    	else if (CommonUtils.mismatch(seedKmerRC, kmer.getKmerString())<=1){	// do not RC kmer, set the string for alignment
	    		kmer.setShift(0);
	    		family.add(kmer);
	    		kmer.setAlignString(kmer.getKmerRC());
	    	}
		}
		return family;
	}
	/**
	 * Get aligned kmers within seed_range using the aligned sequences
	 * @param seqList
	 * @param seed_range
	 * @param kmer_aligned_fraction
	 * @return
	 */
	private ArrayList<Kmer> getAlignedKmers (ArrayList<Sequence> seqList, int seed_range, double kmer_aligned_fraction, ArrayList<Kmer> excludes){

    	/** set kmer consensus position */
		HashMap<Kmer, ArrayList<Integer>> kmer2pos = new HashMap<Kmer, ArrayList<Integer>>();
		for (Sequence s:seqList){
			if (s.pos != UNALIGNED){		// aligned seqs
				HashMap<Kmer, HashSet<Integer>> kmerPos = s.getKmerPos();
				for (Kmer km:kmerPos.keySet()){
					if (excludes.contains(km))
						continue;
//						if (km.getKmerString().equals("CCACGCG")||km.getKmerRC().equals("CCACGCG"))
//							km.getK();;
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
//				if (km.getKmerString().equals("CCACGCG")||km.getKmerRC().equals("CCACGCG"))
//					km.getK();
			ArrayList<Integer> posKmer = kmer2pos.get(km);
			if (posKmer==null || posKmer.isEmpty()){
				km.setAlignString("NOT in 2k region");
				continue;
			}
			if (posKmer.size()<km.getPosHitCount()*kmer_aligned_fraction){			// The kmer hit in 2k region should be at least 1/4 of total hit
				km.setAlignString("Too few hit "+posKmer.size());
				continue;
			}
				
			// find the most frequent kmerPos
			Pair<int[], int[]> sorted = StatUtil.sortByOccurences(posKmer);
			int counts[] = sorted.cdr();
			int posSorted[] = sorted.car();
			int maxCount = counts[counts.length-1];
			if (maxCount<posKmer.size()*0.5 && maxCount<km.getPosHitCount()*0.5)	// for palindromic kmer, posKmer>hitCount
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

	private int buildPWM(ArrayList<Sequence> seqList, KmerCluster cluster, long tic){
		WeightMatrixScorer scorer = null;
		if (cluster.wm!=null && cluster.pwmGoodQuality){
	        scorer = new WeightMatrixScorer(cluster.wm);
		}
    	
		ArrayList<String> alignedSeqs = new ArrayList<String>();
		for (Sequence seq:seqList){
			if (seq.pos==UNALIGNED)
				continue;
			if (scorer!=null){
				Pair<Integer, Double> hit = CommonUtils.scanPWM(seq.getSeq(), cluster.wm.length(), scorer);
				if (hit.cdr()<cluster.pwmThreshold)
					continue;
				if (hit.car()+seq.pos-cluster.pos_pwm_seed!=0)
					continue;
			}
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
			if (end>k_win){
				endPadding = CommonUtils.padding(end-k_win, "N");
				end=k_win;
			}
			if (start>=end)
				continue;
 			String s = startPadding+seq.getSeq().substring(start, end)+endPadding;
 			alignedSeqs.add(s);
    	}
		if (alignedSeqs.size()<seqs.length*motif_hit_factor){
			if (verbose>1)
	    		System.out.println(String.format("%s: Number of sequences %d is too few, stop building PWM.", CommonUtils.timeElapsed(tic), alignedSeqs.size()));
			return -1;
		}
		
		int bestLeft = buildPWMfromAlignedSequences(alignedSeqs, cluster, true);
		
    	return bestLeft;
    }
    
	private String makeTRANSFAC (float[][] pfm, String header){
		// make string in TRANSFAC format
		StringBuilder msb = new StringBuilder();
		msb.append(header);
		for (int p=0;p<pfm.length;p++){
			msb.append(p+1).append(" ");
			int maxBase = 0;
			float maxCount=0;
			for (int b=0;b<LETTERS.length;b++){
				msb.append(String.format("%d ", (int)pfm[p][LETTERS[b]]));
				if (maxCount<pfm[p][LETTERS[b]]){
					maxCount=pfm[p][LETTERS[b]];
					maxBase = b;
				}
			}
			msb.append(LETTERS[maxBase]).append("\n");
		}
		msb.append("XX\n\n");
		return msb.toString();
	}
    
    class PWMHit implements Comparable<PWMHit>{
    	int clusterId;
    	WeightMatrix wm;
    	int seqId;
    	boolean isForward;
    	int start;
    	int end;
    	double score;
    	
		public int compareTo(PWMHit h) {					// descending score
			if(score<h.score){return(1);}
			else if(score>h.score){return(-1);}
			else return(0);
		}
		public int compareByPosition(PWMHit h) {					// ascending score
			if(start<h.start){return(-1);}
			else if(start>h.start){return(1);}
			else return(0);
		}
		boolean overlaps(PWMHit h){
			int minHalf = Math.min(end-start+1, h.end-h.start+1)/2;
			if (start+minHalf>h.end || h.start+minHalf>end)
				return false;
			else
				return true;
		}
		public String toString(){
			return String.format("%d:%s==>%s%d:%d-%d", clusterId, WeightMatrix.getMaxLetters(wm), isForward?"":"-", seqId, start, end);
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
	private void alignSequencesByCoOccurence(ArrayList<Kmer> kmers){
		/* Initialization, setup sequence list, and update kmers */
		// build the kmer search tree
		AhoCorasick oks = new AhoCorasick();
		HashMap<String, Kmer> str2kmer = new HashMap<String, Kmer>();
		for (Kmer km: kmers){
			str2kmer.put(km.getKmerString(), km);
			oks.add(km.getKmerString().getBytes(), km);
	    }
		oks.prepare();
		
		// index kmer->seq, seq->kmer
		ArrayList<Sequence> seqList = new ArrayList<Sequence>();
		HashMap <Kmer, HashSet<Integer>> kmer2seq = new HashMap <Kmer, HashSet<Integer>>();
		for (int i=0;i<seqs.length;i++){
			String seq = seqs[i];
			Sequence s = new Sequence(seq, i+1);
			seqList.add(s);
			HashSet<Kmer> results = KmerEngine.queryTree (seq, oks);
			if (!results.isEmpty()){
				for (Kmer km: results){		
					if (!kmer2seq.containsKey(km)){
						kmer2seq.put(km, new HashSet<Integer>());
					}
					kmer2seq.get(km).add(i);

					// forward strand
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
						s.fPos.get(km).add(getRcCoor(p)+RC);					// 5' end if sequence is RC
					}
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
						s.rPos.get(km).add(getRcCoor(p)+RC);					// 5' end if sequence is RC
					}
					s.setCount();
				}
			}
		}
		// update kmerCount, and hgp()
		ArrayList<Kmer> zeroCount = new ArrayList<Kmer>();
		for (Kmer km:kmers){
			if (kmer2seq.containsKey(km)){
				km.setPosHits(kmer2seq.get(km));
				km.setHgp(computeHGP(km.getPosHitCount(), km.getNegHitCount()));
			}
			else{
				zeroCount.add(km);
			}
		}
		kmers.removeAll(zeroCount);
		
		Collections.sort(seqList);
		
		ArrayList<Sequence> seqAligned = new ArrayList<Sequence>();
		ArrayList<Sequence> seqBadScore = new ArrayList<Sequence>();
		Sequence top = seqList.get(0);
		top.pos = 0;
		seqAligned.add(top);
		seqList.remove(0);
		
		eachSeq: while(!seqList.isEmpty()){
			top = seqList.get(0);
			if (top.id==396)
				top.getKmerPosCount();
			Pair<Double, Integer> max = findBestPosition(top, seqAligned);
			if (max.car()<=0){
				seqBadScore.add(top);
				seqList.remove(0);
				continue eachSeq;
			}
			top.score = max.car();
			if (max.cdr()>RC/2){
				top.pos = max.cdr()-RC;
				top.isForward = false;
			}
			else{
				top.pos = max.cdr();
				top.isForward = true;
			}
			seqAligned.add(top);
			seqList.remove(0);
		}

    	int leftmost = Integer.MAX_VALUE;
		for (int i=0;i<seqAligned.size();i++){
			int pos = seqAligned.get(i).pos;
			if (pos < leftmost )
				leftmost = pos;		
		}

		StringBuilder sb = new StringBuilder();
		for (int i=0;i<seqAligned.size();i++){
			Sequence s = seqAligned.get(i);
			sb.append(String.format("%d\t%.1f\t%d\t%s\t%s%s\n", s.id, s.score, s.pos, s.isForward?"F":"R", CommonUtils.padding(-leftmost+s.pos, '.'), s.getSeq()));
		}
		CommonUtils.writeFile("seqs_aligned.txt", sb.toString());
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
	
	private class KmerCluster implements Comparable<KmerCluster>{
		int clusterId;
		Kmer seedKmer;
		String pfmString="";
		WeightMatrix wm;
		boolean pwmGoodQuality = false;
		double pwmThreshold;
		double pwmThresholdHGP;
		int pwmNegHitCount;
		int pwmPosHitCount;
		int pos_pwm_seed;
		int pos_BS_seed;
		ArrayList<Kmer> alignedKmers;
		MotifThreshold ksmThreshold = new MotifThreshold();
		int total_aligned_seqs;
		HashMap<Integer, PWMHit> seq2hits = null;
		
		KmerCluster(){}	// empty constructor;

		protected KmerCluster clone(){
			KmerCluster cluster = new KmerCluster();
			cluster.clusterId = this.clusterId;
			cluster.seedKmer = this.seedKmer;
			cluster.pfmString = this.pfmString;
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
		ArrayList<Region> negRegions = new ArrayList<Region>();
		negRegions.addAll(neg_region_map.keySet());
		Region.filterOverlapRegions(negRegions, posImpactRegion);	// make sure negative region is not within negRegionDistance of positive regions.
		seqsNegList.clear();
		int negCount = 0;
		for (Region r:negRegions){
			seqsNegList.add(seqsNeg[neg_region_map.get(r)].substring(0, winSize+1));
			negCount++;
			if (negCount==seqs.length)				// limit the neg region count to be same or less than positive region count
				break;
		}
		
		posSeqCount = seqs.length;
	    negSeqCount = seqsNegList.size();
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
	
	public double computeHGP(int posHitCount, int negHitCount){
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
//			km.setStrength(kmerStrength[i]);
	    }
	}
	
	/**
	 * Estimate threshold of a PWM (larger than startingScore) using the positive/negative sequences<br>
	 * Multi-thread to compute HGP<br>
	 * Approximate grid search, to reduce run time
	 */
	public MotifThreshold estimatePwmThreshold(WeightMatrix wm, String outName, boolean printPwmHgp, double startingScore){
		double[] posSeqScores = new double[posSeqCount];
		double[] negSeqScores = new double[negSeqCount];
		for (int i=0;i<posSeqCount;i++){
			posSeqScores[i]=WeightMatrixScorer.getMaxSeqScore(wm, seqs[i]);
		}
		Arrays.sort(posSeqScores);
		int zeroIdx = Arrays.binarySearch(posSeqScores, startingScore);
		if( zeroIdx < 0 ) { zeroIdx = -zeroIdx - 1; }
		
		for (int i=0;i<negSeqCount;i++){
			negSeqScores[i]=WeightMatrixScorer.getMaxSeqScore(wm, seqsNegList.get(i));
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

		ArrayList<Integer> idxs = new ArrayList<Integer>();			// the score ids to compute HGP
		for (int i=posScores_u.length-1;i>=0;i--)
			if (poshits[i]>neghits[i]*2.0*posSeqCount/negSeqCount)	// posHit should be at least 2 fold
				idxs.add(i);
		
		Pair<Double, Integer> best;
		
		if (idxs.size()>100 && use_grid_search){
		
			// coarse search
			int gridStep = (int)Math.ceil(Math.sqrt((double)idxs.size()));
			ArrayList<Integer> idxCoarse = new ArrayList<Integer>();	
			for (int i=0;i<idxs.size();i+=gridStep){
				idxCoarse.add(idxs.get(i));
			}
//			if (idxCoarse.get(idxCoarse.size()-1)!=idxs.get(idxs.size()-1))
//				idxCoarse.add(idxs.get(idxs.size()-1));
			
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

	private Pair<Double, Integer> findBestScore(ArrayList<Integer> idxs, int[] poshits, int[] neghits, double[] hgps){

		int numThread = java.lang.Runtime.getRuntime().availableProcessors();
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
		hgps[0]=0;		// the lowest threshold will match all positive sequences, lead to hgp=0
		
		StringBuilder sb = new StringBuilder();
//		if (printPwmHgp){
//			for (int i=0;i<posScores_u.length;i++)
//				sb.append(String.format("%d\t%.2f\t%d\t%d\t%.1f\n", i, posScores_u[i], poshits[i], neghits[i], hgps[i] ));
//		}
		
		Pair<Double, TreeSet<Integer>> minHgp = StatUtil.findMin(hgps);
		int minIdx = minHgp.cdr().last();
		
		return new Pair<Double, Integer>(minHgp.car(),minIdx);
	}
	
	public class MotifThreshold{
		public double score;
		public int posHit;
		public int negHit;
		public double hgp=0;		
	}

	/**
	 * Estimate threshold of a Kmer Group Score using the positive/negative sequences<br>
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

		int numThread = java.lang.Runtime.getRuntime().availableProcessors();
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
		Kmer.printKmers(kmers, posSeqCount, negSeqCount, outPrefix, false, true);
		
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
			KmerGroup kg = new KmerGroup(result.get(p), p);
			matches[idx]=kg;
			kg.setHgp(computeHGP(kg.getGroupHitCount(), kg.getGroupNegHitCount()));
			idx++;
		}
		return matches;
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
			sb.append(km.getKmerString()).append("\t").append(map.get(km.getKmerString())).append("\n");
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
    		double k = kmers.get(0).getK();
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
		public String toString(){
			return String.format("%s %d: %d+/%d-, hpg=%.2f", getBestKmer().getKmerString(), bs, posHitGroupCount, negHitGroupCount, hgp);
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
	public static void main(String[] args){
		ArrayList<String> pos_seqs = new ArrayList<String>();
		try {	
			BufferedReader bin = new BufferedReader(new InputStreamReader(new FileInputStream(new File(args[0]))));
	        String line;
	        while((line = bin.readLine()) != null) { 
	            line = line.trim();
	            pos_seqs.add(line);
	        }			
	        if (bin != null) {
	            bin.close();
	        }
        } catch (IOException e) {
        	System.err.println("Error when processing "+args[0]);
            e.printStackTrace(System.err);
        }
		ArrayList<String> neg_seqs = new ArrayList<String>();
		try {	
			BufferedReader bin = new BufferedReader(new InputStreamReader(new FileInputStream(new File(args[1]))));
	        String line;
	        while((line = bin.readLine()) != null) { 
	            line = line.trim();
	            neg_seqs.add(line);
	        }			
	        if (bin != null) {
	            bin.close();
	        }
        } catch (IOException e) {
        	System.err.println("Error when processing "+args[1]);
            e.printStackTrace(System.err);
        }
        
        KmerMotifFinder kmf = new KmerMotifFinder();
        kmf.setSequences(pos_seqs, neg_seqs);
        kmf.setParameters(-3, 3, 0.01, 0.05, "Test", false, true, 2, 0.2, 0.5);
        int k = kmf.selectK(6, 12);
        ArrayList<Kmer>kmers = kmf.selectEnrichedKmers(k);
//        kmf.clusterKmers(kmers, k/2, 0.3, false);
        kmf.alignBySimplePWM(kmers, k, 0.5, false, false);
	}
}


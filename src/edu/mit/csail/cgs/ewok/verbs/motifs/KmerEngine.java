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
	private static int ITER_INIT=10;
	private static int ITER_END=20;
	private Genome genome;
	private Organism org;
	private String[] args;
	private String motifString;
	private WeightMatrix motif = null;
	private double motifThreshold;
	
	private double HyperGeomThreshold;
	private int max_kmer_per_seq = 5;
	private int seqLength=-1;
	private int k=10;
	private int round = 3;
	private double smooth = 3;
	private int minHitCount = 3;
	private int numPos;
	
	// input data
	private String[] seqs;		// DNA sequences around binding sites
	private double[] seqProbs;	// Relative strength of that binding site, i.e. ChIP-Seq read count
	private Region[] seqCoors;
	private String[] seqsNeg;		// DNA sequences in negative sets
	
	private ArrayList<Kmer> kmers = new ArrayList<Kmer>();
	private HashMap<String, Kmer> str2kmer = new HashMap<String, Kmer>();
	
	// AhoCorasick algorithm for multi-pattern search
	// Pre-processing is to build the tree with all the patters (kmers)
	// Then each individual search can be done in scan()
	private AhoCorasick tree;
	
	private int maxCount;		// max kmer hit count of whole dataset
	public int getMaxCount() {return maxCount;}
	private int minCount;		// min kmer hit count of whole dataset
	public int getMinCount() {return minCount;}
	
	// The Kmers at the specific position: seq--strand--pos
	// Mapping from sequence to kmer, index as seqId, strand, posId 
	// it can contains null, because the match of kmer and kmer_RC only count once
	private Kmer[][][] seq_kmer_map; 
	
	// Corresponding to seq_kmer_map; The probability of this kmer in the sequence
	private double[][][] probsOfKmerInSeq;	// seq_id, pos_id;
	
	// The average profile/density of kmers along the sequence positions
	private double[] positionProbs;
	
	private ArrayList<Kmer> explainedKmers = new ArrayList<Kmer>();
	private long tic;
	
	public KmerEngine(String[] args){
		this.args = args;
		tic = System.currentTimeMillis();
		
	    ArgParser ap = new ArgParser(args);
	    try {
	      Pair<Organism, Genome> pair = Args.parseGenome(args);
	      if(pair==null){
	        //Make fake genome... chr lengths provided???
	        if(ap.hasKey("geninfo")){
	          genome = new Genome("Genome", new File(ap.getKeyValue("geninfo")));
	            }else{
	              System.err.println("No genome provided; provide a Gifford lab DB genome name or a file containing chromosome name/length pairs.");;System.exit(1);
	            }
	      }else{
	        genome = pair.cdr();
	        org = pair.car();
	      }
	    } catch (NotFoundException e) {
	      e.printStackTrace();
	    }
//	    
//	    // load motif
//	    try {
//	      String motifString = Args.parseString(args, "motif", null);
//	      if (motifString!=null){
//		      String motifVersion = Args.parseString(args, "version", null);
//		      motifThreshold = Args.parseDouble(args, "motifThreshold", 10.0);
//		      if (motifThreshold==10.0)
//		    	  System.err.println("No motif threshold was provided, default=10.0 is used.");
//			  int motif_species_id = Args.parseInteger(args, "motif_species_id", -1);
//			  int wmid = WeightMatrix.getWeightMatrixID(motif_species_id!=-1?motif_species_id:org.getDBID(), motifString, motifVersion);
//			  motif = WeightMatrix.getWeightMatrix(wmid);
//	      }
//	    } 
//	    catch (NotFoundException e) {
//	      e.printStackTrace();
//	    }   

	    k = Args.parseInteger(args, "k", k);
	    round = Args.parseInteger(args, "round", round);
	    smooth = Args.parseDouble(args, "smooth", smooth);
	    String filePos = Args.parseString(args, "seqFile", null);
	    String fileNeg = Args.parseString(args, "negSeqFile", null);
	    
	    max_kmer_per_seq = Args.parseInteger(args, "max_kmer", max_kmer_per_seq);
		if (filePos==null||fileNeg==null)
			System.exit(0);

		loadSeqFile(filePos, fileNeg);
		numPos = seqLength-k+1;
		
	}
	
	public KmerEngine(Genome g, ArrayList<ComponentFeature> events, int winSize, double hgp){
		tic = System.currentTimeMillis();
		genome = g;
		int eventCount = events.size();
		seqLength = winSize+1;
		HyperGeomThreshold = hgp;
		Collections.sort(events);		// sort by location
		seqs = new String[eventCount];
		seqsNeg = new String[eventCount];
		seqProbs = new double[eventCount];
		seqCoors = new Region[eventCount];
		SequenceGenerator<Region> seqgen = new SequenceGenerator<Region>();
		seqgen.useCache(true);

		ArrayList<Region> negRegions = new ArrayList<Region>();
		// prepare for progress reporting
        Vector<Integer> processRegionCount = new Vector<Integer>();		// for counting how many regions are processed by all threads
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
			ComponentFeature f = events.get(i);
			seqProbs[i] = f.getTotalSumResponsibility();
			Region posRegion = f.getPeak().expand(winSize/2);
			seqCoors[i] = posRegion;
			seqs[i] = seqgen.execute(seqCoors[i]).toUpperCase();
			// getting negative sequences
			// exclude negative regions that overlap with positive regions, or exceed start of chrom
			// it is OK if we lose a few sequences here, so some entries of the seqsNeg will be null
			int start = posRegion.getStart()-winSize;
			int end = posRegion.getEnd()-winSize;
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
				System.out.println(trigger+"\t/"+eventCount+"\t"+CommonUtils.timeElapsed(tic));
				reportTriggers.remove(reportTriggers.first());
            }
            seqsNeg[i] = seqgen.execute(negRegion).toUpperCase();
		}
		System.out.println(eventCount+"\t/"+eventCount+"\t"+CommonUtils.timeElapsed(tic));
	}
	
	/*
	 * Find significant Kmers that have high HyperGeometric p-value
	 * and build the kmer AhoCorasick engine
	 */
	public void buildEngine(int k, String outPrefix) {
		this.k = k;
		numPos = seqLength-k+1;
		
		HashMap<String, ArrayList<KmerMatch>> map = new HashMap<String, ArrayList<KmerMatch>>();
		ArrayList<Kmer> kmers = new ArrayList<Kmer>();
		for (int seqId=0;seqId<seqs.length;seqId++){
//			System.out.print(seqId+" ");
			String seq = seqs[seqId];
			for (int i=0;i<numPos;i++){
				if ((i+k)>seq.length()) // endIndex of substring is exclusive
					break;
				String s = seq.substring(i, i+k);
				if (map.containsKey(s)){
					 ArrayList<KmerMatch> hits = map.get(s);
					 hits.add(new KmerMatch(seqId, i, 0));
				}
				else{
					 ArrayList<KmerMatch> hits = new ArrayList<KmerMatch>();
					 hits.add(new KmerMatch(seqId, i, 0));
					 map.put(s, hits);
				}
			}
		}
		System.out.println("\nKmers indexed "+timeElapsed(tic));
		
		seq_kmer_map = new Kmer[seqs.length][2][numPos];
	
		// sort the kmer strings
		//so that we can only use one kmer to represent its reverse compliment (RC)
		ArrayList<String> sortedKeys = new ArrayList<String>();
		sortedKeys.addAll(map.keySet());
		Collections.sort(sortedKeys);
		
		for (String key:sortedKeys){
			if (!map.containsKey(key))		// this kmer has been removed, represented by RC
				continue;

			// consolidate kmer and its reverseComplment kmer
			ArrayList<KmerMatch> hits_all = new ArrayList<KmerMatch>();
			hits_all.addAll(map.get(key));
			String key_rc = SequenceUtils.reverseComplement(key).toUpperCase();
			ArrayList<KmerMatch> hits_rc = map.get(key_rc);
			if (!key_rc.equals(key)){	// if it is not reverse compliment itself
				map.remove(key_rc);		// remove this kmer because it is represented by its RC
				if (hits_rc!=null){
					for (KmerMatch hit: hits_rc){
						hits_all.add(new KmerMatch(hit.seqId, hit.pos, 1));
					}
				}
			}
			
			if (hits_all.size()< minHitCount)
				continue;	// skip low count (<2) kmers, 
			
			// create the kmer object
			Kmer kmer = new Kmer(key, hits_all);
			kmers.add(kmer);
		}
		map=null;
		System.gc();
		System.out.println("Kmers("+kmers.size()+") mapped "+timeElapsed(tic));
		
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
//			System.out.println(seq);
			while (searcher.hasNext()) {
				SearchResult result = (SearchResult) searcher.next();
				kmerHits.addAll(result.getOutputs());
//				System.out.println(result.getOutputs()+"\tat index: " + result.getLastIndex());
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
		System.out.println("Positive sequences: "+n+" \t"+"Negative sequences: "+negSeqCount);
		
		ArrayList<Kmer> toRemove = new ArrayList<Kmer>();
		for (Kmer kmer:kmers){
			if (kmer.seqHitCount<=1){
				toRemove.add(kmer);	
				continue;
			}
			// count Kmers in negative sets
//			int count = 0;
//			String str_rc = SequenceUtils.reverseComplement(kmer.kmerString);
//			//TODO: index the negative set may be more efficient
//			for (String seq: seqsNeg){
//				if (seq.indexOf(kmer.kmerString)!=-1 || seq.indexOf(str_rc)!=-1)
//					count++;
//			}
//			kmer.negCount = count;
//			assert(kmer.negCount==negHitCounts.get(kmer.kmerString) );
			// add one pseudo-count for negative set (zero count in negative set leads to tiny p-value)
			int kmerAllHitCount = kmer.seqHitCount+1;
			if (negHitCounts.containsKey(kmer.kmerString)){
				kmer.negCount = negHitCounts.get(kmer.kmerString);
				kmerAllHitCount += kmer.negCount;
				
			}
			double p = 1-StatUtil.hyperGeometricCDF(kmer.seqHitCount, N, kmerAllHitCount, n);
			kmer.hg = p;
//			System.out.println(String.format("%s\t%d\t%.4f", kmer.kmerString, kmer.seqHitCount, kmer.hg));
			if (kmer.hg>HyperGeomThreshold)
				toRemove.add(kmer);		
		}
		// remove un-enriched kmers		
		kmers.removeAll(toRemove);
		System.out.println("\nKmers selected "+timeElapsed(tic));

		// set Kmers and prepare the search Engine
		loadKmers(kmers, outPrefix);
	}
	
	/** load Kmers and prepare the search Engine
	 * 
	 * @param kmers List of kmers (with kmerString, sequence hit count)
	 */
	public void loadKmers(ArrayList<Kmer> kmers, String outPrefix){
		tic = System.currentTimeMillis();
		this.kmers = kmers;
		Collections.sort(kmers);
		this.maxCount = kmers.get(0).seqHitCount;
		this.minCount = kmers.get(kmers.size()-1).seqHitCount;
		printKmers(kmers, outPrefix);
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
	    System.out.println("Kmers("+kmers.size()+") loaded to the Kmer Engine "+timeElapsed(tic));
	}	
	
	/** 
	 * Search all kmers in the sequence
	 * @param seq
	 * @return a map (pos-->kmer)
	 * if pos is negative, then the start position maches the reverse compliment of kmer string
	 */
	public HashMap<Integer, Kmer> query (String seq){
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
		HashMap<Integer, Kmer> result = new HashMap<Integer, Kmer> ();		
		for (Object o: kmerFound){
			String kmerStr = (String) o;
			ArrayList<Integer> pos = StringUtils.findAllOccurences(seq, kmerStr);
			for (int p: pos){
				result.put(p, str2kmer.get(kmerStr));
			}
			ArrayList<Integer> pos_rc = StringUtils.findAllOccurences(seq, SequenceUtils.reverseComplement(kmerStr));
			for (int p: pos_rc){
				result.put(-p, str2kmer.get(kmerStr));
			}
		}

		return result;
	}
	
	
	// hard assignment by exclusion, 1 kmer will explain the whole sequence
	private void hardExclusion(){
		//sort kmers
		Collections.sort(kmers);
		System.out.println("\nKmers sorted "+timeElapsed(tic));
		
		HashSet<Integer> unexplainedSeq = new HashSet<Integer>();
		for (int i=0;i<seqs.length;i++){
			unexplainedSeq.add(i);
		}
		
		ArrayList<Kmer> temp = new ArrayList<Kmer>();
		for(Kmer kmer:kmers){
			kmer.calcWeight(seqProbs, positionProbs);
			if (kmer.seqHitCount>=minHitCount)
				temp.add(kmer);
		}
		
		while(!unexplainedSeq.isEmpty() && !temp.isEmpty()){
			Kmer top = temp.get(0);
			explainedKmers.add(top);
			temp.remove(0);
			for (Kmer kmer:temp){
				kmer.seqHits.removeAll(top.seqHits);
				kmer.calcWeight(seqProbs, positionProbs);
			}
			unexplainedSeq.removeAll(top.seqHits);
			Collections.sort(temp);
		}
		System.out.println("\nSequences explained "+timeElapsed(tic));
	}
	
	private void softWeighting(){
		for (int r=0;r<round;r++){
			
			//init			
			probsOfKmerInSeq = new double[seq_kmer_map.length][2][numPos];
			
			for (int i=0;i<seq_kmer_map.length;i++){
				double sum=0;
				for (int j=0;j<=1;j++){
					for(int p=0;p<numPos;p++){
						if (seq_kmer_map[i][j][p]!=null){
							probsOfKmerInSeq[i][j][p] = positionProbs[p];
							sum += positionProbs[p];
						}
						else
							probsOfKmerInSeq[i][j][p]=0;
					}
				}
				// normalize kmer prob for each sequece i
				if (sum!=0)
					for (int j=0;j<=1;j++)
						for(int p=0;p<numPos;p++)
							probsOfKmerInSeq[i][j][p]/=sum ;
			}
			
			//iterate
			for (int iter=0;iter<ITER_INIT+numPos-5+ITER_END;iter++){
				// Kmer probs is the sum of the kmer's weight from each sequence
				// i.e. The number of binding strength explained by this kmer over whole genome
				for (Kmer kmer:kmers){
					double kprob=0;
					for (KmerMatch h: kmer.hits){
						kprob+=seqProbs[h.seqId]*probsOfKmerInSeq[h.seqId][h.isMinus][h.pos];
					}
//					kmer.weight = kprob;
					kmer.strength = kprob/kmer.count*kmer.seqHitCount;		// re-scale to disfavor simple repeat pattern
				}
				// re-calculate kmer's contribution inside each sequence
				double[] kmerProbs = new double[numPos*2];
				for (int i=0;i<seq_kmer_map.length;i++){
					double sum = 0;
					for (int j=0;j<=1;j++){
						for(int p=0;p<numPos;p++){
							if (probsOfKmerInSeq[i][j][p]!=0){
								double weight=positionProbs[p]*seq_kmer_map[i][j][p].strength;		// influence of position
								probsOfKmerInSeq[i][j][p] = weight;
								sum += weight;
							}
						}
					}
					// after burn-in period, before Iter_end, eliminate low prob position
					if (iter>ITER_INIT && iter<=ITER_INIT+numPos-5){
						for (int j=0;j<=1;j++)
							for(int p=0;p<numPos;p++)
								kmerProbs[j*numPos+p]=probsOfKmerInSeq[i][j][p];
						Arrays.sort(kmerProbs);
						double cutoff = kmerProbs[numPos*2-max_kmer_per_seq];	// keep 5 kmer per sequence
						for (int j=0;j<=1;j++)
							for(int p=0;p<numPos;p++)
								if (probsOfKmerInSeq[i][j][p]<cutoff){	// if less than at least 5th min
									sum-=probsOfKmerInSeq[i][j][p];
									probsOfKmerInSeq[i][j][p]=0;
								}
					}
					// normalize
					if (sum!=0)
						for (int j=0;j<=1;j++)
							for(int p=0;p<numPos;p++)
								if (probsOfKmerInSeq[i][j][p]!=0)
									probsOfKmerInSeq[i][j][p]/=sum ;
				}
	//			System.out.print(iter+" ");
			}// iterate
			
			// re-estimate position probabilities
			positionProbs = new double[seqLength-k+1];
			for (Kmer m : kmers){
				if (m.strength!=0){
					float[] counts = m.getPositionCounts(seqLength);
					for (int i=0;i<positionProbs.length;i++){
						positionProbs[i] += counts[i];
					}
				}
			}
			positionProbs = StatUtil.gaussianSmoother(positionProbs, smooth);
//			positionProbs=StatUtil.cubicSpline(positionProbs, 3, 3);
			printPositionProbabilities(String.format("kmer_%d_%d_p_pos.txt", k, r+1));
			printKmers(kmers, ""+r+1);
			printSeqStats(r+1);
		} // round
	}
	
	private void printPositionProbabilities(String name){
		StringBuilder sb = new StringBuilder();
		for (int i=0;i<positionProbs.length;i++){
			sb.append(String.format("%d\t%.4f\n",i-seqLength/2, positionProbs[i]));
		}
		CommonUtils.writeFile(name, sb.toString());
	}
	private void printKmers(ArrayList<Kmer> kmers, String outPrefix){
		Collections.sort(kmers);
		
		StringBuilder sb = new StringBuilder();
		sb.append(Kmer.toHeader());
		sb.append("\n");
		for (Kmer kmer:kmers){
			sb.append(kmer.toString()).append("\n");
		}
		CommonUtils.writeFile(String.format("%s_kmer_%d.txt",outPrefix, k), sb.toString());
	}
	
	
	// find all occurrences of kmer from a list of sequences
	private void indexKmers(){
		buildEngine(k, "0");
		
		StringBuilder sb = new StringBuilder();
		for (Kmer kmer:kmers){
			sb.append(String.format("%s\t%d\t%d\t%.4f\n", kmer.kmerString, kmer.seqHitCount, kmer.negCount, kmer.hg));
		}
		CommonUtils.writeFile("Kmer_p_values.txt", sb.toString());
		System.out.println("Enriched Kmer count:\t" + kmers.size());
		System.out.println("\nKmers scored "+timeElapsed(tic));
		
		// map the kmers to position
		for (Kmer kmer:kmers){
			for (KmerMatch hit:kmer.hits){
				seq_kmer_map[hit.seqId][hit.isMinus][hit.pos] = kmer;
			}
		}
		
		// build positional probabilities using the high count kmers
		positionProbs = new double[seqLength-k+1];
		for (Kmer m : kmers){
			if (m.count>3){
				float[] counts = m.getPositionCounts(seqLength);
				for (int i=0;i<positionProbs.length;i++){
					positionProbs[i] += counts[i];
				}
			}
		}
		StatUtil.normalize(positionProbs);
		positionProbs=StatUtil.gaussianSmoother(positionProbs, smooth);
		printPositionProbabilities("kmer_"+k+"_0_p_pos.txt");
		printKmers(kmers, ""+0);
	}
	
	private void printSeqStats(int round){
		StringBuilder sb = new StringBuilder();
		StringBuilder sb2 = new StringBuilder();
		for (int i=0; i<seqs.length; i++){
			sb.append(i).append(" ").append(seqProbs[i]).append("\t\t").append(seqCoors[i].getChrom()+":"+seqCoors[i].getStart()).append("\n");
			ArrayList<Kmer> kmers_seq = new ArrayList<Kmer>();
			Kmer[][] temp = seq_kmer_map[i];
			for (int s=0;s<2;s++)
				for (int p=0;p<temp[0].length;p++)
					if (temp[s][p]!=null)
						kmers_seq.add(seq_kmer_map[i][s][p]);
			Collections.sort(kmers_seq);
			Collections.sort(kmers_seq, new Comparator<Kmer>(){
			    public int compare(Kmer o1, Kmer o2) {
			    		return o1.compareByWeight(o2);
			    }
			});
			char[] letters = seqs[i].toLowerCase().toCharArray();
			char[] symbols = new char[letters.length];
			for (int s=0;s<symbols.length;s++){
				symbols[s]='.';
			}
			for (Kmer kmer : kmers_seq){
				for (KmerMatch hit : kmer.hits){
					if (hit.seqId == i){
						double resp = probsOfKmerInSeq[i][hit.isMinus][hit.pos]*seqProbs[i];
						if (resp>6){	// display only if this kmer can support at least 6 read
							String match = hit.isMinus==0?kmer.kmerString:SequenceUtils.reverseComplement(kmer.kmerString);
							for (int p=0;p<kmer.k;p++){
								letters[hit.pos+p]=match.charAt(p);
								symbols[hit.pos+p]= hit.isMinus==0? '+':'-';
							}
							sb.append(String.format("%s/%s %.2f\t%.1f\t/%.2f\n", kmer.kmerString, 
									SequenceUtils.reverseComplement(kmer.kmerString),
									probsOfKmerInSeq[i][hit.isMinus][hit.pos],
									resp,
									kmer.strength));
							// print the start position of kmers, to be consistent with warp-drive motif track
							sb2.append(seqCoors[i].getChrom()+":"+(seqCoors[i].getStart()+hit.pos)).append("\n"); 
						}
					}
				}
			}
			
			sb.append(letters).append("\n").append(symbols).append("\n").append("\n");
		}
		CommonUtils.writeFile("kmer_"+k+"_"+round+"_seq.txt", sb.toString());
		CommonUtils.writeFile("kmer_"+k+"_"+round+"_coords.txt", sb2.toString());
		System.out.println("print seq stats " + timeElapsed(tic));
	}
	
	private float scorePWM(String seq){
		float score = 0;
//		char[] chars = seq.toCharArray();
//        for (int j = 0; j < motif.length(); j++) {
//            score += motif.matrix[j][chars[j]];
//        }
        return score;
	}
	
	// print top 5 kmers with seq:pos distribution
	public void print5Kmers(){
		for (int i=0;i<5;i++){
			StringBuilder sb = new StringBuilder();
			int index=i*100;
			Kmer kmer = kmers.get(index);
			ArrayList<KmerMatch> hits = kmer.hits;
			for (KmerMatch hit:hits){
				sb.append(hit.seqId).append("\t")
				  .append(hit.pos).append("\t")
				  .append(hit.isMinus).append("\n");
			}
			CommonUtils.writeFile("kmer_"+k+"_"+index+"_"+kmer.kmerString+".txt", sb.toString());
		}
	}
	
	public static void main(String[] args){
		KmerEngine analysis = new KmerEngine(args);
		analysis.indexKmers();
		analysis.softWeighting();
//		analysis.hardExclusion();
//		analysis.print5Kmers();
	}	
	
	// load text file with sequences
	private void loadSeqFile(String posFile, String negFile) {
		BufferedReader bin=null;
		try {
			ArrayList<Region> regions = new ArrayList<Region>();
			ArrayList<Double> probs = new ArrayList<Double>();
			ArrayList<String> ss = new ArrayList<String>();
			bin = new BufferedReader(new FileReader(new File(posFile)));
			String line;
			while((line = bin.readLine()) != null) { 
				if (line.charAt(0)=='>'){
					String[] items = line.split("\t");
					if (seqLength==-1)
						seqLength = Integer.parseInt(items[1]);
					regions.add(Region.fromString(genome, items[0].substring(1)));
					probs.add(Double.parseDouble(items[2]));
				}
				else{
					ss.add(line.trim().toUpperCase());
				}
			}
			
			if (ss.size()!=probs.size()){
				System.err.println("The header line count does not match sequence count!");
				System.exit(0);
			}
				
			seqs = new String[ss.size()];
			seqProbs = new double[ss.size()];
			seqCoors = new Region[ss.size()];
			ss.toArray(seqs);
			for (int i=0;i<ss.size();i++){
				seqProbs[i]=probs.get(i);
				seqCoors[i]=regions.get(i);
			}

			// load negative sets
			bin = new BufferedReader(new FileReader(new File(negFile)));
			ss.clear();
			while((line = bin.readLine()) != null) { 
				line=line.trim();
				if (line.length()==0) continue;
				if (line.charAt(0)=='>'){
					String[] items = line.split("\t");
					regions.add(Region.fromString(genome, items[0].substring(1)));
					probs.add(Double.parseDouble(items[2]));
				}
				else{
					ss.add(line.toUpperCase());
				}
			}
			seqsNeg = new String[ss.size()];
			ss.toArray(seqsNeg);
		}
		catch(IOException ioex) {
			ioex.printStackTrace();
		}
		finally {
			try {
				if (bin != null) {
					bin.close();
				}
			}
			catch(IOException ioex2) {
				System.err.println("Error closing buffered reader");
				ioex2.printStackTrace(System.err);
			}			
		}
	}
	private static String timeElapsed(long tic){
		return timeString(System.currentTimeMillis()-tic);
	}
	private static String timeString(long length){
		float sec = length/1000F;
		return sec>60?
			String.format("%.1f",sec/60)+" min":
			String.format("%.1f",sec)+" sec";
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


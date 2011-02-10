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
	
	private int maxShift;		// max kmer flanking Shift of whole dataset
	public int getMaxShift() {return maxShift;}
	public void setMaxShift(int maxShift) {this.maxShift = maxShift;}
	private int minShift;		// min kmer flanking Shift of whole dataset
	public int getMinShift() {return minShift;}
	public void setMinShift(int minShift) {this.minShift = minShift;}
	
	// The average profile/density of kmers along the sequence positions
	private double[] positionProbs;
	
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
	
	public KmerEngine(Genome g, ArrayList<ComponentFeature> events, int winSize, int winShift, double hgp){
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
		
		HashMap<String, Integer> map = new HashMap<String, Integer>();
		ArrayList<Kmer> kmers = new ArrayList<Kmer>();
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
		System.out.println("\nKmers indexed "+timeElapsed(tic));
	
		// sort the kmer strings
		//so that we can only use one kmer to represent its reverse compliment (RC)
		ArrayList<String> sortedKeys = new ArrayList<String>();
		sortedKeys.addAll(map.keySet());
		Collections.sort(sortedKeys);
		
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
			
			if (map.get(key)< minHitCount)
				continue;	// skip low count (<2) kmers, 
			
			// create the kmer object
			Kmer kmer = new Kmer(key, map.get(key));
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
			Kmer kmer = str2kmer.get(kmerStr);
			ArrayList<Integer> pos = StringUtils.findAllOccurences(seq, kmerStr);
			for (int p: pos){
				result.put(p-kmer.getGlobalShift(), kmer);
			}
			ArrayList<Integer> pos_rc = StringUtils.findAllOccurences(seq, SequenceUtils.reverseComplement(kmerStr));
			for (int p: pos_rc){
				result.put(-(p-kmer.getGlobalShift()), kmer);
			}
		}

		return result;
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
	
	
	private float scorePWM(String seq){
		float score = 0;
//		char[] chars = seq.toCharArray();
//        for (int j = 0; j < motif.length(); j++) {
//            score += motif.matrix[j][chars[j]];
//        }
        return score;
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


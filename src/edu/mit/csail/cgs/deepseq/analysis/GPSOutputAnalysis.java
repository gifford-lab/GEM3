package edu.mit.csail.cgs.deepseq.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.discovery.BindingMixture;
import edu.mit.csail.cgs.deepseq.discovery.kmer.Kmer;
import edu.mit.csail.cgs.deepseq.features.ComponentFeature;
import edu.mit.csail.cgs.deepseq.features.Feature;
import edu.mit.csail.cgs.deepseq.utilities.AnnotationLoader;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSPeak;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public class GPSOutputAnalysis {
  static final int NOHIT_OFFSET = 999;
  static final int maxAnnotDistance=50000;
  private final char[] LETTERS = {'A','C','G','T'};
  
  private int motif_window = 100;  
  private Genome genome;
  private String[] args;
  private WeightMatrix motif = null;
  private double motifThreshold;
  private int extend;				// number of bases to extend from motif hit sequence
  
  private boolean useWeight;
  
  private int[]motif_offsets;
  private double[]motif_scores;
  
  private List<GPSPeak> gpsPeaks;
  private String outputFileName;
  
  // build empirical distribution
  private String chipSeqExpt = null;
  private String chipSeqVersion = null;
  private boolean useMotif = true;
  
  
  /**
   * @param args
   */
  public static void main(String[] args) throws IOException {
    ArgParser ap = new ArgParser(args);
    Set<String> flags = Args.parseFlags(args);
    Organism org=null;
    Genome genome=null;
    
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
    // kmer motif 
    boolean useWeight = flags.contains("use_weight");
    
	// load motif
    String motifString = Args.parseString(args, "motif", null);
    Pair<WeightMatrix, Double> wm = null;
    if (motifString!=null){
		wm = CommonUtils.loadPWM(args, org.getDBID());
    }
    int motif_window = Args.parseInteger(args, "motif_window", 100);
    int extend = Args.parseInteger(args, "extend", 0);
    
    // load GPS results
    String GPSfileName = Args.parseString(args, "GPS", null);
    if (GPSfileName==null){
      System.err.println("GPS file not found!");
      System.exit(0);
    }
    File gpsFile = new File(GPSfileName);
    List<GPSPeak> gpsPeaks = GPSParser.parseGPSOutput(gpsFile.getAbsolutePath(), genome);
//		Collections.sort(gpsPeaks, new Comparator<GPSPeak>(){
//		    public int compare(GPSPeak o1, GPSPeak o2) {
//		        return o1.compareByPValue(o2);
//		    }
//		});
		
		
    GPSOutputAnalysis analysis = new GPSOutputAnalysis(genome, 
    		wm.car(), wm.cdr().doubleValue(), gpsPeaks, GPSfileName, motif_window, extend, useWeight);
	
    int type = Args.parseInteger(args, "type", 0);
    int win = Args.parseInteger(args, "win", 100);
    int top = Args.parseInteger(args, "top", -1);
	switch(type){
	case 0: analysis.jointBindingMotifAnalysis(true);break;
	case 1: analysis.printSequences(win, top);break;
	case 2:	analysis.geneAnnotation();break;
	case 3:	analysis.expressionIntegration();break;
	case 4: analysis.printMotifHitSequence(top);break;
	case 5: analysis.clusterSequences();break;
	default: System.err.println("Unrecognize analysis type: "+type);
	}
  }
  
  public GPSOutputAnalysis(Genome g, WeightMatrix wm, double threshold, 
                           List<GPSPeak> p, String outputFile, int motif_win, 
                           int extend, boolean useWeight) {
	  genome = g;
	  motif = wm;
	  motifThreshold = threshold;
	  gpsPeaks = p;
	  outputFileName = outputFile;
	  motif_window = motif_win;
	  this.extend = extend;
	  this.useWeight = useWeight;
  }
  
  
  public void clusterSequences(){
	ArrayList<GPSPeak> unalignedPeaks = new ArrayList<GPSPeak>();
	for(GPSPeak p : gpsPeaks){
		String kmer = p.getKmer();
		if (kmer!=null && !kmer.equals(""))
			unalignedPeaks.add(p);
	}
	
	// get kmers and their counts
	HashMap<String, Integer> kmer2count = new HashMap<String, Integer>();
	for(GPSPeak p : gpsPeaks){
		String kmer = p.getKmer();
		if (kmer==null || kmer.equals(""))
			continue;
		if (kmer2count.containsKey(kmer))
			kmer2count.put(kmer, kmer2count.get(kmer)+1);
		else
			kmer2count.put(kmer, 1);
	}
	ArrayList<Kmer> kmers = new ArrayList<Kmer>();
	for (String k : kmer2count.keySet()){
		kmers.add(new Kmer(k, kmer2count.get(k)));
	}
	if (kmers.isEmpty())
		return;
	
	Collections.sort(kmers);
	
	// calc background frequencies of bases in all bound sequence
	double gc = 0.21;
	double[] bg = {0.5-gc,gc,gc,0.5-gc};
//	double[] bg = new double[LETTERS.length];
//	for(GPSPeak p : gpsPeaks){	
//		String seq = p.getBoundSequence();
//		if (seq==null || seq.equals(""))
//			continue;
//		for (int i=0;i<seq.length();i++){
//			for (int b=0;b<LETTERS.length;b++){
//				if (seq.charAt(i)==LETTERS[b])
//					bg[b]++;
//			}
//		}
//	}
//	// normalize
//	double sum=0;
//	for (int b=0;b<LETTERS.length;b++){
//		sum += bg[b];
//	}
//	for (int b=0;b<LETTERS.length;b++){
//		bg[b] /=sum;
//	}
//	// consider 2 strands
//	bg[0]=(bg[0]+bg[3])/2;
//	bg[3]=bg[0];
//	bg[1]=(bg[1]+bg[2])/2;
//	bg[2]=bg[1];
	
	ArrayList<Pair<ArrayList<GPSPeak>, ArrayList<Integer>>> clusters = new ArrayList<Pair<ArrayList<GPSPeak>, ArrayList<Integer>>>();
	while(!unalignedPeaks.isEmpty()){
		Pair<ArrayList<GPSPeak>, ArrayList<Integer>> cluster = growSeqCluster(unalignedPeaks, kmers, bg);
		if (cluster!=null)
			clusters.add(cluster);
		else
			break;
	}
	
	// build WeightMatrix from each aligned group of sequences
	System.out.println("Total clusters of motifs: "+clusters.size());
	int count = 0;
	StringBuilder sb = new StringBuilder();
	for (Pair<ArrayList<GPSPeak>, ArrayList<Integer>> cluster : clusters){
		if (cluster.car().size()<=5)
			continue;
		else
			System.out.println("This cluster of motifs is learned from "+cluster.car().size()+" binding events.");
		WeightMatrix wm = makePWM( cluster.car(), cluster.cdr(), bg, false);
//		System.out.println(WeightMatrix.printMatrix(wm));
		System.out.println(WeightMatrix.printMatrixLetters(wm));
		
		sb.append(getPFMString(cluster.car(), cluster.cdr(), wm.length(), count++));
	}
	CommonUtils.writeFile("Ctcf_PFM.txt", sb.toString());
  }
  
  // make PFM from aligned sequences
  private String getPFMString(ArrayList<GPSPeak> alignedPeaks, ArrayList<Integer> peakShifts, int length, int num){
    int MAXLETTERVAL = Math.max(Math.max(Math.max('A','C'),Math.max('T','G')),
            Math.max(Math.max('a','c'),Math.max('t','g'))) + 1;
	boolean useStrength = true;

    String outName = "Ctcf";
	double[][] pfm = new double[length][MAXLETTERVAL];
	for (int i=0;i<alignedPeaks.size();i++){
		GPSPeak peak = alignedPeaks.get(i);
		int shift = peakShifts.get(i);	
		String seq = peak.getBoundSequence().substring(shift, shift+length);		
		for (int p=0;p<length;p++){
			char base = seq.charAt(p);
			double strength = useStrength?alignedPeaks.get(i).getStrength():1;
			pfm[p][base] +=strength;
		}
	}
	
	// make string in TRANSFAC format
	StringBuilder msb = new StringBuilder();
	msb.append(String.format("DE %s_%d_c%d\n", outName, num, alignedPeaks.size()));
	for (int p=0;p<pfm.length;p++){
		msb.append(p+1).append(" ");
		int maxBase = 0;
		double maxCount=0;
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
  
  // greedily grow a cluster from the top count kmer
  private Pair<ArrayList<GPSPeak>, ArrayList<Integer>> growSeqCluster(ArrayList<GPSPeak> unalignedPeaks, ArrayList<Kmer> kmers, double[] bg){
	ArrayList<GPSPeak> alignedPeaks = new ArrayList<GPSPeak>();
	ArrayList<Integer> peakShifts = new ArrayList<Integer>();
	growByKmers(alignedPeaks,peakShifts, kmers, unalignedPeaks);
	if (alignedPeaks.isEmpty())
		return null;
	WeightMatrix wm = makePWM( alignedPeaks, peakShifts, bg, true);

	boolean noMore=false;
	while(!noMore){
		noMore = growByPWM(alignedPeaks, peakShifts, kmers,  unalignedPeaks, wm, bg);
	}
	
	return new Pair<ArrayList<GPSPeak>, ArrayList<Integer>>(alignedPeaks, peakShifts);
  }
	
  //continue to greedily grow a cluster from kmer-aligned peaks by building a PWM
  private boolean growByPWM(ArrayList<GPSPeak> alignedPeaks, ArrayList<Integer> peakShifts, ArrayList<Kmer> kmers, ArrayList<GPSPeak> unalignedPeaks, WeightMatrix wm, double[] bg){
	double wm_factor = 0.3;
	
	// build PWM from Kmer cores
	boolean noMore = true;
	ArrayList<GPSPeak> temp = new ArrayList<GPSPeak>();
	HashSet<String> alignedKmers = new HashSet<String>();
	if (alignedPeaks.size()<=5)
		return true;
	
	double maxWMScore = wm.getMaxScore();	
    WeightMatrixScorer scorer = new WeightMatrixScorer(wm);

	// update the peakShifts for previously identified sequences, w.r.t. PWM
    for (int p=0;p<alignedPeaks.size();p++){
    	GPSPeak peak = alignedPeaks.get(p);
    	String seq = peak.getBoundSequence();
  	  if (seq==null||seq.length()<wm.length()-1){
		  continue;
	  }
        WeightMatrixScoreProfile profiler = scorer.execute(seq);
        double maxSeqScore = Double.NEGATIVE_INFINITY;
        int maxScoringShift = 0;
        char maxScoringStrand = '+';
        for (int i=0;i<profiler.length();i++){
      	  double score = profiler.getMaxScore(i);
      	  if (maxSeqScore<score){
      		  maxSeqScore = score;
      		  maxScoringShift = i;
      		  maxScoringStrand = profiler.getMaxStrand(i);
      	  }
        }
        
  	  	if (maxScoringStrand =='-'){
  		  peak.setBoundSequence(SequenceUtils.reverseComplement(seq));
  		  maxScoringShift = seq.length()-(wm.length())-maxScoringShift;		// WM has a extra column
  	  	}
		peakShifts.set(p, maxScoringShift);		// update
    }

    for (GPSPeak peak: unalignedPeaks){
	  String seq = peak.getBoundSequence();
	  if (seq==null||seq.length()<wm.length()-1){
		  temp.add(peak);
		  continue;
	  }
      WeightMatrixScoreProfile profiler = scorer.execute(seq);
      double maxSeqScore = Double.NEGATIVE_INFINITY;
      int maxScoringShift = 0;
      char maxScoringStrand = '+';
      for (int i=0;i<profiler.length();i++){
    	  double score = profiler.getMaxScore(i);
    	  if (maxSeqScore<score){
    		  maxSeqScore = score;
    		  maxScoringShift = i;
    		  maxScoringStrand = profiler.getMaxStrand(i);
    	  }
      }
      if (maxSeqScore >= maxWMScore * wm_factor){
    	  	if (maxScoringStrand =='-'){
    		  peak.setBoundSequence(SequenceUtils.reverseComplement(seq));
    		  maxScoringShift = seq.length()-wm.length()-maxScoringShift;
    	  	}
			alignedPeaks.add(peak);
			peakShifts.add(maxScoringShift);
			temp.add(peak);
			alignedKmers.add(peak.getKmer());
			noMore = false;
      }
    }

	// sync aligned kmers
	unalignedPeaks.removeAll(temp);
	temp.clear();	
	ArrayList<Kmer> toRemove = new ArrayList<Kmer>();
	for (Kmer kmer: kmers){
		if (alignedKmers.contains(kmer.getKmerString())){
			toRemove.add(kmer);
		}
	}
	kmers.removeAll(toRemove);
	return noMore;
  }
  
  private WeightMatrix makePWM(ArrayList<GPSPeak> alignedPeaks, ArrayList<Integer> peakShifts, double[] bg, boolean isFromKmers){
	double ic_trim = 0.4;
	boolean useStrength = true;
	
	int seqLen = alignedPeaks.get(0).getBoundSequence().length();
	
	// count
	int leftMost = seqLen;
	int shortest = seqLen;	// the shortest length starting from the shift position
	for (int i=0;i<alignedPeaks.size();i++){
		GPSPeak peak = alignedPeaks.get(i);
		String seq = peak.getBoundSequence();
		int shift = peakShifts.get(i);			
		if (isFromKmers){
			int start = seq.indexOf(peak.getKmer());		// start position of kmer in the sequence
			shift += start;
			peakShifts.set(i, shift);
		}
		if (leftMost>shift)
			leftMost = shift;
		if (shortest>seq.length()-shift)
			shortest=seq.length()-shift;
	}

	int MAXLETTERVAL = Math.max(Math.max(Math.max('A','C'),Math.max('T','G')),
            Math.max(Math.max('a','c'),Math.max('t','g'))) + 1;
	double[][] pwm = new double[leftMost+shortest][MAXLETTERVAL];
	for (int i=0;i<alignedPeaks.size();i++){
		int shift = peakShifts.get(i);	
		String seq  = alignedPeaks.get(i).getBoundSequence().substring(shift-leftMost, shift+shortest);
		for (int p=0;p<leftMost+shortest;p++){
			char base = seq.charAt(p);
			double strength = useStrength?alignedPeaks.get(i).getStrength():1;
			pwm[p][base] +=strength;
		}
	}
	
	// normalize and log2
	double[] ic = new double[pwm.length];						// information content
	for (int p=0;p<pwm.length;p++){
		int sum=0;
		for (int b=0;b<LETTERS.length;b++){
			sum += pwm[p][LETTERS[b]];
		}
		for (int b=0;b<LETTERS.length;b++){
			char base = LETTERS[b];
			if (pwm[p][base]==0)
				pwm[p][base]=1;
			double f = pwm[p][base]/sum;						// normalize freq
			pwm[p][base] = Math.log(f/bg[b])/Math.log(2.0);		//log base 2
			ic[p] += f*pwm[p][base];
		}
	}
	// TODO: trim low ic ends
	int leftIdx=0;
	for (int p=0;p<ic.length;p++){
		if (ic[p]>ic_trim){
			leftIdx=p;
			break;
		}
	}
	int rightIdx=ic.length-1;
	for (int p=ic.length-1;p>=0;p--){
		if (ic[p]>ic_trim){
			rightIdx=p;
			break;
		}
	}
	
	// make a WeightMatrix object
	float[][] matrix = new float[rightIdx-leftIdx+1][MAXLETTERVAL];   
	for(int p=leftIdx;p<=rightIdx;p++){
		for (int b=0;b<LETTERS.length;b++){
			matrix[p-leftIdx][LETTERS[b]]=(float) pwm[p][LETTERS[b]];
		}
	}
	WeightMatrix wm = new WeightMatrix(matrix);
//	System.out.println(WeightMatrix.printMatrixLetters(wm));
	return wm;
  }

  //greedily grow a cluster from the top count kmer only by matching kmers
  private void growByKmers(ArrayList<GPSPeak> alignedPeaks, ArrayList<Integer> peakShifts, ArrayList<Kmer> kmers, ArrayList<GPSPeak> unalignedPeaks){
	ArrayList<GPSPeak> temp = new ArrayList<GPSPeak>();
	HashSet<String> alignedKmers = new HashSet<String>();
	String seedKmerStr = kmers.get(0).getKmerString();
	String seedKmerRC = kmers.get(0).getKmerRC();

	// align bound sequences underlying the binding events
	// 1. find perfect kmer matches
	for(GPSPeak p : unalignedPeaks){	
		String kmer = p.getKmer();
		if (kmer==null)
			continue;
		if (kmer.equals(seedKmerStr)){
			alignedPeaks.add(p);
			peakShifts.add(0);
			alignedKmers.add(kmer);
			continue;
		}
		if (kmer.equals(seedKmerRC)){
			alignedPeaks.add(p);
			peakShifts.add(0);
			alignedKmers.add(kmer);
		}
	}
	unalignedPeaks.removeAll(alignedPeaks);
	// 2. with 1 mismatch
	for(GPSPeak p : unalignedPeaks){	
		String kmer = p.getKmer();
		if (kmer==null)
			continue;
		if (mismatch(seedKmerStr, kmer)<=1){
			alignedPeaks.add(p);
			peakShifts.add(0);
			temp.add(p);
			alignedKmers.add(kmer);
		}
		else if(mismatch(seedKmerRC, kmer)<=1){	// if match RC, reset kmer and seq
			p.setBoundSequence(SequenceUtils.reverseComplement(p.getBoundSequence()));
			p.setKmer(SequenceUtils.reverseComplement(p.getKmer()));
			alignedPeaks.add(p);
			peakShifts.add(0);
			temp.add(p);
			alignedKmers.add(kmer);
		}
	}
	unalignedPeaks.removeAll(temp);
	temp.clear();
	// 3. with 1 shift, 0 or 1 mismatch
	for(GPSPeak p : unalignedPeaks){	
		if (p.getKmer()==null)
			continue;
		String kmer = p.getKmer().substring(1);
		String ref = seedKmerStr.substring(0, seedKmerStr.length()-1);
		String refRC = seedKmerRC.substring(1);
		if (mismatch(ref, kmer)<=1){
			alignedPeaks.add(p);
			peakShifts.add(1);
			temp.add(p);
			alignedKmers.add(p.getKmer());
		}
		else if(mismatch(refRC, kmer)<=1){	// if match RC, reset kmer and seq
			p.setBoundSequence(SequenceUtils.reverseComplement(p.getBoundSequence()));
			p.setKmer(SequenceUtils.reverseComplement(p.getKmer()));
			alignedPeaks.add(p);
			peakShifts.add(1);
			temp.add(p);
			alignedKmers.add(p.getKmer());
		}
	}
	unalignedPeaks.removeAll(temp);
	temp.clear();	
	for(GPSPeak p : unalignedPeaks){	
		if (p.getKmer()==null)
			continue;
		String kmer = p.getKmer().substring(0, seedKmerStr.length()-1);
		String ref = seedKmerStr.substring(1);
		String refRC = seedKmerRC.substring(0, seedKmerStr.length()-1);
		if (mismatch(ref, kmer)<=1){
			alignedPeaks.add(p);
			peakShifts.add(-1);
			temp.add(p);
			alignedKmers.add(p.getKmer());
		}
		else if(mismatch(refRC, kmer)<=1){	// if match RC, reset kmer and seq
			p.setBoundSequence(SequenceUtils.reverseComplement(p.getBoundSequence()));
			p.setKmer(SequenceUtils.reverseComplement(p.getKmer()));
			alignedPeaks.add(p);
			peakShifts.add(-1);
			temp.add(p);
			alignedKmers.add(p.getKmer());
		}
	}
	unalignedPeaks.removeAll(temp);
	temp.clear();	
	ArrayList<Kmer> toRemove = new ArrayList<Kmer>();
	for (Kmer kmer: kmers){
		if (alignedKmers.contains(kmer.getKmerString())){
			toRemove.add(kmer);
		}
	}
	kmers.removeAll(toRemove);
  }

private int mismatch(String ref, String seq){
	  if (ref.length()!=seq.length())
		  return -1;
	  int mismatch = 0;
	  for (int i=0;i<ref.length();i++){
		  if (ref.charAt(i)!=seq.charAt(i))
			  mismatch ++;
	  }
	  return mismatch;
  }
  
  private Pair<Integer, Integer> scoreSequence(String ref, String seq){
	int refLen = ref.length();
	int maxScore=0;
	int maxShift=0;
	//try all possible shifts to get the largest score
	byte[] thisBytes = seq.getBytes();
	byte[] refBytes = ref.getBytes();
	for (int s=-refLen;s<=+refLen;s++){
		int largestCount = 0;
		int count = 0;
		for (int i=-refLen;i<refLen+refLen;i++){
			if (i<0 || i>refLen-1 ||i+s<0 || i+s>refLen-1 )
				continue;
			if (refBytes[i]==thisBytes[i+s]){
				count ++;			// increment for match
				if (count>largestCount){
					largestCount = count;
				}
			}else{
				if (count!=0)		// restart when score drops below zero
					count--;		// decrement for mismatch
			}
		}
		
		if (largestCount>maxScore){
			maxScore = largestCount;
			maxShift = s;
		}
	}
	return new Pair<Integer, Integer>(maxScore, maxShift);
  }
  
  
  class KmerGroup{
	  Kmer seed;
	  TreeSet<Kmer> members = new TreeSet<Kmer>();
	  KmerGroup(Kmer seed){
		  this.seed = seed;
	  }
	  void addMember(Kmer member){
		  members.add(member);
	  }
  }
  
  
  /**
   * For each peak find the occurrence of a motif match >= the specified 
   * threshold that is closest to the peak
   * @param analyzeMotif
   */
  public void jointBindingMotifAnalysis(boolean analyzeMotif){
    System.out.println("Analyzing motif matches");
    Collections.sort(gpsPeaks);
    
    if (analyzeMotif)
    	findNearestMotifHit();
    
    int interPeak_distances[] = new int[gpsPeaks.size()];
    interPeak_distances[0]=99999999;
    for (int i=1;i<gpsPeaks.size();i++){
      GPSPeak p = gpsPeaks.get(i);
      GPSPeak prev = gpsPeaks.get(i-1);
      if (p.getChrom().equals(prev.getChrom()))
        interPeak_distances[i]=p.getLocation()-prev.getLocation();
      else
        interPeak_distances[i]=99999999;
    }
     
    //sort by inter-event distance
    int old_index[] = StatUtil.findSort(interPeak_distances);
    
    // output results
    StringBuilder sb = new StringBuilder();
    StringBuilder sb_binary = new StringBuilder();
    // header
    sb.append("# motif: "+motif.getName()+" "+motif.version+" ");
    sb.append("threshold="+motifThreshold+" within "+motif_window+" bps\n");
//    sb.append("GPS Peak\tNearestGene\tDistance\tSample reads\tControl reads\t"+"Motif_score\tMotif_peak\n"); 
    //old style verbose header
    //sb.append("threshold="+motifThreshold+" within "+MOTIF_DISTANCE+" bps\n");
    sb.append(GPSPeak.toGPS_short_Header()+"\tMotif_score\tMotif_peak\tMotif_prevPeak\tPeak_distance\n"); 
    
    // all other peaks
    for (int j=0;j<old_index.length;j++){
    	int i=old_index[j];
    	if (i==0)
    		continue;
//      sb.append(gpsPeaks.get(i).toGPS_motifShort()).append("\t");
      sb.append(gpsPeaks.get(i).toGPS_short()).append("\t");
      sb.append(String.format("%.2f", motif_scores[i])).append("\t");
      sb.append(motif_offsets[i]).append("\t");
      sb.append(motif_offsets[i-1]).append("\t");
      sb.append(interPeak_distances[j]).append("\n");     
      if (motif_offsets[i-1]-motif_offsets[i]!=interPeak_distances[j] 
          && motif_offsets[i]!=NOHIT_OFFSET
          && motif_offsets[i-1]!=NOHIT_OFFSET
          && interPeak_distances[j]<=1000){
    	  sb_binary.append(gpsPeaks.get(i).toGPS_short()).append("\t");
    	  sb_binary.append(String.format("%.2f", motif_scores[i])).append("\t");
    	  sb_binary.append(motif_offsets[i]).append("\t");
    	  sb_binary.append(motif_offsets[i-1]).append("\t");
    	  sb_binary.append(interPeak_distances[j]).append("\n");     
      }
    }
    // first peak
	//  sb.append(gpsPeaks.get(0).toGPS_motifShort()).append("\t");
	  sb.append(gpsPeaks.get(0).toGPS_short()).append("\t");
	  sb.append(String.format("%.2f", motif_scores[0])).append("\t");
	  sb.append(motif_offsets[0]).append("\t");
	  sb.append(NOHIT_OFFSET).append("\t"); //prev motif distance
	  sb.append(99999999).append("\n");  //interpeak distance

  
	  CommonUtils.writeFile(outputFileName+"_JointPeak_Motif.txt", sb.toString());
	  CommonUtils.writeFile(outputFileName+"_BinaryMotifEvents.txt", sb_binary.toString());
  }
  
  /*
   * find nearest motif hit for each peak
   * motif hits are stranded, the offset is the relative position from the motif (positive copy)
   * offset = 999 (NOHIT_OFFSET) if no hit is found
   */
  public void findNearestMotifHit(){
	if (motif_offsets!=null)
		return;
	  
    // motif
    motif_offsets = new int[gpsPeaks.size()];
    motif_scores = new double[gpsPeaks.size()];
      WeightMatrixScorer scorer = new WeightMatrixScorer(motif);
      int halfWidth = motif.length()/2;		// center of motif
      // for each GPS peak, search nearest motif with threshold
      for(int i=0; i<gpsPeaks.size(); i++){
        if (i % 1000 == 0) {
          if (i % 10000 == 0) {
            System.out.println(i);
          }
          else {
            System.out.print(i);
          }
        }
        else if (i % 500==0) {
          System.out.print(".");
        }
        GPSPeak peak = gpsPeaks.get(i);
        Region r= peak.expand(motif_window+halfWidth);
        WeightMatrixScoreProfile profiler = scorer.execute(r);
        //search from BS outwards
        for(int z=0; z<=motif_window; z++){
          double leftScore = Double.NEGATIVE_INFINITY;
          double rightScore = Double.NEGATIVE_INFINITY;
          if ((motif_window+z) < profiler.length()) {
            rightScore= profiler.getMaxScore(motif_window+z);        
            if(rightScore>=motifThreshold){
              motif_offsets[i] = z * (profiler.getMaxStrand(z)=='+'?1:-1);
              motif_scores[i] = rightScore;
              break;
            }
          }
          if ((motif_window-z) >= 0) {
            leftScore= profiler.getMaxScore(motif_window-z);
            if(leftScore>=motifThreshold){
              motif_offsets[i] = -z * (profiler.getMaxStrand(z)=='+'?1:-1);
              motif_scores[i] = leftScore;
              break;
            }
          }
          // if motif score at this position is too small, and reach end of region
          if (z==motif_window){
            motif_offsets[i] = NOHIT_OFFSET;
            motif_scores[i] = Math.max(leftScore, rightScore);
          }
        }
      }
  }
  
  public String printMotifHitList(){
	StringBuilder statStr = new StringBuilder();
	statStr.append(String.format("\n----------------------------------------------\n%s\nTotal Event #:\t%d", 
				outputFileName, gpsPeaks.size()));
	System.out.println(statStr.toString());  
	if (gpsPeaks.size()==0)
		return statStr.toString();
	findNearestMotifHit();
	  
	ArrayList<GPSPeak> hitList = new ArrayList<GPSPeak>();
	ArrayList<GPSPeak> nohitList = new ArrayList<GPSPeak>();
	ArrayList<Integer> offsets_hit= new ArrayList<Integer>();
	for (int i=0;i<gpsPeaks.size();i++){
		  if (motif_offsets[i]!=NOHIT_OFFSET){
			  hitList.add(gpsPeaks.get(i));
			  offsets_hit.add(motif_offsets[i]);
		  }
		  else
			  nohitList.add(gpsPeaks.get(i));
	}
	
	int[]offsets = new int[offsets_hit.size()];
	for (int i=0;i<offsets_hit.size();i++)
		offsets[i]=offsets_hit.get(i);
	double bias = 0;
	// only calc bias when sufficient points, because we assume positive and negative offsets 
	// are balanced when we have sufficient number of points
	if (offsets_hit.size()>200)
		bias= StatUtil.mean(offsets);		// position bias (motif center point may not be binding site)
	for (int i=0;i<offsets_hit.size();i++){
		int unbiasedOffset = offsets_hit.get(i)-(int)bias;
		offsets[i]=Math.abs(unbiasedOffset);
		offsets_hit.set(i, unbiasedOffset);
	}
	double mean = 0;
	if (offsets_hit.size()>0){
		mean = StatUtil.mean(offsets);		// distance mean
		Arrays.sort(offsets);
	}
	
	// motif hit list
	StringBuilder sb = new StringBuilder();
	sb.append("# motif: "+motif.getName()+" "+motif.version+" (" + motif.speciesid+") ");
	sb.append("threshold="+motifThreshold+" within "+motif_window+" bps\n");
	sb.append(GPSPeak.toGPS_short_Header()+"\tOffset\n"); 
	for (int i=0;i<hitList.size();i++){
		sb.append(hitList.get(i).toGPS_short()).append("\t").append(offsets_hit.get(i)).append("\n");
	}
	CommonUtils.writeFile(outputFileName+"_motifHit.txt", sb.toString());
	
	// no hit list
	sb = new StringBuilder();
	sb.append("# motif: "+motif.getName()+" "+motif.version+" ");
	sb.append("threshold="+motifThreshold+" within "+motif_window+" bps\n");
	sb.append(GPSPeak.toGPS_short_Header()+"\n"); 
	for (GPSPeak p : nohitList){
		sb.append(p.toGPS_short()).append("\n");
	}
	CommonUtils.writeFile(outputFileName+"_noHit.txt", sb.toString());
	
	// overall stats
	String msg = String.format("\nWith motif:\t%d (%.1f%%)\nNO motif:\t%d (%.1f%%)\n\n", 
			hitList.size(), (double)hitList.size()/gpsPeaks.size()*100,
			nohitList.size(), (double)nohitList.size()/gpsPeaks.size()*100);
	statStr.append(msg);
	System.out.print(msg);
	msg = String.format("Average SR (spatial resolution) = %.2f bp, bias=%.1f\n"+
			"SR(quantile)  %d(25%%)  %d(50%%)  %d(75%%)\n",
			mean, bias, offsets[offsets.length/4], offsets[offsets.length/2], offsets[offsets.length*3/4]);
	statStr.append(msg);
	System.out.print(msg);
	return statStr.toString();
  }
  
  public void geneAnnotation(){
    boolean annotOverlapOnly=true;	
    
    ArrayList<AnnotationLoader> geneAnnotations = new ArrayList<AnnotationLoader>();
    geneAnnotations.add(new AnnotationLoader(genome, "GPS", "refGene", maxAnnotDistance, annotOverlapOnly));
    StringBuilder sb = new StringBuilder();
    for(AnnotationLoader loader : geneAnnotations){
      for (GPSPeak gpspeak : gpsPeaks) {
        Region coords = gpspeak.expand(1);
        Gene targetGene = null;
        int distToGene = Integer.MIN_VALUE;
                if (annotOverlapOnly) {
                    for(Gene gene : loader.getGenes(coords)){
                        int overlap = gene.getOverlapSize(coords);
                        if (targetGene == null || overlap > distToGene) {
                            targetGene = gene;
                            distToGene = overlap;
                        }
                    }
                    if (targetGene != null) {
                        distToGene = targetGene.distance(coords);
                    }
                } else {
                  distToGene = maxAnnotDistance;
                  Region query = coords.expand(maxAnnotDistance, maxAnnotDistance);
                  for(Gene gene : loader.getGenes(query)){
                        int distance = gpspeak.getLocation() - gene.getFivePrime();
                        if (gene.getStrand()=='-')
                          distance = -distance;
                        if (Math.abs(distance) < Math.abs(distToGene)) {
                            targetGene = gene;
                            distToGene = distance;                            
                        }                            
                    }
                }
                String geneName = targetGene == null ? "NONE" : targetGene.getName();
                String geneID = targetGene == null ? "NONE" : targetGene.getID();
                sb.append(gpspeak.toGPS()).append("\t").append(geneID).append("\t").append(geneName)
                  .append("\t").append(distToGene).append("\n");
            }
      String filename = outputFileName+(annotOverlapOnly?"_overlap_genes.txt":"_nearest_genes.txt");
      CommonUtils.writeFile(filename, sb.toString());  
    }
  }
  /*
   * Output sequences around the binding events
   * --win	window size
   * --top	for top # of events
   */
  public void printSequences(int win, int top){
	StringBuilder sb = new StringBuilder();
	StringBuilder sb2 = new StringBuilder();
	SequenceGenerator<Region> seqgen = new SequenceGenerator<Region>();
	Region peakWin=null;
	top = top==-1? gpsPeaks.size(): Math.min(top, gpsPeaks.size());
	for (int i=0;i<top;i++){
		GPSPeak gpspeak = gpsPeaks.get(i);
        peakWin = gpspeak.expand(win/2);
		sb.append(">"+peakWin.getLocationString()+"\t"+peakWin.getWidth() +"\t"+gpspeak.getStrength() +"\n");
		sb.append(seqgen.execute(peakWin)+"\n");
		peakWin = new Region(genome, peakWin.getChrom(), peakWin.getStart()+1000, peakWin.getEnd()+1000);
		sb2.append(">"+peakWin.getLocationString()+"\t"+peakWin.getWidth() +"\t0\n");
		sb2.append(seqgen.execute(peakWin)+"\n");
	}
	String filename = outputFileName+"_sequence.txt";
	CommonUtils.writeFile(filename, sb.toString());  
	filename = outputFileName+"_sequence_neg.txt";
	CommonUtils.writeFile(filename, sb2.toString());  
  }


  public void printMotifHitSequence(int top){
    WeightMatrixScorer scorer = new WeightMatrixScorer(motif);
    StringBuilder sb = new StringBuilder();
    top = top==-1? gpsPeaks.size(): Math.min(top, gpsPeaks.size());
    for (int i=0;i<top;i++){
	  GPSPeak gps = gpsPeaks.get(i);
      if ( gps.getShape()<-0.3){
          Region r= gps.expand(motif_window);
          String hit = scorer.getMaxScoreSequence(r, motifThreshold, extend);
          if (hit!=null)
        	  sb.append(hit).append("\n");
      }
    }
    CommonUtils.writeFile(outputFileName+"_motifHit.txt", sb.toString());
  }

  public void expressionIntegration(){
    // loading array annotation file
    String arrayFilename = Args.parseString(args, "ArrayProbeMap", null);
    if (arrayFilename==null){
      System.err.println("Array annotation file not found!");
      System.exit(0);
    }
    HashMap<String, ArrayAnnotation> probeMap = new HashMap<String, ArrayAnnotation>();
    try{
      File rFile = new File(arrayFilename);
      if(!rFile.isFile()){System.err.println("Invalid file name");System.exit(1);}
          BufferedReader reader = new BufferedReader(new FileReader(rFile));
          String line;
          while ((line = reader.readLine()) != null) {
              line = line.trim();
              String[] words = line.split("\\t");
              ArrayAnnotation array = new ArrayAnnotation();
              array.probe = words[0];
              array.geneName = words[1];
              array.refSeqId = words[2];
              array.region = new StrandedRegion(genome, words[3], Integer.parseInt(words[4]), 
                  Integer.parseInt(words[5]), words[6].charAt(0));
              if (array.region!=null)
                probeMap.put(array.probe, array);
            }
          reader.close();
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    } catch (IOException e) {
      e.printStackTrace();
    }
    
    // loading differential expression file
    String exprFilename = Args.parseString(args, "DiffExpression", null);
    if (exprFilename==null){
      System.err.println("Differential expression file not found!");
      System.exit(0);
    }
    ArrayList<DiffExpression> diffExpressions = new ArrayList<DiffExpression>();
    try{
      File rFile = new File(exprFilename);
      if(!rFile.isFile()){System.err.println("Invalid file name");System.exit(1);}
          BufferedReader reader = new BufferedReader(new FileReader(rFile));
          String line= reader.readLine(); //skip first header line
          while ((line = reader.readLine()) != null) {
              line = line.trim();
              String[] words = line.split("\\t");
              DiffExpression expr = new DiffExpression();
              expr.probe = words[0];
              expr.geneName = words[1];
              expr.foldChange = Double.parseDouble(words[2]);
              expr.pValue = Double.parseDouble(words[3]);
              diffExpressions.add(expr);
            }
          reader.close();
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    } catch (IOException e) {
      e.printStackTrace();
    }
    
    // for each differentially expressed genes
    HashMap<String, ArrayList<GPSPeak>> gene2peaks = new HashMap<String, ArrayList<GPSPeak>>();
    for (DiffExpression expr: diffExpressions){
      ArrayAnnotation array = probeMap.get(expr.probe);
      if (array==null)
        continue;
      Region probeRegion = array.region;
      ArrayList<GPSPeak> peaks = new ArrayList<GPSPeak>();
      for (GPSPeak gpspeak : gpsPeaks) {
        Region gpsRegion = gpspeak.expand(maxAnnotDistance);
        if (gpsRegion.overlaps(probeRegion))
          peaks.add(gpspeak);
      }
      gene2peaks.put(expr.probe, peaks);
    }
    
    // output
    StringBuilder sb = new StringBuilder();
    sb.append("Probe\tGene\tFoldChange\tProbeRegion\tBindingSites\n");
    for (DiffExpression expr: diffExpressions){
      ArrayList<GPSPeak> peaks = gene2peaks.get(expr.probe);
      if (peaks==null || peaks.size()==0)
        continue;
      sb.append(expr.probe+"\t").append(expr.geneName+"\t").append(String.format("%.2f\t", expr.foldChange));
      sb.append(probeMap.get(expr.probe).region.toString()).append("\t");
      for (GPSPeak peak: peaks){
        sb.append(peak.toString()).append("|");
      }
      sb.append("\n");
    }
    CommonUtils.writeFile(outputFileName+"_GPS_DiffExpr.txt", sb.toString()); 
  }

}
  class DiffExpression {
    String probe;
    String geneName;
    double foldChange;
    double pValue;
  }
  class ArrayAnnotation{
    String probe;
    String geneName;
    String refSeqId;
    StrandedRegion region;
  }

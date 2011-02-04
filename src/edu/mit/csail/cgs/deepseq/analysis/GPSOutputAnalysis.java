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
import edu.mit.csail.cgs.deepseq.features.ComponentFeature;
import edu.mit.csail.cgs.deepseq.features.Feature;
import edu.mit.csail.cgs.deepseq.utilities.AnnotationLoader;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSPeak;
import edu.mit.csail.cgs.ewok.verbs.motifs.Kmer;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public class GPSOutputAnalysis {
  static final int NOHIT_OFFSET = 999;
  static final int maxAnnotDistance=50000;
  
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
	case 5: analysis.kmerGrouping();break;
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
  
  public void kmerGrouping(){
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
	Collections.sort(kmers);
	
	extendSeeds(kmers);
	// generate motifs from groups
	
  }
  private void extendSeeds(ArrayList<Kmer> kmers){
	char[] letters = {'A','C','T','G'};
	int MISS = 1;
	// group kmers by greedy extending method
	Kmer root = kmers.get(0);		
	ArrayList<TreeSet<Kmer>> alignedKmerSets = new ArrayList<TreeSet<Kmer>>();
	TreeSet<Kmer> prevSelected = new TreeSet<Kmer>();
	prevSelected.add(root);
	TreeSet<Kmer> alignedKmerSet = new TreeSet<Kmer>();	// all selected Kmers (aligned (1bp mismatch) or extended (1bp offset))	
	while(!kmers.isEmpty()){
		while (!prevSelected.isEmpty()){					
			TreeSet<Kmer> selectedKmers = new TreeSet<Kmer>();
			
			// select all kmers aligned to or extend from the prevSelected set	
			for(Kmer km: kmers){
				for (Kmer sk: prevSelected){
					km.setReference(sk);	
					if (km.getScore()>=sk.getK()-MISS){		// select only if 1 mismatch or 1 extension
						selectedKmers.add(km);
						break;
					}
				}
			}
			kmers.removeAll(selectedKmers);
			prevSelected = selectedKmers;
			alignedKmerSet.addAll(selectedKmers);
		}
		// no more extension, start a new kmer set
		alignedKmerSets.add(alignedKmerSet);
		if (!kmers.isEmpty()){
			Kmer next = kmers.get(0);
			next.setReference(root);					// adjust orientation of this kmer to better align with root
			prevSelected.add(next);						// add the seed for next group of extension, start with the next largest count kmer
			alignedKmerSet = new TreeSet<Kmer>();		// new set
			alignedKmerSet.addAll(prevSelected);
		}
		
	}
	
	// trace back to set the offset wrt the firstSeed kmer, to aligned globally
	int min = 0;
	int max = 0;
	for (TreeSet<Kmer> ks:alignedKmerSets){
		for (Kmer k: ks){
			Kmer thisKmer = k;
			int offset = thisKmer.getShift();
			while(thisKmer.getRef()!=thisKmer){						// recursively trace back the reference kmer
				offset += thisKmer.getRef().getShift();	
				thisKmer = thisKmer.getRef();
			}
			if (min>offset)	
				min = offset;
			if (max<offset)
				max = offset;
			k.setGlobalShift(offset);
		}
	}
	
	// frequency matrix
	double[][][] pfms = new double[alignedKmerSets.size()][max-min+1][4];
	int k = root.getK();
	for (int i=0;i<alignedKmerSets.size();i++){		
		double[][] pfm = new double[max-min+1][4];
		TreeSet<Kmer> ks = alignedKmerSets.get(i);
		pfms[i] = pfm;
		for (Kmer km: ks){
			int shift = km.getGlobalShift();
			double factor = km.getSeqHitCount();
			if (useWeight)
				factor = km.getWeight();
			for (int p=0;p<max-shift;p++){
				for (int b=0;b<letters.length;b++){
					pfm[p][b] +=0.25*factor;
				}
			}
			for (int p=max-shift;p<max-shift+k;p++){
				char base = km.getKmerString().charAt(p);
				for (int b=0;b<letters.length;b++){
					if (base==letters[b])
						pfm[p][b] +=1*factor;
				}
			}
			for (int p=max-shift+k;p<max-min+1;p++){
				for (int b=0;b<letters.length;b++){
					pfm[p][b] +=0.25*factor;
				}
			}
		}
	}
	// print out motif result
	StringBuilder msb = new StringBuilder();
	for (int i=0;i<alignedKmerSets.size();i++){	
		msb.append("DE "+motif.name+"_"+(i+1)).append("\n");
		double[][] pfm = pfms[i] ;
		for (int p=0;p<pfm.length;p++){
			msb.append(p+1).append(" ");
			int maxBase = 0;
			double maxFreq=0;
			for (int b=0;b<letters.length;b++){
				msb.append(String.format("%.2f ", pfm[p][b]));
				if (maxFreq<pfm[p][b]){
					maxFreq=pfm[p][b];
					maxBase = b;
				}
			}
			msb.append(letters[maxBase]).append("\n");
		}
		msb.append("XX\n");
	}
	CommonUtils.writeFile("Kmer_PFM.txt", msb.toString());	
	
	//print out
	StringBuilder sb = new StringBuilder();
	sb.append("Kmer").append("\t").append("gShift").append("\t").append("refKmer").append("\t")
	  .append("seqCt").append("\t").append("Alignment").append("\n");
	for (TreeSet<Kmer> ks:alignedKmerSets){
		for (Kmer km: ks){
			String shiftedKmer = CommonUtils.padding(max-km.getGlobalShift(), ' ').concat(km.getKmerString());
			sb.append(km.getKmerString()).append("\t").append(km.getGlobalShift()).append("\t").append(km.getRef().getKmerString()).append("\t")
			  .append(km.getSeqHitCount()).append("\t").append(shiftedKmer).append("\n");			
		}
		sb.append("\n");
	}
	sb.append("Total kmer extension set: " + alignedKmerSets.size());
	CommonUtils.writeFile("Extended_Kmers.txt", sb.toString());
  }

  private void matchSeeds(ArrayList<Kmer> kmers){
	// group kmers by greedy growing method
	ArrayList<KmerGroup> groups = new ArrayList<KmerGroup>();
	Kmer firstSeed = kmers.get(0);		
	while(!kmers.isEmpty()){
		Kmer topKmer = kmers.get(0);
		topKmer.setReference(firstSeed);
		KmerGroup group = new KmerGroup(topKmer);
		groups.add(group);
		kmers.remove(topKmer);
		HashSet<Kmer> selected = new HashSet<Kmer>();
		for(Kmer km: kmers){
			km.setReference(topKmer);
			if (km.getScore()>topKmer.getK()*0.8){
				group.addMember(km);
				selected.add(km);
			}
		}
		kmers.removeAll(selected);
	}
	
	// print out kmer groups
	StringBuilder sb = new StringBuilder();
	for (KmerGroup g:groups){
		sb.append(g.seed.getKmerString()).append("\t").append(g.seed.getSeqHitCount()).append("\n");
		for (Kmer km: g.members){
			int shift = km.getShift();
			String str = km.getKmerString();
			if (shift==-99)
				continue;
			else if (shift>=0){
				str = str.substring(shift);
			}
			else{
				str = (shift==-1?" ":"  ")+str.substring(0, str.length()+shift);
			}
			sb.append(str).append("\t").append(km.getSeqHitCount()).append("\n");
		}
		sb.append("\n");
	}
	sb.append("Total kmer groups: " + groups.size());
	CommonUtils.writeFile("Aligned_Kmers.txt", sb.toString());
		
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

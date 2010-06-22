package edu.mit.csail.cgs.deepseq.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.discovery.BindingMixture;
import edu.mit.csail.cgs.deepseq.utilities.AnnotationLoader;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSPeak;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public class GPSOutputAnalysis {
  static final int MOTIF_DISTANCE = 200;  
  static final int maxAnnotDistance=50000;
  
  private Genome genome;
  private Organism org;
  private String[] args;
  private WeightMatrix motif = null;
  double motifThreshold;
  
  List<GPSPeak> gpsPeaks;
  String GPSfileName;
  
  // build empirical distribution
  String chipSeqExpt = null;
  String chipSeqVersion = null;
  boolean useMotif = true;
  
  /**
   * @param args
   */
  public static void main(String[] args) {
    //String ctcf_results = "/home/rca/psrg/projects/GPS/GPS_results/CTCF_Chen08_GPS.peaks.txt";
    String ctcf_results = "/home/rca/psrg/projects/GPS/GPS_results/Sing_CTCF_ES_1_signal.peaks.txt";
    String[] args_ctcf = new String[] {"--motif", "CTCF",
        "--version", "090828",
        "--threshold", "16.08",
        "--GPS", ctcf_results        
    };
    
    //Nanog
//    String nanog_results = "/home/rca/psrg/projects/GPS/GPS_results/Nanog_rep1_YL_0_signal.peaks.txt";
//    String nanog_results = "/home/rca/psrg/projects/GPS/GPS_results/Nanog_YL_2_signal.peaks.txt";
    String nanog_results = "/home/rca/psrg/projects/GPS/GPS_results/Nanog_Chen08_3_signal.peaks.txt";
    
    String[] args_nanog_osnt = new String[] {"--motif", "OSNT",
        "--version", "YoungLab",
        "--threshold", "12.72",
        "--GPS", nanog_results        
    };
    String[] args_nanog_ji = new String[] {"--motif", "NanogSox2",
        "--version", "Ji2008_NanogSox2",
        "--threshold", "8.33",
        "--GPS", nanog_results        
    };
    String[] args_nanog_M01123 = new String[] {"--motif", "Nanog",
        "--version", "TRANSFAC 10.4, M01123",
        "--threshold", "9.13",
        "--GPS", nanog_results        
    };
    String[] args_nanog_M01247 = new String[] {"--motif", "Nanog",
        "--version", "TRANSFAC 2009.4, M01247",
        "--threshold", "9.095",
        "--GPS", nanog_results        
    };
    String[] args_nanog_he = new String[] {"--motif", "Nanog",
        "--version", "He2009",
        "--threshold", "7.73",
        "--GPS", nanog_results        
    };
    String[] args_nanog_mitsui = new String[] {"--motif", "Nanog",
        "--version", "Mitsui2003",
        "--threshold", "7.82",
        "--GPS", nanog_results        
    };
    String[] args_nanog_macisaac_0 = new String[] {"--motif", "Nanog",
        "--version", "MacIsaac2006_Nanog0",
        "--threshold", "10.97",
        "--GPS", nanog_results        
    };
    String[] args_nanog_macisaac_1 = new String[] {"--motif", "Nanog",
        "--version", "MacIsaac2006_Nanog1",
        "--threshold", "11.88",
        "--GPS", nanog_results        
    };
    String[] args_nanog_macisaac_2 = new String[] {"--motif", "Nanog",
        "--version", "MacIsaac2006_Nanog2",
        "--threshold", "11.02",
        "--GPS", nanog_results        
    };
    String[] args_nanog_macisaac_3 = new String[] {"--motif", "Nanog",
        "--version", "MacIsaac2006_Nanog3",
        "--threshold", "11.80",
        "--GPS", nanog_results        
    };
    String[] args_nanog_macisaac_4 = new String[] {"--motif", "Nanog",
        "--version", "MacIsaac2006_Nanog4",
        "--threshold", "11.45",
        "--GPS", nanog_results        
    };
    
    
    //Oct4
//    String oct4_results = "/home/rca/psrg/projects/GPS/GPS_results/Oct4_YL_0_signal.peaks.txt";
//    String oct4_results = "/home/rca/psrg/projects/GPS/GPS_results/YL_Oct4_ES_1_signal.peaks.txt";
    String oct4_results = "/home/rca/psrg/projects/GPS/GPS_results/Sing_Oct4_ES_2_signal.peaks.txt";
    String[] args_oct4_pou5f1 = new String[] {"--motif", "Oct-4 (POU5F1)",
        "--version", "TRANSFAC 10.4, M01124",
        "--threshold", "12.82",
        "--GPS", oct4_results        
    };
    String[] args_oct4_osnt = new String[] {"--motif", "OSNT",
        "--version", "YoungLab",
        "--threshold", "12.72",
        "--GPS", oct4_results        
    };
    
    //Sox2
//    String sox2_results = "/home/rca/psrg/projects/GPS/GPS_results/Sox2_rep1_YL_0_signal.peaks.txt";
//    String sox2_results = "/home/rca/psrg/projects/GPS/GPS_results/YL_Sox2_ES_2_signal.peaks.txt";
    String sox2_results = "/home/rca/psrg/projects/GPS/GPS_results/Sing_Sox2_ES_2_signal.peaks.txt";
    
    String[] args_sox2_osnt = new String[] {"--motif", "OSNT",
        "--version", "YoungLab",
        "--threshold", "12.72",
        "--GPS", sox2_results        
    };
  
  
   // args = args_nanog_macisaac_2;
    
    GPSOutputAnalysis analysis = new GPSOutputAnalysis(args);
//    analysis.buildEmpiricalDistribution();
    analysis.jointBindingMotifAnalysis(true);
//    analysis.geneAnnotation();
//    analysis.expressionIntegration();
  }
  
  public GPSOutputAnalysis(String[] args) {
    this.args = args;
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
//      genome = Organism.findGenome("mm8");
//      org = Organism.getOrganism("Mus musculus");     
    } catch (NotFoundException e) {
      e.printStackTrace();
    }
    // load motif
    try {
      String motifString = Args.parseString(args, "motif", null);
      String motifVersion = Args.parseString(args, "version", null);
      motifThreshold = Args.parseDouble(args, "threshold", 10.0);
      if (motifThreshold==10.0)
        System.err.println("No motif threshold was provided, default=10.0 is used.");
      int wmid = WeightMatrix.getWeightMatrixID(org.getDBID(), motifString, motifVersion);
      motif = WeightMatrix.getWeightMatrix(wmid);
    } 
    catch (NotFoundException e) {
      e.printStackTrace();
    }   
    
    // load GPS results
    GPSfileName = Args.parseString(args, "GPS", null);
    if (GPSfileName==null){
      System.err.println("GPS file not found!");
      System.exit(0);
    }
    File gpsFile = new File(GPSfileName);
    gpsPeaks = GPSParser.parseGPSOutput(gpsFile.getAbsolutePath(), genome);
    
    // parameter for building empirical distribution
//    chipSeqExpt = Args.parseString(args, "chipSeqExpt", null);
//    chipSeqVersion = Args.parseString(args, "chipSeqVersion", null);
//    useMotif = Args.parseInteger(args, "useMotif", 0)==1?true:false;
  }
  
  public void buildEmpiricalDistribution(){
    WeightMatrixScorer scorer=null;
    if (useMotif){
      scorer = new WeightMatrixScorer(motif);
    }
    ArrayList<Point> points = new ArrayList<Point>();
    for (GPSPeak gps: gpsPeaks){
//      if (gps.getMixProb()==1 && gps.getQvalue()>5 && gps.getShape()<1)
      if (gps.getMixProb()==1 && gps.getStrength()>40 && gps.getShape()<1){
        if (useMotif){
          Region r= gps.expand(MOTIF_DISTANCE);
          WeightMatrixScoreProfile profiler = scorer.execute(r);
          int halfWidth = profiler.getMatrix().length()/2;
          //search from BS outwards, to find the nearest strong motif
          for(int z=0; z<=r.getWidth()/2; z++){
            double leftScore= profiler.getMaxScore(MOTIF_DISTANCE-z);
            double rightScore= profiler.getMaxScore(MOTIF_DISTANCE+z);  
            int position = -Integer.MAX_VALUE;  // middle of motif, relative to GPS peak
            if(rightScore>=motifThreshold){
              position = z+halfWidth;
            }
            if(leftScore>=motifThreshold){
              position = -z+halfWidth;
            }
            if(position > -Integer.MAX_VALUE){
              Point motifPos = new Point(genome, gps.getChrom(), gps.getLocation()+position);
              points.add(motifPos);
              break;  // break from the motif search of this peak
            }
          }
        }
        else
          points.add((Point)gps);
      }
    }
    ChipSeqStats.printEmpiricalDistribution(points, chipSeqExpt, chipSeqVersion );
    System.out.println(points.size()+" GPS peaks are used to build empricaldistribution.");
  }

  
  /**
   * For each peak find the occurrence of a motif match >= the specified 
   * threshold that is closest to the peak
   * @param analyzeMotif
   */
  public void jointBindingMotifAnalysis(boolean analyzeMotif){
    System.out.println("Analyzing motif matches");
    Collections.sort(gpsPeaks);
    
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

    // motif
    int[]motif_distances = new int[gpsPeaks.size()];
    double[]motif_scores = new double[gpsPeaks.size()];
    if (analyzeMotif){
      WeightMatrixScorer scorer = new WeightMatrixScorer(motif);
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
        Region r= peak.expand(MOTIF_DISTANCE);
        WeightMatrixScoreProfile profiler = scorer.execute(r);
        int halfWidth = profiler.getMatrix().length()/2;
        //search from BS outwards
        for(int z=0; z<=MOTIF_DISTANCE; z++){
          double leftScore = Double.NEGATIVE_INFINITY;
          double rightScore = Double.NEGATIVE_INFINITY;
          if ((MOTIF_DISTANCE+z) < profiler.length()) {
            rightScore= profiler.getMaxScore(MOTIF_DISTANCE+z);        
            if(rightScore>=motifThreshold){
              motif_distances[i] = z+halfWidth;
              motif_scores[i] = rightScore;
              break;
            }
          }
          if ((MOTIF_DISTANCE-z) >= 0) {
            leftScore= profiler.getMaxScore(MOTIF_DISTANCE-z);
            if(leftScore>=motifThreshold){
              motif_distances[i] = -z+halfWidth;
              motif_scores[i] = leftScore;
              break;
            }
          }
          // if motif score at this position is too small, and reach end of region
          if (z==r.getWidth()/2){
            motif_distances[i] = 999;
            motif_scores[i] = Math.max(leftScore, rightScore);
          }
        }
      }
      System.out.println();
    }
     
    //sort by inter-event distance
    int old_index[] = StatUtil.findSort(interPeak_distances);
    
    // output results
    StringBuilder sb = new StringBuilder();
    StringBuilder sb_binary = new StringBuilder();
    // header
    sb.append("# motif: "+motif.getName()+" "+motif.version+" ");
    sb.append("threshold="+motifThreshold+" within "+MOTIF_DISTANCE+" bps\n");
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
      sb.append(motif_distances[i]).append("\t");
      sb.append(motif_distances[i-1]).append("\t");
      sb.append(interPeak_distances[j]).append("\n");     
      if (motif_distances[i-1]-motif_distances[i]!=interPeak_distances[j] 
          && motif_distances[i]!=999
          && motif_distances[i-1]!=999
          && interPeak_distances[j]<=1000){
    	  sb_binary.append(gpsPeaks.get(i).toGPS_short()).append("\t");
    	  sb_binary.append(String.format("%.2f", motif_scores[i])).append("\t");
    	  sb_binary.append(motif_distances[i]).append("\t");
    	  sb_binary.append(motif_distances[i-1]).append("\t");
    	  sb_binary.append(interPeak_distances[j]).append("\n");     
      }
    }
    // first peak
	//  sb.append(gpsPeaks.get(0).toGPS_motifShort()).append("\t");
	  sb.append(gpsPeaks.get(0).toGPS_short()).append("\t");
	  sb.append(String.format("%.2f", motif_scores[0])).append("\t");
	  sb.append(motif_distances[0]).append("\t");
	  sb.append(999).append("\t"); //prev motif distance
	  sb.append(99999999).append("\n");  //interpeak distance

  
    BindingMixture.writeFile(GPSfileName+"_JointPeak_Motif.txt", sb.toString());
    BindingMixture.writeFile(GPSfileName+"_BinaryMotifEvents.txt", sb_binary.toString());
  }
  
  private void jointEvents(ArrayList<Point> points, int jointCutoff){
	  Collections.sort(points);
	  Point previous = points.get(0);
	  ArrayList<Point> cluster = new ArrayList<Point>();
	  ArrayList<ArrayList<Point>> clusters = new ArrayList<ArrayList<Point>>();
	  clusters.add(cluster);
	  boolean isJoint = false; 	// Previous point is joint event?
	  for (int i=1;i<points.size();i++){
		  Point p = points.get(i);
		  if (p.getChrom().equals(previous.getChrom()) && p.distance(previous)<=jointCutoff){
			  if (isJoint){
				  cluster.add(p);
			  }
			  else{
				  cluster = new ArrayList<Point>();
				  clusters.add(cluster);
				  cluster.add(previous);
				  cluster.add(p);
				  isJoint = true;
			  }
		  }
		  else{
			  isJoint = false;
		  }
	  }
	  StringBuilder sb = new StringBuilder();
	  for (ArrayList<Point> c:clusters){
		  Point p = c.get(0);
		  int start = p.getLocation();
		  int end = c.get(c.size()-1).getLocation();
		  Region r = new Region(p.getGenome(), p.getChrom(), start, end);
		  sb.append(r.toString()).append("\t").append(c.size()).append("\n");
	  }
	  System.out.println(sb.toString());
  }
  
  private void geneAnnotation(){
    boolean annotOverlapOnly=false;	
    
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
      String filename = GPSfileName+(annotOverlapOnly?"_overlap_genes.txt":"_nearest_genes.txt");
      BindingMixture.writeFile(filename, sb.toString());  
    }
  }
  
  private void expressionIntegration(){
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
    BindingMixture.writeFile(GPSfileName+"_GPS_DiffExpr.txt", sb.toString()); 
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

package edu.mit.csail.cgs.deepseq.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
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
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils.SISSRS_Event;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSPeak;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public class JointEventAnalysis {
  static final int JOINT_DISTANCE = 500;
  static final int maxAnnotDistance=50000;
  
  private Genome genome;
  private Organism org;
  private String[] args;
  private WeightMatrix motif = null;
  private double motifThreshold;
  private int motifWindowSize = 20;  
  private int topRank;
  
  ArrayList<Point> events;
  HashMap<Point, Double> eventStrength;
  
  /**
   * @param args
   */
  public static void main(String[] args) throws IOException {
	  JointEventAnalysis analysis = new JointEventAnalysis(args);
	  analysis.jointEvents(JOINT_DISTANCE);
  }
  
  public JointEventAnalysis(String[] args) throws IOException {
    this.args = args;
    ArgParser ap = new ArgParser(args);
    
    genome = CommonUtils.parseGenome(args);
    
    motifWindowSize = Args.parseInteger(args, "windowSize", 50);
	topRank = Args.parseInteger(args, "rank", 0);
	
    // load motif
	Pair<WeightMatrix, Double> wm = CommonUtils.loadPWM(args, org.getDBID());
	if (wm!=null){
		motif = wm.car();
		motifThreshold = wm.cdr(); 
	}
    // load GPS results
    String GPSfileName = Args.parseString(args, "GPS", null);
    String SISSRsFileName = Args.parseString(args, "SISSRs", null);
    if (GPSfileName!=null){
	    File gpsFile = new File(GPSfileName);
	    List<GPSPeak> gpsPeaks = GPSParser.parseGPSOutput(gpsFile.getAbsolutePath(), genome);
	    if (topRank==0)
	    	topRank = gpsPeaks.size();
	    events = new ArrayList<Point>();
	    eventStrength = new HashMap<Point, Double>();
	    int count=0;
	    for (GPSPeak g: gpsPeaks){
	    	events.add(g);
	    	eventStrength.put(g, g.getStrength());
			count++;
			if (count>=topRank)
				break;
	    }
    }
    else if (SISSRsFileName!=null){
    	ArrayList<SISSRS_Event> SISSRs_events = CommonUtils.load_SISSRs_events(genome, SISSRsFileName, true);
    	if (topRank==0)
	    	topRank = SISSRs_events.size();
	    events = new ArrayList<Point>();
    	eventStrength = new HashMap<Point, Double>();
	    int count=0;
	    for (SISSRS_Event se: SISSRs_events){
	    	events.add(se.getPeak());
	    	eventStrength.put(se.getPeak(), se.tags);
			count++;
			if (count>=topRank)
				break;
	    }
    }
    else{
    	events = CommonUtils.loadCgsPointFile(Args.parseString(args, "events", null), genome);
    }
  }
  
	
  private void jointEvents(int jointCutoff){
	  if (events.size()==0){
		  System.err.println("No binding events!!");
		  System.exit(0);
	  }
	  
	  Collections.sort(events);
	  
	  Point previous = events.get(0);
	  ArrayList<Point> cluster = null;
	  ArrayList<ArrayList<Point>> clusters = new ArrayList<ArrayList<Point>>();
	  ArrayList<Double> clusterStrength = new ArrayList<Double>();
	  boolean isJoint = false; 	// Previous point is joint event?
	  for (int i=1;i<events.size();i++){
		  Point p = events.get(i);
		  if (p.getChrom().equals(previous.getChrom()) && p.distance(previous)<=jointCutoff){
			  if (isJoint){				// if in the joint event state, keep adding event
				  cluster.add(p);
				  int clusterIdx = clusterStrength.size()-1;
				  clusterStrength.set(clusterIdx, clusterStrength.get(clusterIdx)+eventStrength.get(p));
			  }
			  else{						// if it is a new joint event, initialized the cluster, add to cluster list
				  cluster = new ArrayList<Point>();
				  clusters.add(cluster);
				  cluster.add(previous);
				  cluster.add(p);
				  clusterStrength.add(eventStrength.get(previous) + eventStrength.get(p));
				  isJoint = true;
			  }
		  }
		  else{
			  isJoint = false;
		  }
		  previous = p;
	  }
	  
	  // sort by total read count in the cluster
	  double[] readCounts = new double[clusters.size()];
	  for (int i=0;i<clusterStrength.size();i++){
		  readCounts[i]=clusterStrength.get(i);
	  }
	  int[] sorted = StatUtil.findSort(readCounts);	                                   
	  
	  StringBuilder sb = new StringBuilder();
	  System.out.println(String.format("Total %d joint event clusters from top %d events.", clusters.size(), topRank));
	  if (motif!=null){
		  sb.append("MidPoint\t      Region       \tEvents\tMotifs\tDistances\tMinDistance");
		  if (eventStrength!=null)
			  sb.append("\tTotalReadCount");
		  sb.append("\n");
	  }
	  else{
		  sb.append("MidPoint\t      Region       \tEvents\tDistances\tMinDistance");
		  if (eventStrength!=null)
			  sb.append("\tTotalReadCount");
		  sb.append("\n");
	  }
	  int eventTotal=0;
	  int motifTotal=0;
	  for (int i=sorted.length-1;i>=0;i--){
		  ArrayList<Point> c = clusters.get(sorted[i]);
		  if (c.size()==0)
			  continue;
		  eventTotal += c.size();
		  Point p = c.get(0);
		  int start = p.getLocation();
		  int end = c.get(c.size()-1).getLocation();
		  Region r = new Region(p.getGenome(), p.getChrom(), start, end);
		  sb.append(r.getMidpoint().toString()).append("\t")
		    .append(r.toString()).append("\t")
		    .append(c.size()).append("\t");
		  if (motif!=null){
			  // if given motif PWM, get motif matches within window size from each event position, keep unique motif only
			  HashSet<Point> uniqueMotifs = new HashSet<Point>();
			  for (Point pp: c){
				  Pair<ArrayList<Point>, ArrayList<Double>> motifs = getAllMotifs(pp.expand(motifWindowSize), motifThreshold);
				  uniqueMotifs.addAll(motifs.car());
			  }
			  sb.append(uniqueMotifs.size()).append("\t");
			  motifTotal += uniqueMotifs.size();
		  }
		  int min = Integer.MAX_VALUE;
		  for (int j=0;j<c.size()-1;j++){
			  int distance = c.get(j+1).distance(c.get(j));
			  sb.append(distance).append(" ");
			  if (distance < min)
				  min = distance;
		  }
		  sb.append("\t").append(min);
		  if (eventStrength!=null)
			  sb.append("\t"+clusterStrength.get(sorted[i]));
		  sb.append("\n");
	  }
	  System.out.println(String.format("Motifs/Events=%d/%d=%.2f", motifTotal, eventTotal, ((double)motifTotal)/eventTotal));
	  String outName = Args.parseString(args, "out", "NoName");
	  CommonUtils.writeFile(String.format("%s_%.2f_%d_jointEvents.txt", outName, motifThreshold, motifWindowSize), sb.toString());
  }
  
	// get the all motif hits in the region (windowSize) 
	// around each peak in the list
	// if a motif have more than one peaks nearby, associate it with the nearest peak
	private Pair<ArrayList<Point>, ArrayList<Double>> getAllMotifs(Region region, double motifScoreThreshold){
		ArrayList<Point> allMotifs = new ArrayList<Point>();
		ArrayList<Double> scores = new ArrayList<Double>();
		WeightMatrixScorer scorer = new WeightMatrixScorer(motif);
		WeightMatrixScoreProfile profiler = scorer.execute(region);
		//search for whole region
		for(int z=0; z<region.getWidth(); z++){		
			double score = profiler.getHigherScore(z);
			if(score >= motifScoreThreshold){
				Point motifPos = new Point(genome, region.getChrom(), region.getStart()+z+profiler.getMatrix().length()/2);
				allMotifs.add(motifPos);
				scores.add(score * (profiler.getHigherScoreStrand(z)=='+'?1:-1)); // positive score if motif match on '+' strand
			}
		}
		return new Pair<ArrayList<Point>, ArrayList<Double>>(allMotifs, scores);
	}
	
}
 

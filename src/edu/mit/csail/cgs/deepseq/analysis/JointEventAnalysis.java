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
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
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
  static final int MOTIF_DISTANCE = 200;  
  static final int JOINT_DISTANCE = 500;
  static final int maxAnnotDistance=50000;
  
  private Genome genome;
  private Organism org;
  private String[] args;
  private WeightMatrix motif = null;
  private double motifThreshold;
  private int topRank;
  
  ArrayList<Point> events;
  
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
    
    topRank = Args.parseInteger(args, "top", -1);
    
    // load motif
	Pair<WeightMatrix, Double> wm = CommonUtils.loadPWM(args, org.getDBID());
	if (wm!=null){
		motif = wm.car();
		motifThreshold = wm.cdr(); 
	}
    // load GPS results
    String GPSfileName = Args.parseString(args, "GPS", null);
    if (GPSfileName!=null){
	    File gpsFile = new File(GPSfileName);
	    List<GPSPeak> gpsPeaks = GPSParser.parseGPSOutput(gpsFile.getAbsolutePath(), genome);
	    if (topRank==-1)
	    	topRank = gpsPeaks.size();
	    events = new ArrayList<Point>();
	    int count=0;
	    for (GPSPeak g: gpsPeaks){
	    	events.add(g);
			count++;
			if (count>=topRank)
				break;
	    }
    }
    else{
    	events = loadCgsPointFile(Args.parseString(args, "events", null));
    }
  }
  
	/** 
	 * load text file in CGS Point format
	 * @param filename
	 * chr:coord, e.g. 1:234234
	 */
	private ArrayList<Point> loadCgsPointFile(String filename) {

		File file = new File(filename);
		FileReader in = null;
		BufferedReader bin = null;
		ArrayList<Point> points = new ArrayList<Point>();
		int count = 0;
		try {
			in = new FileReader(file);
			bin = new BufferedReader(in);
			String line;
			while((line = bin.readLine()) != null) { 
				line = line.trim();
				Region point = Region.fromString(genome, line);
				points.add(new Point(genome, point.getChrom(),point.getStart()));
				count++;
				if (count>=topRank)
					break;
			}
		}
		catch(IOException ioex) {
			//logger.error("Error parsing file", ioex);
		}
		finally {
			try {
				if (bin != null) {
					bin.close();
				}
			}
			catch(IOException ioex2) {
				//nothing left to do here, just log the error
				//logger.error("Error closing buffered reader", ioex2);
			}			
		}
		return points;
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
	  boolean isJoint = false; 	// Previous point is joint event?
	  for (int i=1;i<events.size();i++){
		  Point p = events.get(i);
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
		  previous = p;
	  }
	  
	  StringBuilder sb = new StringBuilder();
	  System.out.println("Total "+clusters.size()+" joint event clusters.");
	  if (motif!=null)
		  sb.append("      Region       \tEvents\tMotifs\tDistances\tMinDistance\n");
	  else
		  sb.append("      Region       \tEvents\tDistances\tMinDistance\n");
	  for (ArrayList<Point> c:clusters){
		  if (c.size()==0)
			  continue;
		  Point p = c.get(0);
		  int start = p.getLocation();
		  int end = c.get(c.size()-1).getLocation();
		  Region r = new Region(p.getGenome(), p.getChrom(), start, end);
		  sb.append(r.toString()).append("\t")
		    .append(c.size()).append("\t");
		  if (motif!=null){
			  Pair<ArrayList<Point>, ArrayList<Double>> motifs = getAllMotifs(r.expand(MOTIF_DISTANCE, MOTIF_DISTANCE), motifThreshold);
			  sb.append(motifs.car().size()).append("\t");
		  }
		  int min = Integer.MAX_VALUE;
		  for (int i=0;i<c.size()-1;i++){
			  int distance = c.get(i+1).distance(c.get(i));
			  sb.append(distance).append(" ");
			  if (distance < min)
				  min = distance;
		  }
		  sb.append("\t").append(min).append("\n");
	  }
//	  System.out.println(sb.toString());
	  String outName = Args.parseString(args, "out", "");
	  CommonUtils.writeFile(outName+"_jointEvents.txt", sb.toString());
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
			double score = profiler.getMaxScore(z);
			if(score >= motifScoreThreshold){
				Point motifPos = new Point(genome, region.getChrom(), region.getStart()+z+profiler.getMatrix().length()/2);
				allMotifs.add(motifPos);
				scores.add(score * (profiler.getMaxStrand(z)=='+'?1:-1)); // positive score if motif match on '+' strand
			}
		}
		return new Pair<ArrayList<Point>, ArrayList<Double>>(allMotifs, scores);
	}
	
}
 

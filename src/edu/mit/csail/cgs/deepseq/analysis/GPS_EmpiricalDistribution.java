package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.BindingModel;
import edu.mit.csail.cgs.deepseq.DeepSeqExpt;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSPeak;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class GPS_EmpiricalDistribution {
	static final int MOTIF_DISTANCE = 50;	
	
	private Genome genome;
	private Organism org;
	private WeightMatrix motif = null;
	private double motifThreshold;
	
	private List<GPSPeak> gpsPeaks;
	private String GPSfileName;
	private int strength=40;
	// build empirical distribution
	private DeepSeqExpt chipSeqExpt = null;
	private int range = 250;
	private boolean useMotif = true;
	private String name = null;
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		GPS_EmpiricalDistribution analysis = new GPS_EmpiricalDistribution(args);
		analysis.buildEmpiricalDistribution();
	}
	
	public GPS_EmpiricalDistribution(String[] args) {
		ArgParser ap = new ArgParser(args);
		
		try {
			Pair<Organism, Genome> pair = Args.parseGenome(args);
			if(pair==null){
				//Make fake genome... chr lengths provided???
				if(ap.hasKey("geninfo")){
					genome = new Genome("Genome", new File(ap.getKeyValue("geninfo")));
	        	}else{
	        		System.err.println("No genome provided."); 
	        		System.exit(1);
	        	}
			}else{
				org = pair.car();
				genome = pair.cdr();
			}
		} catch (NotFoundException e) {
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
		strength = Args.parseInteger(args,"strength",40);
		
		// parameter for building empirical distribution
		String chipSeqFile = Args.parseString(args, "chipSeqFile", null);	
		String fileFormat = Args.parseString(args, "format", "BED");    	
		int readLength = Args.parseInteger(args,"readlen",36);
		List<File> expts = new ArrayList<File>();
		expts.add(new File(chipSeqFile));
		chipSeqExpt = new DeepSeqExpt(genome, expts, false, fileFormat, readLength);
		range = Args.parseInteger(args, "range", range);
		name = Args.parseString(args, "name", "noname");
		
		// load motif
		useMotif = Args.parseInteger(args, "useMotif", 0)==1?true:false;
		if (useMotif){
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
		}
	}
	
	public void buildEmpiricalDistribution(){
		WeightMatrixScorer scorer=null;
		if (useMotif){
			scorer = new WeightMatrixScorer(motif);
		}
		ArrayList<Point> points = new ArrayList<Point>();
		for (GPSPeak gps: gpsPeaks){
//			if (gps.getMixProb()==1 && gps.getQvalue()>5 && gps.getShape()<1)
			if (gps.getMixProb()==1 && gps.getStrength()>strength && gps.getShape()<1){
				if (useMotif){
					Region r= gps.expand(MOTIF_DISTANCE);
					WeightMatrixScoreProfile profiler = scorer.execute(r);
					int halfWidth = profiler.getMatrix().length()/2;
					//search from BS outwards, to find the nearest strong motif
					for(int z=0; z<=r.getWidth()/2; z++){
						double leftScore= profiler.getMaxScore(MOTIF_DISTANCE-z);
						double rightScore= profiler.getMaxScore(MOTIF_DISTANCE+z);	
						int position = -Integer.MAX_VALUE;	// middle of motif, relative to GPS peak
						if(rightScore>=motifThreshold){
							position = z+halfWidth;
						}
						if(leftScore>=motifThreshold){
							position = -z+halfWidth;
						}
						if(position > -Integer.MAX_VALUE){
							Point motifPos = new Point(genome, gps.getChrom(), gps.getLocation()+position);
							points.add(motifPos);
							break; 	// break from the motif search of this peak
						}
					}
				}
				else
					points.add((Point)gps);
			}
		}
		printEmpiricalDistribution(points);
		System.out.println(points.size()+" GPS peaks are used to build emprical distribution.");
	}

	private void printEmpiricalDistribution(ArrayList<Point> points){
		BindingModel model_plus = getStrandedDistribution (points, '+');
		model_plus.printToFile(name+"_plus_Distribution.txt");
		BindingModel model_minus = getStrandedDistribution (points, '-');
		model_minus.printToFile(name+"_minus_Distribution.txt");
		
		double[] prob_plus = model_plus.getProbabilities();
		double[] prob_minus = model_minus.getProbabilities();
		BindingModel.minKL_Shift(prob_plus, prob_minus);
		
		ArrayList<Pair<Integer, Double>> dist = new ArrayList<Pair<Integer, Double>>();
		for(int i=model_plus.getMin();i<=model_plus.getMax();i++){
			int index = i-model_plus.getMin();
			dist.add(new Pair<Integer, Double>(i, prob_plus[index]+prob_minus[index]));
		}
		BindingModel model=new BindingModel(dist);
		model.printToFile("Empirical_Distribution_"+name+".txt");
		System.out.println("Empirical_Distribution_"+name+".txt is written.");
	}	
	
	private BindingModel getStrandedDistribution (ArrayList<Point> points, char strand){
		float[] sum = new float[2*range+1];
		Pair<ArrayList<Integer>,ArrayList<Float>> pair = null;
		for (Point p:points){
			int pos = p.getLocation();
			pair = chipSeqExpt.loadStrandedBaseCounts(p.expand(range), strand);
			for (int i=0;i<pair.car().size();i++){
				sum[pair.car().get(i)-pos+range]+=pair.cdr().get(i);
			}
		}

		ArrayList<Pair<Integer, Double>> dist = new ArrayList<Pair<Integer, Double>>();
		if (strand=='-')
			for(int i=sum.length-1;i>=0;i--){
				int pos = -(i-sum.length/2);
				dist.add(new Pair<Integer, Double>(pos, (double)sum[i]));
			}
		else
			for(int i=0;i<sum.length;i++){
				int pos = i-sum.length/2;
				dist.add(new Pair<Integer, Double>(pos, (double)sum[i]));
			}
		BindingModel model=new BindingModel(dist);
		model.smooth(BindingModel.SMOOTHING_STEPSIZE);
		return model;
	}
 
}

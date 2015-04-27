package edu.mit.csail.cgs.deepseq.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.mit.csail.cgs.datasets.chipseq.ChipSeqLocator;
import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedPoint;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.BindingModel;
import edu.mit.csail.cgs.deepseq.BindingModel2D;
import edu.mit.csail.cgs.deepseq.DeepSeqExpt;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSPeak;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class GPS_ReadDistribution {
	static final int MOTIF_DISTANCE = 50;	
	
	private Genome genome;
	private Organism org;
	
	private ArrayList<Point> points = new ArrayList<Point>(); 
	private String GPSfileName;
	private int strength=40;
	private int mrc = 100;

	// build empirical distribution
	private DeepSeqExpt chipSeqExpt = null;
	private int range = 250;
	private int smooth_step = 10;
	private int top = 1000;
	
	private WeightMatrix motif = null;
	private double motifThreshold;
	
	private String name = null;
	private String motif_strand = null;
	private boolean print_4_orientations;
	private boolean isMirrorSymmetry;
	
	public static void main(String[] args) throws IOException {
		GPS_ReadDistribution analysis = new GPS_ReadDistribution(args);
		analysis.printEmpiricalDistribution(analysis.points);;
	}
	
	public GPS_ReadDistribution(String[] args) throws IOException {
		ArgParser ap = new ArgParser(args);
		Set<String> flags = Args.parseFlags(args);
		print_4_orientations = flags.contains("p4");
		isMirrorSymmetry = flags.contains("mirrow");
		if (isMirrorSymmetry)
			System.out.println("Running MirrorSymmetry mode ... ");
		try {
			Pair<Organism, Genome> pair = Args.parseGenome(args);
			if(pair==null){
				//Make fake genome... chr lengths provided???
				if(ap.hasKey("g")){
					genome = new Genome("Genome", new File(ap.getKeyValue("g")), true);
	        	}else{
	        		System.err.println("No genome information provided."); 
	        		System.exit(1);
	        	}
			}else{
				org = pair.car();
				genome = pair.cdr();
			}
		} catch (NotFoundException e) {
			e.printStackTrace();
		}
		
//		Map<String, Integer> map = genome.getChromLengthMap();
//		for (String chr:map.keySet()){
//			System.out.println(chr+"\t"+map.get(chr));
//		}
		
		// parameter for building empirical distribution
		String chipSeqFile = Args.parseString(args, "chipseq", null);
		if (chipSeqFile!=null){
			String fileFormat = Args.parseString(args, "f", "BED");  
			List<File> expts = new ArrayList<File>();
			expts.add(new File(chipSeqFile));
			chipSeqExpt = new DeepSeqExpt(genome, expts, false, fileFormat, -1);
		}
		else{
			List<ChipSeqLocator> rdbexpts = Args.parseChipSeq(args,"rdb");
			chipSeqExpt = new DeepSeqExpt(genome, rdbexpts, "readdb", -1);
		}
			
		range = Args.parseInteger(args, "range", range);
		top = Args.parseInteger(args, "top", top);			// number of top point positions to use
		mrc = Args.parseInteger(args, "mrc", mrc);			// max read count
		name = Args.parseString(args, "name", "noname");
		motif_strand = Args.parseString(args, "motif_strand", null);
		smooth_step = Args.parseInteger(args, "smooth", smooth_step);
		// load points
		String coordFile = Args.parseString(args, "coords", null);
		ArrayList<Point> coords = null;
		if (coordFile!=null){
			coords = loadCgsPointFile(coordFile, top);
		}
		else{
			// load GPS results
			GPSfileName = Args.parseString(args, "GPS", null);
			if (GPSfileName==null){
				System.err.println("Coordinate file not found!");
				System.exit(0);
			}
			File gpsFile = new File(GPSfileName);
			List<GPSPeak>gpsPeaks = GPSParser.parseGPSOutput(gpsFile.getAbsolutePath(), genome);
			strength = Args.parseInteger(args,"strength",40);
			coords = new ArrayList<Point>();
			for (GPSPeak gps: gpsPeaks){
				if ((!gps.isJointEvent()) && gps.getStrength()>strength )
					coords.add(gps);
			}
		}
		
		// load motif
	    Pair<WeightMatrix, Double> wm = null;
	    WeightMatrixScorer scorer=null;
	    if (Args.parseString(args, "motif", null)!=null||Args.parseString(args, "pfm", null)!=null){
			wm = CommonUtils.loadPWM(args, org.getDBID());
			System.out.println("Using motif "+wm.getFirst().name);
			scorer = new WeightMatrixScorer(wm.getFirst());
			motifThreshold = wm.cdr().doubleValue();
	    }		
		// get points or motif points
		for (Point p: coords){
			if (points.size()>=top)
				break;
			
			if (wm!=null){
				Region r= p.expand(MOTIF_DISTANCE);
				WeightMatrixScoreProfile profiler = scorer.execute(r);
				int halfWidth = profiler.getMatrix().length()/2;
				//search from BS outwards, to find the nearest strong motif
				for(int z=0; z<=r.getWidth()/2; z++){
					if (MOTIF_DISTANCE+z>=profiler.length())
						continue;
					double leftScore= profiler.getHigherScore(MOTIF_DISTANCE-z);
					double rightScore= profiler.getHigherScore(MOTIF_DISTANCE+z);	
					int position = -Integer.MAX_VALUE;	// middle of motif, relative to GPS peak
					char strand = '*';
					if(rightScore>=motifThreshold){
						position = z+halfWidth;
						strand = profiler.getHigherScoreStrand(MOTIF_DISTANCE+z);
					}
					if(leftScore>=motifThreshold && leftScore>rightScore){
						position = -z+halfWidth;
						strand = profiler.getHigherScoreStrand(MOTIF_DISTANCE-z);
					}
					if(position > -Integer.MAX_VALUE){
						Point motifPos = null;
						if (!flags.contains("unstranded"))
							motifPos = new StrandedPoint(genome, p.getChrom(), p.getLocation()+position, strand);
						else if (motif_strand==null || (motif_strand.equalsIgnoreCase("F")&& strand=='+')||(motif_strand.equalsIgnoreCase("R")&& strand=='-'))
							motifPos = new Point(genome, p.getChrom(), p.getLocation()+position);
						else
							break;
						points.add(motifPos);
						break; 	// break from the motif search of this peak
					}
				}
			}
			else
				points.add(p);
		}	// each point
		System.out.println("Computing profiles using "+points.size()+" coordinates ...");
	}

	private void printEmpiricalDistribution(ArrayList<Point> points){
//		BindingModel model_plus = getStrandedDistribution (points, '+');
//		model_plus.printToFile(String.format("Read_Distribution_%s_%d_plus.txt", name, points.size()));
//		BindingModel model_minus = getStrandedDistribution (points, '-');
//		model_minus.printToFile(String.format("Read_Distribution_%s_%d_minus.txt", name, points.size()));
//		
//		double[] prob_plus = model_plus.getProbabilities();
//		double[] prob_minus = model_minus.getProbabilities();
//		BindingModel.minKL_Shift(prob_plus, prob_minus);
//		
//		ArrayList<Pair<Integer, Double>> dist = new ArrayList<Pair<Integer, Double>>();
//		for(int i=model_plus.getMin();i<=model_plus.getMax();i++){
//			int index = i-model_plus.getMin();
//			dist.add(new Pair<Integer, Double>(i, prob_plus[index]+prob_minus[index]));
//		}
		BindingModel model=getDistribution(points);
		String outFile = String.format("Read_Distribution_%s.txt", name);
		model.printToFile(outFile);
		System.out.println(outFile+" has been written.");
	}	
	
	private BindingModel getStrandedDistribution (ArrayList<Point> points, char strand){
		Map<String,Integer> chromLengthMap = genome.getChromLengthMap();
		float[] sum = new float[2*range+1];
		Pair<ArrayList<Integer>,ArrayList<Float>> pair = null;
		for (Point p:points){
			int pos = p.getLocation();
			if (!chromLengthMap.containsKey(p.getChrom()) || pos>chromLengthMap.get(p.getChrom()))
				continue;
			if (p instanceof StrandedPoint){
				char point_strand = ((StrandedPoint) p).getStrand();
				pair = chipSeqExpt.loadStrandedBaseCounts(p.expand(range), point_strand=='+'?strand:(char)(88-strand));
				// convert absolute coordinates to relative offset
				ArrayList<Integer> coords = pair.car();
				for (int i=0;i<coords.size();i++){
					int offset = coords.get(i)-pos;
					if (point_strand=='-')
						offset = -offset;
					coords.set(i, offset);
				}
			}
			else{
				pair = chipSeqExpt.loadStrandedBaseCounts(p.expand(range), strand);
				// convert absolute coordinates to relative offset
				ArrayList<Integer> coords = pair.car();
				for (int i=0;i<coords.size();i++){
					int offset = coords.get(i)-pos;
					coords.set(i, offset);
				}
			}
			for (int i=0;i<pair.car().size();i++){
				sum[pair.car().get(i)+range] += Math.min(pair.cdr().get(i), mrc);
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
		if (smooth_step>0)
			model.smooth(smooth_step, smooth_step);
		return model;
	}
	
	private BindingModel2D getDistribution2D (ArrayList<Point> points) {
	    Map<String,Integer> chromLengthMap = genome.getChromLengthMap();
	    float[][] sum = new float[2*range+1][2*range+1];
	    Pair<Pair<ArrayList<Integer>, ArrayList<ArrayList<Integer>>>, ArrayList<ArrayList<Float>>> pair = null;
	    for (Point p:points) {
	        int pos = p.getLocation();
	        if (!chromLengthMap.containsKey(p.getChrom()) || pos>chromLengthMap.get(p.getChrom()))
                continue;
	        pair = chipSeqExpt.loadStrandedBaseCountsPaired(p.expand(range), '+');
	        Pair<ArrayList<Integer>, ArrayList<ArrayList<Integer>>> coords = pair.car();
	        ArrayList<Integer> plusCoords = coords.car();
	        ArrayList<ArrayList<Integer>> minCoords = coords.cdr();
	        ArrayList<ArrayList<Float>> weight = pair.cdr();
	        for (int i=0; i<minCoords.size(); i++) {
	            int upOffset = plusCoords.get(i)-pos;
	            for (int j=0; j<minCoords.get(i).size(); i++){
	                int downOffset = minCoords.get(i).get(j)-pos;
	                if (downOffset <= range) {
	                    sum[upOffset+range][downOffset+range] += Math.min(weight.get(i).get(j), mrc);
	                }
	            }
	        }
	    }
	    ArrayList<Pair<Integer, List<Pair<Integer, Double>>>> dist = new ArrayList<Pair<Integer, List<Pair<Integer, Double>>>>();
	    return null;
	}
	
	// Either to negate the offset of minus strand depends on the type of data
	// ChIP-Seq data is point-symmetry (i.e. 180 degree rotation, the sequenced ends are from different ends of the pulled down fragment)
	// DNase-Seq data is mirror-symmetry (i.e. the sequenced ends are from the same cut)
	private BindingModel getDistribution (ArrayList<Point> points){
		Map<String,Integer> chromLengthMap = genome.getChromLengthMap();
		int length = range*2+1;
		double[] pp = new double[length];		// motif/point and DNA are all on plus strand
		double[] pm = new double[length];		// motif/point: plus,  DNA: minus
		double[] mp = new double[length];
		double[] mm = new double[length];
		Pair<ArrayList<Integer>,ArrayList<Float>> pair = null;
		for (Point p:points){
			int pos = p.getLocation();
			if (!chromLengthMap.containsKey(p.getChrom()) || pos>chromLengthMap.get(p.getChrom()))
				continue;
			if (p instanceof StrandedPoint){
				char point_strand = ((StrandedPoint) p).getStrand();
				
				pair = chipSeqExpt.loadStrandedBaseCounts(p.expand(range), '+');
				ArrayList<Integer> coords = pair.car();
				for (int i=0;i<coords.size();i++){
					int offset = coords.get(i)-pos;		// convert absolute coordinates to relative offset
					if (point_strand=='-'){
//						if (isDNase)
//							offset = -offset;
						mp[offset+range] += Math.min(pair.cdr().get(i), mrc);
					}
					else{
						pp[offset+range] += Math.min(pair.cdr().get(i), mrc);
					}
				}
				pair = chipSeqExpt.loadStrandedBaseCounts(p.expand(range), '-');
				coords = pair.car();
				for (int i=0;i<coords.size();i++){
					int offset = coords.get(i)-pos;		// convert absolute coordinates to relative offset
					if (point_strand=='-'){
//						offset = -offset;
						mm[offset+range] += Math.min(pair.cdr().get(i), mrc);
					}
					else{
//						if (!isDNase)
//							offset = -offset;
						pm[offset+range] += Math.min(pair.cdr().get(i), mrc);
					}
				}
			}
			else{	// unstranded
				pair = chipSeqExpt.loadStrandedBaseCounts(p.expand(range), '+');
				ArrayList<Integer> coords = pair.car();
				for (int i=0;i<coords.size();i++){
					int offset = coords.get(i)-pos;		// convert absolute coordinates to relative offset
					pp[offset+range] += Math.min(pair.cdr().get(i), mrc);
				}
				pair = chipSeqExpt.loadStrandedBaseCounts(p.expand(range), '-');
				coords = pair.car();
				for (int i=0;i<coords.size();i++){
					int offset = coords.get(i)-pos;		// convert absolute coordinates to relative offset
					mm[offset+range] += Math.min(pair.cdr().get(i), mrc);
				}
			}
		}

		// combine data from plus and minus strand
		ArrayList<Pair<Integer, Double>> dist = new ArrayList<Pair<Integer, Double>>();	
		if (points.get(0) instanceof StrandedPoint){
			if (print_4_orientations){
				ArrayList<Pair<Integer, Double>> one_orientation = new ArrayList<Pair<Integer, Double>>();		// reads on the plus strand 
				for (int i=0;i<length;i++){
					one_orientation.add(new Pair<Integer, Double>(i-range, pp[i]));
				}
				BindingModel model_one_orientation=new BindingModel(one_orientation);
				if (smooth_step>0)
					model_one_orientation.smooth(smooth_step, smooth_step);
				model_one_orientation.printToFile(String.format("Read_Distribution_%s_pp.txt", name));
				
				one_orientation.clear();
				for (int i=0;i<length;i++){
					one_orientation.add(new Pair<Integer, Double>(i-range, mm[i]));
				}
				model_one_orientation=new BindingModel(one_orientation);
				if (smooth_step>0)
					model_one_orientation.smooth(smooth_step, smooth_step);
				model_one_orientation.printToFile(String.format("Read_Distribution_%s_mm.txt", name));
				
				one_orientation.clear();
				for (int i=0;i<length;i++){
					one_orientation.add(new Pair<Integer, Double>(i-range, mp[i]));
				}
				model_one_orientation=new BindingModel(one_orientation);
				if (smooth_step>0)
					model_one_orientation.smooth(smooth_step, smooth_step);
				model_one_orientation.printToFile(String.format("Read_Distribution_%s_mp.txt", name));
				
				one_orientation.clear();
				for (int i=0;i<length;i++){
					one_orientation.add(new Pair<Integer, Double>(i-range, pm[i]));
				}
				model_one_orientation=new BindingModel(one_orientation);
				if (smooth_step>0)
					model_one_orientation.smooth(smooth_step, smooth_step);
				model_one_orientation.printToFile(String.format("Read_Distribution_%s_pm.txt", name));
			}
			
			// reverse the positions to merge
			if (isMirrorSymmetry)
				mp = reverseArray(mp);
			if (!isMirrorSymmetry)
				pm = reverseArray(pm);
			mm = reverseArray(mm);
			
			BindingModel.minKL_Shift(pp, mm, 5);
			BindingModel.minKL_Shift(mp, pm, 5);
			ArrayList<Pair<Integer, Double>> same = new ArrayList<Pair<Integer, Double>>();		// reads on the same strand as motif
			for (int i=0;i<length;i++){
				pp[i]+=mm[i];		// store in pp
				mp[i]+=pm[i];		// store in mp
			}
			for (int i=0;i<length;i++){
				same.add(new Pair<Integer, Double>(i-range, pp[i]));
			}
			BindingModel model_same=new BindingModel(same);
			if (smooth_step>0)
				model_same.smooth(smooth_step, smooth_step);
			model_same.printToFile(String.format("Read_Distribution_%s_same.txt", name));
			ArrayList<Pair<Integer, Double>> diff = new ArrayList<Pair<Integer, Double>>();		// on diff strand
			for (int i=0;i<length;i++){
				diff.add(new Pair<Integer, Double>(i-range, mp[i]));
			}
			BindingModel model_diff=new BindingModel(diff);
			if (smooth_step>0)
				model_diff.smooth(smooth_step, smooth_step);
			model_diff.printToFile(String.format("Read_Distribution_%s_diff.txt", name));
			BindingModel.minKL_Shift(pp, mp, 5);
			for (int i=0;i<length;i++){
				dist.add(new Pair<Integer, Double>(i-range, pp[i]+mp[i]));
			}
		}
		else{
			ArrayList<Pair<Integer, Double>> plus = new ArrayList<Pair<Integer, Double>>();		// reads on the plus strand 
			for (int i=0;i<length;i++){
				plus.add(new Pair<Integer, Double>(i-range, pp[i]));
			}
			BindingModel model_plus=new BindingModel(plus);
			if (smooth_step>0)
				model_plus.smooth(smooth_step, smooth_step);
			model_plus.printToFile(String.format("Read_Distribution_%s_plus.txt", name));
			ArrayList<Pair<Integer, Double>> minus = new ArrayList<Pair<Integer, Double>>();		// on minus strand
			for (int i=0;i<length;i++){
				minus.add(new Pair<Integer, Double>(i-range, mm[i]));
			}
			BindingModel model_minus=new BindingModel(minus);
			if (smooth_step>0)
				model_minus.smooth(smooth_step, smooth_step);
			model_minus.printToFile(String.format("Read_Distribution_%s_minus.txt", name));

			BindingModel.minKL_Shift(pp, mm, 5);
			for (int i=0;i<length;i++){
				dist.add(new Pair<Integer, Double>(i-range, pp[i]+mm[i]));
			}
		}

		BindingModel model=new BindingModel(dist);
		if (smooth_step>0)
			model.smooth(smooth_step, smooth_step);
		return model;
	}
	
	private double[] reverseArray(double[] a){
		double[] b = new double[a.length];
		for (int i=0;i<a.length;i++)
			b[i]=a[a.length-1-i];
		return b;
	}
	private ArrayList<Point> loadCgsPointFile(String filename, int ptCount) {

		File file = new File(filename);
		if(!file.isFile()){
			System.err.println("\nCannot find coordinate file!");
			System.exit(1);
		}
		FileReader in = null;
		BufferedReader bin = null;
		ArrayList<Point> points = new ArrayList<Point>();
		try {
			in = new FileReader(file);
			bin = new BufferedReader(in);
			String line;
			while((line = bin.readLine()) != null) { 
				if (points.size()>=ptCount)
					break;
				line = line.trim();
				Region point = Region.fromString(genome, line);
				if (point!=null){
					if (point instanceof StrandedRegion){
						StrandedRegion sr = (StrandedRegion)point;
						points.add(new StrandedPoint(genome, point.getChrom(),point.getStart(),sr.getStrand()));
					}
					else
						points.add(new Point(genome, point.getChrom(),point.getStart()));
				}
			}
		}
		catch(IOException ioex) {
			System.err.println("Error when parsing coordinate file! ");
			ioex.printStackTrace(System.err);
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
}

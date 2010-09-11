package edu.mit.csail.cgs.deepseq.analysis;

import java.awt.Color;
import java.io.*;
import java.util.*;

import edu.mit.csail.cgs.datasets.chipseq.ChipSeqLocator;
import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.BindingModel;
import edu.mit.csail.cgs.ewok.verbs.chipseq.ChipSeqExpander;
import edu.mit.csail.cgs.ewok.verbs.chipseq.MACSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.MACSPeakRegion;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.metagenes.BinningParameters;
import edu.mit.csail.cgs.metagenes.ChipSeq5PrimeProfiler;
import edu.mit.csail.cgs.metagenes.MetaNonFrame;
import edu.mit.csail.cgs.metagenes.PointProfile;
import edu.mit.csail.cgs.metagenes.MetaProfile;
import edu.mit.csail.cgs.metagenes.MetaProfileHandler;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class ChipSeqStats {

	static final double default_p_value = 90;
	static final double weak_p_value = 50;
	static final int WINDOW = 600;
	static final int BIN_NUM = 600;	
	static final int MOTIF_DISTANCE = 150;	
	static final int SAMPLE_NUM = 100;	
	static final double MOTIF_THRESH=14;
	
	static final File dir = new File("C:\\Data\\ES");
	static final String[] allFactors = {"Oct4", "Sox2", "Nanog", "Tcf3",  
		"smad1", "stat3", "esrrb", "tcfcp2I1", "e2f1", "zfx", 
		"klf4", "c-myc", "n-myc", /*"p300",*/ "suz12", "ctcf", 
		"oct4S", "sox2S", "nanogS"};

	static Genome g;
	static Organism org;
	static{
		try {
			g = Organism.findGenome("mm8");
			org = Organism.getOrganism("Mus musculus");
		} catch (NotFoundException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
//		peakMotifAnalysis();
//		printPeakReadCount("Oct4", "YL_Oct4_ES", "bowtie_unique");
//		printMotifHits();
//		refinePeakModel("Oct4", "YoungLab_Solexa_Oct4", 20000);
//		printNarrowPeaks("Oct4", 20000);
//		readStartDistance();
//		sortBEDFile();
	}
	private static void sortBEDFile(){
		Genome g=null;
		try {
			g = Organism.findGenome("mm8");
		} catch (NotFoundException e) {
			e.printStackTrace();
		}
		File file = new File("C:\\Data\\ES\\data\\mES_Oct4.bed");
		FileReader in = null;
		BufferedReader bin = null;
		ArrayList<BEDRecord> reads = new ArrayList<BEDRecord>();
		try {
			in = new FileReader(file);
			bin = new BufferedReader(in);
			String line;
			while((line = bin.readLine()) != null) { 
				line = line.trim();
	            String[] f=line.split("\t");
	            BEDRecord bed = new BEDRecord(g, f[0], Integer.parseInt(f[1]), Integer.parseInt(f[2]), f[5].charAt(0), f[3]);
	            reads.add(bed);
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
		Collections.sort(reads);
		
		StringBuffer sb = new StringBuffer();
		try{
			File outFile = new File("C:\\Data\\ES\\data\\mES_Oct4_sorted.bed");
			PrintStream ps = new PrintStream(new FileOutputStream(outFile));
			int i=100000;
			for(BEDRecord r: reads){
				sb.append(r.bedString());
				if (i==0){
					ps.print(sb.toString());
					sb =  new StringBuffer();
					i=100000;
				}
				i--;
			}
			ps.print(sb.toString());
			ps.close();
			System.out.println("File is written to "+outFile.getAbsolutePath());
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	private static void getBEDbyRegion(){
		
	}
	public static void printNarrowPeaks(String tf, int number){
		File file = new File(dir, tf+".tab");
		List<MACSPeakRegion> peaks = MACSPeakRegion.filterByPValue(MACSParser.parseMACSOutput(file.getAbsolutePath(), g),  default_p_value, Double.MAX_VALUE);  
		peaks = MACSPeakRegion.filterByTags(peaks, 50, 250);
		//sort by tags/width
		Collections.sort(peaks, new Comparator<MACSPeakRegion>() {
	        public int compare(MACSPeakRegion p1, MACSPeakRegion p2) {
	        	double ratio1 = (double)p1.getTags()/p1.getWidth();
	        	double ratio2 = (double)p2.getTags()/p2.getWidth();
	        	return (ratio1==ratio2)?0:((ratio1>ratio2)?-1:1);
	        };});
		
		int outputCount = Math.min(peaks.size(), number);
		StringBuilder sb=new StringBuilder();
		for (int i=0; i<outputCount;i++){
			MACSPeakRegion p = peaks.get(i);
			sb.append(p.getPeak().getLocationString()).append("\t");
			sb.append(p.getWidth()).append("\t");
			sb.append(p.getTags()).append("\t").append("\n");
		}
		writeToFile(tf+"_narrowPeak_"+outputCount+".txt", sb.toString());
	}
	
	public static void refinePeakModel(String tf, String expt, String version, int number){
		File file = new File(dir, tf+".tab");
		List<MACSPeakRegion> all_peaks = MACSPeakRegion.filterByPValue(MACSParser.parseMACSOutput(file.getAbsolutePath(), g),  default_p_value, Double.MAX_VALUE);  
		List<MACSPeakRegion> peaks = MACSPeakRegion.filterByTags(all_peaks, 50, 250);
//		List<MACSPeakRegion> peaks = all_peaks;
		//sort by tags/width
		Collections.sort(peaks, new Comparator<MACSPeakRegion>() {
	        public int compare(MACSPeakRegion p1, MACSPeakRegion p2) {
	        	double ratio1 = (double)p1.getTags()/p1.getWidth();
	        	double ratio2 = (double)p2.getTags()/p2.getWidth();
	        	return (ratio1==ratio2)?0:((ratio1>ratio2)?-1:1);
	        };});
		
		int outputCount = Math.min(peaks.size(), number);
		ArrayList<Point> points = new ArrayList<Point>();
		HashMap<Point, MACSPeakRegion> point2macs=new HashMap<Point, MACSPeakRegion>();
		for (int i=0; i<outputCount;i++){
			MACSPeakRegion p=peaks.get(i);
			points.add(p.getPeak());
			point2macs.put(p.getPeak(), p);
		}
		// build meta peak
		MetaProfile metaProfile = buildMetaPeak (points, expt, version, '+');
		ArrayList<PointProfile> profiles = new ArrayList<PointProfile>();
		for (int i=0;i<metaProfile.getNumProfiles();i++){
			profiles.add((PointProfile)metaProfile.profile(i));
		}
//		metaProfile.saveToFile(tf+"_0_metaProfile.txt");

		MetaProfile newMetaProfile = refineMetaProfile(profiles, metaProfile, outputCount/2);
//		newMetaProfile.saveToFile(tf+"_1_metaProfile.txt");
		for(int i=2;i<6;i++){
			newMetaProfile = refineMetaProfile(profiles, newMetaProfile, outputCount/2);
//			newMetaProfile.saveToFile(tf+"_"+i+"_metaProfile.txt");
		}

		StringBuilder sb=new StringBuilder();
		final Map<Integer, Double> quantiles = getEqualSpaces(newMetaProfile, SAMPLE_NUM);
//		for (int b:quantiles.keySet()){
//			System.out.println(b+"\t"+quantiles.get(b)+"\n");
//		}
		//for (int i=0;i<newMetaProfile.getNumProfiles();i++){
//			PointProfile p = (PointProfile)newMetaProfile.profile(i);
//			sb.append(p.getPoint().getLocationString()).append("\t");
//			sb.append((int)(getDistributionDistance(p, quantiles)*100000)+"\t");
//			for (int j=0;j<p.length();j++){
//				sb.append((int)p.value(j)).append("\t");
//			}
//			sb.append("\n");
//		}
		WeightMatrix motif = null;
		try {
			int wmid = WeightMatrix.getWeightMatrixID(org.getDBID(), "Oct-4 (POU5F1)", "TRANSFAC 10.4, M01124");
			motif = WeightMatrix.getWeightMatrix(wmid);
		} catch (NotFoundException e) {
			e.printStackTrace();
		}
		WeightMatrixScorer scorer = new WeightMatrixScorer(motif);
		for (int i=0;i<profiles.size();i++){
			PointProfile p = profiles.get(i);
			sb.append(p.getPoint().getLocationString()).append("\t");
			Region r= p.getPoint().expand(MOTIF_DISTANCE);
			WeightMatrixScoreProfile profiler = scorer.execute(r);
			double bestScore=0;
			for(int z=0; z<r.getWidth(); z++){
				bestScore= Math.max(bestScore, profiler.getMaxScore(z));
			}
			sb.append(bestScore+"\t");
			sb.append(point2macs.get(p.getPoint()).getPvalue()+"\t");
			sb.append(point2macs.get(p.getPoint()).getTags()+"\t");
			sb.append((int)(getDistributionDistance(p, quantiles)*100000)+"\t");
			for (int j=0;j<p.length();j++){
				sb.append((int)p.value(j)).append("\t");
			}
			sb.append("\n");
		}
		writeToFile(tf+"_plus_best_profiles.txt", sb.toString());
	}
	public static void printPeakReadCount(ArrayList<Point> points, String expt, String version){

		// build meta peak
		MetaProfile metaProfile = buildMetaPeak (points, expt, version, '+');
		ArrayList<PointProfile> profiles = new ArrayList<PointProfile>();
		for (int i=0;i<metaProfile.getNumProfiles();i++){
			profiles.add((PointProfile)metaProfile.profile(i));
		}
		
		MetaProfile minus_Profile = buildMetaPeak (points, expt, version, '-');
		
		StringBuilder sb=new StringBuilder();

		for (int i=0;i<profiles.size();i++){
			PointProfile p = profiles.get(i);
			sb.append(p.getPoint().getLocationString()).append("\t");
			for (int j=0;j<p.length();j++){
				sb.append((int)p.value(j)).append("\t");
			}
			PointProfile p_minus = (PointProfile)minus_Profile.profile(i);
			for (int j=0;j<p_minus.length();j++){
				sb.append((int)p_minus.value(j)).append("\t");
			}
			sb.append("\n");
		}
		writeToFile(expt+"_peak_readCount_profiles.txt", sb.toString());
	}
	public static void printEmpiricalDistribution(ArrayList<Point> points, String expt, String version){
		BindingModel model_plus = getStrandedDistribution (points, expt, version, '+');
//		model_plus.printToFile(expt+"_plus_Distribution.txt");
		BindingModel model_minus = getStrandedDistribution (points, expt, version, '-');
//		model_minus.printToFile(expt+"_minus_Distribution.txt");
		
		double[] prob_plus = model_plus.getProbabilities();
		double[] prob_minus = model_minus.getProbabilities();
		BindingModel.minKL_Shift(prob_plus, prob_minus);
		
		ArrayList<Pair<Integer, Double>> dist = new ArrayList<Pair<Integer, Double>>();
		for(int i=model_plus.getMin();i<=model_plus.getMax();i++){
			int index = i-model_plus.getMin();
			dist.add(new Pair<Integer, Double>(i, prob_plus[index]+prob_minus[index]));
		}
		BindingModel model=new BindingModel(dist);
		model.printToFile("Empirical_Distribution_"+expt+".txt");
	}	
	
	private static BindingModel getStrandedDistribution (ArrayList<Point> points, String expt, String version, char strand){
		MetaProfile metaProfile = buildMetaPeak (points, expt, version, strand);
		int[] sum = new int[metaProfile.profile(0).length()];
		for (int i=0;i<metaProfile.getNumProfiles();i++){
			PointProfile p = (PointProfile)metaProfile.profile(i);
			for (int j=0;j<p.length();j++){
				sum[j]+=p.value(j);
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
		model.smooth(BindingModel.SMOOTHING_STEPSIZE, BindingModel.SMOOTHING_AVG_PTS);
		return model;
	}
	
	private static MetaProfile buildMetaPeak (List<Point> points, String expt, String version, char strand){
		BinningParameters params = new BinningParameters(WINDOW, BIN_NUM);
		MetaProfile metaProfile=null;
		ArrayList<ChipSeqLocator> exptlocs = new ArrayList<ChipSeqLocator>();
		exptlocs.add(new ChipSeqLocator(expt, version));
		try{
			ArrayList<ChipSeqExpander> exptexps = new ArrayList<ChipSeqExpander>();
			for(ChipSeqLocator loc : exptlocs){
				exptexps.add(new ChipSeqExpander(loc));
			}
			System.out.println("Loading data...");
			ChipSeq5PrimeProfiler profiler = new ChipSeq5PrimeProfiler(params, exptexps, strand);
			MetaNonFrame nonframe = new MetaNonFrame(g, params, profiler, false);
			nonframe.setColor(Color.blue);
			MetaProfileHandler handler = nonframe.getHandler();
			handler.addPoints(points);
			while(handler.addingPoints()){}
			metaProfile = handler.getProfile();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return metaProfile;
	}
	
	private static Map<Integer, Double> getQuantiles(MetaProfile profile, int num){
		Map<Integer, Double> map = new TreeMap<Integer, Double> ();
		int count = profile.length();
		double total=0;
		for (int i=0; i<count; i++){
			total += profile.value(i);
		}
		
		double avg = total/num;
		double binSum=0;
		for (int i=0; i<count; i++){
			binSum += profile.value(i);
			if (binSum>=avg){
				map.put(i, binSum/total);
				binSum=0;
			}
		}
		if (binSum!=0){
			map.put(count-1, binSum/total);
		}
		return map;
	}

	private static Map<Integer, Double> getEqualSpaces(MetaProfile profile, int num){
		Map<Integer, Double> map = new TreeMap<Integer, Double> ();
		int count = profile.length();
		double total=0;
		for (int i=0; i<count; i++){
			total += profile.value(i);
		}
		
		int next = count/num;
		double binSum=0;
		for (int i=0; i<count; i++){
			binSum += profile.value(i);
			if (i>=next){
				map.put(i, binSum/total);
				binSum=0;
				next+=count/num;
			}
		}
		if (binSum!=0){
			map.put(count-1, binSum/total);
		}
		return map;
	}
	private static double getDistributionDistance(PointProfile profile, Map<Integer, Double> distribution){
		double total = profile.total();
		int start=0;
		double distance = 0;
		for (int end:distribution.keySet()){
			double binSum=0;
			for (int i=start; i<=end; i++){
				binSum += profile.value(i);
			}
//			distance+=Math.abs((distribution.get(end)-binSum/total));
			distance+=(distribution.get(end)-binSum/total)*(distribution.get(end)-binSum/total);
			start = end+1;
		}
		return distance;
	}
	
	private static MetaProfile refineMetaProfile(ArrayList<PointProfile> profiles, MetaProfile metaProfile, int num){
		// use meta peak to rank most similar peaks
		System.out.println("\n"+num+"\t"+profiles.size());
		final Map<Integer, Double> quantiles = getEqualSpaces(metaProfile, SAMPLE_NUM);
		Collections.sort(profiles, new Comparator<PointProfile>() {
	        public int compare(PointProfile p1, PointProfile p2) {
	        	double dist1 = getDistributionDistance(p1, quantiles);
	        	double dist2 = getDistributionDistance(p2, quantiles);
	        	return (dist1==dist2)?0:((dist1<dist2)?-1:1);
	        };});
		// use top most similar peaks to build new metaProfile
		MetaProfile newMetaProfile = new MetaProfile("MetaProfile", metaProfile.getBinningParameters());
		for (int i=0;i<num;i++){
			newMetaProfile.addProfile(profiles.get(i));
			System.out.print((int)(getDistributionDistance(profiles.get(i), quantiles)*100000)+"  ");
		}
		return newMetaProfile;
	}
	
	public static void printMotifHits(){
		File file = new File(dir, "Oct4.tab");
		List<MACSPeakRegion> peaks = MACSPeakRegion.filterByPValue(MACSParser.parseMACSOutput(file.getAbsolutePath(), g),  default_p_value, Double.MAX_VALUE);  

		WeightMatrix motif = null;
		try {
			//Load motif
			int wmid = WeightMatrix.getWeightMatrixID(org.getDBID(), "Oct-4 (POU5F1)", "TRANSFAC 10.4, M01124");
			motif = WeightMatrix.getWeightMatrix(wmid);
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		int numHits=0, peaksWithHits=0, totalPeaks=0;
		WeightMatrixScorer scorer = new WeightMatrixScorer(motif);
		ArrayList<Point> hitPoints = new ArrayList<Point>();
		Map<Point, Point> motif2peak = new HashMap<Point, Point>();
		for(int i=0; i<peaks.size(); i++){
			totalPeaks++;
			MACSPeakRegion peak = peaks.get(i);
			Region r= peak.getPeak().expand(MOTIF_DISTANCE);
			WeightMatrixScoreProfile profiler = scorer.execute(r);
			boolean goodMotif=false;
			for(int z=0; z<r.getWidth(); z++){
				double currScore= profiler.getMaxScore(z);				
				if(currScore>=MOTIF_THRESH){
//					System.out.println(currScore);
					numHits++;
					goodMotif=true;
					Point point=new Point(g, r.getChrom(), r.getStart()+z+motif.length()/2);
					hitPoints.add(point);
					motif2peak.put(point, peak.getPeak());
				}
			}
			if(goodMotif){
				peaksWithHits++;
			}
		}
		double perc = (double)peaksWithHits/(double)totalPeaks;
		StringBuilder sb = new StringBuilder();
		sb.append((motif.name+" hits: "+numHits+" hits in "+peaksWithHits+" peaks from "+totalPeaks+" total peaks ("+perc+").\n"));
		for(int h=0; h<hitPoints.size(); h++){
			sb.append(hitPoints.get(h)+"\n");
		}
		writeToFile("Oct4_motif_"+MOTIF_DISTANCE+"_"+MOTIF_THRESH+".txt", sb.toString());
		sb = new StringBuilder();
		for (Point p:hitPoints){
			sb.append(motif2peak.get(p)+"\n");
		}
		writeToFile("Oct4_peak_"+MOTIF_DISTANCE+"_"+MOTIF_THRESH+".txt", sb.toString());
		
	}
	private static void peakMotifAnalysis(){
		File file;
		
//		file = new File("Oct4.MACS.sortPvalue.txt");
//		List<MACSPeakRegion> macsPeaks = MACSParser.parseMACSOutput(file.getAbsolutePath(), g);
//		ArrayList<Point> peaks_macs = new ArrayList<Point>();
//		for (MACSPeakRegion p: macsPeaks){
//			peaks_macs.add(p.getPeak());
//		}
		file = new File("CTCF_MM_1.txt");
		ArrayList<Point> peaks_mm = new ArrayList<Point>();
		FileReader in = null;
		BufferedReader bin = null;
		try {
			in = new FileReader(file);
			bin = new BufferedReader(in);
			String line;
			while((line = bin.readLine()) != null) { 
				line = line.trim();
	            String[] f=line.split("\t");
	            Region r = Region.fromString(g, f[0]);
	            peaks_mm.add(new Point(g, r.getChrom(), r.getStart()));
			}
		}
		catch(IOException ioex) {
			//logger.error("Error parsing file", ioex);
		}
		bindingPeakMotifOccurence(peaks_mm);
	}
	
	private static void bindingPeakMotifOccurence(ArrayList<Point> peaks){
		WeightMatrix motif = null;
		try {
			//Load motif
//			int wmid = WeightMatrix.getWeightMatrixID(org.getDBID(), "Oct-4 (POU5F1)", "TRANSFAC 10.4, M01124");
			int wmid = WeightMatrix.getWeightMatrixID(org.getDBID(), "CTCF", "Shaun");
			motif = WeightMatrix.getWeightMatrix(wmid);
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		WeightMatrixScorer scorer = new WeightMatrixScorer(motif);
		
		double[] motifThresholds = new double[]{0,4,8,12,16,20,24,28};
		// for each peak and motif threshold, the shortest distance between BS and motif
		int[][]distance = new int[motifThresholds.length][peaks.size()];
		for(int i=0; i<peaks.size(); i++){
			if (i % 100==0)
				System.out.println(i);
			Point peak = peaks.get(i);
			Region r= peak.expand(MOTIF_DISTANCE);
			WeightMatrixScoreProfile profiler = scorer.execute(r);
			for(int j=0;j<motifThresholds.length;j++){
				double threshold = motifThresholds[j];
				//search from BS outwards
				for(int z=0; z<=r.getWidth()/2; z++){
					double leftScore= profiler.getMaxScore(MOTIF_DISTANCE-z);
					double rightScore= profiler.getMaxScore(MOTIF_DISTANCE+z);				
					if(rightScore>=threshold){
						distance[j][i] = z;
						break;
					}
					if(leftScore>=threshold){
						distance[j][i] = -z;
						break;
					}
					// if motif score at this position is too small, and reach end of region
					if (z==r.getWidth()/2){
						distance[j][i] = 999;
					}
				}
			}
		}
		
		// output results
		StringBuilder sb = new StringBuilder();
		sb.append("#"+motif.name+" "+ motif.version+"\n").append("#Motif Score");
		for(int j=0;j<motifThresholds.length;j++){
			sb.append("\t"+motifThresholds[j]);
		}
		sb.append("\n");
		for(int i=0; i<peaks.size(); i++){
			sb.append(peaks.get(i).getLocationString());
			for(int j=0;j<motifThresholds.length;j++){
				sb.append("\t"+distance[j][i]);
			}
			sb.append("\n");
		}
		writeToFile("CTCF_MM_motif_distance_byScore_"+MOTIF_DISTANCE+".txt", sb.toString());
	}	

	// Given a sorted list of binding events, 
	// report percentage of top peaks that have a motif with some threshold at some bp away
	private static void bindEvent_motifOccurence(ArrayList<Point> peaks){
		WeightMatrix motif = null;
		try {
			//Load motif
//			int wmid = WeightMatrix.getWeightMatrixID(org.getDBID(), "Oct-4 (POU5F1)", "TRANSFAC 10.4, M01124");
			int wmid = WeightMatrix.getWeightMatrixID(org.getDBID(), "CTCF", "Shaun");
			motif = WeightMatrix.getWeightMatrix(wmid);
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		WeightMatrixScorer scorer = new WeightMatrixScorer(motif);
		
		double[] motifThresholds = new double[]{0,4,8,12,16,20,24,28};
		// for each peak and motif threshold, the shortest distance between BS and motif
		int[][]distance = new int[motifThresholds.length][peaks.size()];
		for(int i=0; i<peaks.size(); i++){
//			if (i % 100==0)
//				System.out.println(i);
			Point peak = peaks.get(i);
			Region r= peak.expand(MOTIF_DISTANCE);
			WeightMatrixScoreProfile profiler = scorer.execute(r);
			for(int j=0;j<motifThresholds.length;j++){
				double threshold = motifThresholds[j];
				//search from BS outwards, find the closest match
				for(int z=0; z<=r.getWidth()/2; z++){
					double leftScore= profiler.getMaxScore(MOTIF_DISTANCE-z);
					double rightScore= profiler.getMaxScore(MOTIF_DISTANCE+z);				
					if(rightScore>=threshold){
						distance[j][i] = z;
						break;
					}
					if(leftScore>=threshold){
						distance[j][i] = -z;
						break;
					}
					// if motif score at this position is too small, and reach end of region
					if (z==r.getWidth()/2){
						distance[j][i] = 999;
					}
				}
			}
		}
		
		int distanceCutoff[] = {10, 30, 50, 100, 150};
		
		// output results
		StringBuilder sb = new StringBuilder();
		sb.append("#"+motif.name+" "+ motif.version+"\n").append("#Motif Score");
		for(int j=0;j<motifThresholds.length;j++){
			sb.append("\t"+motifThresholds[j]);
		}
		sb.append("\n");
		for(int i=0; i<peaks.size(); i++){
			sb.append(peaks.get(i).getLocationString());
			for(int j=0;j<motifThresholds.length;j++){
				sb.append("\t"+distance[j][i]);
			}
			sb.append("\n");
		}
		writeToFile("CTCF_MM_motif_distance_byScore_"+MOTIF_DISTANCE+".txt", sb.toString());
	}	

	private static void readStartDistance(){
		Genome g=null;
		try {
			g = Organism.findGenome("mm8");
		} catch (NotFoundException e) {
			e.printStackTrace();
		}
		File file = new File("C:\\Data\\ES\\data\\mES_Oct4.bed");
		FileReader in = null;
		BufferedReader bin = null;
		List<StrandedRegion> reads = new ArrayList<StrandedRegion>();
		try {
			in = new FileReader(file);
			bin = new BufferedReader(in);
			String line;
			while((line = bin.readLine()) != null) { 
				line = line.trim();
	            String[] f=line.split("\t");
	            StrandedRegion sr = new StrandedRegion(g, f[0], Integer.parseInt(f[1]), Integer.parseInt(f[2]), f[5].charAt(0));
	            reads.add(sr);
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
		
		HashMap<String, List<Integer>> plus = new HashMap<String, List<Integer>>();
		HashMap<String, List<Integer>> minus = new HashMap<String, List<Integer>>();
		for (StrandedRegion sr: reads){
			String chr = sr.getChrom();
			if (sr.getStrand()=='+'){
				if (plus.containsKey(chr)){
					plus.get(chr).add(sr.getStart());
				}
				else{
					plus.put(chr, new ArrayList<Integer>());
				}
			}
			else{
				if (minus.containsKey(chr)){
					minus.get(chr).add(sr.getEnd());
				}
				else{
					minus.put(chr, new ArrayList<Integer>());
				}
			}
		}
		StringBuffer sb = new StringBuffer();
		for(String chr: plus.keySet()){
			List<Integer> starts=plus.get(chr);
			Collections.sort(starts);
			for (int i=1; i<starts.size();i++){
				sb.append(starts.get(i)-starts.get(i-1)).append("\t").append(chr).append("\t").append(starts.get(i-1)).append("\n");
			}
		}
		writeToFile("plus_start_distance.txt", sb.toString());
	}
	
    private static void writeToFile(String filename, String text) {
    	try {
    		File outFile = new File(filename);
    		PrintStream ps = new PrintStream(new FileOutputStream(outFile));
			ps.print(text);
			ps.flush();
			ps.close();
			System.out.println("File is written to "+outFile.getAbsolutePath());
		} catch (Exception e) {
			e.printStackTrace();
		}
    }   
    

}
class BEDRecord extends StrandedRegion{
	String uCode;
	public BEDRecord(Genome g, String chr, int start, int end, char strand, String uCode){
		super(g, chr, start, end, strand);
		this.uCode = uCode;
	}
	public String bedString(){
		return String.format("%s\t%d\t%d\t%s\t0\t%s\n", getChrom(), getStart(), getEnd(), uCode, getStrand());
	}
}


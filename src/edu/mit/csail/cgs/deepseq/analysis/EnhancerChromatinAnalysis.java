package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import edu.mit.csail.cgs.datasets.chipseq.ChipSeqLocator;
import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.BindingModel;
import edu.mit.csail.cgs.deepseq.DeepSeqExpt;
import edu.mit.csail.cgs.deepseq.StrandedBase;
import edu.mit.csail.cgs.deepseq.features.ComponentFeature;
import edu.mit.csail.cgs.deepseq.features.Feature;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.deepseq.utilities.ReadCache;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils.SISSRS_Event;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSPeak;
import edu.mit.csail.cgs.ewok.verbs.chipseq.MACSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.MACSPeakRegion;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.SetTools;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public class EnhancerChromatinAnalysis {
	static String h3k27ac_esc = "Wysocka hESC H3K27ac H9;prealigned_unique";
	static String h3k27me3_esc = "Wysocka hESC H3K27me3 H9;prealigned_unique";
	static String h3k4me3_esc = "Wysocka hESC H3K4me3 H9;prealigned_unique";
	static String h3k4me1_esc = "Wysocka hESC H3K4me1 H9;prealigned_unique";
	static String input_esc = "Wysocka hESC input H9;prealigned_unique";
	static String h3k27ac_nec = "Wysocka hNEC H3K27ac H9;prealigned_unique";
	static String h3k27me3_nec = "Wysocka hNEC H3K27me3 H9;prealigned_unique";
	static String h3k4me3_nec = "Wysocka hNEC H3K4me3 H9;prealigned_unique";
	static String h3k4me1_nec = "Wysocka hNEC H3K4me1 H9;prealigned_unique";	
	static String input_nec = "Wysocka hNEC input H9;prealigned_unique";
	
	private final int NOHIT_OFFSET = 999;
	
	private boolean isPreSorted = false;	
	private int windowSize = 1000;	
	private int rank = 1000;
	private int smooth_step = 30;
	private boolean sortByStrength = false;
	
	private Genome genome;
	private Organism org;
	private String[] args;
	private String outName="out";
	
	// each element in the list is for one ChIP-Seq method
	private ArrayList<String> markNames = new ArrayList<String>();
	private HashMap<String, Integer> markIDs = new HashMap<String, Integer>();
	private ArrayList<ReadCache> caches = null;
	private ArrayList<Point> events = new ArrayList<Point>();
	private double[][] profiles_I;
	private double[][] profiles_II;

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {
		EnhancerChromatinAnalysis analysis = new EnhancerChromatinAnalysis(args);

		int type = Args.parseInteger(args, "analysisType", 0);
		switch(type){
		case 0:analysis.findEnhancers();break;
		case 1:analysis.computeLikelihoodRatio();break;
		default: System.err.println("Unrecognize analysis type: "+type);
		}
	}
	
	EnhancerChromatinAnalysis(String[] args){
		this.args = args;
		ArgParser ap = new ArgParser(args);
		Set<String> flags = Args.parseFlags(args);		
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
	    	    
		// some parameters
		windowSize = Args.parseInteger(args, "windowSize", windowSize);
		rank = Args.parseInteger(args, "rank", rank);
		smooth_step = Args.parseInteger(args, "smooth", smooth_step);
		outName = Args.parseString(args,"out",outName);
		sortByStrength = flags.contains("ss");	// only for GPS
		isPreSorted = flags.contains("sorted");
		
		markNames.add("H3K4me1");
		markNames.add("H3K4me3");
		markNames.add("H3K27me3");
		markNames.add("H3K27ac");
		markNames.add("input");
		for (int i=0;i<markNames.size();i++){
			markIDs.put(markNames.get(i), i);
		}
		try{
			readPeakLists();
		}
		catch (IOException e){
			e.printStackTrace(System.err);
			System.exit(-1);
		}

		loadChIPSeqData();
	}
	/**
	 * Find 2 classes of enhancers using P300 binding events and chromatin marks
	 * @throws IOException
	 */
	// This method print 2 text files: class I and II enhancer coordinates
	private void findEnhancers() throws IOException {
		long tic = System.currentTimeMillis();
		
		ArrayList<Point> classI = new ArrayList<Point>();
		ArrayList<Point> classII = new ArrayList<Point>();
		for (Point p: events){
			Region r = p.expand(windowSize);
			if (isEnriched(r, "H3K4me1", "input", 10) && isEnriched(r, "H3K27ac", "input", 10) && (!isEnriched(r, "H3K4me3", "input", 30)))
				classI.add(p);
			if (isEnriched(r, "H3K4me1", "input", 10) && isEnriched(r, "H3K27me3", "input", 8) && (!isEnriched(r, "H3K27ac", "input", 10)) && (!isEnriched(r, "H3K4me3", "input", 30)))
				classII.add(p);
		}
		System.out.println(outName+"_enhancer_I:\t"+classI.size());
		StringBuilder sb = new StringBuilder();
		for (Point p:classI)
			sb.append(p.toString()).append("\n");
		CommonUtils.writeFile(outName+"_enhancer_I.txt", sb.toString());
		
		generateProfiles("I", classI, profiles_I, windowSize, windowSize);
		
		System.out.println(outName+"_enhancer_II:\t"+classII.size());
		sb = new StringBuilder();
		for (Point p:classII)
			sb.append(p.toString()).append("\n");
		CommonUtils.writeFile(outName+"_enhancer_II.txt", sb.toString());
		
		generateProfiles("II", classII, profiles_II, windowSize, windowSize);

		System.out.println("Done! " + CommonUtils.timeElapsed(tic));
	}

	
	private boolean isEnriched(Region region, String ip, String ctrl, double threshold){
		float ip_count = caches.get(markIDs.get(ip)).countHits(region);
		float ctrl_count = caches.get(markIDs.get(ctrl)).countHits(region);
		if (ctrl_count==0)
			ctrl_count = 1;
		return ip_count>threshold && ip_count/ctrl_count>2.5;
	}

	private ArrayList<Point> readPeakLists() throws IOException {
        String peakTag=null;
        for(String s : args)
        	if(s.contains("peakCall"))
        		peakTag=s;
    	
        // each tag represents a peak caller output file
    	String name="";
    	if(peakTag.startsWith("--peakCall")){
    		name = peakTag.replaceFirst("--peakCall", ""); 
    	}
    	List<File> files = Args.parseFileHandles(args, "peakCall"+name);
    	File file= null;
    	if (!files.isEmpty()){
    		file = files.get(0);
    	}

    	String filePath = file.getAbsolutePath();
		ArrayList<Point> peakPoints = new ArrayList<Point>(); 
    	if (name.contains("GPS")||name.contains("GEM")){
			List<GPSPeak> gpsPeaks = GPSParser.parseGPSOutput(filePath, genome);
			// sort by descending pValue (-log P-value)
			if (!isPreSorted){
				if (sortByStrength){
					Collections.sort(gpsPeaks, new Comparator<GPSPeak>(){
					    public int compare(GPSPeak o1, GPSPeak o2) {
					        return o1.compareByIPStrength(o2);
					    }
					});	
				}
				else{
					Collections.sort(gpsPeaks, new Comparator<GPSPeak>(){
					    public int compare(GPSPeak o1, GPSPeak o2) {
					        return o1.compareByPValue(o2);
					    }
					});
				}
			}
			for (GPSPeak p: gpsPeaks){
				if (p.getStrength()>30)				// p300 strength > 30
					peakPoints.add(p);
			}
    	}
    	else{
			peakPoints = CommonUtils.loadCgsPointFile(filePath, genome);
    	}  
    	System.out.println(peakPoints.size()+" p300 events are loaded.");
    	return peakPoints;
	}
	
	/*****************************************************
	 * Load ChIP-Seq data (adapted from KPPMixture.java)
	 * ***************************************************/
	private void loadChIPSeqData(){
		boolean fromReadDB = true;
		final int MAXREAD = 1000000;
		String cellStr = "hESC";
		
		this.caches = new ArrayList<ReadCache>();
		for (int i=0;i<markNames.size();i++){
			String expt = "Wysocka "+cellStr+" "+markNames.get(i)+" H9";
			String align = "prealigned_unique";
			List<ChipSeqLocator> locators = new ArrayList<ChipSeqLocator>();
			locators.add(new ChipSeqLocator(expt, align));
			DeepSeqExpt ip = new DeepSeqExpt(genome, locators, "readdb", -1);
			ReadCache ipCache = new ReadCache(genome, markNames.get(i));
			caches.add(ipCache);

			// cache sorted start positions and counts of all positions
			long tic = System.currentTimeMillis();
			System.out.print("\nLoading  "+ipCache.getName()+" data from ReadDB ... \t");
			for (String chrom: genome.getChromList()){
				// load  data for this chromosome.
				int length = genome.getChromLength(chrom);
				Region wholeChrom = new Region(genome, chrom, 0, length-1);
				int count = ip.countHits(wholeChrom);
				ArrayList<Region> chunks = new ArrayList<Region>();
				// if there are too many reads in a chrom, read smaller chunks
				if (count>MAXREAD){
					int chunkNum = count/MAXREAD*2+1;
					int chunkLength = length/chunkNum;
					int start = 0;
					while (start<=length){
						int end = Math.min(length, start+chunkLength-1);
						Region r = new Region(genome, chrom, start, end);
						start = end+1;
						chunks.add(r);
					}
				}else
					chunks.add(wholeChrom);

				for (Region chunk: chunks){
					Pair<ArrayList<Integer>,ArrayList<Float>> hits = ip.loadStrandedBaseCounts(chunk, '+');
					ipCache.addHits(chrom, '+', hits.car(), hits.cdr());
					hits = ip.loadStrandedBaseCounts(chunk, '-');
					ipCache.addHits(chrom, '-', hits.car(), hits.cdr());
				}
			} // for each chrom

			ipCache.populateArrays(true);
			ip.closeLoaders();
			ip=null;
			System.gc();

			if (fromReadDB){
				System.out.println(CommonUtils.timeElapsed(tic));
			}
		}
	}

	private void generateProfiles(String classType, ArrayList<Point> coords, double[][] profiles, int left, int right){
		int width = left+right+1;
        // data for all conditions
		int numMarks = 4;
        double[][] profile_plus=new double[numMarks][width];
        double[][] profile_minus=new double[numMarks][width];
		// sum the read profiles from all qualified binding events for updating model
		for (int c=0;c<numMarks;c++){
			ReadCache ip = caches.get(c);
			for (int i=0;i<Math.min(rank, coords.size());i++){
				Point p = coords.get(i);
				Region region = p.expand(0).expand(left, right);
				List<StrandedBase> bases = ip.getStrandedBases(region, '+');
				for (StrandedBase base:bases){
					profile_plus[c][base.getCoordinate()-region.getStart()]+=base.getCount();
				}
				region = p.expand(0).expand(right, left);
				bases = ip.getStrandedBases(region, '-');
				for (StrandedBase base:bases){
					profile_minus[c][region.getEnd()-base.getCoordinate()]+=base.getCount();
				}
			}

			// smooth the model profile
			if (smooth_step>0){
				profile_plus[c] = StatUtil.cubicSpline(profile_plus[c], smooth_step, smooth_step);
				profile_minus[c] = StatUtil.cubicSpline(profile_minus[c], smooth_step, smooth_step);
			}

			//compare models from 2 strands, shift binding position if needed
			BindingModel.minKL_Shift(profile_plus[c], profile_minus[c] );
		}
		
		// output the read density
		profiles=new double[numMarks][width];
		double[] scaleFactors = new double[numMarks];
		for (int c=0;c<numMarks;c++){
			scaleFactors[c] = caches.get(c).getHitCount();
		}
		StatUtil.mutate_normalize(scaleFactors);
		
		StringBuilder sb = new StringBuilder();
		for (int i=0;i<width;i++){
			sb.append(i-left).append("\t");
			for (int c=0;c<numMarks;c++){
				double sum = (profile_plus[c][i]+profile_minus[c][i])/4.0/scaleFactors[c];
				sum = sum>=0?sum:2.0E-300;
				profiles[c][i] = sum;
				sb.append(sum).append("\t");
			}
			sb.append("\n");
		}
		
		// normalize profiles across all marks
		double total = 0;
		for (int i=0;i<width;i++){
			for (int c=0;c<numMarks;c++){
				total += profiles[c][i];
			}
		}	
		if (total>0){
			for (int i=0;i<width;i++){
				for (int c=0;c<numMarks;c++){
					profiles[c][i] /= total;
				}
			}
		}
		CommonUtils.writeFile(outName+"_"+classType+"_profiles.txt", sb.toString());
	}
	
	// Compute the likelihood ratio of a p300 event region fitting to class I and II enhancer profiles
	private void computeLikelihoodRatio() throws IOException {
		long tic = System.currentTimeMillis();
		readPeakLists();
		
		BindingModel h3k4me1_I = getProfile("Read_Distribution_h3k4me1_esc_I_1000.txt");
		BindingModel h3k4me3_I = getProfile("Read_Distribution_h3k4me3_esc_I_1000.txt");
		BindingModel h3k27me3_I = getProfile("Read_Distribution_h3k27me3_esc_I_1000.txt");
		BindingModel h3k27ac_I = getProfile("Read_Distribution_h3k27ac_esc_I_1000.txt");

		BindingModel h3k4me1_II = getProfile("Read_Distribution_h3k4me1_esc_II_1000.txt");
		BindingModel h3k4me3_II = getProfile("Read_Distribution_h3k4me3_esc_II_1000.txt");
		BindingModel h3k27me3_II = getProfile("Read_Distribution_h3k27me3_esc_II_1000.txt");
		BindingModel h3k27ac_II = getProfile("Read_Distribution_h3k27ac_esc_II_1000.txt");

		StringBuilder sb = new StringBuilder();
		sb.append("Event\t\tk4me1I\tk4me1II\tk4me3I\tk4me3II\tk27me3I\tk27m3II\tk27acI\tk27acII\tOverall\tK27m3ac\n");
//		for (Point p: events){
//			sb.append(p.toString()+"\t");
//			Pair<Double, Double> pair1 = getLikelihood(p, h3k4me1, h3k4me1_I, h3k4me1_II);
//			sb.append(String.format("%.1f\t%.1f\t", pair1.car(), pair1.cdr()));
//			Pair<Double, Double> pair2 = getLikelihood(p, h3k4me3, h3k4me3_I, h3k4me3_II);
//			sb.append(String.format("%.1f\t%.1f\t", pair2.car(), pair2.cdr()));
//			Pair<Double, Double> pair3 = getLikelihood(p, h3k27me3, h3k27me3_I, h3k27me3_II);
//			sb.append(String.format("%.1f\t%.1f\t", pair3.car(), pair3.cdr()));
//			Pair<Double, Double> pair4 = getLikelihood(p, h3k27ac, h3k27ac_I, h3k27ac_II);
//			sb.append(String.format("%.1f\t%.1f\t", pair4.car(), pair4.cdr()));
//			double total = pair1.car()-pair1.cdr()+pair2.car()-pair2.cdr()+pair3.car()-pair3.cdr()+pair4.car()-pair4.cdr();
//			double k27 = pair3.car()-pair3.cdr()+pair4.car()-pair4.cdr();
//			sb.append(String.format("%.1f\t%.1f\n", total, k27));
//		}
		CommonUtils.writeFile("enhancer_all_likelihood_ratio.txt", sb.toString());
		
		System.out.println("Done! " + CommonUtils.timeElapsed(tic));
	}
	private BindingModel getProfile(String modelFile){
		//Load Binding Model
		File pFile = new File(modelFile);
		if(!pFile.isFile()){
			System.err.println("\nCannot find read distribution file!");
			System.exit(1);
		}
        return new BindingModel(pFile);
	}
	private Pair<Double, Double> getLikelihood(Point p, DeepSeqExpt ip, BindingModel classI, BindingModel classII){
		double ll_I = 0;
		double ll_II = 0;
		int pos = p.getLocation();
		Pair<ArrayList<Integer>,ArrayList<Float>> pair = ip.loadStrandedBaseCounts(p.expand(1500), '+');
		for (int i=0;i<pair.car().size();i++){
			ll_I += Math.log10(classI.probability(pair.car().get(i)-pos))*pair.cdr().get(i);
			ll_II += Math.log10(classII.probability(pair.car().get(i)-pos))*pair.cdr().get(i);
		}
		return new Pair<Double, Double>(ll_I, ll_II);
	}
	


}

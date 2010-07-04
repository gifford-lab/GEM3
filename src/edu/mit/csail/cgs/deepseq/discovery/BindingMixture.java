package edu.mit.csail.cgs.deepseq.discovery;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.Iterator;
import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.BindingModel;
import edu.mit.csail.cgs.deepseq.DeepSeqExpt;
import edu.mit.csail.cgs.deepseq.StrandedBase;
import edu.mit.csail.cgs.deepseq.features.*;
import edu.mit.csail.cgs.deepseq.multicond.MultiIndependentMixtureCounts;
import edu.mit.csail.cgs.deepseq.utilities.ReadCache;
import edu.mit.csail.cgs.ewok.verbs.chipseq.MACSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.MACSPeakRegion;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.Utils;
import edu.mit.csail.cgs.utils.probability.NormalDistribution;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public class BindingMixture extends MultiConditionFeatureFinder{
	/****************************
	 * Constants
	 ***************************/
	private final static boolean LOG_ALL=false;
	
	// width for smoothing a read (used as a stddev for creating the Gaussian kernel of probability)
	public final static int READ_KERNEL_ESTIMATOR_WIDTH = 25;
	// the range to scan a peak if we know position from EM result
	protected final static int SCAN_RANGE = 50;
	
	//Maximum region (in bp) considered for running EM
	protected final int MAX_REGION_LEN=200000;
	//IP/Control Fold change threshold
	protected final int IP_CTRL_FOLD = 2;
	
	//Maximum number of components that EM can handle efficiently
	protected final int MAX_NUM_COMPONENTS=1000;
	protected final int OPTIMAL_NUM_COMPONENTS=100;
	protected final int INIT_RESOLUTION=5;
	// Maximum number of EM iterations
	protected final int MAX_EM_ITER=10000;
	protected final int MAX_EM_ML_ITER=1;
	//Run EM up until <tt>ML_ITER</tt> without using sparse prior
	protected final int ML_ITER=10;
	//Run EM up until <tt>ANNEALING_ITER</tt> with smaller alpha based on the current iteration
	protected final int ANNEALING_ITER=200;

	//EM convergence between the likelihood of the current and the previous step
	protected final double EM_CONVERGENCE = 1e-5;
	protected final double EM_ML_CONVERGENCE = 1e-8;
	
	protected final int num_top_mfold_feats = 1000;
	
	protected final double cpenalty = 1;
	protected final double lambda = 0;
	
	// Max number of reads to load from DB
	private final static int MAXREAD = 1000000;
    private final static int WINDOW_SIZE_FACTOR = 3;
    // true:  eliminated component in batch, as long as matching criteria in EM derivation
    // false: eliminate only the worse case components, re-distribute the reads of eliminated component
	private final static boolean BATCH_ELIMINATION = false;	
	private final static boolean SPLINE_SMOOTH = true;
    private final static boolean SMART_SPACING = true;		// dynamically determine init comopent spacing
    private final static boolean MAKE_HARD_ASSIGNMENT=false;
	private StringBuilder pi_sb = new StringBuilder();
	
	private boolean development_mode = false;
	private boolean do_model_selection=true;
	private boolean linear_model_expansion=false;
    private boolean print_mixing_probabilities=false;
    private boolean use_multi_event = false;
	private boolean TF_binding = true;
	private boolean pre_artifact_filter=false;	
	private boolean post_artifact_filter=false;	
    private int max_hit_per_bp = -1;
    private double q_value_threshold = 2.0;
    private double alpha_factor = 3.0;
    private int top_event_percentile = 50;
    private double mappable_genome_length = 2.08E9; // mouse genome 
    private double sparseness=6.0;
    private int first_lambda_region_width  =  1000;
    private int second_lambda_region_width =  5000;
    private int third_lambda_region_width  = 10000;
    private boolean use_dynamic_sparseness = true;
    private boolean use_internal_em_train  = true;
    private boolean use_scanPeak  = true;
    private boolean needToCleanBases = false;
    private int bmverbose=1;		// BindingMixture verbose mode
	private int needle_height_factor = 2;	// threshold to call a base as needle, (height/expected height) 
	private double needle_hitCount_fraction = 0.1;	//threshold to call a region as tower, (hit_needles / hit_total)
	
	//Binding model representing the empirical distribution of a binding event 
	private BindingModel model;
	private int modelRange;			// max length of one side, 250bp
	private int modelWidth; 		// width of empirical read distribution, 500bp
	private int windowSize;
	
	// Max number of reads on a base position, to filter out towers and needles.
	private int max_HitCount_per_base = 3; 
	
	//Gaussian kernel for density estimation
	protected double[] gaussian;
	// peak strength --> <mean, std>
	private HashMap<Integer, Pair<Double, Double>> shapeParameters;

	protected boolean doScanning=true;
	//Spacing between the components
	protected int componentSpacing=1;
	
	/****************
	 * Data
	 ****************/
	private boolean wholeGenomeDataLoaded = false;
	// memory cache to store all the read data, loaded from DB or file
	protected ArrayList<Pair<ReadCache, ReadCache>> caches;
	protected ArrayList<String> conditionNames = new ArrayList<String>();
	// Do we have matched control data?
	protected boolean controlDataExist=false;
	// Ratio of each pair of IP/Ctrl for all conditions
	protected double[] ratio_total;
	protected double[] ratio_non_specific_total;
	protected double shift_size;
	/****************
	 * Prediction
	 ****************/	
	// List of all EM predictions 
	protected ArrayList<ComponentFeature> allFeatures;
	protected ArrayList<Feature> insignificantFeatures;	// does not pass statistical significance
	protected List<Feature>[] condSignalFeats; //Signal channel features in each condition
	
	//Regions to be evaluated by EM, estimated whole genome, or specified by a file
	protected ArrayList<Region> restrictRegions = new ArrayList<Region>();
	// Regions that have towers (many reads stack on same base positions)
	protected ArrayList<Region> towerRegions = new ArrayList<Region>();
	protected ArrayList<Float> towerStrength = new ArrayList<Float>();
	/**
	 * <tt>HashMap</tt> containing the single event regions. <br>
	 * Each <tt>key</tt> contains the (single event) region of interest and the corresponding <tt>value</tt> 
	 * the smaller sub-region (within the <tt>key</tt> region) containing the binding location
	 * The purpose is that if we re-evaluate with new empirical distribution, we do not re-run EM for unary event.
	 * The assumption is that changing empirical distribution will not have changed the existence of events, 
	 * it will only modify the position of the event slightly, so we can just do the scanning.
	 */
	protected HashMap<Region, Point> singleEventRegions=new HashMap<Region, Point>();

	//Number of reads for a specified region of the IP experiment for each condition
	protected double[] sigHitCounts;
	//(Total) number of reads for a specified region of the IP experiment (across all conditions)
	protected double totalSigCount=0;
	//Components representing the binding events in the analyzed window
	protected ArrayList<BindingComponent> components = new ArrayList<BindingComponent>();
	//Number of non zero components of specified region
	protected double nonZeroComponentNum=0;
	// Maximum number of components determined by the data
	protected int componentMax;
	
	/** Stores the reads as <tt>StrandedBase</tt> structures.
	 * key:chrom
	 * value:{first dimension: IP or CTRL channel, second:condition, third:reads as StrandedBases ('+' and '-' are mixed)}
	 */
	private Map<String,List<ArrayList<List<StrandedBase>>>> chrom_signals = new HashMap<String,List<ArrayList<List<StrandedBase>>>>();
	
	/** key:chrom, value:{first dimension: IP or CTRL channel, second: condition, 
	 * third:strand, fourth:five primes of this chromosome for that condition and channel} */
	private Map<String, int[][][][]> chromFivePrimes; 
	
	/** Total IP counts for each chromosome (summed over all conditions) */
	private Map<String, Integer> totalIPCounts;
	
	/** Read counts for each chromosome for each condition. 
	 * List 0 => IP channel, List 1 => CTRL channel */
	private Map<String, List<List<Integer>>> condHitCounts;

	
    private StringBuilder config = new StringBuilder();
	private StringBuilder log_all_msg = new StringBuilder();
	private FileWriter logFileWriter;
	private boolean reportProgress = true;
	private StringBuilder needleReport = new StringBuilder();

	/**
	 * This constructor will get most parameters from property file 
	 * It is usually called from ChipSeqAnalyzer.
	 */
	public BindingMixture(Genome g, ArrayList<Pair<DeepSeqExpt,DeepSeqExpt>> expts,
						  ArrayList<String> conditionNames, 
						  String[] args){
		super (args, g, expts);
		try{
			logFileWriter = new FileWriter("GPS_Log.txt", true); //append
			logFileWriter.write("\n==============================================\n");
			logFileWriter.write(getDateTime());
			logFileWriter.flush();
		} 
		catch (IOException e) {
			e.printStackTrace();
		}
		
		/* ***************************************************
		 * Load parameters and properties
		 * ***************************************************/
		StringBuffer sb = new StringBuffer();
		sb.append("\nParameters>> ");
		for (String arg:args){
			sb.append(arg).append(" ");
		}
		log(1, sb.toString());
		
		/* ***************************************************
		 * Load Binding Model, empirical distribution
		 * ***************************************************/
		String modelFile = Args.parseString(args, "d", null);	// read distribution file

		commonInit(modelFile);

        if (SPLINE_SMOOTH)
        	model.smooth(BindingModel.SMOOTHING_STEPSIZE);
		model.printToFile(outName+"_0_Read_distribution.txt");
		
		// Required input parameter
    	mappable_genome_length = Args.parseDouble(args, "s", 2.08E9);	// size of mappable genome
    	
    	// Optional input parameter
    	q_value_threshold = Args.parseDouble(args, "q", 2.0);	// q-value
    	sparseness = Args.parseDouble(args, "a", 6.0);	// alpha parameter for sparse prior
    	alpha_factor = Args.parseDouble(args, "alpha_factor", 3.0); // denominator in calculating alpha value
    	max_hit_per_bp = Args.parseInteger(args, "max_hit_per_bp", -1);
    	top_event_percentile = Args.parseInteger(args, "top_event_percentile", 50);
    	needle_height_factor = Args.parseInteger(args, "needle_height_factor", 2);
    	needle_hitCount_fraction = Args.parseDouble(args, "needle_hitCount_fraction", 0.1);
    	
    	// flags
    	Set<String> flags = Args.parseFlags(args);
    	// default as false, need the flag to turn it on
    	print_mixing_probabilities = flags.contains("print_mixing_probabilities");
    	linear_model_expansion = flags.contains( "linear_model_expansion");
    	use_multi_event = flags.contains("refine_using_multi_event");
    	needToCleanBases = flags.contains("needToCleanBases");
    	pre_artifact_filter =  flags.contains("pre_artifact_filter");
    	post_artifact_filter = flags.contains("post_artifact_filter");
    	development_mode = flags.contains("development_mode");
    	// default as true, need the opposite flag to turn it off
    	use_dynamic_sparseness = ! flags.contains( "fa"); // fix alpha parameter
    	TF_binding = ! flags.contains("non_punctate_binding");
    	reportProgress =! flags.contains("no_report_progress");
    	use_internal_em_train = ! flags.contains( "use_multi_condition_em_train");
    	use_scanPeak = ! flags.contains( "do_not_scanPeak");
    	boolean loadWholeGenome = ! flags.contains( "loadRegionOnly");
    	do_model_selection = !flags.contains( "no_model_selection");    	
		/* **************************************************
		 * Determine the focus regions to run EM
 		 * It can be specified as a file from command line.
		 * If no file pre-specified, estimate the enrichedRegions after loading the data.
		 * We do not want to run EM on regions with little number of reads
		 * **************************************************/	
    	String focusFormat = Args.parseString(args, "focusFormat", null);
    	String focusFile = Args.parseString(args, "focusFile", null);
    	
    	if (focusFormat!=null){
	    	// take candidate regions from MACS output, expand point to regions, merge overlaps
	    	if (focusFormat.equals("MACS")){
	    		File peakFlie = new File(focusFile);
	    		List<MACSPeakRegion> peaks = MACSParser.parseMACSOutput(peakFlie.getAbsolutePath(), gen);
	    		setRegions(peaks);
	    	}
	    	// take candidate regions from output of statistical peak finder by Shaun, 
	    	// will expand point to regions, merge overlaps
	    	else if (focusFormat.equals("StatPeak")){
	    		setRegions(focusFile, false);    		
	    	}
	    	// take candidate regions from previous analysis of mixture model
	    	// Do not expand the regions, precisely defined already
	    	else if (focusFormat.equals("Regions")){
	        	setRegions(focusFile, true);    		
	    	}
	    	else if (focusFormat.equals("RegionsToMerge")){
	        	setRegions(focusFile, false);    		
	    	}
    	}
    	
		log(1, "\nSorting reads and selecting regions for analysis ...");

		/* ***************************************************
		 * ChIP-Seq data
		 * ***************************************************/
		this.caches = new ArrayList<Pair<ReadCache, ReadCache>>();
		this.numConditions = experiments.size();
		this.conditionNames = conditionNames;
		ComponentFeature.setConditionNames(conditionNames);
		condSignalFeats = new ArrayList[this.numConditions];
		for(int c = 0; c < this.numConditions; c++) { condSignalFeats[c] = new ArrayList<Feature>(); }
		long tic = System.currentTimeMillis();
		for (int i=0;i<experiments.size();i++){
			Pair<DeepSeqExpt, DeepSeqExpt> pair = experiments.get(i);
			DeepSeqExpt ip = pair.car();
			DeepSeqExpt ctrl = pair.cdr();
			if(ctrl.getHitCount()>0){
	    		controlDataExist = true;
	    	}
			ReadCache ipCache = new ReadCache(gen, conditionNames.get(i)+"_IP  ");
			ReadCache ctrlCache = null;
			if (controlDataExist)
				ctrlCache = new ReadCache(gen, conditionNames.get(i)+"_CTRL");
			this.caches.add(new Pair<ReadCache, ReadCache>(ipCache, ctrlCache));
			
			// cache sorted start positions and counts of all positions
			if (ip.isFromReadDB()){		// load from read DB
				System.out.println("Loading data from ReadDB ... \n");
				if (loadWholeGenome){		// if whole genome
					
					for (String chrom: gen.getChromList()){									
						// load  data for this chromosome.
						int length = gen.getChromLength(chrom);
						Region wholeChrom = new Region(gen, chrom, 0, length-1);
						int count = Math.max(ip.countHits(wholeChrom), ctrl.countHits(wholeChrom));
						ArrayList<Region> chunks = new ArrayList<Region>();
						// if there are too many reads in a chrom, read smaller chunks
						if (count>MAXREAD){
							int chunkNum = count/MAXREAD+1;
							int chunkLength = length/chunkNum;
							int start = 0;
							while (start<=length){
								int end = Math.min(length, start+chunkLength-1);
								Region r = new Region(gen, chrom, start, end);
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
							if (controlDataExist){
								hits = ctrl.loadStrandedBaseCounts(chunk, '+');
								ctrlCache.addHits(chrom, '+', hits.car(), hits.cdr());
								hits = ctrl.loadStrandedBaseCounts(chunk, '-');
								ctrlCache.addHits(chrom, '-', hits.car(), hits.cdr());
							}
						}
					} // for each chrom
					wholeGenomeDataLoaded = true;
				} //if (loadWholeGenome==1){
				else{	// if loadWholeGenome==0, only load these reads
					for (Region region: restrictRegions){
						Pair<ArrayList<Integer>,ArrayList<Float>> hits = ip.loadStrandedBaseCounts(region, '+');
						ipCache.addHits(region.getChrom(), '+', hits.car(), hits.cdr());
						hits = ip.loadStrandedBaseCounts(region, '-');
						ipCache.addHits(region.getChrom(), '-', hits.car(), hits.cdr());
						if (controlDataExist){
							hits = ctrl.loadStrandedBaseCounts(region, '+');
							ctrlCache.addHits(region.getChrom(), '+', hits.car(), hits.cdr());
							hits = ctrl.loadStrandedBaseCounts(region, '-');
							ctrlCache.addHits(region.getChrom(), '-', hits.car(), hits.cdr());
						}
					}					
				}
				ipCache.populateArrays();
				ipCache.displayStats();
//				ipCache.printBin500Counts();
//				ipCache.printBinCounts();
				if (controlDataExist){
					ctrlCache.populateArrays();
					ctrlCache.displayStats();
//					ctrlCache.printBinCounts();
//					ctrlCache.printBin500Counts();
				}
				
			}
			else if (ip.isFromFile()){		// load from File
				ipCache.addAllFivePrimes(ip.getAllStarts(), ip.getReadLen());
				ipCache.populateArrays();
				if (controlDataExist){
					ctrlCache.addAllFivePrimes(ctrl.getAllStarts(), ctrl.getReadLen());
					ctrlCache.populateArrays();
				}
				wholeGenomeDataLoaded = true;
			}
//			max_HitCount_per_base = Math.max(max_HitCount_per_base, ipCache.getMaxHitPerBP(DUPLICATE_READ_FRACTION));
//			System.out.println("max_HitCount_per_base = "+max_HitCount_per_base);
			System.gc();
		} //end for each condition
		
		// If data was loaded from ReadDB, clean up the database connection 
		if(experiments.get(0).car().isFromReadDB()) {
			for(Pair<DeepSeqExpt,DeepSeqExpt> e : experiments){
				e.car().closeLoaders();
				e.cdr().closeLoaders();
			}
			System.out.println("Finish loading data from ReadDB, " + timeElapsed(tic));			
		}
		experiments = null;
		System.gc();
		
    	ratio_total=new double[numConditions];
    	ratio_non_specific_total = new double[numConditions];
        sigHitCounts=new double[numConditions];
        seqwin=100;
        if (focusFormat == null){		// estimate some parameters if whole genome data
	        for(int i=0; i<caches.size(); i++){
				Pair<ReadCache,ReadCache> e = caches.get(i);
				double ipCount = e.car().getHitCount();
				// estimate max hit count per BP
				// if user supply using max_hit_per_bp, use it
				// if not supplied, take the max of default (3) and possion expected count
				if (max_hit_per_bp==-1){
					int maxPerBP = calcHitCount_per_BP(ipCount, 1e-9);
					max_HitCount_per_base = Math.max(max_HitCount_per_base, maxPerBP);
				}
				else
					max_HitCount_per_base = max_hit_per_bp;
				
				// estimate IP/Ctrl ratio from all reads
				double ctrlCount = 0;
				if (e.cdr()!=null){
					ctrlCount = e.cdr().getHitCount();
					ratio_total[i]=ipCount/ctrlCount;
					System.out.println(String.format(conditionNames.get(i)+"\tIP: %.0f\tCtrl: %.0f\t IP/Ctrl: %.2f", ipCount, ctrlCount, ratio_total[i]));
				}
	        }
        }
        else{	// want to analyze only specified regions, set default
        	if (max_hit_per_bp!=-1)
				max_HitCount_per_base = max_hit_per_bp;
        	for(int i=0; i<caches.size(); i++){
        		ratio_total[i]=1;
        		ratio_non_specific_total[i]=1;
        	}
        }
        if (development_mode)
        	log(1, "\nmax_HitCount_per_base = "+max_HitCount_per_base);
        
    	// if no focus list, directly estimate candidate regions from data
		if (focusFormat==null){
    		setRegions(selectEnrichedRegions());    		
		}
		if (development_mode)
			printNoneZeroRegions(true);
		log(1, "\n"+restrictRegions.size()+" regions loaded for analysis.");
		
		log(2, "BindingMixture initialized. "+numConditions+" conditions.");
	}

	/**
	 * This constructor is only called from Robustness analysis
	 * The read data of each sub-sampling will be suppled with simpleExecute() call
	 */
	public BindingMixture(String[] args, boolean doScanning){
		super (args);

		// ChIP-Seq data
		this.doScanning = doScanning;
	    controlDataExist = false;
	    numConditions = 1;
	    
	    String modelFile = Args.parseString(args, "read_distribution", null);
		commonInit(modelFile);
	           
        // the configuration of mixture model
        addConfigString("Do ML scanning\t", doScanning);
        addConfigString("Batch elimination", BATCH_ELIMINATION);
        addConfigString("Smart inital event spacing", SMART_SPACING);
//        addConfigString("Make hard assignment", properties.make_hard_assignment);
        System.out.println(getConfigString());
	}
	
	private void commonInit(String modelFile){
		//Load Binding Model
		File pFile = new File(modelFile);
		if(!pFile.isFile()){System.err.println("Cannot find binding model file");System.exit(1);}
        model = new BindingModel(pFile);
        modelWidth = model.getWidth();
        modelRange = model.getRange();
        windowSize = modelWidth * WINDOW_SIZE_FACTOR;
        
        //init the Guassian kernel prob. for smoothing the read profile of called events
		gaussian = new double[modelWidth]; 
		NormalDistribution gaussianDist = new NormalDistribution(0, READ_KERNEL_ESTIMATOR_WIDTH*READ_KERNEL_ESTIMATOR_WIDTH);
		for (int i=0;i<gaussian.length;i++)
			gaussian[i]=gaussianDist.calcProbability((double)i);
	}

	/**
	 * Calls EM training on the Mixture Model.
	 * In this implementation, we have chosen to run EM independently on each testRegion.
	 */
	public List<Feature> execute() {
		long tic = System.currentTimeMillis();
		// if we have predicted events from previous round, setup the restrictRegions
		// because when we update the model, the model range might changed.
		if (allFeatures!=null && (!allFeatures.isEmpty())){
			Collections.sort(allFeatures);
			ArrayList<Region> potentialRegions = new ArrayList<Region>();
			for (ComponentFeature cf:allFeatures){
				potentialRegions.add(cf.getPeak().expand(0));
			}
			// expand with modelRange, and merge overlapped regions ==> Refined enriched regions
			this.restrictRegions=mergeRegions(potentialRegions, true);
		}
		
        signalFeatures.clear();	
		ArrayList<ComponentFeature> compFeatures = new ArrayList<ComponentFeature>();

		int totalRegionCount = restrictRegions.size();
		System.out.println("\nRunning EM for each of "+totalRegionCount+" regions, please wait ...");
		int displayStep = 5000;
		if (totalRegionCount<10000)
			displayStep = 1000;
		if (totalRegionCount<1000)
			displayStep = 100;
		if (totalRegionCount<100)
			displayStep = 10;	
		if (totalRegionCount<20)
			displayStep = 2;	
		System.out.println("(Progress will be reported in steps of "+displayStep+" regions).\n");
		
		//for each test region
		for (int j=0;j<restrictRegions.size();j++) {
			Region rr = restrictRegions.get(j);
			// Cut long regions into windowSize(4kb) sliding window (500bp overlap) to analyze
			ArrayList<Region> windows = new ArrayList<Region>();
			if (rr.getWidth()<=windowSize)
				windows.add(rr);
			else{
				int lastEnd = rr.getStart()+windowSize-1;
				windows.add(new Region(rr.getGenome(), rr.getChrom(), rr.getStart(), lastEnd));
				while(lastEnd<rr.getEnd()){
					int newStart = lastEnd+1-modelWidth;	//overlap length = modelWidth 
					lastEnd = Math.min(newStart+windowSize-1, rr.getEnd());
					windows.add(new Region(rr.getGenome(), rr.getChrom(), newStart, lastEnd));
				}
			}
			
			// process for each window
			ArrayList<BindingComponent> comps= new ArrayList<BindingComponent>();
			for (Region w : windows){
				ArrayList<BindingComponent> result = analyzeWindow(w);
				if (result!=null){
					comps.addAll(result);
				}
			} 
			
			/* ****************************************************************
			 * fix sliding window boundary effect
			 */						
			if (windows.size()>1){
				Collections.sort(comps);
				ArrayList<Region> rs = new ArrayList<Region>();
				for (BindingComponent m:comps){
					rs.add(m.getLocation().expand(0));
				}
				// expand with modelRange padding, and merge overlapped regions 
				// ==> Refined enriched regions
				rs=mergeRegions(rs, true); 
				for (Region r:rs){					
					if (r.getWidth()>2*modelRange+1){ // for non-unary events (for unary event, one of the windows should contain it in full, no need to evaluate again)
						// the regions in rs includes the influence paddings, remove it here
						Region tightRegion = new Region(r.getGenome(), r.getChrom(), r.getStart()+modelRange, r.getEnd()-modelRange); 
						for (int i=0;i<windows.size()-1;i++){
//							int end = windows.get(i).getEnd();
							Region boundary = windows.get(i).getOverlap(windows.get(i+1));
//							Region boundary = new Region(r.getGenome(), r.getChrom(), end-modelWidth, end);
							// if the predicted component overlaps with boundary of sliding window
							if (boundary.overlaps(tightRegion)){
								// remove old 
								ArrayList<BindingComponent> toRemove = new ArrayList<BindingComponent>();
								for (BindingComponent m:comps){
									if (tightRegion.overlaps(m.getLocation().expand(0)))
										toRemove.add(m);
								}
								comps.removeAll(toRemove);
								// reprocess the refined region
								ArrayList<BindingComponent> result = analyzeWindow(r);
								if (result!=null){
									comps.addAll(result);
								}
								break;	// break from testing more boundary
							}
						}
					}
				}
			}
			
			/* ****************************************************************
			 * refine unary events and collect all the events as features
			 * this is last step because fixing boundary may result in some new unary events
			 */			
			singleEventRegions.clear();
			Collections.sort(comps);
			ArrayList<Region> rs = new ArrayList<Region>();
			for (BindingComponent m:comps){
				rs.add(m.getLocation().expand(0));
			}
			// expand with modelRange padding, and merge overlapped regions 
			// ==> Refined enriched regions
			rs=mergeRegions(rs, true); 
			for (Region r:rs){					
				ArrayList<BindingComponent> bs = new ArrayList<BindingComponent>();
				if (r.getWidth()==2*modelRange+1){ // refine unary event
					BindingComponent b;
					// Scan Peak unary region
					if(use_scanPeak) {
						Point p = r.getMidpoint();
						b = scanPeak(p);
						if (b==null){
							//						System.err.println("No read to scan! "+r.toString());
							continue;
						}
						b.setEMPosition(p);
						for (BindingComponent m:comps){
							if (p.getLocation()==m.getLocation().getLocation()) {
								b.setAlpha(m.getAlpha());	// inherit alpha value from previous EM component

								if(!use_internal_em_train) {
									for(int c = 0; c < numConditions; c++) { b.setConditionBeta(c, m.getConditionBeta(c)); }
								}

								break;   // Once you found a previous component on the same location as b, there is no need to search all the others
							}
						}
					}
					// Run EM again on the unary region and take the component with the maximum strength
					else {
						ArrayList<BindingComponent> bl = analyzeWindow(r);
						if(bl == null || bl.size() == 0)      { continue; }
						else if(bl.size() == 1) { b = bl.get(0); }
						else {
							double[] compStrength = new double[bl.size()];
							for(int k = 0; k < compStrength.length; k++) { compStrength[k] = bl.get(k).getMixProb(); }
							Pair<Double, TreeSet<Integer>> max_maxCompIdx = StatUtil.findMax(compStrength);
							b = bl.get(max_maxCompIdx.cdr().first());
						}
						
					}
					
					bs.add(b);
					compFeatures.addAll(callFeatures(bs));
					singleEventRegions.put(r, b.getLocation());
				}
				else{	// for joint events
					for (BindingComponent m:comps){
						if (r.overlaps(m.getLocation().expand(0)))
							bs.add(m);
					}
					compFeatures.addAll(callFeatures(bs));
				}		
			}
			
			if ((j+1) % displayStep==0 && reportProgress)
				System.out.println((j+1)+"\t/"+totalRegionCount+"\t"+timeElapsed(tic));

		}// end of for (Region rr : restrictRegions) 
		System.out.println(totalRegionCount+"\t/"+totalRegionCount+"\t"+timeElapsed(tic));
		
		/* ********************************************************
		 * merge nearby tower regions, filter events in or at the edge of tower regions 
		 */
		towerRegions = mergeRegions(towerRegions, false);
		towerStrength.clear();
		for (Region tower: towerRegions){
			float allCount=0;
			for(Pair<ReadCache,ReadCache> e : caches){
				List<StrandedBase> bases_p= e.car().getStrandedBases(tower, '+');  // reads of the current region - IP channel
				List<StrandedBase> bases_m= e.car().getStrandedBases(tower, '-');  // reads of the current region - IP channel
				allCount += StrandedBase.countBaseHits(bases_p)+ StrandedBase.countBaseHits(bases_m);
			}
			towerStrength.add(allCount);
		}
		
		
		ArrayList<ComponentFeature> toBeRemoved=new ArrayList<ComponentFeature>();
//		for (ComponentFeature cf:compFeatures){
//			for (Region tower: towerRegions){
//				if (tower.overlaps(cf.getPeak().expand(modelRange))) { toBeRemoved.add(cf); break; }
//			}
//		}
		
		Region[] compfeatRegions = new Region[compFeatures.size()];
		for(int i = 0; i < compfeatRegions.length; i++) { compfeatRegions[i] = compFeatures.get(i).getPeak().expand(modelRange); }
		boolean[] isTower = Region.overlap(compfeatRegions, towerRegions.toArray(new Region[0]));
		for(int i = 0; i < isTower.length; i++)
			if(isTower[i]) { toBeRemoved.add(compFeatures.get(i)); }
		
		
		for (ComponentFeature cf:toBeRemoved){			
			while(compFeatures.remove(cf));  // remove all the features that are towers
			System.out.println(cf.getPeak().toString()+"\tPrediction in tower region.");
		}
		
		/* ********************************************************
		 * refine the specific regions that contains binding events
		 */
		Collections.sort(compFeatures);
		ArrayList<Region> refinedRegions = new ArrayList<Region>();
		for (ComponentFeature cf:compFeatures){
			refinedRegions.add(cf.getPeak().expand(0));
		}
		// expand with modelRange, and merge overlapped regions ==> Refined enriched regions
		this.restrictRegions=mergeRegions(refinedRegions, true); 
		
		if (development_mode){
			printNoneZeroRegions(false);		// print refined regions
			printTowerRegions();
			writeFile(outName+"_needles.txt", needleReport.toString());
		}
		log(1, "\nFinish calling events: "+timeElapsed(tic)+"\n");
		
		// post processing
		postEMProcessing(compFeatures);
		for(int c = 0; c < numConditions; c++) { condSignalFeats[c] = condPostFiltering(signalFeatures, c); }
		
		return signalFeatures;
	}// end of execute method
	
	
	private ArrayList<BindingComponent> analyzeWindow(Region w){
		
		ArrayList<List<StrandedBase>> signals = loadBasesInWindow(w, "IP", this.pre_artifact_filter);		
		if (signals==null || signals.isEmpty())
			return null;
		
		//initialize count arrays
		totalSigCount=0;
		for(int c=0; c<signals.size(); c++){
			sigHitCounts[c]=StrandedBase.countBaseHits(signals.get(c));
			totalSigCount+=sigHitCounts[c];
		}
		// if less than significant number of reads
		if (totalSigCount<sparseness)
			return null;
		
		// We want to run EM only for potential overlapping regions
		// If after first round, we are sure the region contains unary event, we will just scan for peak
		if (!singleEventRegions.containsKey(w)) {
			log(3, "\nWindow:\t"+w.getLocationString()+" ("+w.getWidth()+" bp, "+(int)totalSigCount+" reads)");

			// dynamically determine an alpha value for this sliding window
			double alpha = sparseness;
			if (use_dynamic_sparseness)
				alpha = Math.max(estimateAlpha2(w, signals), sparseness);
			
			HashMap<Integer, double[][]> responsibilities = null;
			
			// Choose either Yuchun's internal EM or Temporal Coupling
			if(use_internal_em_train) {
				//Run EM and increase resolution
				initializeComponents(w, numConditions);
				int lastResolution;
		
				while(nonZeroComponentNum>0){
					lastResolution = componentSpacing;
					// EM learning, components list will only contains non-zero components
					responsibilities = EMTrain(signals, alpha);
					// log(4, componentSpacing+" bp\t"+(int)nonZeroComponents+" components.");

					// increase resolution
					updateComponentResolution(w, numConditions, lastResolution);
					if(componentSpacing==lastResolution)
						break;
				} 	// end of while (resolution)

			}
			else {
				// Run MultiIndependentMixture (Temporal Coupling)
				
				// set up params for TC to be consistent with those defined in BingingMixture
				MultiIndependentMixtureCounts mim = new MultiIndependentMixtureCounts();
				mim.set_trainVersion(1);
				mim.set_alpha(alpha);
				mim.set_maxIters(MAX_EM_ITER);
				mim.set_ML_maxIters(ML_ITER);
//				mim.set_anneal_maxIters(ANNEALING_ITER);
//				mim.set_convergence_thres(EM_CONVERGENCE);
				mim.set_anneal_maxIters(120);
				mim.set_convergence_thres(1e-8);
				mim.set_decreasing_thres(-1e-6);
				mim.set_check_increased(true);
				mim.set_isBatchElimOn(BATCH_ELIMINATION);
				
				int[][] pos     = new int[numConditions][];  // position is relative to region's start
				int[][] count   = new int[numConditions][];
				char[][] strand = new char[numConditions][];
				for(int t = 0; t < numConditions; t++) {
					pos[t]    = new int[signals.get(t).size()];
					count[t]  = new int[signals.get(t).size()];
					strand[t] = new char[signals.get(t).size()];
					
					for(int k = 0; k < signals.get(t).size(); k++) {
						pos[t][k]    = signals.get(t).get(k).getCoordinate() - w.getStart();  // position is relative to region's start
						count[t][k]  = (int) signals.get(t).get(k).getCount();
						strand[t][k] = signals.get(t).get(k).getStrand();
					}
				}
				
				int[] compPos;
				if(w.getWidth() > 2*modelRange+1) { // non-unary events
					componentSpacing = 3;
					compPos = setComponentPositions(w.getWidth(), pos, modelRange, componentSpacing);
				}
				else { // unary events
					componentSpacing = 2;
					compPos = setComponentPositions(w.getWidth(), pos, modelRange, componentSpacing);
				}
				
				int M  = compPos.length;           

				mim.exec(pos, count, strand, compPos, model);
				
				double[]   glob_pi = mim.get_glob_prior_weight();
				double[][] pi = mim.get_cond_prior_weight();
				responsibilities = (HashMap<Integer, double[][]>) mim.get_resp();
				
				boolean atLeast_one_nonZero_comp = false;
				for(int j = 0; j < M; j++) {
					if(glob_pi[j] > 0.0) { 
						if(numConditions > 1) {
							for(int t = 0; t < numConditions; t++)
								if(pi[t][j] > 0.0) { atLeast_one_nonZero_comp = true; break; }
						}
						else { atLeast_one_nonZero_comp = true; break; }
						 
					}
					if(atLeast_one_nonZero_comp) { break; }
				}
				
				if(!atLeast_one_nonZero_comp) { components.clear(); }
				
				List<BindingComponent> nonZeroComponents = new ArrayList<BindingComponent>(); 
				for(int j = 0; j < M; j++) {
					if(glob_pi[j] > 0.0) {
						// Multi-condition experiment				
						if(numConditions > 1) {
							BindingComponent comp = null;
							int numNonZeroCondEvents = 0;
							for(int t = 0; t < numConditions; t++) {
								if(pi[t][j] > 0.0) {
									numNonZeroCondEvents++;
									// if there is at least one non zero condition-event initialize the component
									if(numNonZeroCondEvents == 1) {
										comp = new BindingComponent(model, new Point(w.getGenome(), w.getChrom(), w.getStart() + compPos[j]), numConditions);
										comp.setMixProb(glob_pi[j]);
										comp.setAlpha(alpha);
										comp.setOld_index(j);
									}
								}
								if(comp != null) { comp.setConditionBeta(t, pi[t][j]); }
							}
							
							if(comp != null) { nonZeroComponents.add(comp); }
						}
						// Single-condition experiment
						else { 
							BindingComponent comp = new BindingComponent(model, new Point(w.getGenome(), w.getChrom(), w.getStart() + compPos[j]), numConditions);
							comp.setMixProb(glob_pi[j]);
							comp.setAlpha(alpha);
							comp.setOld_index(j);
							nonZeroComponents.add(comp);
						}
					}
				}
				
				components = (ArrayList<BindingComponent>) nonZeroComponents;	
				nonZeroComponentNum  = components.size();
			}// end of else condition (running it with TC)
			
			if (nonZeroComponentNum==0)	return null;
			setComponentResponsibilities(signals, responsibilities);
		} 	// if not single event region

		else{// if single event region, just scan it
			BindingComponent peak = scanPeak(singleEventRegions.get(w));
			components = new ArrayList<BindingComponent>();
			components.add(peak);
		}
		return components;
	}
	
	private ArrayList<List<StrandedBase>> loadBasesInWindow(Region w){
		return loadBasesInWindow(w, "IP", true);
	}
	
	// Load the reads from all conditions in the region
	// filter out needles, but still return reads for EM processing
	// if many needles, record as a tower
	private ArrayList<List<StrandedBase>> loadBasesInWindow(Region w, String channel, boolean needToFilterBases){
		boolean isTowerRegion = false;
//		float allCount = 0;
		ArrayList<List<StrandedBase>> signals = new ArrayList<List<StrandedBase>>();
		if(channel.equalsIgnoreCase("IP")) {
			//Load each condition's read hits
			for(Pair<ReadCache,ReadCache> e : caches){
				List<StrandedBase> bases_p= e.car().getStrandedBases(w, '+');  // reads of the current region - IP channel
				List<StrandedBase> bases_m= e.car().getStrandedBases(w, '-');  // reads of the current region - IP channel
				
				if(needToCleanBases) {
					cleanBases(bases_p);
					cleanBases(bases_m);
				}
				
				if(needToFilterBases) {
					boolean isTower_p = filterBases(bases_p, e.car(), w, '+');
					boolean isTower_m = filterBases(bases_m, e.car(), w, '-');
					// if it is tower
					if (isTower_p || isTower_m){
						isTowerRegion = true;
					}					
				}
				bases_p.addAll(bases_m);
				signals.add(bases_p);
			} // for loop
			
			if (isTowerRegion){
//				System.err.println("Tower region!\t" + w.getLocationString()+"\tHit#:"+ allCount);
				towerRegions.add(w);
//				towerStrength.add(allCount);
				return null;
			}			
		}
		else if(channel.equalsIgnoreCase("CTRL")) {
			//Load each condition's read hits
			for(Pair<ReadCache,ReadCache> e : caches){
				List<StrandedBase> bases_p= e.cdr().getStrandedBases(w, '+');  // reads of the current region - Ctrl channel
				List<StrandedBase> bases_m= e.cdr().getStrandedBases(w, '-');  // reads of the current region - Ctrl channel
				bases_p.addAll(bases_m);
				signals.add(bases_p);
			} // for loop

		}
		else
			throw new IllegalArgumentException("The only valid values for channel is either IP or CTRL.");

		return signals;
	}
	
	/** Sets the count of a base to the maximum hit count per base (<tt>max_HitCount_per_base</tt>) */
	private void cleanBases(List<StrandedBase> bases) {
		for (StrandedBase base:bases){
			if(base.getCount() > max_HitCount_per_base)
				base.setCount(max_HitCount_per_base);
		}
	}//end of cleanBases method
	
	// if only few needles, reset the baseCount for needles
	// if many needles, report as tower
	private boolean filterBases(List<StrandedBase> bases, ReadCache cache, Region region, char strand){
		if (bases.isEmpty())
			return false;
		int summit = model.getSummit();
		int min = model.getMin();
		int max = model.getMax();
		float topProb = (float) model.probability(summit);
		boolean isTower = false;
		boolean[] isNeedles = new boolean[bases.size()];
		float needlesHitCount = 0;
		
		for (int i=0;i<bases.size();i++){
			StrandedBase base = bases.get(i);
			float count = base.getCount();
			// only check bases that pass the threshold
			if (count > max_HitCount_per_base){
				// assume this position corresponds to summit of empirical distribution.
				Region r = null;
				int b = base.getCoordinate();
				if (strand=='+')
					r = new Region(gen, region.getChrom(), b-(summit-min), b+(max-summit));
				else
					r = new Region(gen, region.getChrom(), b-(max-summit), b+(summit-min));
				List<StrandedBase> neighbours = cache.getStrandedBases(r, strand);
				// excluding this base for calculation
				float totalCount = (StrandedBase.countBaseHits(neighbours)-count) / (1-topProb);
				float expectedCount = totalCount * topProb;
				// if count of this base is higher than needle_height_factor (4) times expected count from emp.dist.
				if (count> expectedCount * needle_height_factor){
					isNeedles[i]=true;
					needlesHitCount+= expectedCount;
					base.setCount(expectedCount);
				}
//				// find the most likely binding event that contains this base position 
//				// using empirical distribution, only based on this condition.
//				int b = base.getCoordinate();
//				Region scanRegion = new Point(gen, chrom, b).expand(modelRange);
//				Region dataRegion = scanRegion.expand(modelRange, modelRange);
//				ArrayList<List<StrandedBase>> data = new ArrayList<List<StrandedBase>>();
//				data.add(cache.getUnstrandedBases(dataRegion));
//				BindingComponent event = scanPeak(data, scanRegion);
//				if (event!=null){
//					Region eventRegion = event.getLocation().expand(modelRange);
//					List<StrandedBase> neighbours = cache.getStrandedBases(eventRegion, strand);
//					float totalCount = StrandedBase.countBaseHits(neighbours);
//					// if count of this base is higher than 2 times expected from emp.dist.
//					if (count> totalCount*model.probability(b-event.getLocation().getLocation())*2){
//						isNeedles[i]=true;
//						base.setCount(max_HitCount_per_base);
//					}
//				}
//				else{
////					System.out.println("event==null: "+dataRegion.toString());
//					isNeedles[i]=true;
//					base.setCount(max_HitCount_per_base);
//				}
			}
		}
		int needleCount = 0;
		for (boolean isNeedle:isNeedles)
			if (isNeedle)
				needleCount ++;
		int width = bases.get(bases.size()-1).getCoordinate()-bases.get(0).getCoordinate();
		float correctedRatio = needlesHitCount/StrandedBase.countBaseHits(bases);
		if (correctedRatio > needle_hitCount_fraction)
			isTower = true;
		// more than 10 needles --> tower 
		// allow 5 needles per 500bp, 1 needle per 8 positions
//		if (needleCount>=10 ||(needleCount>5.0*width/modelWidth && needleCount>bases.size()/8.0))
//			isTower = true;
		if (needleCount>0 && cache.getName().endsWith(":IP  "))
			needleReport.append(region.toString()+"\t"+strand+"\t"+
					needleCount+"\t"+width+"\t"+bases.size()+"\t"+isTower+"\t"+correctedRatio+"\n");
		return isTower;
	}
	// evaluate if a predicted peak is a tower
	

	private void setComponentResponsibilities(ArrayList<List<StrandedBase>> signals, HashMap<Integer, double[][]> responsibilities) {
		// Set responsibility profile for each component (for kernel density and KL calculation)
		for(int j=0;j<components.size();j++){
			BindingComponent comp = components.get(j);
			int jr = comp.getOld_index();
			for(int c=0; c<numConditions; c++){
				List<StrandedBase> bases = signals.get(c);
				double[][] rc = responsibilities.get(c);

				// store binding profile (read responsibilities in c condition) of this component
				double[] profile_plus = new double[modelWidth];
				double[] profile_minus = new double[modelWidth];
				for(int i=0;i<bases.size();i++){
					StrandedBase base = bases.get(i);
					if (rc[i][jr]>0){
						try{
						comp.setResponsibility(c, base, rc[i][jr]);
						if (base.getStrand()=='+')
							profile_plus[base.getCoordinate()-comp.getLocation().getLocation()-model.getMin()]=rc[i][jr]*base.getCount();
						else{
							profile_minus[comp.getLocation().getLocation()-base.getCoordinate()-model.getMin()]=rc[i][jr]*base.getCount();
						}
						}
						catch (Exception e){
							System.err.println(comp.toString()+"\t"+base.getStrand()+"\t"+base.getCoordinate());
							e.printStackTrace(System.err);
						}
					}
				}

				comp.setReadProfile(c, profile_plus, '+');
				comp.setReadProfile(c, profile_minus, '-');
			}
		}
	}

	// dynamically determine an alpha value for this sliding window based on IP reads
	// find a 500bp(modelWidth) region with maximum read counts, set alpha=sqrt(maxCount)/3  

	private double estimateAlpha2(Region window, ArrayList<List<StrandedBase>> signals){
		float maxCount = 0;
		int left = Math.abs(model.getMin());
		int right = Math.abs(model.getMax());
		int edge = Math.min(left, right);
		for (int i=window.getStart()+edge;i<=window.getEnd()-edge; i++){	
			// i is each possible event position
			float count=0;
			for (int c=0;c<signals.size();c++){
				List<StrandedBase> bases = signals.get(c);
				for (StrandedBase b: bases){
					int coor = b.getCoordinate();
					if (b.getStrand()=='+'){
						if (coor>i-left && coor<i+right)
							count += b.getCount();
					}
					else{
						if (coor>i-right && coor<i+left)
							count += b.getCount();
					}
				}
			}
			if (maxCount<count)
				maxCount = count;
		}
//		double alphaEstimate = maxCount/ALPHA_FACTOR;
		double alphaEstimate = Math.sqrt(maxCount)/alpha_factor;
		return alphaEstimate;
	}

	/**
	 * Calls EM training on the Mixture Model.
	 * This is called by SimulatedReadAnalysis for batch analysis in a region
	 */
	public ArrayList<ComponentFeature> simpleExecute(List<StrandedBase> bases, Region r) {
		components = new ArrayList<BindingComponent>();
		totalSigCount=StrandedBase.countBaseHits(bases);
		ArrayList<List<StrandedBase>> signals = new ArrayList<List<StrandedBase>>();
		signals.add(bases);

		HashMap<Integer, double[][]> responsibilities = null;
		if (!doScanning){
			//Run EM and increase resolution
			initializeComponents(r, numConditions );
			while(nonZeroComponentNum>0){
				double alpha = Math.max(Math.sqrt(StrandedBase.countBaseHits(bases))/alpha_factor, sparseness);
				responsibilities = EMTrain(signals, alpha);
				if(componentSpacing==1)
					break;
				updateComponentResolution(r, 1, componentSpacing);
			} 	
			setComponentResponsibilities(signals, responsibilities);
		}
		else{		// scan it
			BindingComponent peak = scanPeak(signals, r);
			components.add(peak);
		}
	
		return callFeatures(components);
	}// end of simpleExecute method
	
	private List<Feature> condPostFiltering(List<Feature> signifFeats, int cond){
		List<Feature> condFeats = new ArrayList<Feature>();
		for (Feature sf:signifFeats){
			boolean significant = false;
			// If TF binding keep the condition-wise component that exceeds the threshold
			if (TF_binding) {
				if(((ComponentFeature)sf).getQValueLog10(cond)>q_value_threshold &&
						((ComponentFeature)sf).getEventReadCounts(cond)>sparseness)
					significant = true;
			}
				
			// if not TF binding, it is Chromatin data, keep all components
			else { significant = true; }
			
			if (significant && ((ComponentFeature)sf).getCondBetas()[cond] > 0.0 )
				condFeats.add((ComponentFeature)sf);
			else
				significant = false;
		}
		
		return condFeats;
	}//end of condPostFiltering method
	
	/**
	 * It performs statistical tests for determining significant peaks         <br>
	 * If a control is used, the current method is used. Otherwise, MACS proposed method is used.
	 * @param compFeatures
	 */
	private void postEMProcessing(ArrayList<ComponentFeature> compFeatures){
		// use the refined regions to count non-specific reads
		countNonSpecificReads(compFeatures);
				
		// calculate p-values with or without control
		evaluateConfidence(compFeatures);
		
		// sort features for final output, by location
		Collections.sort(compFeatures);
		allFeatures = compFeatures;
		
		insignificantFeatures = new ArrayList<Feature>();
		for (ComponentFeature cf: compFeatures){
			boolean significant = false;
			// for multi-condition, at least be significant in one condition
			// The read count test is for each condition
			if (TF_binding)
				for (int cond=0; cond<conditionNames.size(); cond++){
					if(cf.getQValueLog10(cond)>q_value_threshold &&
							cf.getEventReadCounts(cond)>sparseness)
						significant = true;
				}
			else	// if not TF binding, it is Chromatin data, keep all components
				significant = true;
			
			if (significant)
				signalFeatures.add(cf);
			else
				insignificantFeatures.add(cf);				
		}
		
		log(1, "Significant events: "+signalFeatures.size()+"\tInsignificant: "+insignificantFeatures.size()+"\n");
	}

	
	/**
	 * EM training
	 * 
	 * Returns the responsibility of EM training
	 * 
	 * purely matrix/array operations
	 * After EM training, components list will only contains non-zero components
	 */
	protected HashMap<Integer, double[][]>  EMTrain(ArrayList<List<StrandedBase>> signals, double alpha){
		int numComp = components.size();

		// Number of reads in different conditions may be different,
		// use HashMap, instead of array to store the values
		HashMap<Integer, double[]>   counts= new HashMap<Integer, double[]>();	// Hit Count
		HashMap<Integer, double[][]> h= new HashMap<Integer, double[][]>(); 	// H function
		HashMap<Integer, double[][]> r= new HashMap<Integer, double[][]>();		// Responsibility
		HashMap<Integer, double[]>   b= new HashMap<Integer, double[]>();		// Beta
		double[] pi = new double[numComp];										// Pi
		
		// nz_start, nz_end to mark the non-zero (for H) range of i for each j. condition specific.
		// used to reduce number of unnecessary computation
		// assume the reads are sorted by genome location
		HashMap<Integer, int[]> i_start= new HashMap<Integer, int[]>();		// start read index for NZ Hij
		HashMap<Integer, int[]> i_end = new HashMap<Integer, int[]>();			// end read index for NZ Hij

		StringBuilder sb = new StringBuilder();
		for(int j=0;j<numComp;j++){
			BindingComponent comp = components.get(j);
			pi[j]= comp.getMixProb();
			if(print_mixing_probabilities)
				sb.append(comp.getLocation().getLocation()+"\t");
		}
		if(print_mixing_probabilities)
			pi_sb.append(sb.toString().trim()+"\n");
		
		// create an expanded read distribution model by extending the tails with min and max prob.
		// The purpose is to have a wide model to cover every data point. Thus every data point 
		// will at least get some minimal probability, so that we can compare likelihood when we 
		// eliminating components (model selection)
		
		// figure out the largest distance between bases and components
		int minusEnd = 0;
		int plusEnd = 0;
		for(int c=0; c<numConditions; c++){
			List<StrandedBase> bases = signals.get(c);
			int numBases = bases.size();	
			for(int i=0;i<numBases;i++){
				StrandedBase base = bases.get(i);
				for(int j=0;j<numComp;j++){
					BindingComponent comp = components.get(j);
					int dist = base.getStrand()=='+' ? base.getCoordinate()-comp.getLocation().getLocation(): comp.getLocation().getLocation()-base.getCoordinate();
					if (dist<0)
						minusEnd = Math.min(dist, minusEnd);
					else
						plusEnd = Math.max(dist, plusEnd);
				}
			}
		}
		// expand the model to cover the ends we found
		BindingModel expandedModel;
		if (!linear_model_expansion)
			expandedModel = model.getExponentialExpandedModel(model.getMin()-minusEnd, 
				plusEnd - model.getMax());
		else
			expandedModel = model.getLinearExpandedModel(model.getMin()-minusEnd, 
					plusEnd - model.getMax());
		//expandedModel.printToFile("Expanded_Read_Distribution.txt");

		for(int c=0; c<numConditions; c++){
			List<StrandedBase> bases = signals.get(c);
			int numBases = bases.size();
			double[] bc= new double[numComp];
			for(int j=0;j<numComp;j++)
				bc[j]=1.0/numConditions;
			b.put(c, bc);
			
			double[] countc= new double[numBases];
			for(int i=0;i<numBases;i++)
				countc[i]=bases.get(i).getCount();
			counts.put(c, countc);
			
			double[][] hc= new double[numBases][numComp];
			for(int i=0;i<numBases;i++){
				StrandedBase base = bases.get(i);
				for(int j=0;j<numComp;j++){
					BindingComponent comp = components.get(j);
					int dist = base.getStrand()=='+' ? base.getCoordinate()-comp.getLocation().getLocation(): comp.getLocation().getLocation()-base.getCoordinate();
					hc[i][j]=expandedModel.probability(dist);
				}
			}
			h.put(c, hc);			

			//TODO: read has + and - strand, H prob. is not completely sorted
			int[] c_i_start=new int[numComp];
			int[] c_i_end=new int[numComp];
			for(int j=0;j<numComp;j++){
				for(int i=0;i<numBases;i++){
					if (hc[i][j]>0){
						c_i_start[j]=i;
						break;
					}
				}
				for(int k=numBases-1;k>=0;k--){
					if (hc[k][j]>0){
						c_i_end[j]=k;
						break;
					}
				}
			}
			i_start.put(c, c_i_start);
			i_end.put(c, c_i_end);

			//Initial Semi-E-step: initialize unnormalized responsibilities
			double[][] rc= new double[numBases][numComp];
			for(int j=0;j<numComp;j++){
				for(int i=c_i_start[j];i<c_i_end[j];i++){
					rc[i][j] = hc[i][j]*pi[j]*bc[j];
				}	
			}
			r.put(c,rc);
		}
		log(5, "\n"+componentSpacing+" bp:\t");
		//////////
		// Run EM steps
		//////////
		EM_MAP(counts, h, r, b, pi, i_start, i_end, alpha);
	
		//////////
		// re-assign EM result back to the objects
		//////////
		if (this.nonZeroComponentNum==0){
//			for(int j=0;j<numComp;j++)
//				components.get(j).setMixProb(0);
			components.clear();
			return r;
		}
		ArrayList<BindingComponent> nonZeroComponents = new ArrayList<BindingComponent>();
		for(int j=0;j<numComp;j++){
			if (pi[j]>0){		// non-zero components
				BindingComponent comp = components.get(j);
				comp.setMixProb(pi[j]);
				comp.setAlpha(alpha);
				for(int c=0; c<numConditions; c++){
					double[] bc = b.get(c);
					comp.setConditionBeta(c, bc[j]);
					comp.setOld_index(j);
				}
				nonZeroComponents.add(comp);
			}
		}
		components = nonZeroComponents;
		
		// limit the range of the responsibilities, and re-normalize it
		// this is to reverse the the effect of using a expandedModel
		for(int c=0; c<numConditions; c++){
			double[][]rc = r.get(c);
			List<StrandedBase> bases = signals.get(c);
			int numBases = bases.size();	
			for(int i=0;i<numBases;i++){
				StrandedBase base = bases.get(i);
				double respSum = 0;
				for(int j=0;j<components.size();j++){
					BindingComponent comp = components.get(j);
					int oldIndex = comp.getOld_index();
					int dist = base.getStrand()=='+' ? base.getCoordinate()-comp.getLocation().getLocation(): comp.getLocation().getLocation()-base.getCoordinate();
					if (dist<model.getMin() || dist>model.getMax())
						rc[i][oldIndex] = 0;		// limit the range of the responsibilities
					respSum += rc[i][oldIndex];
				}
				if (respSum!=0)
					for(int j=0;j<numComp;j++){	// re-normalize
						rc[i][j] /= respSum ;
					}
			}
		}

		return r;
	}//end of EMTrain method 


	/** 
	 * core EM steps, sparse prior (component elimination), multi-condition 
	 */
	private void EM_MAP(  	HashMap<Integer, double[]>   counts,
							HashMap<Integer, double[][]> h, 
							HashMap<Integer, double[][]> r, 
							HashMap<Integer, double[]>   b,
							double[] pi,
							HashMap<Integer, int[]> i_start,
							HashMap<Integer, int[]> i_end,
							double alpha)
	{
		int numComp = pi.length;
		
		ArrayList<EM_State> models = new  ArrayList<EM_State> ();
		
		// print initial PI values
		if(print_mixing_probabilities){
			StringBuilder sb = new StringBuilder();
			for (int i=0; i<pi.length; i++){
				sb.append(pi[i]+"\t");
			}
			pi_sb.append(sb.toString().trim()+"\n");
			appendEMMixProbFile(componentSpacing);
		}
		
		double lastLAP=0, LAP=0; // log posterior prob
		int t=0;
		double currAlpha = alpha/2;
		log(5, (int)nonZeroComponentNum+" ");
		//Run EM while not converged
		for(t=0; t<MAX_EM_ITER ; t++){
			
			lastLAP=LAP;

			//////////
			//Semi-E-step (presumes unnormalized responsibilities have been calculated)
			//Simply normalize responsibilities here
			//////////
			for(int c=0; c<numConditions; c++){
				double[][] rc = r.get(c);
				int numReads = rc.length;
				double totalResp[] = new double[numReads];
				int c_i_start[]=i_start.get(c);
				int c_i_end[] = i_end.get(c);
				// init
				for(int i=0;i<numReads;i++){
					totalResp[i] = 0;
				}
				// sum
				for(int j=0;j<numComp;j++)
					if (pi[j]!=0)
						for(int i=c_i_start[j];i<c_i_end[j];i++)
							totalResp[i] += rc[i][j];
				// normalize
				for(int j=0;j<numComp;j++)
					if (pi[j]!=0)
						for(int i=c_i_start[j];i<c_i_end[j];i++)
							if (totalResp[i]>0)
								rc[i][j] = rc[i][j]/totalResp[i];
			}


			//////////
			//M-step
			//////////

			//Pi
			nonZeroComponentNum=0;
			// ML: standard EM 
			if (t<=ML_ITER){
            	for(int j=0;j<numComp;j++){
                    double r_sum=0;
                    for(int c=0; c<numConditions; c++){
        				double[][] rc = r.get(c);
        				double[]counts_c = counts.get(c);
        				int c_i_start[]=i_start.get(c);
        				int c_i_end[] = i_end.get(c);
        				for(int i=c_i_start[j];i<c_i_end[j];i++){
                    		r_sum += rc[i][j]*counts_c[i];
                    	}
                    }
                    pi[j]=r_sum;    // standard EM ML
            	}
            }
			// MAP: EM with sparce prior
            else {
            	if (BATCH_ELIMINATION){
                	if (t >ML_ITER && t <= ANNEALING_ITER)
                		currAlpha = alpha * (t-ML_ITER)/(ANNEALING_ITER-ML_ITER);
                    else
                    	currAlpha = alpha;
                	// adjust alpha for coarse resolution
//                	if (currAlpha>1 && INITIAL_SPARCENESS)
//                		currAlpha = Math.max(1, currAlpha/componentSpacing); 
                	
            		// Batch elimination, all the component with support less than
                	// currAlpha will be eliminated altogether      
                	// Risk: all the components might be eliminated pre-maturely
                	// 		 at the initial EM runs (has not reached proper assignment)
                	// 		 when pi[j] is too small (too many components)
                	// 		 That is the reason of having the annealing cycles.
	                for(int j=0;j<numComp;j++){
	                	if (pi[j]!=0){
		                	double r_sum=0;
		                    for(int c=0; c<numConditions; c++){
		        				double[][] rc = r.get(c);
		        				double[]counts_c = counts.get(c);
		        				int c_i_start[]=i_start.get(c);
		        				int c_i_end[] = i_end.get(c);
		        				for(int i=c_i_start[j];i<c_i_end[j];i++){
		                    		r_sum += rc[i][j]*counts_c[i];
		                    	}
		                    }
	
		                    // component elimination
		                    pi[j]=Math.max(0,r_sum-currAlpha);
	
	                		// if component prob becomes 0, clear responsibility
		                    if (pi[j]==0){
		                		for(int c=0; c<numConditions; c++){
		            				double[][] rc = r.get(c);
		            				int c_i_start[]=i_start.get(c);
		            				int c_i_end[] = i_end.get(c);
		            				for(int i=c_i_start[j];i<c_i_end[j];i++){
		                        		rc[i][j] = 0;
		                        	}
		            				c_i_start[j]=0;
		            				c_i_end[j]=0;
		                        }
		                	}
	                	}
	                }
            	}
            	else{ 	//batch_elimination is false
            		// eliminate only the worst cases
            		// redistribute to boost neighbor component by next round of EM
                   	// in this case, we do not need annealing schedule for alpha
                   	
            		double r_sum[]=new double[numComp];
                    for(int j=0;j<numComp;j++){
	                	if (pi[j]!=0){
		                	for(int c=0; c<numConditions; c++){
		        				double[][] rc = r.get(c);
		        				double[]counts_c = counts.get(c);
		        				int c_i_start[]=i_start.get(c);
		        				int c_i_end[] = i_end.get(c);
		        				for(int i=c_i_start[j];i<c_i_end[j];i++)
		                    		r_sum[j] += rc[i][j]*counts_c[i];
		                    }
	                	}
	                	else
	                		// set a large value to prevent this component being selected in findMin()
	                		r_sum[j]=9999; 
                	}
                    // find the worst component
                    Pair<Double, TreeSet<Integer>> min = StatUtil.findMin(r_sum);
                    if (min.car()-currAlpha>0){ 
                    	// no component to be eliminated, update pi(j)
                    	 for(int j=0;j<numComp;j++)
                    		 if (pi[j]!=0)
                    			 pi[j]=r_sum[j]-currAlpha;

//                     	System.out.println(t+":\t"+currAlpha+"\t iterating");
                    }
                    else{	
                    	// eliminate worst case components, could be 1 or multiple components
                    	// redistribute responsibilities in next E step
                    	for (int j: min.cdr()){
                        	pi[j]=0;
	                		// clear responsibility
	                		for(int c=0; c<numConditions; c++){
	            				double[][] rc = r.get(c);
	            				int c_i_start[]=i_start.get(c);
	            				int c_i_end[] = i_end.get(c);
	            				for(int i=c_i_start[j];i<c_i_end[j];i++){
	                        		rc[i][j] = 0;
	                        	}
	            				c_i_start[j]=0;
	            				c_i_end[j]=0;
	                        }
                    	}
                    	// keep iterating on this Alpha value, until converge, then we raise it up to eliminate next one
                    	// give EM a time to stabilize before eliminating the next components
//                    	System.out.println(t+":\t"+currAlpha+"\t elimination");
                    	currAlpha = Math.max(min.car(), alpha/2);
                    }
            	}
            }

            // update component count, normalize pi
            double totalPi=0;
            double count = 0;
            for(int j=0;j<numComp;j++){
            	totalPi+=pi[j];
            	if (pi[j]!=0)
            		count++;
            }
            nonZeroComponentNum = count;
            if (totalPi!=0){
	            for(int j=0;j<numComp;j++){
	            	pi[j]=pi[j]/totalPi;
	            }
            }
			if(print_mixing_probabilities){
				StringBuilder sb = new StringBuilder();
				for (int i=0; i<pi.length; i++){
					sb.append(pi[i]+"\t");
				}
				pi_sb.append(sb.toString().trim()+"\n");
				appendEMMixProbFile(componentSpacing);
			}

			//Beta parameters
			if(numConditions>1){
				for(int j=0;j<numComp;j++){
	                if (pi[j]==0)
	                	for(int c=0; c<numConditions; c++){
	                		double[]bc = b.get(c);
	                        bc[j] = 1.0/numConditions;
	                	}
	                else{
	                    double b_sum=0;
		                for(int c=0; c<numConditions; c++){
	        				double[][] rc = r.get(c);
	        				double[]bc = b.get(c);
	        				double[]counts_c = counts.get(c);
	        				double sum_i=0;
	        				int c_i_start[]=i_start.get(c);
	        				int c_i_end[] = i_end.get(c);
	        				for(int i=c_i_start[j];i<c_i_end[j];i++)
	                    		sum_i += rc[i][j]*counts_c[i];	//TODO: verify
	        				
	                    	bc[j]=sum_i;
	                    	b_sum+=sum_i;
	                    }

	                    // normalize across conditions
		                if (b_sum!=0){
			                for(int c=0; c<numConditions; c++){
			                	double[]bc = b.get(c);
			                	bc[j] = bc[j]/b_sum;
			                }
		                }
	                } // else
				} // for
				// Beta clustering would go here.
			}

			//Semi-E-step:Calculate next un-normalized responsibilities
			for(int j=0;j<numComp;j++){
				if (pi[j]!=0){
					for(int c=0; c<numConditions; c++){
						double[][] rc = r.get(c);
						double[][] hc = h.get(c);
						double[] bc = b.get(c);
						int c_i_start[]=i_start.get(c);
        				int c_i_end[] = i_end.get(c);
        				for(int i=c_i_start[j];i<c_i_end[j];i++){
							rc[i][j] = pi[j]*bc[j]*hc[i][j];
						}
					}
				}
			}	

			//Log-likelihood calculation 
			double LL =0;
			for(int c=0; c<numConditions; c++){
				double[][] rc = r.get(c);
				double[]counts_c = counts.get(c);
				int numReads = rc.length;
				for(int i=0;i<numReads;i++){
					double sum_j=0;
					for(int j=0;j<numComp;j++)
						sum_j += rc[i][j]*counts_c[i];

					if (sum_j!=0)
						LL+=Math.log(sum_j);
				}
			}
			// log prior
			double LP=0;
			for(int j=0;j<numComp;j++)
				if(pi[j]>0)
					LP+=Math.log(pi[j]);

			LP= -currAlpha*LP;
			LAP = LL+LP;
			
			log(5, (int)nonZeroComponentNum+" ");
			
//			log(5, "\nEM: "+t+"\t"+currAlpha+"\t"+LAP+"\t"+lastLAP+"\t("+(int)nonZeroComponents+" non-zero components).");
			if (t<2 || Math.abs(LAP-lastLAP)>EM_CONVERGENCE){
				continue;
			}
			else{	
				// if converge on a smaller alpha value
				if (currAlpha<alpha){
					currAlpha = alpha;		// raise alpha, it may come down again after eliminating the next comp
//					System.out.println(t+"::\t"+currAlpha);
					continue;
				}
				// else: converge on full alpha value
				if (componentSpacing==1 && do_model_selection && (!BATCH_ELIMINATION)){
				// BIC model selection to decide which solution is better
				// 1. save the state
					EM_State state = new EM_State(this.numConditions);
					state.LL = LL;
					state.numComponent = nonZeroComponentNum;
					for (int c: r.keySet()){
						double[][]rc = r.get(c);
						double[][]copy = new double[rc.length][];
						for (int i=0;i<rc.length;i++){
							copy[i] = rc[i].clone();
						}
						state.resp[c]= copy;
					}
					for (int c: b.keySet()){
						state.beta[c]=b.get(c).clone();
					}
					state.pi = pi.clone();
					state.alpha = currAlpha;

					models.add(state);
				//2. more than one component left, raise alpha, continue EM
					if (nonZeroComponentNum>1){
						currAlpha *= 2;
						continue;
					}
					else	// single component left, done
						break;
				}
				
				else // if at 5bp resolution step, done with EM
					break;
			}
		} //LOOP: Run EM while not converged

		if (componentSpacing==1 && do_model_selection && models.size()>1 && (!BATCH_ELIMINATION)){
			// BIC model selection to decide which solution is better
			// 3. find best solution
			double totalCount = 0;
			for(int c=0; c<numConditions; c++){
				double[]counts_c = counts.get(c);
				for (double cc:counts_c){
					totalCount +=cc;
				}
			}
			int best = 0;
			double bestBIC = models.get(0).BIC(totalCount);
			for (int i=1;i<models.size();i++){
				double bic = models.get(i).BIC(totalCount);
//				System.out.println(String.format("%.3f\t%.3f\t", bestBIC, bic)+models.get(i).toString());
				if (bestBIC <= bic){
					best = i;
					bestBIC = bic;
				}
			}
			// copy the best model state back to memory
//			System.out.println("************BEST MODEL***************");
			EM_State bestModel = models.get(best);
//			for (int c:r.keySet()){
//				double[][] rc = r.get(c);
//				System.out.print(arrayToString(rc));
//			}
//			System.out.println("************copy-back best model***************");
			r.clear();
			for (int c=0;c<numConditions; c++){
				r.put(c, bestModel.resp[c]);
			}
//			for (int c:r.keySet()){
//				double[][] rc = r.get(c);
//				System.out.print(arrayToString(rc));
//			}
			b.clear();
			for (int c=0;c<numConditions; c++){
				b.put(c, bestModel.beta[c]);
			}			
			for (int i=0;i<pi.length;i++){
				pi[i] = bestModel.pi[i];
//				System.out.print(String.format("%.2f ", pi[i]));				
			}	
//			System.out.println();			
		}
		// normalize responsibilities here
		for(int c=0; c<numConditions; c++){
			double[][] rc = r.get(c);
			int numBases = rc.length;
			for(int i=0;i<numBases;i++){
				double totalResp = 0;
				for(int j=0;j<numComp;j++){
					totalResp += rc[i][j];
				}
				if (totalResp!=0){
					for(int j=0;j<numComp;j++){
						rc[i][j] = rc[i][j]/totalResp;
					}
				}
			}
		}

		// make hard assignment
		if (componentSpacing==1 && MAKE_HARD_ASSIGNMENT){
			for(int c=0; c<numConditions; c++){
				double[][] rc = r.get(c);
				int numReads = rc.length;
				for(int i=0;i<numReads;i++){
					Pair<Double,  TreeSet<Integer>> max = StatUtil.findMax(rc[i]);
					if (max.car()==0)	// max is 0, this read is not assigned to any event
						continue;
					int index = max.cdr().first();
					for(int j=0;j<numComp;j++){
						if (j==index)
							rc[i][j]=1;
						else
							rc[i][j]=0;
					}
				}
			}
		}
		
//		log(4, "EM_MAP(): "+timeElapsed(tic)+"\tt="+t+"\t"+
//				String.format("%.6f",LAP)+"\t("+(int)nonZeroComponents+" events)");

	}//end of EM_Steps method
	
	/**
	 * Loads a set of regions from a MACS peak file. <br>
	 * For the proper function of the method, regions contained in the MACS file should be sorted,
	 *  something done inside the method's body.
	 * @param peaks The <tt>List</tt> created by a MACS file
	 * @return
	 */
	public void setRegions(List<MACSPeakRegion> peaks){
		ArrayList<Region> regions= new ArrayList<Region>();
		for (MACSPeakRegion p: peaks){
			// filter towers
			if (p.getTags()>1000 && p.getFold_enrichment()<5){
				//TODO: should compare with WCE to decide it is really a tower
				continue;
			}
			regions.add(p);
		}
		restrictRegions = mergeRegions(regions, true);
	}
	
	
	public void setRegions(ArrayList<Region> regions){
		restrictRegions = regions;
	}
	

	/**
	 * Loads a set of regions (Format: Chrom:start-end). <br>
	 * For the proper function of the method, regions contained in the file should
	 * be sorted something done inside the method's body.
	 * @param fname the file containing the regions
	 * @return
	 */
	public void setRegions(String fname, boolean expandedRegion) {
		ArrayList<Region> rset = new ArrayList<Region>();
		try{
			File rFile = new File(fname);
			if(!rFile.isFile()){System.err.println("Invalid file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(rFile));
	        String line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
            	Region r = Region.fromString(gen, words[0]);
            	if (r!=null)
            		rset.add(r);
            }
	        reader.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

		if (expandedRegion){			// if regions are from previous round of analysis
			restrictRegions = mergeRegions(rset, false);
		}
		else{							// if regions are from other sources (i.e. StatPeakFinder)
			restrictRegions = mergeRegions(rset, true);	
		}
	}
	
	// merge the overlapped regions 
	// if "toExpandRegion"=true, expand each region on both side to leave enough space, 
	// to include every potential reads, then merge
	private ArrayList<Region> mergeRegions(ArrayList<Region> regions, boolean toExpandRegion){
		ArrayList<Region> mergedRegions = new ArrayList<Region>();
		if (regions.isEmpty())
			return mergedRegions;
		Collections.sort(regions);
		Region previous = regions.get(0);
		if (toExpandRegion)
			previous = previous.expand(modelRange, modelRange);
		mergedRegions.add(previous);
		
		for (Region region: regions){
			if (toExpandRegion)
				region=region.expand(modelRange, modelRange);
			// if overlaps with previous region, combine the regions
			if (previous.overlaps(region)){
				mergedRegions.remove(previous);
				previous = previous.combine(region);
			}
			else{
				previous = region;
			}
			mergedRegions.add(previous);
		}
		return mergedRegions;
	}//end of mergeRegions method

	// get the list of enriched regions directly from data
	private ArrayList<Region> selectEnrichedRegions(){
		long tic = System.currentTimeMillis();
		ArrayList<Region> regions = new ArrayList<Region>();
		for (String chrom: gen.getChromList()){
			int length = gen.getChromLength(chrom);
			Region wholeChrom = new Region(gen, chrom, 0, length-1);
			List<Region> rs = new ArrayList<Region>();
			List<StrandedBase> allBases = new ArrayList<StrandedBase>();
			for (Pair<ReadCache, ReadCache> pair: caches){
				ReadCache ip = pair.car();
				List<StrandedBase> bases = ip.getUnstrandedBases(wholeChrom); 
				if (bases==null || bases.size()==0){
					continue;
				}
				allBases.addAll(bases);	// pool all conditions
			}
			Collections.sort(allBases); 
			int start=0;
			for (int i=1;i<allBases.size();i++){
				int distance = allBases.get(i).getCoordinate()-allBases.get(i-1).getCoordinate();
				if (distance > modelWidth){ // a large enough gap to cut
					// only select region with read count larger than minimum count
					float count = 0;
					for(int m=start;m<=i-1;m++){
						count += allBases.get(m).getCount();
					}
					if (count >= sparseness){							
						Region r = new Region(gen, chrom, allBases.get(start).getCoordinate(), allBases.get(i-1).getCoordinate());
						rs.add(r);						
					}
					start = i;
				}
			}
			// check last continuous region
			float count = 0;
			for(int m=start;m<allBases.size();m++){
				count += allBases.get(m).getCount();
			}
			if (count>=sparseness){							
				Region r = new Region(gen, chrom, allBases.get(start).getCoordinate(), allBases.get(allBases.size()-1).getCoordinate());
				rs.add(r);						
			}
			
			// check regions, exclude un-enriched regions based on control counts
			ArrayList<Region> toRemove = new ArrayList<Region>();
			for (Region r: rs){
				if (r.getWidth()<=100){
					toRemove.add(r);
				}
				// for regions <= 500bp, most likely single event, can be compared to control
				if (r.getWidth()<=modelWidth){
					boolean enriched = false;
					for (int c=0;c<numConditions;c++){
						if (countIpReads(r,c)/countCtrlReads(r,c)/this.ratio_total[c]>IP_CTRL_FOLD){
							enriched = true;
							break;
						}
					}
					if (!enriched)	// remove this region if it is not enriched in any condition
						toRemove.add(r);
				}
			}
			rs.removeAll(toRemove);
			
			if (!rs.isEmpty()){
					regions.addAll(rs);
			}
		} // each chrom
		
		log(3, "selectEnrichedRegions(): "+timeElapsed(tic));
		return regions;
	}
	
	/**
	 * Initializes the components. Spaces them evenly along the region for now.
	 * 
	 * @param currReg
	 * @param signals
	 * @param weighted (uniform or weighted initialization?)
	 */
	protected void initializeComponents(Region currReg,int numCond){
		components = new ArrayList<BindingComponent>();

		//decide how many components
		if (currReg.getWidth()/50 > MAX_NUM_COMPONENTS)
			System.err.println("Very large region, "+currReg.toString());
		componentMax = Math.max(OPTIMAL_NUM_COMPONENTS, currReg.getWidth()/50);
		if(componentMax>0){
			if (SMART_SPACING){				
				int spacing = componentMax>=currReg.getWidth() ? 1 : Math.max(2, (currReg.getWidth()/componentMax)+1);
				componentSpacing = Math.max(currReg.getWidth()/MAX_NUM_COMPONENTS, Math.min(spacing, INIT_RESOLUTION));
			}
			else
				componentSpacing = Math.max(currReg.getWidth()/MAX_NUM_COMPONENTS, INIT_RESOLUTION);
			
			int numComponents = currReg.getWidth()/componentSpacing;

			//Set up the components
			double totalMP =0;
			for(int i=0; i<numComponents; i++){
				Point pos = new Point(gen, currReg.getChrom(), currReg.getStart()+(i*componentSpacing));
				BindingComponent currComp = new BindingComponent(model, pos, numCond);
				currComp.setMixProb(1);
				totalMP+=currComp.getMixProb();
				components.add(currComp);						
			}
			//Normalize mixProbs
			for(BindingComponent b : components)
				b.setMixProb(b.getMixProb()/totalMP);
		}
		nonZeroComponentNum=components.size();		
	}

	//Update the resolution of the components 
	//Add new components in around non-zero components	
	protected void updateComponentResolution(Region currReg, int numCond, int lastResolution){
		double extend=2;
		ArrayList<BindingComponent> newComponents = new ArrayList<BindingComponent>();

		//First map the valid bases around non-zero components
		int numValidBases=0;
		int[] valid = new int[currReg.getWidth()];
//		for(int v=0; v<currReg.getWidth(); v++){valid[v]=0;}
		for(BindingComponent b : components){
			for(int v=(int)Math.max(0, (b.getLocation().getLocation()-lastResolution*extend)-currReg.getStart()); 
			    v<=Math.min(currReg.getWidth()-1, (b.getLocation().getLocation()+lastResolution*extend)-currReg.getStart()); 
			    v++){
				if(valid[v]==0){
					valid[v]=1;
					numValidBases++;
				}	
			}	
		}	

		//Reset the componentMax to reflect the number of reads in the valid areas
		//Calc new resolution
		componentSpacing = componentMax>=numValidBases ? 1 : Math.max(2, (numValidBases/componentMax)+1);
		if(componentSpacing<lastResolution){
			//Set up the components
			double totalMP =0;
			for(BindingComponent b: components){
				for(int v=(int)Math.max(0, (b.getLocation().getLocation()-lastResolution*extend+1)-currReg.getStart()); 
				    v<=Math.min(currReg.getWidth()-1, (b.getLocation().getLocation()+lastResolution*extend-1)-currReg.getStart()); 
				    v+=componentSpacing){
					Point pos = new Point(gen, currReg.getChrom(), currReg.getStart()+v);
					//Reusing the valid array for the sake of it
					if(valid[v]!=-1){
						BindingComponent currComp = new BindingComponent(model, pos, numCond);
						currComp.setMixProb(1.0);
						totalMP+=currComp.getMixProb();
						newComponents.add(currComp);
						valid[v]=-1;
			}	}	}	

			//Normalize mixProbs
			for(BindingComponent b : newComponents){
				b.setMixProb(b.getMixProb()/totalMP);				
			}
			components=newComponents;
			nonZeroComponentNum = components.size();
		}
	}

	/* Given the bases data, and the position of events
	 * run standard EM (no component elimination) to assign bases to events
	 */
	private double[][] assignBasesToEvents(List<StrandedBase> bases, ArrayList<BindingComponent> comps){
		// assuming components only contains non-zero events
		int numComp = comps.size();
		int numBases = bases.size();

		double[][] r;		//[base][components]
		if (numComp==1){
			r= new double[numBases][1];
			BindingComponent comp = comps.get(0);			
			for(int i=0;i<numBases;i++){
				StrandedBase base = bases.get(i);
				if (comp.scoreBase(base)>0)
					r[i][0]=1;
			}	
		}
		else{
			// pi initialize with IP event strength
			double[] pi = new double[numComp];				// Pi
			for(int j=0;j<numComp;j++){
				pi[j]= comps.get(j).getTotalSumResponsibility();
			}
			StatUtil.mutate_normalize(pi);
			
			double[][] h= new double[numBases][numComp]; 	// H function
			for(int i=0;i<numBases;i++){
				StrandedBase base = bases.get(i);
				for(int j=0;j<numComp;j++){
					BindingComponent comp = comps.get(j);
					h[i][j]=comp.scoreBase(base);
				}
			}
			
			r= new double[numBases][numComp];	// Responsibility
			for(int j=0;j<numComp;j++){
				for(int i=0;i<numBases;i++){
					r[i][j] = h[i][j]*pi[j];
				}	
			}
			
			double[] counts= new double[numBases];
			for(int i=0;i<numBases;i++)
				counts[i]=bases.get(i).getCount();
			
			// run EM
			EM_ML( counts, h, r, pi);
		}
		
		return r;
	}
	/** 
	 * core EM steps, standard ML, single condition 
	 * this is used to assign control reads to corresponding events 
	 */
	private void EM_ML(  double[] counts, double[][] h, double[][] r, double[] pi) {
//		long tic = System.currentTimeMillis();
		int numComp = pi.length;		
		int numBases = r.length;
		double totalResp[] = new double[numBases];
		
		//Run EM while not converged
		double lastLL=0, LL=0;
		for(int t=0; t<MAX_EM_ML_ITER ; t++){
			lastLL=LL;

			//////////
			//Semi-E-step (presumes unnormalized responsibilities have been calculated)
			//Simply normalize responsibilities here
			//////////	
			for(int i=0;i<numBases;i++){
				totalResp[i] = 0;
				for(int j=0;j<numComp;j++){
					totalResp[i] += r[i][j];
				}
				for(int j=0;j<numComp;j++){
					if (totalResp[i]>0)
						r[i][j] = r[i][j]/totalResp[i];
				}
			}

			//////////
			//M-step : Pi
			//////////
			for(int j=0;j<numComp;j++){
                double r_sum=0;
                for(int i=0;i<numBases;i++)
                	r_sum += r[i][j]*counts[i];
                pi[j]=r_sum;    // standard EM ML
            }
            StatUtil.mutate_normalize(pi);

			//Semi-E-step:Calculate next un-normalized responsibilities
    		for(int j=0;j<numComp;j++){
    			for(int i=0;i<numBases;i++){
    				r[i][j] = h[i][j]*pi[j];
    			}	
    		}
			//Log-likelihood calculation 
			LL =0;
			for(int i=0;i<numBases;i++){
				double sum_j=0;
				for(int j=0;j<numComp;j++){
					sum_j += r[i][j]*counts[i];
				}

				if (sum_j!=0){
					LL+=Math.log(sum_j);
				}
			}

//			log(4, "\nEM_ML: "+t+"\t"+LL+"\t"+lastLL+".");
			if (Math.abs(LL-lastLL)>EM_ML_CONVERGENCE ){
				continue;
			}
			else{
				break;
			}

		} //LOOP: Run EM while not converged

		// normalize responsibilities here
		for(int i=0;i<numBases;i++){
			totalResp[i] = 0;
			for(int j=0;j<numComp;j++){
				totalResp[i] += r[i][j];
			}
			for(int j=0;j<numComp;j++){
				if (totalResp[i]>0)
					r[i][j] = r[i][j]/totalResp[i];
			}
		}
		
		// make hard assignment
		if (MAKE_HARD_ASSIGNMENT){
			for(int i=0;i<numBases;i++){
				Pair<Double,  TreeSet<Integer>> max = StatUtil.findMax(r[i]);
				if (max.car()==0)	// max is 0, this read is not assigned to any event
					continue;
				int index = max.cdr().first();
				for(int j=0;j<numComp;j++){
					if (j==index)
						r[i][j]=1;
					else
						r[i][j]=0;
				}
			}
		}
//		log(4, "EM time, EM_ML(): "+timeElapsed(tic));

	}//end of EM_ML method
	
	/**
	 * Returns the component positions based on the data (unique) positions, the window chosen,
	 * and the spacing between the components.
	 * @param width width of the region whose components will be determined
	 * @param pos unique positions of the reads
	 * @param window expansion window
	 * @param spacing spacing between the component
	 * @return
	 */
	private int[] setComponentPositions(int width, int[][] pos, int window, int spacing) {
		int[] compPos;
		Set<Integer> initialSet = new TreeSet<Integer>();
		for(int t = 0; t < pos.length; t++) {
			for(int k = 0; k < pos[t].length; k++) {
				int start = Math.max(      0, pos[t][k] - window);
				int   end = Math.min(width-1, pos[t][k] + window);
				for(int i = start; i <= end; i+=spacing) { initialSet.add(i); }
			}
		}
		
		int prev = -1;
		Iterator<Integer> itr = initialSet.iterator();
		Set<Integer> finalSet   = new TreeSet<Integer>();
		while(itr.hasNext()) {
			int curr = itr.next();			 
			if(prev != -1 && 0 < curr - prev && curr - prev < spacing) {
				finalSet.add(Math.min(width-1, prev+spacing));
				prev += spacing;
			}
			else if(prev == -1 || curr - prev >= spacing) {
				finalSet.add(curr);
				prev = curr;
			}
		}
				
		compPos = Utils.ref2prim(finalSet.toArray(new Integer[0]));
		return compPos;
	}//end of setComponentPositions method

	protected ArrayList<ComponentFeature> callFeatures(ArrayList<BindingComponent> comps){
		if (post_artifact_filter){
			comps = postArtifactFilter(comps);
		}		
		
		ArrayList<ComponentFeature> features = new ArrayList<ComponentFeature>();
		
		if (comps.size()==0)
			return features;
		
		//Make feature calls (all non-zero components)
		for(BindingComponent b : comps){
			if(b.getMixProb()>0){
				features.add(callFeature(b, use_internal_em_train));
			}
		}

		// assign control reads to IP events for each pair of expt/ctrl (condition)
		// if no control data, skip
		if (controlDataExist){		
			// expand to the region covering all events
			Collections.sort(comps);
			Region r = comps.get(0).getLocation().expand(modelRange);
			if (comps.size()>1){
				r=new Region(r.getGenome(), r.getChrom(), r.getStart(),comps.get(comps.size()-1).getLocation().getLocation()+modelRange);
			}
				
			for(int c=0; c<caches.size(); c++){	// each condition
				Pair<ReadCache,ReadCache> exptPair = caches.get(c);
				ReadCache control = exptPair.cdr();
				// read count from Control data in specified region 
				List<StrandedBase> bases= control.getStrandedBases(r, '+');  // plus reads of the current region - CTRL channel
				List<StrandedBase> bases_m= control.getStrandedBases(r, '-');  // minus reads of the current region - CTRL channel
				
				if(needToCleanBases) {
					cleanBases(bases);
					cleanBases(bases_m);
				}
				
				if (pre_artifact_filter){
					filterBases(bases, control, r, '+');
					filterBases(bases_m, control, r, '-');
				}
				
				if (post_artifact_filter){
					// filter bases using probability landscape (knowledge of predicted events and read distribution)
					double landscape_p[] = new double[r.getWidth()];
					double landscape_m[] = new double[r.getWidth()];
					double readDensity[] = model.getProbabilities();
					boolean skip = false;
					for (BindingComponent b: comps){
						for (int i=0; i<readDensity.length;i++){
							int index = b.getLocation().getLocation()-r.getStart()+model.getMin()+i;
							// if the position is on the edge of chrom, the index will be out of bound, skip filtering
							if (index<0 || index>=landscape_p.length){
								//System.out.println("postArtifactFilter\t"+b.getLocation().getLocationString()+"\tindex="+index);
								skip=true;
							}
							else
								landscape_p[index] += b.getTotalSumResponsibility()*readDensity[i];
							int index2 = b.getLocation().getLocation()-r.getStart()-model.getMin()-i;
							if (index2<0 || index2>=landscape_p.length){
								//System.out.println("postArtifactFilter\t"+b.getLocation().getLocationString()+"\tindex="+index2);
								skip=true;
							}
							else
								landscape_m[index2] += b.getTotalSumResponsibility()*readDensity[i];
						}
					}
					if (!skip){
						StatUtil.normalize(landscape_p);
						StatUtil.normalize(landscape_m);
						
						filterBases( bases, r.getStart(), landscape_p);
						filterBases( bases_m, r.getStart(), landscape_m);
					}
				}
				bases.addAll(bases_m);

				double[][] assignment = assignBasesToEvents(bases, comps);
				
				for(int j=0;j<features.size();j++){
					ComponentFeature comp = features.get(j);
					int pos = comp.getPeak().getLocation();
					// set control reads count and logKL
					double total=0;
					double[] profile_plus = new double[modelWidth];
					double[] profile_minus = new double[modelWidth];
					for(int i=0;i<bases.size();i++){
						StrandedBase base = bases.get(i);
						if (assignment[i][j]>0){
							if (base.getStrand()=='+'){
								if (base.getCoordinate()-pos-model.getMin()>modelWidth)
									System.err.println(pos+"\t+\t"+base.getCoordinate());
								profile_plus[base.getCoordinate()-pos-model.getMin()]=assignment[i][j]*base.getCount();
							}
							else{
								if (pos-base.getCoordinate()-model.getMin()>modelWidth)
									System.err.println(pos+"\t-\t"+base.getCoordinate());
								profile_minus[pos-base.getCoordinate()-model.getMin()]=assignment[i][j]*base.getCount();
							}
							total += assignment[i][j]*base.getCount();
						}
					}
					comp.setControlReadCounts(total, c);
					double logKL_plus  = StatUtil.log_KL_Divergence(model.getProbabilities(), StatUtil.symmetricKernelSmoother(profile_plus, gaussian));
					double logKL_minus = StatUtil.log_KL_Divergence(model.getProbabilities(), StatUtil.symmetricKernelSmoother(profile_minus, gaussian));
					comp.setCtrlProfileLogKL(c, logKL_plus, logKL_minus);
				}
			}
		}
		
		return(features);
	}

	// make componentFeature based on the bindingComponent
	private ComponentFeature callFeature(BindingComponent b, boolean use_internal_em_train){
		ComponentFeature cf = new ComponentFeature(b, use_internal_em_train);
//		double mfold = eval_mfold(b);
//		cf.set_mfold(mfold);
//		
//		
//		mm
		
		// KL
		double[] logKL_plus = new double[numConditions];
		double[] logKL_minus = new double[numConditions];
		for (int c=0;c<numConditions;c++){
			logKL_plus[c]  = StatUtil.log_KL_Divergence(model.getProbabilities(), StatUtil.symmetricKernelSmoother(b.getReadProfile_plus(c), gaussian));
			logKL_minus[c] = StatUtil.log_KL_Divergence(model.getProbabilities(), StatUtil.symmetricKernelSmoother(b.getReadProfile_minus(c), gaussian));
		}
		cf.setProfileLogKL(logKL_plus, logKL_minus);

		return cf;
	}

	private ArrayList<BindingComponent> postArtifactFilter(ArrayList<BindingComponent> comps){
		ArrayList<BindingComponent> filteredComponents = new ArrayList<BindingComponent>();
		
		// expand to the region covering all events
		Collections.sort(comps);
		Region r = comps.get(0).getLocation().expand(modelRange);
		if (comps.size()>1){
			r=new Region(r.getGenome(), r.getChrom(), r.getStart(),comps.get(comps.size()-1).getLocation().getLocation()+modelRange);
		}
		// filter bases using probability landscape (knowledge of predicted events and read distribution)
		double landscape_p[] = new double[r.getWidth()];
		double landscape_m[] = new double[r.getWidth()];
		double readDensity[] = model.getProbabilities();
		for (BindingComponent b: comps){
			for (int i=0; i<readDensity.length;i++){
				int index = b.getLocation().getLocation()-r.getStart()+model.getMin()+i;
				// if the position is on the edge of chrom, the index will be out of bound, skip it
				if (index<0 || index>=landscape_p.length){
					//System.out.println("postArtifactFilter\t"+b.getLocation().getLocationString()+"\tindex="+index);
					return comps;
				}
				landscape_p[index] += b.getTotalSumResponsibility()*readDensity[i];
				int index2 = b.getLocation().getLocation()-r.getStart()-model.getMin()-i;
				if (index2<0 || index2>=landscape_p.length){
					//System.out.println("postArtifactFilter\t"+b.getLocation().getLocationString()+"\tindex="+index2);
					return comps;
				}
				landscape_m[index2] += b.getTotalSumResponsibility()*readDensity[i];
			}
		}
		StatUtil.normalize(landscape_p);
		StatUtil.normalize(landscape_m);
		
		ArrayList<List<StrandedBase>> signals = new ArrayList<List<StrandedBase>>();
		//Load and filter each condition's read hits
		for(Pair<ReadCache,ReadCache> e : caches){
			List<StrandedBase> bases_p= e.car().getStrandedBases(r, '+');  // reads of the current region - IP channel
			List<StrandedBase> bases_m= e.car().getStrandedBases(r, '-');  // reads of the current region - IP channel
			// filter
			filterBases( bases_p, r.getStart(), landscape_p);
			filterBases( bases_m, r.getStart(), landscape_m);
			
			bases_p.addAll(bases_m);
			signals.add(bases_p);
		}
		
		// We want to run EM only for joint events
		if (comps.size()>1) {
			// dynamically determine an alpha value for this sliding window
			double alpha = sparseness;
			if (use_dynamic_sparseness)
				alpha = Math.max(estimateAlpha2(r, signals), sparseness);
			
			HashMap<Integer, double[][]> responsibilities = null;
			
			// Multi-condition??
			if(use_internal_em_train) {
				//Run EM and increase resolution
				initializeComponents(r, numConditions);
				int lastResolution;
		
				while(nonZeroComponentNum>0){
					lastResolution = componentSpacing;
					// EM learning, components list will only contains non-zero components
					responsibilities = EMTrain(signals, alpha);
					// log(4, componentSpacing+" bp\t"+(int)nonZeroComponents+" components.");

					// increase resolution
					updateComponentResolution(r, numConditions, lastResolution);
					if(componentSpacing==lastResolution)
						break;
				} 	// end of while (resolution)

			}
			else {
				//TODO: Multi-condition??
			}
			
			if (nonZeroComponentNum==0)	return filteredComponents;
			setComponentResponsibilities(signals, responsibilities);
			filteredComponents = this.components;
		} 	// if not single event region
		else{// if single event region, just scan it
			BindingComponent peak = scanPeak(signals, r);
			filteredComponents = new ArrayList<BindingComponent>();
			filteredComponents.add(peak);
		}
		
		return filteredComponents;
	}
			
	void filterBases(List<StrandedBase> bases, int startCoor, double[] landscape){
		float totalCount = StrandedBase.countBaseHits(bases);
		for (int i=0;i<bases.size();i++){
			StrandedBase base = bases.get(i);
			float count = base.getCount();
			// only check bases that pass the threshold
			if (count > max_HitCount_per_base){
				// assume this position corresponds to summit of empirical distribution.
				double prob = landscape[base.getCoordinate()-startCoor];
				// excluding this base for calculation
				double expectedCount = (totalCount-count) / (1-prob) * prob* needle_height_factor;
				// if count of this base is higher than needle_height_factor (2) times expected count from emp.dist.
				if (count> expectedCount ){
					base.setCount((float)Math.max(1, expectedCount));
				}
			}
			//System.out.println(base.getCoordinate()+"\t"+count+"\t"+base.getCount());
		}
	}
	
	
//	private double eval_mfold(BindingComponent b) {
//		Region r = b.getRegion();
//		ArrayList<List<StrandedBase>> signals = loadBasesInWindow(r);
//		float[] sigHitCounts = new float[signals.size()];
//		int totalSigCount=0;
//		for(int c=0; c<signals.size(); c++){
//			sigHitCounts[c]=StrandedBase.countBaseHits(signals.get(c));
//			totalSigCount += sigHitCounts[c];
//		}		
//		
//	}
	
	/**
	 * Finds the exact peak in the region by scanning for maximum-likelihood using a <tt>BindingModel</tt>
	 * assuming that the region contains only a single event
	 * @param signals <tt>ReadHits</tt> of the IP channel corresponding to region <tt>region</tt>
	 * @param region Regions to be scanned
	 * @return the peak as a <tt>BindingComponent</tt>
	 */
	private BindingComponent scanPeak(Point p){
		Region scanRegion = p.expand(SCAN_RANGE);
		Region plusRegion = scanRegion.expand(-model.getMin(), model.getMax());
		Region minusRegion = scanRegion.expand(model.getMax(), -model.getMin());
		ArrayList<List<StrandedBase>> signals = new ArrayList<List<StrandedBase>>();
		//Load each condition's read hits
		int totalHitCounts = 0;
		for(int c=0;c<numConditions;c++){
			Pair<ReadCache,ReadCache> e = caches.get(c);
			List<StrandedBase> bases_p= e.car().getStrandedBases(plusRegion, '+');  // reads of the current region - IP channel
			List<StrandedBase> bases_m= e.car().getStrandedBases(minusRegion, '-');  // reads of the current region - IP channel
			
			if(needToCleanBases) {
				cleanBases(bases_p);
				cleanBases(bases_m);
			}
			
			if (pre_artifact_filter){
				filterBases(bases_p, e.car(), plusRegion, '+');
				filterBases(bases_m, e.car(), minusRegion, '-');
			}
			
			bases_p.addAll(bases_m);
			totalHitCounts += StrandedBase.countBaseHits(bases_p);
			signals.add(bases_p);
		} // for loop
		
		if (signals==null||totalHitCounts==0) // check for empty read region
			return null;

		return scanPeak( signals, scanRegion);
	}
	
	private BindingComponent scanPeak(ArrayList<List<StrandedBase>> signals, Region scanRegion){
		ArrayList<StrandedBase> allBases = new ArrayList<StrandedBase>();
		// for multi-condition data, pool all reads
		for (int c=0;c<numConditions;c++){
			allBases.addAll(signals.get(c));
		}

		int newModelRange = scanRegion.getWidth()-1+modelRange;
		BindingModel newModel = model.getExpandedModel(newModelRange-Math.abs(model.getMin()), 
				newModelRange-Math.abs(model.getMax()));
		
		BindingComponent maxComp=null; 
		double maxScore=Double.NEGATIVE_INFINITY;
		// scan in the given region the find the component that 
		// gives the max likelihood from all reads
		for(int i=scanRegion.getStart(); i<=scanRegion.getEnd(); i++){
			BindingComponent component = new BindingComponent(newModel, new Point(gen, scanRegion.getChrom(), i), numConditions);
			double logLikelihood=0;
			for(StrandedBase base : allBases){
				double prob = component.scoreBase(base);
				// single peak assumption, pi=1, prob*pi = prob
				logLikelihood += Math.log(prob*base.getCount());
			}
			if(logLikelihood>maxScore){
				maxComp=component;
				maxScore = logLikelihood;
			}
//			System.out.println(i-region.getStart()+"\t"
//					+String.format("%.2f", logLikelihood)+"\t"
//					+String.format("%.2f", maxScore));
		}
		
		// make a new binding component with original binding model
		if (maxComp==null)
			maxComp = new BindingComponent(model, scanRegion.getMidpoint(), numConditions);
		else
			maxComp = new BindingComponent(model, maxComp.getLocation(), numConditions);
		// set pi, beta, responsibility, note it is a single event
		maxComp.setMixProb(1);
		
		float[] hitCounts = new float[numConditions];
		int totalHitCounts = 0;
		for(int c=0;c<numConditions;c++){
			List<StrandedBase> bases= signals.get(c);
			hitCounts[c] = StrandedBase.countBaseHits(bases);
			totalHitCounts += hitCounts[c];
		} 
		
		if(use_internal_em_train) {
			for (int c=0;c<signals.size();c++)
				maxComp.setConditionBeta(c, hitCounts[c]/totalHitCounts);  //beta is being evaluated for one point across all conditions
		}
		else {
			for (int c=0;c<signals.size();c++)
				maxComp.setConditionBeta(c, 1.0);   //beta is being evaluated for one condition separately
		}
		
		for(StrandedBase base : allBases){
			for(int c = 0; c < numConditions; c++) {
				if (maxComp.scoreBase(base)>0 && maxComp.getConditionBeta(c) > 0.0)
					maxComp.setResponsibility(c, base, 1.0);				
			}
		}

		int coord = maxComp.getLocation().getLocation();
		for(int c=0; c<numConditions; c++){
			List<StrandedBase> bases = signals.get(c);

			// store binding profile (read responsibilities in c condition) of this component
			double[] profile_plus = new double[modelWidth];
			double[] profile_minus = new double[modelWidth];
			for(StrandedBase base : bases){
				int fivePrime = base.getCoordinate();
				if (base.getStrand()=='+'){
					int index = fivePrime-coord-model.getMin();
					// if the read is in the influence range of binding event
					if (index>=0 && index<modelWidth)
						profile_plus[index] += base.getCount();
				}
				else{
					int index = coord-fivePrime-model.getMin();
					if (index>=0 && index<modelWidth)
						profile_minus[index] += base.getCount();
				}
			}
			maxComp.setReadProfile(c, profile_plus, '+');
			maxComp.setReadProfile(c, profile_minus, '-');
		}		
		return maxComp;		
	}//end of scanPeak method

	
	public double updateBindingModel(int left, int right){
		if (signalFeatures.size()==0)
			return -100;
		double[] strengths = new double[signalFeatures.size()];
		double[] shapes = new double[signalFeatures.size()];
		for (int i=0; i<signalFeatures.size();i++){
			ComponentFeature cf = (ComponentFeature)signalFeatures.get(i);	
			strengths[i] = cf.getTotalSumResponsibility();
			shapes[i] = cf.getAverageLogKL();
		}
		Arrays.sort(strengths);
		Arrays.sort(shapes);
		double strengthThreshold = strengths[strengths.length*(100-top_event_percentile)/100];
		double shapeThreshold = shapes[shapes.length*top_event_percentile/100];
		
		int eventCountForModelUpdating=0;
        // data for all conditions
        double[][] newModel_plus=new double[numConditions][modelWidth];
        double[][] newModel_minus=new double[numConditions][modelWidth];    	
		// sum the read profiles from all qualified binding events for updating model			
		for (int c=0;c<numConditions;c++){
			ReadCache ip = caches.get(c).car();
			for (Feature f: signalFeatures){
				ComponentFeature cf = (ComponentFeature)f;	
				// the events that are used to refine read distribution should be
				// having strength and shape at the upper half ranking
				if ((use_multi_event || cf.getMixProb()==1)
						&& cf.getAverageLogKL()<=shapeThreshold 
						&& cf.getTotalSumResponsibility()>=strengthThreshold){
					if (cf.getQValueLog10(c)>q_value_threshold){
						Region region = cf.getPeak().expand(0).expand(left, right);
						List<StrandedBase> bases = ip.getStrandedBases(region, '+');
						for (StrandedBase base:bases){
							newModel_plus[c][base.getCoordinate()-region.getStart()]+=base.getCount();
						}
						region = cf.getPeak().expand(0).expand(right, left);
						bases = ip.getStrandedBases(region, '-');
						for (StrandedBase base:bases){
							newModel_minus[c][region.getEnd()-base.getCoordinate()]+=base.getCount();
						}
						eventCountForModelUpdating++;
					}
				}
			}
		}

		// we have collected data for both strands, all conditions 
		// but here we only use single model for all conditions
		double[] model_plus = new double[modelWidth];
		double[] model_minus = new double[modelWidth];
		for (int i=0;i<modelWidth;i++){
			for (int c=0; c<numConditions; c++){
				model_plus[i]+= newModel_plus[c][i];
				model_minus[i]+= newModel_minus[c][i];
			}
		}

		// smooth the model profile
		if (SPLINE_SMOOTH){
			model_plus=StatUtil.cubicSpline(model_plus, BindingModel.SMOOTHING_STEPSIZE, 3);
			model_minus=StatUtil.cubicSpline(model_minus, BindingModel.SMOOTHING_STEPSIZE, 3);
		}
		StatUtil.mutate_normalize(model_plus);
		StatUtil.mutate_normalize(model_minus);

		//compare models from 2 strands, shift binding position if needed
		int shift = BindingModel.minKL_Shift(model_plus, model_minus);
		
		//use single model for now
		List<Pair<Integer, Double>> dist = new ArrayList<Pair<Integer, Double>>();
		for (int i=0;i<modelWidth;i++){
			double sum = model_plus[i]+model_minus[i];
			sum = sum>=0?sum:2.0E-300;	
			Pair<Integer, Double> p = new Pair<Integer, Double>(i-left, sum);
			dist.add(p);
		}

		double[] oldModel = model.getProbabilities();
		String oldName = model.getFileName();
		model = new BindingModel(dist);
		model.setFileName(oldName);
		model.smooth(BindingModel.SMOOTHING_STEPSIZE);
		if (development_mode)
			model.printToFile(outName+"_-"+left+"_"+right+"_Read_distribution.txt");
		else
			model.printToFile(outName+"_Read_distribution.txt");
		modelRange = model.getRange();
		modelWidth = model.getWidth();

		double logKL = StatUtil.log_KL_Divergence(oldModel, model.getProbabilities());
		String details = "(Strength>"+String.format("%.1f", strengthThreshold) +
			", Shape<"+String.format("%.2f", shapeThreshold) +
			"). \nlogKL=" + String.format("%.2f",logKL) + 
			", +/- strand shift "+shift+" bp.\n";
		log(1, "Refine read distribution from " + eventCountForModelUpdating +" binding events. "+
				(development_mode?details:"\n"));
		
		return logKL;
	}
	
	public BindingModel getModel(){
		return model;
	}
	
	// determine the max hit count by threshold of Poisson p-value
	// set average read as lambda parameter for Poisson distribution
	protected int calcHitCount_per_BP(double totalReads, double threshold){
		int countThres=0;
		DRand re = new DRand();
		Poisson P = new Poisson(0, re);
		double lambda = totalReads/mappable_genome_length; 
		P.setMean(lambda);
		double l=1;
		for(int b=1; l>threshold; b++){
			l=1-P.cdf(b);	//p-value as the tail of Poisson
			countThres=b;
		}
		return Math.max(1,countThres);
	}
	
	// evaluate confidence of each called events
	// calculate p-value from binomial distribution, and peak shape parameter
	private void evaluateConfidence(ArrayList<ComponentFeature> compFeatures){	
		if(controlDataExist) {
			for (ComponentFeature cf: compFeatures){
				for(int cond=0; cond<caches.size(); cond++){
					// scale control read count by non-specific read count ratio
					double controlCount = cf.getScaledControlCounts(cond);
					double ipCount = cf.getEventReadCounts(cond);
					double pValue=1;
					if ((controlCount+ipCount)>=0){
						try{
							pValue = StatUtil.binomialPValue(controlCount, 
									controlCount+ipCount);
						}
						catch(Exception err){
							err.printStackTrace();
							System.err.println(cf.toString());
						}
					}
					else{
						System.err.println(cf.toString());
					}
					cf.setPValue(pValue, cond);

					//				// peak shape parameter
					//				double logKL = cf.getLogKL(cond).car();
					//				int strength = (int)Math.round(cf.getEventReadCounts(cond));
					//				if (strength>99) strength = 99;
					//				if (strength<5) strength = 5;
					//				// get mean and std of the normal distribution (of logKL values)
					//				Pair<Double, Double> params = shapeParameters.get(strength);
					//				if (params==null){
					//					System.err.println(cf.getPeak().toString()+" "+strength);
					//					continue;
					//				}
					//				double z_score = (logKL-params.car())/params.cdr();
					//				cf.setShapeZScore(z_score, cond);
					//				double shapeParam = z_score*Math.log(strength);
					//				cf.setShapeParameter(shapeParam, cond);
				}
			}
		}
		else {
			shift_size = eval_avg_shift_size(compFeatures, num_top_mfold_feats);
			for (ComponentFeature cf: compFeatures){
				for(int cond=0; cond<caches.size(); cond++){
					double pValue_wo_ctrl = evalFeatureSignificance(cf, cond);
					cf.setPValue_wo_ctrl(pValue_wo_ctrl, cond);
				}
			}			
		}
		
		// calculate q-values, correction for multiple testing
		benjaminiHochbergCorrection(compFeatures);
	}

	//Multiple hypothesis testing correction 
	// -- assumes peaks ordered according to p-value
	// this method will set log10(q-value) of component features
	private void benjaminiHochbergCorrection(ArrayList<ComponentFeature> compFeatures){
		if(controlDataExist) {
			for (int cond=0; cond<conditionNames.size(); cond++){
				ComponentFeature.setSortingCondition(cond);
				Collections.sort(compFeatures, new Comparator<ComponentFeature>(){
					public int compare(ComponentFeature o1, ComponentFeature o2) {
						return o1.compareByPValue(o2);
					}
				});

				double total = compFeatures.size();
				int rank =1;
				for(ComponentFeature cf : compFeatures){
					cf.setQValueLog10( -Math.log10(cf.getPValue(cond)*(total/(double)rank)), cond);
					cf.addRank_sum(rank);
					rank++;
				}
			}
		}
		else {
			for (int cond=0; cond<conditionNames.size(); cond++){
				ComponentFeature.setSortingCondition(cond);
				Collections.sort(compFeatures, new Comparator<ComponentFeature>(){
					public int compare(ComponentFeature o1, ComponentFeature o2) {
						return o1.compareByPValue_wo_ctrl(o2);
					}
				});

				double total = compFeatures.size();
				int rank = 1;
				for(ComponentFeature cf : compFeatures){
					cf.setQValueLog10( -Math.log10(cf.getPValue_wo_ctrl(cond)*(total/(double)rank)), cond);
					cf.addRank_sum(rank);
					rank++;
				}
			}			
		}
	}//end of benjaminiHochbergCorrection method
	
	private double eval_avg_shift_size(ArrayList<ComponentFeature> compFeatures, int num_top_mfold_feats) {
		double avg_shift_size;
		ComponentFeature[] top_mfold_Feats; 
		
		/******************************************
		 *****  Create chromosome statistics  *****
		 *****************************************/
		totalIPCounts = new HashMap<String, Integer>();
		chromFivePrimes = new HashMap<String, int[][][][]>();
		condHitCounts = new HashMap<String, List<List<Integer>>>();
		// In this loop we will create chromosome statistics.
		// That is, we will store the total IP counts for each chromosome (for all conditions): totalIPCounts
		// as well as the counts for each channel (IP, CTRL) and for each condition: condHitCounts
		for(String chrom:gen.getChromList()){
			int chromLen       = gen.getChromLength(chrom);
			Region chromRegion = new Region(gen, chrom, 0, chromLen-1);
			List<ArrayList<List<StrandedBase>>> curr_chrom_signals = new ArrayList<ArrayList<List<StrandedBase>>>();
			
			ArrayList<List<StrandedBase>> ip_chrom_signals = loadBasesInWindow(chromRegion, "IP", false);
			ArrayList<List<StrandedBase>> ctrl_chrom_signals = new ArrayList<List<StrandedBase>>();
			if(controlDataExist)
				ctrl_chrom_signals = loadBasesInWindow(chromRegion, "CTRL", false);
			else
				for(int t = 0; t < numConditions; t++)
					ctrl_chrom_signals.add(new ArrayList<StrandedBase>());
				
			curr_chrom_signals.add(ip_chrom_signals);
			curr_chrom_signals.add(ctrl_chrom_signals);
			chrom_signals.put(chrom, curr_chrom_signals);
			
			// first dimension: IP or CTRL channel
			// second: condition, 
			// third:strand of this chromosome for this channel and condition,
			// Each list contains five primes of this chromosome for that condition, channel and strand
			List<Integer>[][][] currChromFivePrimesList = new ArrayList[curr_chrom_signals.size()][numConditions][2];
			for(int channel = 0; channel < currChromFivePrimesList.length; channel++)
				for(int t = 0; t < numConditions; t++)
					for(int k = 0; k <= 1; k++)  // '+' or '-' strand
						currChromFivePrimesList[channel][t][k] = new ArrayList<Integer>();
			
			for(int channel = 0; channel < currChromFivePrimesList.length; channel++) {
				for(int t = 0; t < numConditions; t++) {
					// We store the bases of each strand separately
					List<StrandedBase>[] curr_channel_cond_signals = new ArrayList[2];
					curr_channel_cond_signals[0] = new ArrayList<StrandedBase>(); // '+' strand
					curr_channel_cond_signals[1] = new ArrayList<StrandedBase>(); // '-' strand
					
					int used_channel = -1;
					if(channel==0 || (channel==1 && controlDataExist))  // IP or CTRL
						used_channel = channel;
					else if(channel==1 && !controlDataExist)            // There is no CTRL => Use IP instead
						used_channel = channel-1;
					
					for(StrandedBase sb:curr_chrom_signals.get(used_channel).get(t)) {
						if(sb.getStrand() == '+')
							curr_channel_cond_signals[0].add(sb);
						else
							curr_channel_cond_signals[1].add(sb);
					}
					
					for(int k = 0; k <= 1; k++) {
						for(StrandedBase sb:curr_channel_cond_signals[k])
							currChromFivePrimesList[channel][t][k].add(sb.getCoordinate());
						Collections.sort(currChromFivePrimesList[channel][t][k]);
					}
				}
			}
			
			// dim 1: channel, dim 2: condition, dim 3: strand, dim 4: five primes
			int[][][][] currChromFivePrimes = new int[currChromFivePrimesList.length][numConditions][2][];
			for(int channel = 0; channel < currChromFivePrimesList.length; channel++)
				for(int t = 0; t < numConditions; t++)
					for(int k = 0; k <= 1; k++)
						currChromFivePrimes[channel][t][k] = Utils.ref2prim(currChromFivePrimesList[channel][t][k].toArray(new Integer[0]));
			
			chromFivePrimes.put(chrom, currChromFivePrimes);
			int counts = 0;
			// read counts. List 0 for IP, List 1 for CTRL
			// For each of the above lists the read counts for each condition are stored
			List<List<Integer>> currChromCondCounts = new ArrayList<List<Integer>>();
			currChromCondCounts.add(new ArrayList<Integer>()); currChromCondCounts.add(new ArrayList<Integer>());
			for(List<StrandedBase> ip_chrom_signal_cond:curr_chrom_signals.get(0)) {
				int currCondHitCounts = (int)StrandedBase.countBaseHits(ip_chrom_signal_cond); 
				currChromCondCounts.get(0).add(currCondHitCounts);
				counts += currCondHitCounts;
			}
						
			if(controlDataExist) {
				for(List<StrandedBase> ctrl_chrom_signal_cond:curr_chrom_signals.get(1)) {
					int currCondHitCounts = (int)StrandedBase.countBaseHits(ctrl_chrom_signal_cond); 
					currChromCondCounts.get(1).add(currCondHitCounts);
				}
			}
			
			totalIPCounts.put(chrom, counts);
			condHitCounts.put(chrom, currChromCondCounts);
		}//end of Create chromosome statistics
		
		for(ComponentFeature cf:compFeatures) {
			String chrom     = cf.getPosition().getChrom();
			double lambda_bg = (double)totalIPCounts.get(chrom)/gen.getChromLength(chrom); // IP mapped bases
			double mfold = cf.getTotalSumResponsibility()/(lambda_bg*cf.getPeak().expand(modelRange).getWidth());
			cf.set_mfold(mfold);
		}
		
		Collections.sort(compFeatures, new Comparator<ComponentFeature>(){
			public int compare(ComponentFeature o1, ComponentFeature o2) {
				return o1.compareByMfold(o2);
			}
		});
		
		if(num_top_mfold_feats > compFeatures.size()) { num_top_mfold_feats = compFeatures.size(); }
		
		top_mfold_Feats = new ComponentFeature[num_top_mfold_feats];
		int count = 0;
		for(int i = compFeatures.size()-1; i > compFeatures.size() - num_top_mfold_feats-1; i--)
			top_mfold_Feats[count++] = compFeatures.get(i);
		
		double[] shift_size = new double[top_mfold_Feats.length];
		for(int i = 0; i < top_mfold_Feats.length; i++)
			shift_size[i] = eval_shift_size(top_mfold_Feats[i]);
		
		double sum = 0.0;
		int countNonNan = 0;
		for(int i = 0; i < shift_size.length; i++) 
			if(!Double.isNaN(shift_size[i])) { sum += shift_size[i]; countNonNan++;}
		
		avg_shift_size = Math.round(sum/countNonNan);
		return avg_shift_size;
	}//end of eval_avg_shift_size method
	
	private double eval_shift_size(ComponentFeature cf) {
		double shift_size = Double.NaN;
		int plus_mode_pos, minus_mode_pos;
		
		Map<String, StrandedBase> bases_p_map = new HashMap<String, StrandedBase>();
		Map<String, StrandedBase> bases_m_map = new HashMap<String, StrandedBase>();
		StrandedBase[] bases_p;
		StrandedBase[] bases_m;
		
		Region r = cf.getPosition().expand(modelRange);
		ArrayList<List<StrandedBase>> ip_signals = loadBasesInWindow(r);
		
		if(ip_signals != null) {
			for(List<StrandedBase> cond_bases:ip_signals) {
				for(StrandedBase base:cond_bases) {
					String baseKey = String.format("%d\t%c", base.getCoordinate(), base.getStrand());
					if(base.getStrand() == '+') {
						if(!bases_p_map.containsKey(baseKey)) { bases_p_map.put(baseKey, new StrandedBase(base.getStrand(), base.getCoordinate(), 0)); }
						float totalCount = base.getCount() + bases_p_map.get(baseKey).getCount();
						bases_p_map.put(baseKey, new StrandedBase(base.getStrand(), base.getCoordinate(), totalCount));
					}
					else {
						if(!bases_m_map.containsKey(baseKey)) { bases_m_map.put(baseKey, new StrandedBase(base.getStrand(), base.getCoordinate(), 0)); }
						float totalCount = base.getCount() + bases_m_map.get(baseKey).getCount();
						bases_m_map.put(baseKey, new StrandedBase(base.getStrand(), base.getCoordinate(), totalCount));					
					}
				}
			}
			
			bases_p = bases_p_map.values().toArray(new StrandedBase[0]);
			bases_m = bases_m_map.values().toArray(new StrandedBase[0]);
			
			double[] counts_p = new double[bases_p.length];
			for(int i = 0; i < counts_p.length; i++) { counts_p[i] = bases_p[i].getCount(); }
			
			double[] counts_m = new double[bases_m.length];
			for(int i = 0; i < counts_m.length; i++) { counts_m[i] = bases_m[i].getCount(); }
			
			if(counts_p.length != 0 && counts_m.length != 0) {
				plus_mode_pos  = StatUtil.findMax(counts_p).cdr().first();
				minus_mode_pos = StatUtil.findMax(counts_m).cdr().first();
				shift_size = Math.max(0, bases_m[minus_mode_pos].getCoordinate() - bases_p[plus_mode_pos].getCoordinate())/2.0;				
			}
		}
		
		return shift_size;
	}//end eval_shift_size method
	
	
	private double evalFeatureSignificance(ComponentFeature cf, int c) {
		double pVal;
		double local_lambda, lambda_bg, lambda_peak_ip, first_lambda_ip, second_lambda_ip, third_lambda_ip, lambda_peak_ctrl, first_lambda_ctrl, second_lambda_ctrl, third_lambda_ctrl;
		
		String chrom       = cf.getPosition().getChrom();
		int chromLen       = gen.getChromLength(chrom);
		Region peakRegion  = cf.getPeak().expand(modelRange);
		
		List<StrandedBase> ip_signals   = chrom_signals.get(chrom).get(0).get(c);
		List<StrandedBase> ctrl_signals = new ArrayList<StrandedBase>();
		if(controlDataExist)
			ctrl_signals = chrom_signals.get(chrom).get(1).get(c);
		else
			ctrl_signals.addAll(ip_signals);
		
		Collections.sort(ip_signals);
		Collections.sort(ctrl_signals);
				
		double ipCounts = condHitCounts.get(chrom).get(0).get(c);
		double ctrlCounts; 
		if(controlDataExist)		
			ctrlCounts = condHitCounts.get(chrom).get(1).get(c);
		else
			ctrlCounts = ipCounts;
				
		int left_peak           = (int)Math.min(chromLen-1, peakRegion.getStart() + shift_size);
		int left_third_region   = (int)Math.max(         0, cf.getPosition().getLocation() - third_lambda_region_width/2);
		int left_second_region  = (int)Math.max(         0, cf.getPosition().getLocation() - second_lambda_region_width/2);
		int left_first_region   = (int)Math.max(         0, cf.getPosition().getLocation() - first_lambda_region_width/2);
		int right_peak          = (int)Math.max(         0, peakRegion.getEnd()   - shift_size);
		int right_third_region  = (int)Math.min(chromLen-1, cf.getPosition().getLocation() + third_lambda_region_width/2);
		int right_second_region = (int)Math.min(chromLen-1, cf.getPosition().getLocation() + second_lambda_region_width/2);
		int right_first_region  = (int)Math.min(chromLen-1, cf.getPosition().getLocation() + first_lambda_region_width/2);
		
		int num_peak_ip = 0, num_first_lambda_ip = 0, num_second_lambda_ip = 0, num_third_lambda_ip = 0;
		int num_peak_ctrl = 0, num_first_lambda_ctrl = 0, num_second_lambda_ctrl = 0, num_third_lambda_ctrl = 0;
		
		int smallestLeft = (int) Math.min(left_peak, Math.min(left_third_region, Math.min(left_second_region, left_first_region)));
		int largestRight = (int) Math.max(right_peak, Math.max(right_third_region, Math.max(right_second_region, right_first_region)));
		
		// k = 0: '+' strand, k = 1: '-' strand
		for(int k = 0; k <= 1; k++) {
			// IP Channel
			int[] ipFivePrimes = chromFivePrimes.get(chrom)[0][c][k];		
			int s_idx = Arrays.binarySearch(ipFivePrimes, smallestLeft);
			int e_idx = Arrays.binarySearch(ipFivePrimes, largestRight);
			if(s_idx < 0) { s_idx = -s_idx-1; }
			if(e_idx < 0) { e_idx = -e_idx-1; }
			s_idx = StatUtil.searchFrom(ipFivePrimes, ">=", smallestLeft, s_idx);
			e_idx = StatUtil.searchFrom(ipFivePrimes, "<=", largestRight, e_idx);
			for(int i = s_idx; i <= e_idx; i++) {
				if(left_peak <= ip_signals.get(i).getCoordinate() && ip_signals.get(i).getCoordinate() <= right_peak)
					num_peak_ip += ip_signals.get(i).getCount();
				
				if(left_first_region <= ip_signals.get(i).getCoordinate() && ip_signals.get(i).getCoordinate() <= right_first_region) {
					num_first_lambda_ip  += ip_signals.get(i).getCount();
					num_second_lambda_ip += ip_signals.get(i).getCount();
					num_third_lambda_ip  += ip_signals.get(i).getCount();
				}
				else if(left_second_region <= ip_signals.get(i).getCoordinate() && ip_signals.get(i).getCoordinate() <= right_second_region) {
					num_second_lambda_ip += ip_signals.get(i).getCount();
					num_third_lambda_ip  += ip_signals.get(i).getCount();
				}
				else if(left_third_region <= ip_signals.get(i).getCoordinate() && ip_signals.get(i).getCoordinate() <= right_third_region) {
					num_third_lambda_ip  += ip_signals.get(i).getCount();
				}	
			}
			
			// Ctrl Channel
			int[] ctrlFivePrimes = chromFivePrimes.get(chrom)[1][c][k];		
			s_idx = Arrays.binarySearch(ctrlFivePrimes, smallestLeft);
			e_idx = Arrays.binarySearch(ctrlFivePrimes, largestRight);
			if(s_idx < 0) { s_idx = -s_idx-1; }
			if(e_idx < 0) { e_idx = -e_idx-1; }
			s_idx = StatUtil.searchFrom(ctrlFivePrimes, ">=", smallestLeft, s_idx);
			e_idx = StatUtil.searchFrom(ctrlFivePrimes, "<=", largestRight, e_idx);
			for(int i = s_idx; i <= e_idx; i++) {
				if(left_peak <= ctrl_signals.get(i).getCoordinate() && ctrl_signals.get(i).getCoordinate() <= right_peak)
					num_peak_ctrl += ctrl_signals.get(i).getCount();
				
				if(left_first_region <= ctrl_signals.get(i).getCoordinate() && ctrl_signals.get(i).getCoordinate() <= right_first_region) {
					num_first_lambda_ctrl  += ctrl_signals.get(i).getCount();
					num_second_lambda_ctrl += ctrl_signals.get(i).getCount();
					num_third_lambda_ctrl  += ctrl_signals.get(i).getCount();
				}
				else if(left_second_region <= ctrl_signals.get(i).getCoordinate() && ctrl_signals.get(i).getCoordinate() <= right_second_region) {
					num_second_lambda_ctrl += ctrl_signals.get(i).getCount();
					num_third_lambda_ctrl  += ctrl_signals.get(i).getCount();
				}
				else if(left_third_region <= ctrl_signals.get(i).getCoordinate() && ctrl_signals.get(i).getCoordinate() <= right_third_region) {
					num_third_lambda_ctrl  += ctrl_signals.get(i).getCount();
				}	
			}	
		}
		
		double ip2ctrl_ratio = ipCounts/ctrlCounts;
		lambda_bg            = ctrlCounts/chromLen; //numMappedBases.get(chrom).get(1);

		lambda_peak_ip   = 1.0*num_peak_ip;
		first_lambda_ip  = 1.0*num_first_lambda_ip*((double)peakRegion.getWidth()/(double)first_lambda_region_width);
		second_lambda_ip = 1.0*num_second_lambda_ip*((double)peakRegion.getWidth()/(double)second_lambda_region_width);
		third_lambda_ip  = 1.0*num_third_lambda_ip*((double)peakRegion.getWidth()/(double)third_lambda_region_width);
		
		lambda_peak_ctrl   = 1.0*num_peak_ctrl*ip2ctrl_ratio;
		first_lambda_ctrl  = 1.0*num_first_lambda_ctrl*ip2ctrl_ratio*((double)peakRegion.getWidth()/(double)first_lambda_region_width);
		second_lambda_ctrl = 1.0*num_second_lambda_ctrl*ip2ctrl_ratio*((double)peakRegion.getWidth()/(double)second_lambda_region_width);
		third_lambda_ctrl  = 1.0*num_third_lambda_ctrl*ip2ctrl_ratio*((double)peakRegion.getWidth()/(double)third_lambda_region_width);
		lambda_bg          = 1.0*lambda_bg*ip2ctrl_ratio*(double)peakRegion.getWidth();

		if(controlDataExist)
			local_lambda = Math.max(lambda_bg, Math.max(first_lambda_ip, Math.max(second_lambda_ip, Math.max(third_lambda_ip, Math.max(lambda_peak_ctrl, Math.max(first_lambda_ctrl, Math.max(second_lambda_ctrl, third_lambda_ctrl)))))));
		else //skip first lambda regions (1K)
			local_lambda = Math.max(lambda_bg, Math.max(second_lambda_ip, Math.max(third_lambda_ip, Math.max(second_lambda_ctrl, third_lambda_ctrl))));
		
		Poisson poisson = new Poisson(local_lambda, new DRand());
		pVal = 1 - poisson.cdf((int)lambda_peak_ip);
		return pVal;
	}//end of evalFeatureSignificance method
	
	/**
	 * Count non-specific reads  (all the reads not in the refinedRegions)
	 * for scaling experiment and control reads.
	 * Because the signal (specific peak region) can vary a lot in different conditions,
	 * the non-specific regions are a fair estimate for noise.
	 * We assume the noise in expt and control should be comparable, can be used to scale reads.
	 */
	public void countNonSpecificReads(ArrayList<ComponentFeature> compFeatures){
		if (!wholeGenomeDataLoaded){	// default ratio = 1;
			ComponentFeature.setNon_specific_ratio(ratio_non_specific_total);
			return;
		}
//		System.out.println("\nCounting non-specific read numbers ...\n");
		
		//Count the specific read numbers
		int expt_test_region_total[]=new int[numConditions];
		int crtl_test_region_total[]=new int[numConditions];
		int expt_non_specific_total[]=new int[numConditions];
		int crtl_non_specific_total[]=new int[numConditions];
		int totalLength=0;	// total length of non-overlapping peak regions
		for(Region r : restrictRegions){
			totalLength += r.getWidth();
			for(int i=0; i<caches.size(); i++){
				Pair<ReadCache,ReadCache> e = caches.get(i);
				expt_test_region_total[i] += e.car().countHits(r);
				if(controlDataExist)
					crtl_test_region_total[i] += e.cdr().countHits(r);
			}
		}
		log(1, "\nPeak regions total length: "+totalLength);
		
		// non-specific = total - specific
		for(int i=0; i<caches.size(); i++){
			Pair<ReadCache,ReadCache> e = caches.get(i);
			expt_non_specific_total[i]=(int)e.car().getHitCount()-expt_test_region_total[i];
			if(controlDataExist) {
				crtl_non_specific_total[i]=(int)e.cdr().getHitCount()-crtl_test_region_total[i];
				ratio_non_specific_total[i] = (double)expt_non_specific_total[i]/crtl_non_specific_total[i];
			}
	    	
	    	double noiseReadNum_per_kbp = (double) expt_non_specific_total[i] * 1000
	    		/ (mappable_genome_length - totalLength);
	    	
	    	log(1, conditionNames.get(i)+" read stats:");
	    	log(1, "Signal total\t\t" +(int)e.car().getHitCount()+
	    	              "\nControl total\t\t" + (controlDataExist ? (int)e.cdr().getHitCount() : 0 )+
	    	              "\nRatio total\t\t" +String.format("%.3f", (controlDataExist ? e.car().getHitCount()/e.cdr().getHitCount() : 0 ))+
	    	              "\nSignal non-specific\t" +expt_non_specific_total[i]+
	    	              "\nControl non-specific\t" +crtl_non_specific_total[i]+
	    	              "\nRatio non-specific\t" +String.format("%.3f",ratio_non_specific_total[i])+
	    	              "\nNoise reads per 1000bp\t" +String.format("%.2f",noiseReadNum_per_kbp)+
	    	              "\n");
		}
		ComponentFeature.setNon_specific_ratio(ratio_non_specific_total);

//		log(1, "countNonSpecificReads(): "+timeElapsed(tic));		
	}
	
	private void loadPeakShapeParameters(String peak_shape_distribution) {
		try{
			BufferedReader file = new BufferedReader(new FileReader(new File(peak_shape_distribution)));
			String[] readNums = file.readLine().trim().split("\t");
			String[] means = file.readLine().trim().split("\t");
			String[] stds = file.readLine().trim().split("\t");
			shapeParameters = new HashMap<Integer, Pair<Double, Double>>();
			for (int i=0;i<readNums.length;i++){
				shapeParameters.put(Integer.parseInt(readNums[i]), 
						new Pair<Double, Double>(Double.parseDouble(means[i]), 
								Double.parseDouble(stds[i])));
			}
			file.close();
		}catch(IOException e){
			System.err.println(e.toString());
		}
		
	}
	private void addConfigString(String name, boolean value){
		config.append(name).append("\t").append(new Boolean(value).toString()).append("\n");
	}
	public String getConfigString(){
		return "BindingMixture Configurations\n"+config.toString();
	}
		
	// rank for peak-shape
	private void rankPeakShape(ArrayList<ComponentFeature> compFeatures){
		for (int cond=0; cond<conditionNames.size(); cond++){
			ComponentFeature.setSortingCondition(cond);
			Collections.sort(compFeatures, new Comparator<ComponentFeature>(){
			    public int compare(ComponentFeature o1, ComponentFeature o2) {
			        return o1.compareByPeakShape(o2);
			    }
			});
	
			int rank =1;
			for(ComponentFeature cf : compFeatures){
				cf.addRank_sum(rank);
				rank++;
			}
		}
	}
	
	public void addAnnotations(){
		//Add closest genes
		System.out.println("Adding gene annotations");
		addClosestGenes(signalFeatures);

		//Add other annotations
		System.out.println("Adding other annotations");
		addRegionAnnotations(signalFeatures);
	}
	
	public void printFeatures(){
		ArrayList<ComponentFeature> fs = new ArrayList<ComponentFeature>();
		for(Feature f : signalFeatures){
			fs.add((ComponentFeature)f);
		}
		String fname = outName+"_GPS_events.txt"; 
		printFeatures(fname, fs);		
	}
	
	public void printInsignificantFeatures(){
		ArrayList<ComponentFeature> fs = new ArrayList<ComponentFeature>();
		for(Feature f : insignificantFeatures){
			fs.add((ComponentFeature)f);
		}
		String fname = outName+"_Insignificant_events.txt"; 
		printFeatures(fname, fs);
	}
	
	public void printPsortedCondFeatures(){
		for(int c = 0; c < numConditions; c++) {
			ComponentFeature.setSortingCondition(c);
			ArrayList<ComponentFeature> fs = new ArrayList<ComponentFeature>();
			for (Feature f:condSignalFeats[c])
				fs.add((ComponentFeature)f);
				
			if(controlDataExist) {
				Collections.sort(fs, new Comparator<ComponentFeature>(){
				    public int compare(ComponentFeature o1, ComponentFeature o2) {
				        return o1.compareByPValue(o2);
				    }
				});				
			}
			else {
				Collections.sort(fs, new Comparator<ComponentFeature>(){
				    public int compare(ComponentFeature o1, ComponentFeature o2) {
				        return o1.compareByPValue_wo_ctrl(o2);
				    }
				});
			}

			String fname = outName + "_" + conditionNames.get(c) + "_GPS_events_P.txt"; 
			printFeatures(fname, fs);
		}
	}

	public void printPsortedFeatures(){
		ComponentFeature.setSortingCondition(0);
		ArrayList<ComponentFeature> fs = new ArrayList<ComponentFeature>();
		for (Feature f:signalFeatures){
			fs.add((ComponentFeature)f);
		}
		if(controlDataExist) {
			Collections.sort(fs, new Comparator<ComponentFeature>(){
			    public int compare(ComponentFeature o1, ComponentFeature o2) {
			        return o1.compareByPValue(o2);
			    }
			});				
		}
		else {
			Collections.sort(fs, new Comparator<ComponentFeature>(){
			    public int compare(ComponentFeature o1, ComponentFeature o2) {
			        return o1.compareByPValue_wo_ctrl(o2);
			    }
			});
		}

		String fname = outName+"_GPS_events_P.txt"; 
		printFeatures(fname, fs);
	}	
	
	private void printFeatures(String fname, ArrayList<ComponentFeature> fs){
		try{
			FileWriter fw = new FileWriter(fname);
			boolean first=true;
			for(ComponentFeature f : fs){
				if(first){
					fw.write(development_mode?f.headString():f.headString_v1()); 
					first=false;
				}
				fw.write(development_mode?f.toString():f.toString_v1());
			}
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}		
	}
	public void printExpandedPeaks(int range){
		try{
			String fname = outName+"_peaks_"+range+"bp.txt"; 
			FileWriter fw = new FileWriter(fname);
			for(Feature f : signalFeatures){
				fw.write(f.getPeak().expand(range/2).toString()+"\n");
			}
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	private void printNoneZeroRegions(boolean initial){
		// save the list of regions to file
		try{
			StringBuilder fileName = new StringBuilder(outName).append("_");
			if (initial)
				fileName.append("Init_");
			else
				fileName.append("Refined_");
			
			fileName.append("Regions.txt");
			FileWriter fw = new FileWriter(fileName.toString());
			
			StringBuilder txt = new StringBuilder();
			for (Region r:restrictRegions){
				txt.append(r.toString()).append("\t").append(r.getWidth()).append("\t");
				txt.append(countIpReads(r)).append("\n");
			}
			fw.write(txt.toString());
			fw.close();
		} 
		catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private float countIpReads(Region r){
		float count=0;
		for(Pair<ReadCache,ReadCache> e : caches){
			List<StrandedBase> bases_p= e.car().getStrandedBases(r, '+');  
			List<StrandedBase> bases_m= e.car().getStrandedBases(r, '-');  
			count+=StrandedBase.countBaseHits(bases_p)+StrandedBase.countBaseHits(bases_m);
		}
		return count;
	}
	private float countIpReads(Region r, int cond){
		List<StrandedBase> bases_p= caches.get(cond).car().getStrandedBases(r, '+');  
		List<StrandedBase> bases_m= caches.get(cond).car().getStrandedBases(r, '-');  
		return StrandedBase.countBaseHits(bases_p)+StrandedBase.countBaseHits(bases_m);
	}
	private float countCtrlReads(Region r){
		float count=0;
		for(Pair<ReadCache,ReadCache> e : caches){
			List<StrandedBase> bases_p= e.cdr().getStrandedBases(r, '+');  // reads of the current region - IP channel
			List<StrandedBase> bases_m= e.cdr().getStrandedBases(r, '-');  // reads of the current region - IP channel
			count+=StrandedBase.countBaseHits(bases_p)+StrandedBase.countBaseHits(bases_m);
		}
		return count;
	}	
	private float countCtrlReads(Region r, int cond){
		List<StrandedBase> bases_p= caches.get(cond).cdr().getStrandedBases(r, '+');  
		List<StrandedBase> bases_m= caches.get(cond).cdr().getStrandedBases(r, '-');  
		return StrandedBase.countBaseHits(bases_p)+StrandedBase.countBaseHits(bases_m);
	}
	private void printTowerRegions(){
		// save the list of regions to file
		try{
			StringBuilder fileName = new StringBuilder(outName).append("_");
//			for (String cond: conditionNames){
//				fileName.append(cond).append("_");
//			}
			fileName.append("TowerRegions.txt");
			FileWriter fw = new FileWriter(fileName.toString());
			
			StringBuilder txt = new StringBuilder();
			for (int i=0;i<towerRegions.size();i++){
				txt.append(towerRegions.get(i).toString()).append("\t");
				txt.append(towerStrength.get(i)).append("\n");
			}
			fw.write(txt.toString());
			fw.close();
		} 
		catch (IOException e) {
			e.printStackTrace();
		}
	}	
	public static String timeElapsed(long tic){
		return timeString(System.currentTimeMillis()-tic);
	}
	private static String timeString(long length){
		float sec = length/1000F;
		return sec>60? 
			String.format("%.1f",sec/60)+" min":
			String.format("%.1f",sec)+" sec";
	}
	public static String getDateTime() {
        DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH-mm-ss");
        Date date = new Date();
        return dateFormat.format(date);
    } 
	public static void printArray(double[]array, String msgBefore, String msgAfter){
		System.out.print(msgBefore);
		System.out.print(arrayToString(array));
		System.out.print(msgAfter);
	}
	public static void printArray(int[]array, String msgBefore, String msgAfter){
		System.out.print(msgBefore);
		System.out.print(arrayToString(array));
		System.out.print(msgAfter);
	}
	public static String arrayToString(int[] array){
        StringBuilder output = new StringBuilder();
        for (int i=0;i<array.length-1;i++){
        	output.append(String.format("%d\t",array[i]));
        }
        output.append(String.format("%d",array[array.length-1]));
        return output.toString();
	}
	public static String arrayToString(double[] array){
        StringBuilder output = new StringBuilder();
        for (int i=0;i<array.length-1;i++){
        	output.append(String.format("%.2f\t",array[i]));
        }
        output.append(String.format("%.2f",array[array.length-1]));
        return output.toString();
	}
	
	public void writeDebugFile(){
		try{
			FileWriter fw = new FileWriter("Binding Mixture Debug.txt", false); //new file
			fw.write(log_all_msg.toString());
			fw.close();
		} 
		catch (IOException e) {
			e.printStackTrace();
		}
	}	
	
	// show the message in STDOUT and log file
	// depending on "verbose" setting in property file
	private void log(int mode, String msg){
		if (LOG_ALL)
			log_all_msg.append(msg);
		if ( bmverbose>=mode){
	    	System.out.println(msg);
			try{
				logFileWriter.write(msg+"\n");
				logFileWriter.flush();
			} 
			catch (IOException e) {
				e.printStackTrace();
			}
			if (LOG_ALL)
				log_all_msg.append("\n");
		}
	}

	public void appendEMMixProbFile(int resolution){
		try{
			FileWriter fw = new FileWriter(this.outName+"_PI Dump_"+resolution+"bp.txt", true); //new file
			fw.write(pi_sb.toString());
			fw.close();
			pi_sb = new StringBuilder();
		} 
		catch (IOException e) {
			e.printStackTrace();
		}
	}	
	public void closeLogFile(){
		try{
			logFileWriter.close();
		} 
		catch (IOException e) {
			e.printStackTrace();
		}
	}
	public static void writeFile(String fileName, String text){
		try{
			FileWriter fw = new FileWriter(fileName, false); //new file
			fw.write(text);
			fw.close();
		} 
		catch (IOException e) {
			e.printStackTrace();
		}
//		System.out.println("File was written to "+fileName);
	}	
	public void printError() {
		
	}
	class EM_State{
		double[][][] resp;
		double[][]beta;
		double[] pi;
		double alpha;
		
		double LL;
		double numComponent;
		EM_State(int numCond){
			resp=new double[numCond][][];
			beta=new double[numCond][];
		}
		// BIC=LL-#param/2*ln(n)
		// # param: Each component has 2 parameters, mixing prob and position, thus "*2";
		// "-1" comes from the fact that total mix prob sum to 1.
		// n: is the number of data point, i.e. the count of reads summing over all base positions.
		double BIC(double n){
			return LL - (numComponent*2-1)/2*Math.log(n);
		}
		public String toString(){
			return String.format("%.3f\t%.0f", LL, numComponent);
		}
	}
}


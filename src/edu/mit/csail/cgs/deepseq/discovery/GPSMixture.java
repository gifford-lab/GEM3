package edu.mit.csail.cgs.deepseq.discovery;

import java.io.*;
import java.util.*;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import javax.imageio.ImageIO;

import cern.jet.random.Poisson;
import cern.jet.random.Binomial;
import cern.jet.random.engine.DRand;

import edu.mit.csail.cgs.datasets.chipseq.ChipSeqLocator;
import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.*;
import edu.mit.csail.cgs.deepseq.features.*;
import edu.mit.csail.cgs.deepseq.multicond.MultiIndependentMixtureCounts;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.deepseq.utilities.ReadCache;
import edu.mit.csail.cgs.ewok.verbs.chipseq.MACSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.MACSPeakRegion;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Utils;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.models.data.DataFrame;
import edu.mit.csail.cgs.utils.models.data.DataRegression;
import edu.mit.csail.cgs.utils.probability.NormalDistribution;
import edu.mit.csail.cgs.utils.stats.StatUtil;


class GPSMixture extends MultiConditionFeatureFinder {

    private GPSConfig config;
    private GPSConstants constants;

	// Binding model representing the empirical distribution of a binding event
    // only updated in updateBindingModel, which is outside the EM loop
	private BindingModel model;
	private int modelRange;			// max length of one side
	private int modelWidth; 		// width of empirical read distribution, 500bp

	private double[] profile_plus_sum;
	private double[] profile_minus_sum;
	
	private HashMap<String, BindingModel> allModels = new HashMap<String, BindingModel>();
	
	// Max number of reads on a base position, to filter out towers and needles.
	private int max_HitCount_per_base = 3;

	//Gaussian kernel for density estimation
	private double[] gaussian;
	private boolean doScanning=true;

	/****************
	 * Data
	 ****************/
	private boolean wholeGenomeDataLoaded = false;
	// memory cache to store all the read data, loaded from DB or file
	private ArrayList<Pair<ReadCache, ReadCache>> caches;
	private ArrayList<String> conditionNames = new ArrayList<String>();
	// Do we have matched control data?
	private boolean controlDataExist=false;
	// Ratio of each pair of IP/Ctrl for all conditions
	private double[] ratio_total;
	private double[] ratio_non_specific_total;
	private double shift_size;
	/****************
	 * Prediction
	 ****************/
	// List of all EM predictions
	private List<ComponentFeature> allFeatures;
	private ArrayList<Feature> insignificantFeatures;	// does not pass statistical significance
	private ArrayList<Feature> filteredFeatures;	// filtered artifacts (bad shape or low fold enrichment)
	private List<Feature>[] condSignalFeats; //Signal channel features in each condition

	//Regions to be evaluated by EM, estimated whole genome, or specified by a file
	private ArrayList<Region> restrictRegions = new ArrayList<Region>();
	//Regions to be excluded, supplied by user
	private ArrayList<Region> excludedRegions = new ArrayList<Region>();

	// Regions that have towers (many reads stack on same base positions)
	private ArrayList<Region> towerRegions = new ArrayList<Region>();
	private ArrayList<Float> towerStrength = new ArrayList<Float>();

	//Number of reads for a specified region of the IP experiment for each condition
	private double[] sigHitCounts;
	//(Total) number of reads for a specified region of the IP experiment (across all conditions)
	private double totalSigCount=0;

	/** reads of this chromosome (position and count) as <tt>StrandedBases</tt> for this channel, condition and strand
	 * key:chrom
	 * value:{first dimension: IP or CTRL channel, second: condition, third:strand}  */
    //	private Map<String, List<StrandedBase>[][][]> chromUniqueFivePrimes;

	/** key:chrom, value:{first dimension: IP or CTRL channel, second: condition,
	 * third:strand, fourth:unique positions of reads of this chromosome for this channel, condition and strand} */
    //	private Map<String, int[][][][]> chromUniqueFivePrimePos;

	/** Total IP counts for each chromosome (summed over all conditions) */
	private Map<String, Integer> totalIPCounts;

	/** For each chromosome, the counts of reads for each channel and condition are stored.     <br>
        /** key:chrom, value:{first list: IP or CTRL channel, second list: condition */
	private Map<String, List<List<Integer>>> condHitCounts;

	/** Contains the StrandedBases of the current chromosome for the IP chanel <br>
	 *  1x2 array list. 1st element: a list for strand '+', 2nd element: a list for strand '-' */
	private List<StrandedBase>[] ipStrandFivePrimes;
	
	/** A 2-D matrix containing the five prime positions of IP channel the current chromosome.
	 * Row 1: five primes for strand '+', Row 2: five primes for strand '-' */
	private int[][] ipStrandFivePrimePos;
	
	/** Contains the StrandedBases of the current chromosome for the CTRL chanel <br>
	 *  1x2 array list. 1st element: a list for strand '+', 2nd element: a list for strand '-' */
	private List<StrandedBase>[] ctrlStrandFivePrimes;
	
	/** A 2-D matrix containing the five prime positions of CTRL channel the current chromosome.
	 * Row 1: five primes for strand '+', Row 2: five primes for strand '-' */
	private int[][] ctrlStrandFivePrimePos;

    private StringBuilder configsb = new StringBuilder();
	private StringBuilder log_all_msg = new StringBuilder();
	private FileWriter logFileWriter;

	public GPSMixture(Genome g, 
                      ArrayList<Pair<DeepSeqExpt,DeepSeqExpt>> expts,
                      ArrayList<String> conditionNames,
                      String[] args) {
		super (args, g, expts);

		try{
			logFileWriter = new FileWriter("GPS_Log.txt", true); //append
			logFileWriter.write("\n==============================================\n");
			logFileWriter.write(CommonUtils.getDateTime());
			logFileWriter.flush();
		} catch (IOException e) {
			e.printStackTrace();
		}

		/* ***************************************************
		 * Load Binding Model, empirical distribution
		 * ***************************************************/
		String modelFile = Args.parseString(args, "d", null);	// read distribution file

		commonInit(modelFile);
		model.printToFile(outName+"_0_Read_distribution.txt");
		allModels.put(outName+"_0", model);

		/* ***************************************************
		 * Load parameters and properties
		 * ***************************************************/
		StringBuffer sb = new StringBuffer();
		sb.append("\nOptions:\n");
		for (String arg:args){
			if (arg.trim().indexOf(" ")!=-1)
				sb.append("\"").append(arg).append("\" ");
			else
				sb.append(arg).append(" ");
		}
		log(1, sb.toString());
		
    	/* *********************************
    	 * Flags
    	 ***********************************/
        
        config.parseArgs(args);             	
    	if(config.second_lambda_region_width < config.first_lambda_region_width) {
    		System.err.println("\nThe first control region width (w2) has to be more than " + config.first_lambda_region_width + " bp.");
    		System.exit(-1);
    	}    	
    	if(config.third_lambda_region_width < config.second_lambda_region_width) {
    		System.err.println("\nThe second control region width (w3) has to be more than " + config.second_lambda_region_width + " bp.");
    		System.exit(-1);
    	}
   
		/* **************************************************
		 * Determine the subset of regions to run EM.
 		 * It can be specified as a file from command line.
		 * If no file pre-specified, estimate the enrichedRegions after loading the data.
		 * We do not want to run EM on regions with little number of reads
		 * **************************************************/
    	String subset_str = Args.parseString(args, "subs", null);
    	double subsetRatio = Args.parseDouble(args, "subr", -1);		// user input ratio for IP/control, because if the region is too small, the non-specific definition is not applied
     	String subsetFormat = Args.parseString(args, "subFormat", "");
     	boolean toSegmentRegions = !subsetFormat.equals("Regions");
    	if(subset_str == null) { subset_str = Args.parseString(args, "subf", null); }
    	ArrayList<Region> subsetRegions = new ArrayList<Region>();
     	if (subset_str!=null && Args.parseString(args, "subf", null) != null) {
	    	if (subsetFormat.equals("Points")){	// To expand
	    		subsetRegions = CommonUtils.loadRegionFile(subset_str, gen);
	    		subsetRegions = mergeRegions(subsetRegions, true);
	    	}
	    	else {
	    		subsetRegions = CommonUtils.loadRegionFile(subset_str, gen);
	    		subsetRegions = mergeRegions(subsetRegions, false);
	    	}
        } else if (subset_str!=null && Args.parseString(args, "subs", null) != null) {
     		String COMPLETE_REGION_REG_EX = "^\\s*([\\w\\d]+):\\s*([,\\d]+[mMkK]?)\\s*-\\s*([,\\d]+[mMkK]?)\\s*";
     		String CHROMOSOME_REG_EX = "^\\s*([\\w\\d]+)\\s*$";
     		String[] reg_toks = subset_str.split("\\s");
     		for(String regionStr:reg_toks) {
     			
     			// Region already presented in the form x:xxxx-xxxxx
     			if(regionStr.matches(COMPLETE_REGION_REG_EX)) {
     				subsetRegions.add(Region.fromString(gen, regionStr));
     			}
     			
     			// Check if it is about a whole chromosome
     			else if(regionStr.matches(CHROMOSOME_REG_EX)) {
     				regionStr = regionStr.replaceFirst("^chromosome", "");
     				regionStr = regionStr.replaceFirst("^chrom", "");
     				regionStr = regionStr.replaceFirst("^chr", "");
     					
     	     		for(String chrom:gen.getChromList()) {
     	     			String chromNumber = chrom;
     	     			chromNumber = chromNumber.replaceFirst("^chromosome", "");
     	     			chromNumber = chromNumber.replaceFirst("^chrom", "");
     	     			chromNumber = chromNumber.replaceFirst("^chr", "");

     	     			if(regionStr.equalsIgnoreCase(chromNumber)) {
     	     				subsetRegions.add(new Region(gen, chrom, 0, gen.getChromLength(chrom)-1));
     	     				break;
     	     			}
     	     		}
     			}
     		}//end of for(String regionStr:reg_toks) LOOP
     		subsetRegions = mergeRegions(subsetRegions, false);
     	}
     	
		/* ***************************************************
		 * ChIP-Seq data
		 * ***************************************************/
		this.caches = new ArrayList<Pair<ReadCache, ReadCache>>();
		this.numConditions = experiments.size();
		this.conditionNames = conditionNames;
		boolean fromReadDB = experiments.get(0).car().isFromReadDB();
		ComponentFeature.setConditionNames(conditionNames);
		condSignalFeats = new ArrayList[numConditions];
		for(int c = 0; c < numConditions; c++) { condSignalFeats[c] = new ArrayList<Feature>(); }
		long tic = System.currentTimeMillis();
		System.out.println("\nGetting 5' positions of all reads...");

		for (int i=0;i<numConditions;i++){
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
			if (fromReadDB){		// load from read DB
				System.out.println("\nLoading data from ReadDB ...");
				if (subsetRegions.isEmpty()){		// if whole genome

					for (String chrom: gen.getChromList()){
						// load  data for this chromosome.
						int length = gen.getChromLength(chrom);
						Region wholeChrom = new Region(gen, chrom, 0, length-1);
						int count = Math.max(ip.countHits(wholeChrom), ctrl.countHits(wholeChrom));
						ArrayList<Region> chunks = new ArrayList<Region>();
						// if there are too many reads in a chrom, read smaller chunks
						if (count>constants.MAXREAD){
							int chunkNum = count/constants.MAXREAD+1;
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
				} //if (subsetRegions.isEmpty()){
				else{	// if subsetRegions notEmpty, only load these regions
					for (Region region: subsetRegions){
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
                //				ipCache.displayStats();
				if (controlDataExist){
					ctrlCache.populateArrays();
                    //					ctrlCache.displayStats();
				}
			}
			else if (!fromReadDB){		// load from File
				ipCache.addAllFivePrimes(ip.getAllStarts());
				ipCache.populateArrays();
				if (controlDataExist){
					ctrlCache.addAllFivePrimes(ctrl.getAllStarts());
					ctrlCache.populateArrays();
				}
				wholeGenomeDataLoaded = true;
			}
			// cleanup			// clean up connection if it is readDB, cleanup data object if it is file
			ip.closeLoaders();
			ctrl.closeLoaders();
			ip=null;
			ctrl=null;
			System.gc();
		} // for each condition

		if (fromReadDB){
			System.out.println("Finish loading data from ReadDB, " + CommonUtils.timeElapsed(tic));
			System.out.println();
		}
		
		// exclude some regions
     	String excludedName = Args.parseString(args, "ex", "yes");
     	if (!excludedName.equals("yes")){
     		excludedRegions = mergeRegions(CommonUtils.loadRegionFile(excludedName, gen), false);
     		log(1, "Exclude " + excludedRegions.size() + " regions.");
			for(int c = 0; c < numConditions; c++) {
				caches.get(c).car().excludeRegions(excludedRegions);
				if(controlDataExist) {
					caches.get(c).cdr().excludeRegions(excludedRegions);
				}
			}     	
     	}

		// Filtering/reset bases 
		doBaseFiltering();			 		
		
		// print resulting dataset counts
		for(int c = 0; c < numConditions; c++) {
			caches.get(c).car().displayStats();
			if(controlDataExist) {
				caches.get(c).cdr().displayStats();
			}
		}
		
		// Normalize conditions
		normExpts(caches);
				
		// estimate ratio and segment genome into regions
    	ratio_total=new double[numConditions];
    	ratio_non_specific_total = new double[numConditions];
        sigHitCounts=new double[numConditions];
        seqwin=100;
        if (wholeGenomeDataLoaded){		// estimate some parameters if whole genome data
	        // estimate IP/Ctrl ratio from all reads
        	for(int i=0; i<numConditions; i++){
				Pair<ReadCache,ReadCache> e = caches.get(i);
				double ipCount = e.car().getHitCount();				
				double ctrlCount = 0;
				if (e.cdr()!=null){
					ctrlCount = e.cdr().getHitCount();
					ratio_total[i]=ipCount/ctrlCount;
					System.out.println(String.format("\n%s\tIP: %.0f\tCtrl: %.0f\t IP/Ctrl: %.2f", conditionNames.get(i), ipCount, ctrlCount, ratio_total[i]));
				}
	        }
        }
        else{	// want to analyze only specified regions, set default
        	for(int i=0; i<numConditions; i++){
        		if (subsetRatio!=-1){
        			ratio_total[i]=subsetRatio;
        			ratio_non_specific_total[i]=subsetRatio;
        		}
        		else {
	        		Pair<ReadCache,ReadCache> e = caches.get(i);
	        		if (e.cdr()!=null)
						ratio_total[i]=e.car().getHitCount()/e.cdr().getHitCount();
	        		else
	        			ratio_total[i]=1;
	        		ratio_non_specific_total[i]=1;
        		}
        	}
        }        	
        System.out.println("\nSorting reads and selecting enriched regions ...");
    	// if no focus list, directly estimate candidate regions from data
		if (wholeGenomeDataLoaded || toSegmentRegions){
     		// ip/ctrl ratio by regression on non-enriched regions
			if (subsetRatio==-1){
     			setRegions(selectEnrichedRegions(subsetRegions, false));
    			calcIpCtrlRatio(restrictRegions);
    			if(controlDataExist) {
    				for(int t = 0; t < numConditions; t++)
    					System.out.println(String.format("For condition %s, IP/Control = %.2f", conditionNames.get(t), ratio_non_specific_total[t]));
    			}
    		}
			if (config.subtract_for_segmentation)
				setRegions(selectEnrichedRegions(subsetRegions, true));
			else
				setRegions(selectEnrichedRegions(subsetRegions, false));
		}
		else{
			setRegions(subsetRegions);
		}
		
		log(1, "\nThe genome is segmented into "+restrictRegions.size()+" regions for analysis.");

        
        
        // exclude regions?
        //        if (excludedRegions!=null){
        //	        	ArrayList<Region> toRemove = new ArrayList<Region>();
        //	        ArrayList<Region> toAdd = new ArrayList<Region>();
        //	
        //	    	for (Region ex: excludedRegions){
        //	    		for (Region r:restrictRegions){
        //	        		if (r.overlaps(ex)){
        //	        			toRemove.add(r);
        //	        			if (r.getStart()<ex.getStart())
        //	        				toAdd.add(new Region(gen, r.getChrom(), r.getStart(), ex.getStart()-1));
        //	        			if (r.getEnd()>ex.getEnd())
        //	        				toAdd.add(new Region(gen, r.getChrom(), ex.getEnd()+1, r.getEnd()));
        //	        		}
        //	        	}
        //	        }
        //	        restrictRegions.addAll(toAdd);
        //	        restrictRegions.removeAll(toRemove);
        //        }
        
		log(2, "BindingMixture initialized. "+numConditions+" conditions.");
		experiments = null;
		System.gc();
	}//end of BindingMixture constructor

	/**
	 * This constructor is only called from Robustness analysis
	 * The read data of each sub-sampling will be suppled with simpleExecute() call
	 */
	public GPSMixture(String[] args, boolean doScanning){
		super (args);

		// ChIP-Seq data
		this.doScanning = doScanning;
	    controlDataExist = false;
	    numConditions = 1;

	    String modelFile = Args.parseString(args, "d", null);
		commonInit(modelFile);

        // the configuration of mixture model
        addConfigString("Do ML scanning\t", doScanning);
        addConfigString("Batch elimination", constants.BATCH_ELIMINATION);
        addConfigString("Smart inital event spacing", constants.SMART_SPACING);
        //        addConfigString("Make hard assignment", properties.make_hard_assignment);
        System.out.println(getConfigString());
	}

	private void commonInit(String modelFile){
        constants = new GPSConstants();
        config = new GPSConfig();

		//Load Binding Model
		File pFile = new File(modelFile);
		if(!pFile.isFile()){
			System.err.println("\nCannot find read distribution file!");
			System.exit(1);
		}
        model = new BindingModel(pFile);
        modelWidth = model.getWidth();
        modelRange = model.getRange();
        config.windowSize = modelWidth * config.window_size_factor;

        //init the Guassian kernel prob. for smoothing the read profile of called events
		gaussian = new double[modelWidth];
		NormalDistribution gaussianDist = new NormalDistribution(0, constants.READ_KERNEL_ESTIMATOR_WIDTH*constants.READ_KERNEL_ESTIMATOR_WIDTH);
		for (int i=0;i<gaussian.length;i++)
			gaussian[i]=gaussianDist.calcProbability((double)i);
	}

	/**
	 * Calls EM training on the Mixture Model.
	 * In this implementation, we have chosen to run EM independently on each testRegion.
	 */
	public List<Feature> execute() {
		long tic = System.currentTimeMillis();
		
		int totalRegionCount = restrictRegions.size();
		if (totalRegionCount==0)
			return new ArrayList<Feature>();
		
		// refresh the total profile sums every round
		profile_plus_sum = new double[modelWidth];
		profile_minus_sum = new double[modelWidth];
		
		// if we have predicted events from previous round, setup the restrictRegions
		// because when we update the model, the model range might changed.
		if (config.refine_regions){
			if (allFeatures!=null && (!allFeatures.isEmpty())){
				Collections.sort(allFeatures);
				ArrayList<Region> potentialRegions = new ArrayList<Region>();
				for (ComponentFeature cf:allFeatures){
					potentialRegions.add(cf.getPosition().expand(0));
				}
				// expand with modelRange, and merge overlapped regions ==> Refined enriched regions
				this.restrictRegions=mergeRegions(potentialRegions, true);
			}
		}
        signalFeatures.clear();
        Vector<ComponentFeature> compFeatures = new Vector<ComponentFeature>();
        Vector<Integer> processRegionCount = new Vector<Integer>();		// for counting how many regions are processed by all threads
		int displayStep = (int) Math.pow(10, (int) (Math.log10(totalRegionCount)));
		TreeSet<Integer> reportTriggers = new TreeSet<Integer>();
		for (int i=1;i<=totalRegionCount/displayStep; i++){
			reportTriggers.add(i*displayStep);
		}
		reportTriggers.add(100);
		reportTriggers.add(1000);
		reportTriggers.add(10000);
		
        Thread[] threads = new Thread[maxThreads];
        log(1,String.format("Creating %d threads", maxThreads));
        for (int i = 0 ; i < threads.length; i++) {
            ArrayList<Region> threadRegions = new ArrayList<Region>();
            for (int j = i; j < restrictRegions.size(); j += maxThreads) {
                threadRegions.add(restrictRegions.get(j));
            }
            Thread t = new Thread(new GPSThread(threadRegions,
            									processRegionCount,
                                                compFeatures,
                                                this,
                                                constants,
                                                config));
            t.start();
            threads[i] = t;
        }
        boolean anyrunning = true;
        while (anyrunning) {
            anyrunning = false;
            try {
                Thread.sleep(1000);
            } catch (InterruptedException e) { }
            for (int i = 0; i < threads.length; i++) {
                if (threads[i].isAlive()) {
                    anyrunning = true;
                    break;
                }
            }    
            int count = processRegionCount.size();
            int trigger = totalRegionCount;
            if (!reportTriggers.isEmpty())
            	trigger = reportTriggers.first();
            if (count>trigger){
				System.out.println(trigger+"\t/"+totalRegionCount+"\t"+CommonUtils.timeElapsed(tic));
				reportTriggers.remove(reportTriggers.first());
            }
        }
        System.out.println(totalRegionCount+"\t/"+totalRegionCount+"\t"+CommonUtils.timeElapsed(tic));
        processRegionCount.clear();
        
        log(1,String.format("%d threads have finished running", maxThreads));

		for (int j=0;j<restrictRegions.size();j++) {
			Region rr = restrictRegions.get(j);
		}// end of for (Region rr : restrictRegions)
		
		/* ********************************************************
		 * refine the specific regions that contain binding events
		 * ********************************************************/
		Collections.sort(compFeatures);
		ArrayList<Region> refinedRegions = new ArrayList<Region>();
		for (ComponentFeature cf:compFeatures){
			refinedRegions.add(cf.getPosition().expand(0));
		}
		// expand with modelRange, and merge overlapped regions ==> Refined enriched regions
		if (config.refine_regions)
			this.restrictRegions=mergeRegions(refinedRegions, true);

		log(2, "Finish predicting events: "+CommonUtils.timeElapsed(tic)+"\n");

		// post processing
		postEMProcessing(compFeatures);
		
		for(int c = 0; c < numConditions; c++)
			condSignalFeats[c] = condPostFiltering(signalFeatures, c);

		log(1, "Finish prediction: "+CommonUtils.timeElapsed(tic)+"\n");

		return signalFeatures;
	}// end of execute method

	// split a region into smaller windows if the width is larger than windowSize
	private ArrayList<Region> splitWindows(Region r, int windowSize, int overlapSize){
		ArrayList<Region> windows=new ArrayList<Region>();
		if (r.getWidth()<=windowSize)
			return windows;
		List<StrandedBase> allBases = new ArrayList<StrandedBase>();
		for (int c=0; c<caches.size(); c++){
			ReadCache ip = caches.get(c).car();
			List<StrandedBase> bases = ip.getUnstrandedBases(r);
			if (bases==null || bases.size()==0)
				continue;
			allBases.addAll(bases); // pool all conditions
		}
		Collections.sort(allBases);
		
		int[] distances = new int[allBases.size()-1];
		int[]coords = new int[allBases.size()];
		for (int i=0;i<distances.length;i++){
			distances[i] = allBases.get(i+1).getCoordinate()-allBases.get(i).getCoordinate();
			coords[i]=allBases.get(i).getCoordinate();
		}
		coords[coords.length-1] = allBases.get(coords.length-1).getCoordinate();
		TreeSet<Integer> breakPoints = new TreeSet<Integer>();
		split(distances, coords, breakPoints, windowSize);

		breakPoints.add(coords[coords.length-1]);
		
		int start = r.getStart();
		for (int p:breakPoints){
			windows.add(new Region(r.getGenome(), r.getChrom(), 
                                   Math.max(r.getStart(),start-overlapSize), 
                                   Math.min(r.getEnd(), p+overlapSize)));
			start = p;
		}
		return windows;
	}
	
	private void split(int[]distances, int coords[], TreeSet<Integer> breakPoints, int windowSize){
		int[] idx = StatUtil.findSort(distances.clone());
		int breakIdx = idx[idx.length-1];
		int breakCoor = (coords[breakIdx]+coords[breakIdx+1])/2;
		breakPoints.add(breakCoor);
		// if width of region is still large, recursively split
		if (coords[breakIdx]-coords[0]>windowSize){
			int[] left_distances = new int[breakIdx];
			int[] left_coords = new int[breakIdx+1];
			System.arraycopy(distances, 0, left_distances, 0, breakIdx);
			System.arraycopy(coords, 0, left_coords, 0, breakIdx+1);
			split(left_distances, left_coords, breakPoints, windowSize);
		}
		// skip distances[breakIdx], b/c it is selected
		if (coords[coords.length-1]-coords[breakIdx+1]>windowSize){
			int[] right_distances = new int[distances.length-breakIdx-1];
			int[] right_coords = new int[coords.length-breakIdx-1];
			System.arraycopy(distances, breakIdx+1, right_distances, 0, distances.length-breakIdx-1);
			System.arraycopy(coords, breakIdx+1, right_coords, 0, coords.length-breakIdx-1);
			split(right_distances, right_coords, breakPoints, windowSize);
		}
	}

	private ArrayList<List<StrandedBase>> loadData_checkEnrichment(Region w){
		ArrayList<List<StrandedBase>> signals = loadBasesInWindow(w, "IP");
		if (signals==null || signals.isEmpty())
			return null;

		//initialize count arrays
		totalSigCount=0;
		for (int c=0; c<numConditions; c++){
			sigHitCounts[c]=StrandedBase.countBaseHits(signals.get(c));
			totalSigCount+=sigHitCounts[c];
		}
		// if less than significant number of reads
		if (totalSigCount<config.sparseness)
			return null;
		
		// if IP/Control enrichment ratios are lower than cutoff for all 500bp sliding windows in all conditions, skip
		if (controlDataExist && config.exclude_unenriched){
			boolean enriched = false;
			for (int s=w.getStart(); s<w.getEnd();s+=modelWidth/2){
				int startPos = s;
				int endPos = s+modelWidth;
				if (endPos>w.getEnd()){
					endPos = w.getEnd();
					startPos = endPos - modelWidth;
				}
				Region sw = new Region(gen, w.getChrom(), startPos, endPos);		//sliding window
				for (int c=0; c<numConditions; c++){
					float ip = countIpReads(sw, c);
					float ctrl = countCtrlReads(sw, c);
					if (ctrl==0)
						enriched = true;
					else
						if (ip/ctrl/ratio_non_specific_total[c] >= config.fold)
							enriched = true;
					if (enriched)
						break;
				}
				if (enriched)
					break;
			}
			if (!enriched){
				return null;
			}
		}
		return signals;
	}

	/*
	 * return the bottom elements of x 
	 * so that their sum is less than maxSum
	 * alpha is passed in as maxSum, this essentially ensure that 
	 * we do not eliminate too much, but as quick as possible to reduce run time
	 */
	private Pair<Double, TreeSet<Integer>> findSmallestCases(double[] x, double maxSum){
		TreeSet<Integer> cases = new TreeSet<Integer>();
		double maxSmallest = 0; 
		double sum = 0;
		int[] oldIdx = StatUtil.findSort(x);
		for (int i=0;i<x.length;i++){
			sum += x[i];
			if (sum < maxSum){
				cases.add(oldIdx[i]);
				maxSmallest = x[i];
			}
			else
				break;
		}
		if (cases.isEmpty())
			maxSmallest = x[0];
		Pair<Double, TreeSet<Integer>> result = new Pair<Double, TreeSet<Integer>>(maxSmallest, cases);
		return result;
	}
	/*
	 * return the bottom number of x that are less than max
	 */
	private Pair<Double, TreeSet<Integer>> findSmallestCases_old(double[] x, double number, double max){
		TreeSet<Integer> cases = new TreeSet<Integer>();
		double maxSmallest = 0; 
		number = Math.min(number, x.length);
		int[] oldIdx = StatUtil.findSort(x);
		for (int i=0;i<number;i++){
			if (x[i] < max){
				cases.add(oldIdx[i]);
				maxSmallest = x[i];
			}
			else
				break;
		}
		if (cases.isEmpty())
			maxSmallest = x[0];
		Pair<Double, TreeSet<Integer>> result = new Pair<Double, TreeSet<Integer>>(maxSmallest, cases);
		return result;
	}
	
	private List<BindingComponent> determineNonZeroComps(Region w, MultiIndependentMixtureCounts mim, int[] compPos, double alpha) {
		List<BindingComponent> nonZeroComponents = new ArrayList<BindingComponent>();
		double[]   glob_pi = mim.get_glob_prior_weight();
		double[][] pi      = mim.get_cond_prior_weight();
		int M              = compPos.length;

		for(int j = 0; j < M; j++) {
			if(glob_pi[j] > 0.0) {
				BindingComponent comp = new BindingComponent(model, new Point(w.getGenome(), w.getChrom(), w.getStart() + compPos[j]), numConditions);
				comp.setMixProb(glob_pi[j]);
				comp.setAlpha(alpha);
				comp.setOld_index(j);
				
				for(int t = 0; t < numConditions; t++) {
					comp.setConditionBeta(t, pi[t][j]);
					comp.setCondSumResponsibility(t, mim.get_sum_resp()[t][j]);
				}
					
				// If there is only one condition, pi[0][j], that is, 
				// the condition-specific pi will be set by default to 0.
				// So, we set it to 1.
				if(numConditions == 1)
					comp.setConditionBeta(0, 1.0);
					
				nonZeroComponents.add(comp);
			}
		}
		
		return nonZeroComponents;
	}// end of determineNonZeroComps method

	private boolean checkCompExactMatching(List<BindingComponent> comp_group1, List<BindingComponent> comp_group2) {
		if(comp_group1.size() != comp_group2.size())
			return false;
		else if(comp_group1.size() == 0)
			return true;

		Set<String> loc_comp_group1 = new LinkedHashSet<String>();
		for(BindingComponent b:comp_group1)
			loc_comp_group1.add(b.getLocation().toString());

		for(BindingComponent b:comp_group2)
			loc_comp_group1.remove(b.getLocation().toString());

		if(loc_comp_group1.size() == 0)
			return true;
		else
			return false;
	}//end of checkCompExactMatching method

	private void doBaseFiltering(){
		//  set max read count for bases to filter PCR artifacts (optional)
		if (config.max_hit_per_bp==0){
			// Mandatory base reset for extremely high read counts
			for(int i=0; i<numConditions; i++){
				Pair<ReadCache,ReadCache> e = caches.get(i);
				ArrayList<Pair<Point, Float>> f = e.car().resetHugeBases(config.base_reset_threshold);
				for (Pair<Point, Float> p:f)
					System.err.printf("Warning: %s IP read counts %.0f, reset to 1\n",p.car().toString(), p.cdr());
				if (controlDataExist){
					f = e.cdr().resetHugeBases(config.base_reset_threshold);
					for (Pair<Point, Float> p:f)
						System.err.printf("Warning: %s Ctrl read counts %.0f, reset to 1\n",p.car().toString(), p.cdr());
				}
			}
			return;
		}
		
        for(int i=0; i<numConditions; i++){
			Pair<ReadCache,ReadCache> e = caches.get(i);
			double ipCount = e.car().getHitCount();
			// estimate max hit count per BP
			// if user supply using config.max_hit_per_bp, use it
			// if not supplied, take the max of default (3) and possion expected count
			if (config.max_hit_per_bp==-1){
				int maxPerBP = calcHitCount_per_BP(ipCount, 1e-9);
				max_HitCount_per_base = Math.max(max_HitCount_per_base, maxPerBP);
			}
			else
				max_HitCount_per_base = config.max_hit_per_bp;
        }

        log(1, "\nFiltering duplicate reads, max read count on each base = "+max_HitCount_per_base+"\n");
        
		// Re-scale counts by multiplying each read (in each condition) with the corresponding ratio
		for(int t = 0; t < numConditions; t++) {
			ReadCache ipCache = caches.get(t).car();
			ipCache.filterAllBases(max_HitCount_per_base);
			if(controlDataExist) {
				ReadCache ctrlCache = caches.get(t).cdr();
				ctrlCache.filterAllBases(max_HitCount_per_base);
			}
		}
	}
	
	// determine the max hit count by threshold of Poisson p-value
	// set average read as lambda parameter for Poisson distribution
	private int calcHitCount_per_BP(double totalReads, double threshold){
		int countThres=0;
		DRand re = new DRand();
		Poisson P = new Poisson(0, re);
		double lambda = totalReads/config.mappable_genome_length;
		P.setMean(lambda);
		double l=1;
		for(int b=1; l>threshold; b++){
			l=1-P.cdf(b);	//p-value as the tail of Poisson
			countThres=b;
		}
		return Math.max(1,countThres);
	}

	// Load the reads from all conditions in the region
	private ArrayList<List<StrandedBase>> loadBasesInWindow(Region w, String channel){
		ArrayList<List<StrandedBase>> signals = new ArrayList<List<StrandedBase>>();
		if(channel.equalsIgnoreCase("IP")) {
			//Load each condition's read hits
			for(Pair<ReadCache,ReadCache> e : caches){
				List<StrandedBase> bases_p= e.car().getStrandedBases(w, '+');  // reads of the current region - IP channel
				List<StrandedBase> bases_m= e.car().getStrandedBases(w, '-');  // reads of the current region - IP channel
				bases_p.addAll(bases_m);
				signals.add(bases_p);
			} // for loop
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

	// dynamically determine an alpha value for this sliding window based on IP reads
	// find a 500bp(modelWidth) region with maximum read counts, set alpha=sqrt(maxCount)/3

	private double estimateAlpha(Region window, ArrayList<List<StrandedBase>> signals){
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
		double alphaEstimate = Math.sqrt(maxCount)/config.alpha_factor;
		return alphaEstimate;
	}

	/**
	 * Calls EM training on the Mixture Model.
	 * This is called by SimulatedReadAnalysis for batch analysis in a region
	 */
	public ArrayList<ComponentFeature> simpleExecute(List<StrandedBase> bases, Region r) {
		totalSigCount=StrandedBase.countBaseHits(bases);
        ArrayList<ComponentFeature> out = new ArrayList<ComponentFeature>();
        GPSThread t = new GPSThread(null, null, out, this, constants,config);
        t.simpleRun(bases,r);
        return out;
	}// end of simpleExecute method

	/*
	 * TODO: purpose of this method is to set significance for each condition
	 */
	private List<Feature> condPostFiltering(List<Feature> signifFeats, int cond){
		if (numConditions==1)
			return signifFeats;
		
		List<Feature> condFeats = new ArrayList<Feature>();
		for (Feature sf:signifFeats) {
			boolean significant = false;
			// If TF binding keep the condition-wise component that exceeds the threshold
			if (config.TF_binding) {
				if(((ComponentFeature)sf).getQValueLog10(cond)>config.q_value_threshold &&
                   ((ComponentFeature)sf).getEventReadCounts(cond)>config.sparseness)
					significant = true;
			}

			// if not TF binding, it is Chromatin data, keep all components
			else { significant = true; }

            //			if (significant && ((ComponentFeature)sf).getCondBetas()[cond] > 0.0 ) {
			if (significant)  {
				((ComponentFeature)sf).setCondSignificance(cond, significant);
				condFeats.add((ComponentFeature)sf);
			}
		}

		return condFeats;
	}//end of condPostFiltering method

	/**
	 * It performs statistical tests for determining significant peaks         <br>
	 * If a control is used, the current method is used. Otherwise, MACS proposed method is used.
	 * @param compFeatures
	 */
	private void postEMProcessing(List<ComponentFeature> compFeatures) {
		// use the refined regions to count non-specific reads
		countNonSpecificReads(compFeatures);
		
		// collect enriched regions to exclude to define non-specific region
		Collections.sort(compFeatures, new Comparator<ComponentFeature>() {
                public int compare(ComponentFeature o1, ComponentFeature o2) {
                    return o1.compareByTotalResponsibility(o2);
                }
            });
		// get median event strength
		double medianStrength = compFeatures.get(compFeatures.size()/2).getTotalSumResponsibility();
		System.out.println("Median event strength = "+medianStrength);
		
		// Determine how many of the enriched regions will be included during the regression evaluation
		int numIncludedCandRegs = (int)Math.round(config.pcr*compFeatures.size());	
		ArrayList<Region> exRegions = new ArrayList<Region>();
		for(int i = 0; i < compFeatures.size()-numIncludedCandRegs; i++) {
			exRegions.add(compFeatures.get(i).getPosition().expand(modelRange));
		}
		exRegions.addAll(excludedRegions);	// also excluding the excluded regions (user specified + un-enriched)
		
		calcIpCtrlRatio(exRegions);
		if(controlDataExist) {
			for(int c = 0; c < numConditions; c++)
				System.out.println(String.format("\nScaling condition %s, IP/Control = %.2f", conditionNames.get(c), ratio_non_specific_total[c]));
			System.out.println();
		}
		
		// calculate p-values with or without control
		evaluateSignificance(compFeatures);

		// sort features for final output, by location
		Collections.sort(compFeatures);
		allFeatures = compFeatures;

		insignificantFeatures = new ArrayList<Feature>();
		filteredFeatures = new ArrayList<Feature>();
		for (ComponentFeature cf:compFeatures){
			boolean significant = false;
			// for multi-condition, at least be significant in one condition
			// The read count test is for each condition
		
			for (int cond=0; cond<numConditions; cond++){
				if(cf.getEventReadCounts(cond)>config.sparseness){	// first pass read count test
					if (config.TF_binding){	// single event IP/ctrf only applies to TF
						if (cf.getQValueLog10(cond)>config.q_value_threshold){
							significant = true;
							break;
						}
					}
					else		//TODO: if histone data, need to test as a set of events
						significant = true;
				}
			}

			boolean notFiltered = false;
			// only filter high read count events (more likely to be artifacts) 
			if (config.filterEvents && cf.getTotalSumResponsibility()>medianStrength){
				for (int cond=0; cond<numConditions; cond++){
					// if one condition is good event, this position is GOOD
					// logKL of event <= 2.5, and IP/control >= 4 --> good (true)
					if (cf.getShapeDeviation(cond)<=config.shapeDeviation){
						if (!controlDataExist){
							notFiltered = true;
							break;
						}
						else{
							double ratio = cf.getEventReadCounts(cond)/cf.getScaledControlCounts(cond);
//							if ((ratio>=config.fold && cf.getAverageIpCtrlLogKL()>config.kl_ic) || (cf.getAverageIpCtrlLogKL()<config.kl_ic && ratio>=config.fold*2)){
								if (ratio>=config.fold && cf.getAverageIpCtrlLogKL()>config.kl_ic){
								notFiltered = true;
								break;
							}
						}						
					}
				}
			}
			else
				notFiltered = true;

			if (significant){
				if (notFiltered)
					signalFeatures.add(cf);
				else
					filteredFeatures.add(cf);
			}
			else
				insignificantFeatures.add(cf);
		}
		
		// set joint event flag in the events
		if (signalFeatures.size()>=2){
			for (int i=0;i<signalFeatures.size();i++){
				ComponentFeature cf = (ComponentFeature)signalFeatures.get(i);
				cf.setJointEvent(false);		// clean the mark first
			}
			for (int i=0;i<signalFeatures.size()-1;i++){
				ComponentFeature cf = (ComponentFeature)signalFeatures.get(i);
				ComponentFeature cf2 = (ComponentFeature)signalFeatures.get(i+1);
				if (cf.onSameChrom(cf2)){
					if (cf.getPosition().distance(cf2.getPosition())<=config.joint_event_distance){
						cf.setJointEvent(true);
						cf2.setJointEvent(true);
					}
				}
			}
		}

		log(1, "Events discovered \nSignificant:\t"+signalFeatures.size()+
            "\nInsignificant:\t"+insignificantFeatures.size()+
            "\nFiltered:\t"+filteredFeatures.size()+"\n");
	}//end of post EM Processing


	/* 
	 * Calc the ratio of IP vs Control channel, exlcluding the specific regions
	 */
	private void calcIpCtrlRatio(ArrayList<Region> specificRegions) {
		// linear regression to get the IP/control ratio
		// now we do not require whole genome data, because partial data could be run on 1 chrom, still enough data to esitmates
		if(controlDataExist) {
			ratio_non_specific_total = new double[numConditions];
			for(int t = 0; t < numConditions; t++){
				ratio_non_specific_total[t] = getSlope(t, t, "IP/CTRL", specificRegions);
			}
			ComponentFeature.setNon_specific_ratio(ratio_non_specific_total);
		}
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
	public void setRegions(String fname, boolean toExpandRegion) {
		ArrayList<Region> rset = CommonUtils.loadRegionFile(fname, gen);
		restrictRegions = mergeRegions(rset, toExpandRegion);
	}

	// merge the overlapped regions
	// if "toExpandRegion"=true, expand each region on both side to leave enough space,
	// to include every potential reads, then merge
	private ArrayList<Region> mergeRegions(ArrayList<Region> regions, 
                                             boolean toExpandRegion){
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

	/**
	 * Select the regions where the method will run on.    <br>
	 * It will be either whole genome, chromosome(s) or extended region(s).
	 * @param focusRegion
	 * @return
	 */
	private ArrayList<Region> selectEnrichedRegions(List<Region> focusRegions, boolean useSubtractedData){
		long tic = System.currentTimeMillis();
		ArrayList<Region> regions = new ArrayList<Region>();
		Map<String, List<Region>> chr_focusReg_pair = new HashMap<String, List<Region>>();
		
		if(focusRegions.size() > 0) {
			for(Region focusRegion:focusRegions) {
				if(!chr_focusReg_pair.containsKey(focusRegion.getChrom()))
					chr_focusReg_pair.put(focusRegion.getChrom(), new ArrayList<Region>());
				chr_focusReg_pair.get(focusRegion.getChrom()).add(focusRegion);
			}			
		}
		else {
			for (String chrom:gen.getChromList()){
				Region wholeChrom = new Region(gen, chrom, 0, gen.getChromLength(chrom)-1);
				chr_focusReg_pair.put(chrom, new ArrayList<Region>());
				chr_focusReg_pair.get(chrom).add(wholeChrom);
			}
		}
		
		for(String chrom:chr_focusReg_pair.keySet()) {
			for(Region focusRegion:chr_focusReg_pair.get(chrom)) {
				List<Region> rs = new ArrayList<Region>();
				List<StrandedBase> allBases = new ArrayList<StrandedBase>();
				for (int c=0; c<caches.size(); c++){
					ReadCache ip = caches.get(c).car();
					List<StrandedBase> bases = null;
					if (this.controlDataExist && useSubtractedData)
						bases = ip.getSubtractedBases(focusRegion, caches.get(c).cdr(), ratio_non_specific_total[c]);
					else
						bases = ip.getUnstrandedBases(focusRegion);
					if (bases==null || bases.size()==0){
						continue;
					}
					allBases.addAll(bases); // pool all conditions
				}
				Collections.sort(allBases);
				int start=0;
				for (int i=2;i<allBases.size();i++){
					int distance = allBases.get(i).getCoordinate()-allBases.get(i-2).getCoordinate();
					if (distance > modelWidth){ // a large enough gap (with 1 base in middle) to cut
						// only select region with read count larger than minimum count
						int breakPoint = -1;
						if (allBases.get(i-1).getStrand()=='+')
							breakPoint = i-2;
						else
							breakPoint = i-1;
						float count = 0;
						for(int m=start;m<=breakPoint;m++){
							count += allBases.get(m).getCount();
						}
						if (count >= config.sparseness){
							Region r = new Region(gen, chrom, allBases.get(start).getCoordinate(), allBases.get(breakPoint).getCoordinate());
							// if the average read count per modelWidth is less than config.sparseness/2, find sparse point to further split
							rs.add(r);
						}
						start = breakPoint+1;
					}
				}
				// check last continuous region
				float count = 0;
				for(int m=start;m<allBases.size();m++){
					count += allBases.get(m).getCount();
				}
				if (count>=config.sparseness){
					Region r = new Region(gen, chrom, allBases.get(start).getCoordinate(), allBases.get(allBases.size()-1).getCoordinate());
					rs.add(r);
				}

				// check regions, exclude un-enriched regions based on control counts
				ArrayList<Region> toRemove = new ArrayList<Region>();
				for (Region r: rs){
					if (r.getWidth()<=config.min_region_width){
						toRemove.add(r);
						continue;
					}
					// for regions <= modelWidth, most likely single event, can be compared to control
					if (this.controlDataExist)
						if (r.getWidth()<=modelWidth){
							boolean enriched = false;
							for (int c=0;c<numConditions;c++){
								if (countIpReads(r,c)/countCtrlReads(r,c)/this.ratio_total[c]>config.fold){
									enriched = true;
									break;
								}
							}
							if (!enriched){	// remove this region if it is not enriched in any condition
								toRemove.add(r);
								continue;
							}
						}
				}
				rs.removeAll(toRemove);
				if (!rs.isEmpty())
					regions.addAll(rs);
			}
		}
		
		log(3, "selectEnrichedRegions(): " + CommonUtils.timeElapsed(tic));
		return regions;
	}//end of selectEnrichedRegions method

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
		for(int t=0; t<constants.MAX_EM_ML_ITER ; t++){
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

			//Semi-E-step:Calculate new un-normalized responsibilities
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
					sum_j += r[i][j];
				}

				if (sum_j!=0){
					LL+=Math.log(sum_j)*counts[i];
				}
			}

            //			log(4, "\nEM_ML: "+t+"\t"+LL+"\t"+lastLL+".");
			if (Math.abs(LL-lastLL)>constants.EM_ML_CONVERGENCE ){
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
		if (constants.MAKE_HARD_ASSIGNMENT){
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
	 * @param spacing spacing between the components
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

	private int[] setComponentPositions(int width, int spacing) {
		int[] compPos;
		List<Integer> compPosList = new ArrayList<Integer>();
		for(int j = 0; j < width-1; j += spacing)
			compPosList.add(j);

		compPosList.add(width-1);
		compPos = Utils.ref2prim(compPosList.toArray(new Integer[0]));
		return compPos;
	}//end of setComponentPositions method

	private ArrayList<ComponentFeature> callFeatures(ArrayList<BindingComponent> comps){
        //		if (config.post_artifact_filter){
        //			comps = postArtifactFilter(comps);
        //		}

		ArrayList<ComponentFeature> features = new ArrayList<ComponentFeature>();

		if (comps.size()==0)
			return features;

		//Make feature calls (all non-zero components)
		// remove events with IP < alpha (this happens because EM exit before convergence)
		for(BindingComponent b : comps){
			if(b.getMixProb()>0 && b.getTotalSumResponsibility()>=config.sparseness){
				features.add(callFeature(b));
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

				if (config.post_artifact_filter){
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

					}
				}
				bases.addAll(bases_m);

				double[][] assignment = assignBasesToEvents(bases, comps);

				for(int j=0;j<features.size();j++){
					ComponentFeature comp = features.get(j);
					int pos = comp.getPosition().getLocation();
					// set control reads count and logKL
					double total=0;
					double[] ctrl_profile_plus = new double[modelWidth];
					double[] ctrl_profile_minus = new double[modelWidth];
					for(int i=0;i<bases.size();i++){
						StrandedBase base = bases.get(i);
						if (assignment[i][j]>0){
							if (base.getStrand()=='+'){
								if (base.getCoordinate()-pos-model.getMin()>modelWidth)
									System.err.println(pos+"\t+\t"+base.getCoordinate());
								ctrl_profile_plus[base.getCoordinate()-pos-model.getMin()]=assignment[i][j]*base.getCount();
							}
							else{
								if (pos-base.getCoordinate()-model.getMin()>modelWidth)
									System.err.println(pos+"\t-\t"+base.getCoordinate());
								ctrl_profile_minus[pos-base.getCoordinate()-model.getMin()]=assignment[i][j]*base.getCount();
							}
							total += assignment[i][j]*base.getCount();
						}
					}
					comp.setControlReadCounts(total, c);
					BindingComponent bb = null;
					for (BindingComponent b:comps){
						if (b.getLocation().getLocation()==pos)
							bb = b;
					}
					double logKL_plus=99;		// default to a large value if no control data
					double logKL_minus=99; 
					if (total!=0){
						logKL_plus  = StatUtil.log_KL_Divergence(StatUtil.symmetricKernelSmoother(bb.getReadProfile_plus(c), gaussian), StatUtil.symmetricKernelSmoother(ctrl_profile_plus, gaussian));
						logKL_minus = StatUtil.log_KL_Divergence(StatUtil.symmetricKernelSmoother(bb.getReadProfile_minus(c), gaussian), StatUtil.symmetricKernelSmoother(ctrl_profile_minus, gaussian));
					}
					comp.setIpCtrlLogKL(c, logKL_plus, logKL_minus);
				}
			}
		}

		return(features);
	}

	// make componentFeature based on the bindingComponent
	private ComponentFeature callFeature(BindingComponent b){
		ComponentFeature cf = new ComponentFeature(b);

		// KL
		double[] logKL_plus = new double[numConditions];
		double[] logKL_minus = new double[numConditions];
		double[] shapeDeviation = new double[numConditions];
		for (int c=0;c<numConditions;c++){
			double[] profile_plus = b.getReadProfile_plus(c);
			double[] profile_minus = b.getReadProfile_minus(c);
			if (config.kl_count_adjusted){
				logKL_plus[c]  = StatUtil.log_KL_Divergence(model.getProbabilities(), StatUtil.symmetricKernelSmoother(b.getReadProfile_plus(c), gaussian));
				logKL_minus[c] = StatUtil.log_KL_Divergence(model.getProbabilities(), StatUtil.symmetricKernelSmoother(b.getReadProfile_minus(c), gaussian));
			}
			else{
				double count = 0.0;
				for(int n = 0; n < profile_plus.length; n++) { 
					count += profile_plus[n]+profile_minus[n]; 
				}
				double width;
				if (count>6){
					width = 50/Math.log(count);	// just some empirical formula, no special theory here
				}
				else{
					width=28;
				}
				logKL_plus[c]  = logKL_profile(profile_plus, width);
				logKL_minus[c] = logKL_profile(profile_minus, width);
			}

			if (config.KL_smooth_width!=0)
				shapeDeviation[c]  = calc2StrandSmoothKL(profile_plus,profile_minus);
			else
				shapeDeviation[c]  = calc2StrandNzKL(profile_plus,profile_minus);

            //			System.err.println(String.format("%.2f\t%.2f\t%.2f\t%s", 
            //					shapeDeviation[c],
            //					calcAbsDiff(profile_plus), 
            //					calcAbsDiff(profile_minus),
            //					b.getLocation().toString()));
			
			// sum the read profiles (optional)
			// at this point, we do not have p-value, etc, this is our best guess
			if (config.use_joint_event && 
                shapeDeviation[c]<config.shapeDeviation &&
                (logKL_plus[c]<0 || logKL_minus[c]<0)){
				for (int i=0;i<profile_plus.length;i++){
					profile_plus_sum[i] += profile_plus[i];
					profile_minus_sum[i] += profile_minus[i];
				}
			}
		}
		cf.setProfileLogKL(logKL_plus, logKL_minus);
		cf.setShapeDeviation(shapeDeviation);
		return cf;
	}
	/*
	 * Calculate logKL for read profile around an event
	 * Gaussian Width is given
	 */
	private double logKL_profile(double[]profile, double width){
		double[] gaus = new double[modelWidth];
		NormalDistribution gaussianDist = new NormalDistribution(0, width*width);
		for (int i=0;i<gaus.length;i++)
			gaus[i]=gaussianDist.calcProbability((double)i);
		return StatUtil.log_KL_Divergence(model.getProbabilities(), StatUtil.symmetricKernelSmoother(profile, gaus));
	}

	private double calcAbsDiff(double[]profile){
		double absDiff = 0;
		double[] m = model.getProbabilities();
		double totalCount = 0;
		double nzProbSum = 0;
		for (int i=0;i<profile.length;i++){
			if (profile[i]!=0){
				totalCount+=profile[i];
				nzProbSum += m[i];
			}
		}
		if (totalCount==0)
			return config.shapeDeviation*2*0.75;
		for (int i=0;i<profile.length;i++){
			if (profile[i]!=0){
				double expected = totalCount*m[i]/nzProbSum;
				absDiff += Math.abs(profile[i]-expected);
			}
		}
		return absDiff/totalCount;
	}

	/*
	 *   logKL of non-zero discrete profile, concat 2 strands
	 *   no gaussian density, use only non-zero read counts
	 */
	private double calc2StrandNzKL(double[]profile_p, double[]profile_m){
		double asym_p = profile_p.length>2?0:0.2;	// penalty of having too few reads on one strand
		double asym_m = profile_m.length>2?0:0.2;
		
		double[] profile = new double[profile_p.length+profile_m.length];
		System.arraycopy(profile_p, 0, profile, 0, profile_p.length); 
		System.arraycopy(profile_m, 0, profile, profile_p.length, profile_m.length); 
		
		double[] m = model.getProbabilities();
		double[] m2 = new double[profile.length];
		System.arraycopy(m, 0, m2, 0, m.length); 
		System.arraycopy(m, 0, m2, m.length, m.length); 

		ArrayList<Integer> nzPos = new ArrayList<Integer>();
		for (int i=0;i<profile.length;i++){
			if (profile[i]!=0){
				nzPos.add(i);
			}
		}
		double[] m_nz = new double[nzPos.size()];
		double[] p_nz = new double[nzPos.size()];
		for (int i=0;i<nzPos.size();i++){
			m_nz[i] = m2[nzPos.get(i)];
			p_nz[i] = profile[nzPos.get(i)];
		}
		
		double log10KL = StatUtil.log10_KL_Divergence(m_nz, p_nz)+asym_p+asym_m;
		return Math.max(log10KL, -4);
	}	

	/*
	 *   logKL of smooth gaussian profile, concat 2 strands
	 */
	private double calc2StrandSmoothKL(double[]profile_p, double[]profile_m){
		double asym_p = profile_p.length>2?0:0.2;	// penalty of having too few reads on one strand
		double asym_m = profile_m.length>2?0:0.2;
		
		profile_p = StatUtil.gaussianSmoother(profile_p, config.KL_smooth_width);
		profile_m = StatUtil.gaussianSmoother(profile_m, config.KL_smooth_width);
		
		double[] profile = new double[profile_p.length+profile_m.length];
		System.arraycopy(profile_p, 0, profile, 0, profile_p.length); 
		System.arraycopy(profile_m, 0, profile, profile_p.length, profile_m.length); 
		
		double[] m = model.getProbabilities();
		double[] m2 = new double[profile.length];
		System.arraycopy(m, 0, m2, 0, m.length); 
		System.arraycopy(m, 0, m2, m.length, m.length); 
		
		double log10KL = StatUtil.log10_KL_Divergence(m2, profile);
		return Math.max(log10KL, -4);
	}	


	/**
	 * Finds the exact peak in the region by scanning for maximum-likelihood using a <tt>BindingModel</tt>
	 * assuming that the region contains only a single event
	 * @param signals <tt>ReadHits</tt> of the IP channel corresponding to region <tt>region</tt>
	 * @param region Regions to be scanned
	 * @return the peak as a <tt>BindingComponent</tt>
	 */
	private BindingComponent scanPeak(Point p){
		Region scanRegion = p.expand(config.SCAN_RANGE);
		Region plusRegion = scanRegion.expand(-model.getMin(), model.getMax());
		Region minusRegion = scanRegion.expand(model.getMax(), -model.getMin());
		ArrayList<List<StrandedBase>> signals = new ArrayList<List<StrandedBase>>();
		//Load each condition's read hits
		int totalHitCounts = 0;
		for(int c=0;c<numConditions;c++){
			Pair<ReadCache,ReadCache> e = caches.get(c);
			List<StrandedBase> bases_p= e.car().getStrandedBases(plusRegion, '+');  // reads of the current region - IP channel
			List<StrandedBase> bases_m= e.car().getStrandedBases(minusRegion, '-');  // reads of the current region - IP channel

			bases_p.addAll(bases_m);
			totalHitCounts += StrandedBase.countBaseHits(bases_p);
			signals.add(bases_p);
		} // for loop

		if (signals==null||totalHitCounts==0) // check for empty read region
			return null;

		return scanPeak(signals, scanRegion);
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

		double[] likelihoods = scoreLikelihoods(newModel, allBases, scanRegion, 1);
		BindingComponent maxComp=null;
		int maxIdx = 0;
		double maxScore=Double.NEGATIVE_INFINITY;
		for (int i=0;i<likelihoods.length;i++){
			if(likelihoods[i]!=0 && likelihoods[i]>maxScore ){
				maxIdx = i;
				maxScore = likelihoods[i];
			}
		}
		// make a new binding component with original binding model
		maxComp = new BindingComponent(model, new Point(gen, scanRegion.getChrom(), scanRegion.getStart()+maxIdx), numConditions);
		// set pi, beta, responsibility, note it is a single event
		maxComp.setMixProb(1);
	
		// Since we expanded the binding model, all the reads in a condition are assigned to the unique component
		float[] hitCounts = new float[numConditions];
		int totalHitCounts = 0;
		for(int c = 0; c < numConditions; c++) {
			List<StrandedBase> bases= signals.get(c);
			hitCounts[c]    = StrandedBase.countBaseHits(bases);
			totalHitCounts += hitCounts[c];
			maxComp.setCondSumResponsibility(c, hitCounts[c]);
		}
		
		if(config.use_betaEM) {
			for (int c=0;c<numConditions;c++)
				maxComp.setConditionBeta(c, hitCounts[c]/totalHitCounts);  //beta is being evaluated for one point across all conditions
		}
		else {
			for(int c=0; c<numConditions; c++)
				if(maxComp.getSumResponsibility(c) > 0.0)
					maxComp.setConditionBeta(c, 1.0);   //beta is being evaluated for each condition separately
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

	private double[]scoreLikelihoods(BindingModel distr, ArrayList<StrandedBase> bases, Region scanRegion, int step){
		double[] likelihoods = new double[scanRegion.getWidth()];
		
		Region dataRegion = scanRegion.expand(modelRange, modelRange);
		ArrayList<StrandedBase> data = new ArrayList<StrandedBase>();
		for (StrandedBase b: bases){
			if (b.getCoordinate()>=dataRegion.getStart() && b.getCoordinate()<=dataRegion.getEnd())
				data.add(b);
		}
		// score LL for all the data when we have an event on each position in the region
		for(int i=scanRegion.getStart(); i<=scanRegion.getEnd(); i+=step){
			double logLikelihood=0;
			for(StrandedBase b : data){
				int dist = b.getStrand()=='+' ? b.getCoordinate()-i:i-b.getCoordinate();
				logLikelihood += Math.log(distr.probability(dist))*b.getCount();
			}
			likelihoods[i-scanRegion.getStart()] = logLikelihood;
		}
		return likelihoods;
	}
	
	/**
	 * This method normalizes conditions as follows:							<br>
	 * It evaluates the slope between each pair of an IP condition and a reference one.
	 * Here, chosen to be condition 1. It does the same for the Ctrl channel.    <br>
	 * Then, it re-scales all reads for each condition by multipying with the corresponding slopes.  <br>
	 * Lastly, it evaluates the ratio between the original and the 'inflated' read counts
	 * and divides every read weight with this amount to ensure that the total read counts remain constant.
	 * @param caches
	 */
	private void normExpts(ArrayList<Pair<ReadCache, ReadCache>> caches) {
		int C = caches.size();   // # conds
		if(C <= 1) return;
		
		System.out.println("\nNormalizing read counts across multiple conditions ...");
		
		// Get the total read counts in IP and Control
		double totIPCounts   = 0.0;
		double totCtrlCounts = 0.0;
		for(int t = 0; t < C; t++) {
			ReadCache ipCache = caches.get(t).car();
			totIPCounts += ipCache.getHitCount();
			if(controlDataExist) {
				ReadCache ctrlCache = caches.get(t).cdr();
				totCtrlCounts += ctrlCache.getHitCount();
			}
		}
		
		// Get the ratios between IP pairs and Ctrl pairs separately based on linear regression
		ArrayList<Region> emptyReg = new ArrayList<Region>();
		double[] ip_ratio   = new double[C];;
		double[] ctrl_ratio = new double[C];;
		ip_ratio[0] = 1.0;
		for(int t = 1; t < C; t++)
			ip_ratio[t] = getSlope(0, t, "IP", emptyReg);
		if(controlDataExist) {
			ctrl_ratio = new double[C];
			ctrl_ratio[0] = 1.0;
			for(int t = 1; t < C; t++)
				ctrl_ratio[t] = getSlope(0, t, "CTRL", emptyReg);
		}
		
		// Re-scale counts by multiplying each read (in each condition) with the corresponding ratio
		for(int t = 1; t < C; t++) {
			ReadCache ipCache = caches.get(t).car();
			ipCache.normalizeCounts(ip_ratio[t]);
			if(controlDataExist) {
				ReadCache ctrlCache = caches.get(t).cdr();
				ctrlCache.normalizeCounts(ctrl_ratio[t]);
			}
		}

		// The total read counts now in IP and Ctrl are different. Calculate them.
		double new_totIPCounts   = 0.0;
		double new_totCtrlCounts = 0.0;
		for(int t = 0; t < C; t++) {
			ReadCache ipCache = caches.get(t).car();
			new_totIPCounts += ipCache.getHitCount();
			if(controlDataExist) {
				ReadCache ctrlCache = caches.get(t).cdr();
				new_totCtrlCounts += ctrlCache.getHitCount();
			}
		}

		// Find how many times the new read counts are larger than the old ones
		double ip_norm_factor   = 0.0;
		double ctrl_norm_factor = 0.0;
		ip_norm_factor = new_totIPCounts/totIPCounts;
		if(controlDataExist)
			ctrl_norm_factor = new_totCtrlCounts/totCtrlCounts;
		
		// Re-scale counts by dividing each read (in each condition) with the 
		// corresponding ratio so that the total read counts remain constant
		for(int t = 0; t < C; t++) {
			ReadCache ipCache = caches.get(t).car();
			ipCache.normalizeCounts(1.0/ip_norm_factor);
			if(controlDataExist) {
				ReadCache ctrlCache = caches.get(t).cdr();
				ctrlCache.normalizeCounts(1.0/ctrl_norm_factor);
			}
		}
		
		
	}//end of normExpts method
	

	/**
	 * This method will run regression on two dataset (can be IP/IP, or IP/Ctrl, or Ctrl/Ctrl), 
	 * using the non-specific region (all the genome regions, excluding the enriched (specific) regions.
	 * @param condX_idx index of condition (channel) where it will be used a dependent variable 
	 * @param condY_idx index of condition (channel) where it will be used an independent variable
	 * @param flag <tt>IP</tt> for pairs of IP conditions, <tt>CTRL</tt> for pairs of control conditions, <tt>IP/CTRL</tt> for ip/control pairs 
	 * @param specificRegions regions to exclude for data for regression
	 * @return slope, as the ratio of ( 1st dataset / 2nd dataset )
	 */
	private double getSlope(int condX_idx, int condY_idx, String flag, ArrayList<Region> specificRegions) {
		if(condX_idx < 0 || condX_idx >= numConditions)
			throw new IllegalArgumentException(String.format("cond1_idx, cond2_idx have to be numbers between 0 and %d.", numConditions-1));
		if(condY_idx < 0 || condY_idx >= numConditions)
			throw new IllegalArgumentException(String.format("cond1_idx, cond2_idx have to be numbers between 0 and %d.", numConditions-1));
		if(!(flag.equalsIgnoreCase("IP") || flag.equalsIgnoreCase("CTRL") || flag.equalsIgnoreCase("IP/CTRL")))
			throw new IllegalArgumentException("The valid arguments for the flag argument is: IP, CTRL, IP/CTRL.");
		if(flag.equalsIgnoreCase("IP/CTRL") && (condX_idx != condY_idx))
			throw new IllegalArgumentException("When you submit IP/CTRL pairs the condition has to be the same.");
			
		double slope;
		List<PairedCountData> scalePairs = new ArrayList<PairedCountData>();
		int non_specific_reg_len = 10000;	
		Map<String, ArrayList<Region>> chrom_regions_map = new HashMap<String, ArrayList<Region>>();
		// group regions by chrom
		for(Region r:specificRegions) {
			String chrom = r.getChrom();
			if(!chrom_regions_map.containsKey(chrom))
				chrom_regions_map.put(chrom, new ArrayList<Region>());
			chrom_regions_map.get(chrom).add(r);
		}

		for(String chrom:gen.getChromList()) {
			List<Region> chrom_non_specific_regs = new ArrayList<Region>();
			int chromLen = gen.getChromLength(chrom);
			// Get the candidate (enriched) regions of this chrom sorted by location
			List<Region> chr_enriched_regs = new ArrayList<Region>();
			if(chrom_regions_map.containsKey(chrom)) {
				chr_enriched_regs = chrom_regions_map.get(chrom);
				Collections.sort(chr_enriched_regs);
			}
			else{	// We only estimate using the chrom that contains enriched data regions
				continue;
			}
				
			// Construct the non-specific regions and check for overlapping with the candidate (enriched) regions
			int start = 0;
			int prev_reg_idx = 0;
			int curr_reg_idx = 0;
			while(start < chromLen) {
				Region non_specific_reg = new Region(gen, chrom, start, Math.min(start + non_specific_reg_len -1, chromLen-1));
				
				if(chr_enriched_regs.size() > 0) {
					if(!(non_specific_reg.overlaps(chr_enriched_regs.get(prev_reg_idx)) || non_specific_reg.overlaps(chr_enriched_regs.get(curr_reg_idx)))) {
						chrom_non_specific_regs.add(non_specific_reg);
					}
					else {
						while(curr_reg_idx < chr_enriched_regs.size() && non_specific_reg.overlaps(chr_enriched_regs.get(curr_reg_idx)))
							curr_reg_idx++;
						curr_reg_idx = Math.min(chr_enriched_regs.size()-1, curr_reg_idx);
						prev_reg_idx = Math.max(prev_reg_idx, curr_reg_idx-1);
					}
				}
				else {
					chrom_non_specific_regs.add(non_specific_reg);
				}
				
				start += non_specific_reg_len;
			}//end of while LOOP
			
			if(flag.equalsIgnoreCase("IP")) {
				// Estimate the (ipCounts_cond1, ipCounts_cond2) pairs
				for(Region r:chrom_non_specific_regs) {				
					double ipCounts_condX = countIpReads(r, condX_idx);
					double ipCounts_condY = countIpReads(r, condY_idx);
					if (ipCounts_condX!=0 && ipCounts_condY!=0)	// only for non-zero counts, as PeakSeq Paper
						scalePairs.add(new PairedCountData(ipCounts_condY, ipCounts_condX));  // we want to see how many time ipCounts_condY are larger from ipCounts_condX
				}
			}
			else if(flag.equalsIgnoreCase("CTRL")) {
				// Estimate the (ctrlCounts_cond1, ctrlCounts_cond2) pairs
				for(Region r:chrom_non_specific_regs) {				
					double ctrlCounts_condX = countCtrlReads(r, condX_idx);
					double ctrlCounts_condY = countCtrlReads(r, condY_idx);
					if (ctrlCounts_condX!=0 && ctrlCounts_condY!=0)
						scalePairs.add(new PairedCountData(ctrlCounts_condY, ctrlCounts_condX));  // we want to see how many time ctrlCounts_condY are larger from ctrlCounts_condX
				}
			}
			else {
				// Estimate the (ipCount, ctrlCount) pairs
				for(Region r:chrom_non_specific_regs) {
					double ctrlCounts_X = countCtrlReads(r, condX_idx);
					double ipCounts_Y   = countIpReads(r, condY_idx);
					if (ctrlCounts_X!=0 && ipCounts_Y!=0)
						scalePairs.add(new PairedCountData(ctrlCounts_X, ipCounts_Y));  // we want to see how many time ipCounts are larger from ctrlCounts
				}
			}
		}//end of for(String chrom:gen.getChromList()) LOOP
			
		// Calculate the slope for this condition
		slope = calcSlope(scalePairs);
		scalePairs.clear();
		return slope;
	}//end of getSlope method
		
	private static double calcSlope(List<PairedCountData> scalePairs) {
		double slope;
        if(scalePairs==null || scalePairs.size()==0) { return 1.0; }
        DataFrame df = new DataFrame(edu.mit.csail.cgs.deepseq.PairedCountData.class, scalePairs.iterator());
        DataRegression r = new DataRegression(df, "y~x - 1");
        r.calculate();
        Map<String, Double> map = r.collectCoefficients();
        slope = map.get("x");
        return slope;
    }//end of calcSlope method

	
	public double updateBindingModel(int left, int right){
		if (signalFeatures.size()<config.min_event_count){
			System.err.println("Warning: The read distribution is not updated due to too few ("+signalFeatures.size()+"<"+config.min_event_count+") significant events.");
			return -100;
		}
		int width = left+right+1;
		
		ArrayList<ComponentFeature> cfs = new ArrayList<ComponentFeature>();
		for (Feature f: signalFeatures){
			ComponentFeature cf = (ComponentFeature)f;
			// the events that are used to refine read distribution should be
			// having strength and shape at the upper half ranking
			if (config.use_joint_event || !cf.isJointEvent() )
				cfs.add(cf);
		}
		Collections.sort(cfs, new Comparator<ComponentFeature>(){
                public int compare(ComponentFeature o1, ComponentFeature o2) {
                    return o1.compareByTotalResponsibility(o2);
                }
            });
		
		int eventCounter=0;
        // data for all conditions
        double[][] newModel_plus=new double[numConditions][width];
        double[][] newModel_minus=new double[numConditions][width];
		// sum the read profiles from all qualified binding events for updating model
        double strengthThreshold=-1;
		for (int c=0;c<numConditions;c++){
			ReadCache ip = caches.get(c).car();
			eventCounter=0;
			for (ComponentFeature cf: cfs){
				if (cf.getQValueLog10(c)>config.q_value_threshold){
					Region region = cf.getPosition().expand(0).expand(left, right);
					List<StrandedBase> bases = ip.getStrandedBases(region, '+');
					for (StrandedBase base:bases){
						newModel_plus[c][base.getCoordinate()-region.getStart()]+=base.getCount();
					}
					region = cf.getPosition().expand(0).expand(right, left);
					bases = ip.getStrandedBases(region, '-');
					for (StrandedBase base:bases){
						newModel_minus[c][region.getEnd()-base.getCoordinate()]+=base.getCount();
					}
					eventCounter++;
					if (eventCounter>config.top_events-1){	// reach the top counts
						strengthThreshold = cf.getTotalSumResponsibility();
						break;
					}
				}
			}
		}
		if (strengthThreshold==-1)		// if not set, then we are using all events
			strengthThreshold=cfs.get(cfs.size()-1).getTotalSumResponsibility();

		// we have collected data for both strands, all conditions
		// but here we only use single model for all conditions
		double[] model_plus = new double[width];
		double[] model_minus = new double[width];
		if (config.use_joint_event){
			model_plus = profile_plus_sum;
			model_minus = profile_minus_sum;
		}
		else{
			for (int i=0;i<width;i++){
				for (int c=0; c<numConditions; c++){
					model_plus[i]+= newModel_plus[c][i];
					model_minus[i]+= newModel_minus[c][i];
				}
			}
		}
		// smooth the model profile
		if (config.smooth_step>0){
			model_plus = StatUtil.cubicSpline(model_plus, config.smooth_step, config.smooth_step);
			model_minus = StatUtil.cubicSpline(model_minus, config.smooth_step, config.smooth_step);
		}
		StatUtil.mutate_normalize(model_plus);
		StatUtil.mutate_normalize(model_minus);

		//compare models from 2 strands, shift binding position if needed
		int shift = BindingModel.minKL_Shift(model_plus, model_minus);

		//use single model for now
		List<Pair<Integer, Double>> dist = new ArrayList<Pair<Integer, Double>>();
		for (int i=0;i<width;i++){
			double sum = model_plus[i]+model_minus[i];
			sum = sum>=0?sum:2.0E-300;
			Pair<Integer, Double> p = new Pair<Integer, Double>(i-left, sum);
			dist.add(p);
		}

		double[] oldModel = model.getProbabilities();
		String oldName = model.getFileName();
		model = new BindingModel(dist);
		model.setFileName(oldName);
		model.printToFile(outName+"_Read_distribution.txt");
		modelRange = model.getRange();
		modelWidth = model.getWidth();
		allModels.put(outName, model);

		double logKL = 0;
		if (oldModel.length==modelWidth)
			logKL = StatUtil.log_KL_Divergence(oldModel, model.getProbabilities());
		log(1, "Refine read distribution from " + eventCounter +" binding events. \n");

		return logKL;
	}

	public BindingModel getModel(){
		return model;
	}

	public void plotAllReadDistributions(){
		Color[] colors = {Color.black, Color.red, Color.blue, Color.green, Color.cyan, Color.orange};
		String filename = outName.substring(0, outName.length()-2) + "_All_Read_Distributions.png";
		File f = new File(filename);
		int w = 1000;
		int h = 600;
		int margin= 50;
	    BufferedImage im = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
	    Graphics g = im.getGraphics();
	    Graphics2D g2 = (Graphics2D)g;
	    g2.setRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
	    g2.setColor(Color.white);
	    g2.fillRect(0, 0, w, h);	
	    g2.setColor(Color.gray);
	    g2.drawLine(20, h-margin, w-20, h-margin);		// x-axis
	    g2.drawLine(w/2, margin, w/2, h-margin);	// y-axis    
	    g.setFont(new Font("Arial",Font.PLAIN,16));
	    for (int p=-2;p<=2;p++){
	    	g2.drawLine(w/2+p*200, h-margin-10, w/2+p*200, h-margin);	// tick  
	    	g2.drawString(""+p*200, w/2+p*200-5, h-margin+22);			// tick label
	    }
	    
	    double maxProb = 0;	    
	    ArrayList<String> rounds = new ArrayList<String>();
	    rounds.addAll(allModels.keySet());
	    Collections.sort(rounds);
	    for (String key:rounds){
	    	int summit = allModels.get(key).getSummit();
	    	maxProb = Math.max(maxProb, allModels.get(key).probability(summit));
	    }
	    
	    for (int i=0;i<rounds.size();i++){
	    	BindingModel m = allModels.get(rounds.get(i));
	    	List<Pair<Integer, Double>> points = m.getEmpiricalDistribution();
		    g2.setColor(colors[i % colors.length]);
		    g2.setStroke(new BasicStroke(4));
		    for (int p=0;p<points.size()-1;p++){
		    	int x1=points.get(p).car()+w/2;
		    	int y1=(int) (h-points.get(p).cdr()/maxProb*(h-margin*2)*0.8)-margin;
		    	int x2=points.get(p+1).car()+w/2;
		    	int y2=(int) (h-points.get(p+1).cdr()/maxProb*(h-margin*2)*0.8)-margin;
		    	g2.drawLine(x1, y1, x2, y2);	    
		    }
		    g.setFont(new Font("Arial",Font.PLAIN,20));
		    g2.drawString(rounds.get(i), w-300, i*25+margin+25);
	    }

	    try{
	    	ImageIO.write(im, "png", f);
	    }
	    catch(IOException e){
	    	System.err.println("Error in printing file "+filename);
	    }
	}


	// evaluate significance of each called events
	// calculate p-value from binomial distribution, and peak shape parameter
	private void evaluateSignificance(List<ComponentFeature> compFeatures) {
        Binomial binomial = new Binomial(100, .5, new DRand());
		Poisson poisson = new Poisson(1, new DRand());
        double totalIPCount[] = new double[caches.size()];
        double totalControlCount[] = new double[caches.size()];
        for (int i = 0; i < caches.size(); i++) {
            totalIPCount[i] += caches.get(i).car().getHitCount();
        }
		if(controlDataExist) {
            for (int i = 0; i < caches.size(); i++) {
                totalControlCount[i] += caches.get(i).cdr().getHitCount();
            }
			for (ComponentFeature cf: compFeatures){
				for(int cond=0; cond<caches.size(); cond++){
					// scale control read count by non-specific read count ratio
					double controlCount = cf.getUnscaledControlCounts()[cond];
                    double scaledControlCount = cf.getScaledControlCounts(cond);
					double pValueControl = 1, pValueUniform = 1, pValueBalance = 1, pValuePoisson = 1;
                    int ipCount = (int)Math.ceil(cf.getEventReadCounts(cond));
                    if (ipCount==0){			// if one of the conditino does not have reads, set p-value=1
                    	cf.setPValue(1, cond);
                    	continue;
                    }
                    try{
                        assert (totalIPCount[cond] > 0);
                        assert (ipCount <= totalIPCount[cond]);
                        double p = controlCount / totalControlCount[cond];
                        if (p <= 0) {
                            p = 1.0/totalControlCount[cond];
                        } else if (p >= 1) {
                            System.err.println(String.format("p>=1 at evaluateConfidence from %f/%f", controlCount, totalControlCount[cond]));
                            p = 1.0 - 1.0/totalControlCount[cond];
                        } 
                        binomial.setNandP((int)totalIPCount[cond],p);
                        pValueControl = 1 - binomial.cdf(ipCount) + binomial.pdf(ipCount);
                        p = modelWidth / config.mappable_genome_length;
                        binomial.setNandP((int)totalIPCount[cond],p);
                        pValueUniform = 1 - binomial.cdf(ipCount) + binomial.pdf(ipCount);
                        binomial.setNandP((int)Math.ceil(ipCount + scaledControlCount), .5);
                        pValueBalance = 1 - binomial.cdf(ipCount) + binomial.pdf(ipCount);
                        poisson.setMean(Math.max(scaledControlCount, totalIPCount[cond] * modelWidth / config.mappable_genome_length  ));
                        pValuePoisson = 1 - poisson.cdf(ipCount) + poisson.pdf(ipCount);
                    } catch(Exception err){
                        err.printStackTrace();
                        System.err.println(cf.toString());
                        throw new RuntimeException(err.toString(), err);
                    }
                    if (config.testPValues)
                    	cf.setPValue(Math.max(Math.max(pValuePoisson,pValueBalance),Math.max(pValueControl,pValueUniform)), cond);
                    else
                    	cf.setPValue(StatUtil.binomialPValue(scaledControlCount, scaledControlCount+ipCount), cond);
				}
			}
		} else {
            //			I commented the evaluation of the average shift size
            //			We won't need it cause in our data, we do not shift reads.
            //			Besides, we will consider a fixed region of modelRange bp as the peak region.
			createChromStats(compFeatures);
            //			shift_size = eval_avg_shift_size(compFeatures, constants.num_top_mfold_feats);
			
			Map<String, ArrayList<Integer>> chrom_comp_pair = new HashMap<String, ArrayList<Integer>>();
			for(int i = 0; i < compFeatures.size(); i++) {
				String chrom = compFeatures.get(i).getPosition().getChrom();
				if(!chrom_comp_pair.containsKey(chrom))
					chrom_comp_pair.put(chrom, new ArrayList<Integer>());
				
				chrom_comp_pair.get(chrom).add(i);
			}
			
			ipStrandFivePrimes   = new ArrayList[2];
			Arrays.fill(ipStrandFivePrimes, new ArrayList<StrandedBase>());
			ctrlStrandFivePrimes = new ArrayList[2];
			Arrays.fill(ctrlStrandFivePrimes, new ArrayList<StrandedBase>());
			ipStrandFivePrimePos   = new int[2][];
			ctrlStrandFivePrimePos = new int[2][];
			
			for(String chrom:chrom_comp_pair.keySet()) {
				int chromLen = gen.getChromLength(chrom);
				for(int c = 0; c < numConditions; c++) {
					
					for(int i:chrom_comp_pair.get(chrom)) {
						ComponentFeature cf = compFeatures.get(i);
						Region expandedRegion = new Region(gen, chrom, Math.max(0, cf.getPosition().getLocation()-config.third_lambda_region_width), Math.min(chromLen-1, cf.getPosition().getLocation()+config.third_lambda_region_width));
						ipStrandFivePrimes[0] = caches.get(c).car().getStrandedBases(expandedRegion, '+');
						ipStrandFivePrimes[1] = caches.get(c).car().getStrandedBases(expandedRegion, '-');
						
                        ctrlStrandFivePrimes = ipStrandFivePrimes.clone();
						
						for(int k = 0; k < ipStrandFivePrimes.length; k++) {
							ipStrandFivePrimePos[k]   = new int[ipStrandFivePrimes[k].size()];
							ctrlStrandFivePrimePos[k] = new int[ctrlStrandFivePrimes[k].size()];
							for(int v = 0; v < ipStrandFivePrimePos[k].length; v++)
								ipStrandFivePrimePos[k][v] = ipStrandFivePrimes[k].get(v).getCoordinate();
							for(int v = 0; v < ctrlStrandFivePrimePos[k].length; v++)
								ctrlStrandFivePrimePos[k][v] = ctrlStrandFivePrimes[k].get(v).getCoordinate();
						}

						double local_lambda = estimateLocalLambda(cf, c);
						cf.setControlReadCounts(local_lambda, c);                        
						if (config.testPValues)
							poisson.setMean(Math.max(local_lambda, totalIPCount[c] * modelWidth / config.mappable_genome_length));
						else
							poisson.setMean(local_lambda);

                        int count = (int)Math.ceil(cf.getEventReadCounts(c));
                        double pValue = 1 - poisson.cdf(count) + poisson.pdf(count);
						cf.setPValue_wo_ctrl(pValue, c);						
						for(int k = 0; k < ipStrandFivePrimes.length; k++) {
							ipStrandFivePrimes[k].clear();
							ctrlStrandFivePrimes[k].clear();
						}
					}//end of for(int i:chrom_comp_pair.get(chrom)) LOOP	
				}				
			}			
		}
		// calculate q-values, correction for multiple testing
		benjaminiHochbergCorrection(compFeatures);
	}//end of evaluateConfidence method

	//Multiple hypothesis testing correction
	// sort peaks by p-value
	// this method will set -log10(q-value) of component features
	private void benjaminiHochbergCorrection(List<ComponentFeature> compFeatures){
		for (int cond=0; cond<conditionNames.size(); cond++){
			ComponentFeature.setSortingCondition(cond);
			Collections.sort(compFeatures, new Comparator<ComponentFeature>(){
                    public int compare(ComponentFeature o1, ComponentFeature o2) {
                        if(controlDataExist) 
                            return o1.compareByPValue(o2);
                        else
                            return o1.compareByPValue_wo_ctrl(o2);
                    }
                });
			double total = compFeatures.size();
			int rank =1;
			for(ComponentFeature cf : compFeatures){
				if(controlDataExist) 
					cf.setQValueLog10( -Math.log10(cf.getPValue(cond)*(total/(double)rank)), cond);
				else
					cf.setQValueLog10( -Math.log10(cf.getPValue_wo_ctrl(cond)*(total/(double)rank)), cond);
				
				rank++;
			}
		}
	}//end of benjaminiHochbergCorrection method

	private void createChromStats(List<ComponentFeature> compFeatures) {
		/******************************************
		 *****  Create chromosome statistics  *****
		 *****************************************/
		totalIPCounts           = new HashMap<String, Integer>();
		condHitCounts           = new HashMap<String, List<List<Integer>>>();
		
		// In this loop we will store the total IP counts for each chromosome (for all conditions): totalIPCounts
		// as well as the counts for each channel (IP, CTRL) and for each condition: condHitCounts
		for(String chrom:gen.getChromList()){

			Region chromRegion = new Region(gen, chrom, 0, gen.getChromLength(chrom)-1);			
			List<List<StrandedBase>> ip_chrom_signals = loadBasesInWindow(chromRegion, "IP");
			List<List<StrandedBase>> ctrl_chrom_signals = new ArrayList<List<StrandedBase>>();
			if(controlDataExist) {
				ctrl_chrom_signals = loadBasesInWindow(chromRegion, "CTRL");
			}
			else {
				for(int t = 0; t < numConditions; t++)
					ctrl_chrom_signals.add(new ArrayList<StrandedBase>());
			}
			
			int counts = 0;
			// Read counts. List 0 for IP, List 1 for CTRL.
			// Each List contains the read counts for each condition
			List<List<Integer>> currChromCondCounts = new ArrayList<List<Integer>>();
			currChromCondCounts.add(new ArrayList<Integer>()); currChromCondCounts.add(new ArrayList<Integer>());
			for(List<StrandedBase> ip_chrom_signal_cond:ip_chrom_signals) {
				int currCondHitCounts = (int)StrandedBase.countBaseHits(ip_chrom_signal_cond);
				currChromCondCounts.get(0).add(currCondHitCounts);
				counts += currCondHitCounts;
			}

			if(controlDataExist) {
				for(List<StrandedBase> ctrl_chrom_signal_cond:ctrl_chrom_signals) {
					int currCondHitCounts = (int)StrandedBase.countBaseHits(ctrl_chrom_signal_cond);
					currChromCondCounts.get(1).add(currCondHitCounts);
				}
			}

			totalIPCounts.put(chrom, counts);
			condHitCounts.put(chrom, currChromCondCounts);
		}//end of for(String chrom:gen.getChromList()) loop

	}//end of createChromStats method


	private double eval_shift_size(ComponentFeature cf) {
		double shift_size = Double.NaN;
		int plus_mode_pos, minus_mode_pos;

		Map<String, StrandedBase> bases_p_map = new HashMap<String, StrandedBase>();
		Map<String, StrandedBase> bases_m_map = new HashMap<String, StrandedBase>();
		StrandedBase[] bases_p;
		StrandedBase[] bases_m;

		Region r = cf.getPosition().expand(modelRange);
		ArrayList<List<StrandedBase>> ip_signals = loadBasesInWindow(r, "IP");

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

	/*
	 * evalute significance when there is no control.  Returns the estimate of the 
     * local read density 
	 */
	private double estimateLocalLambda(ComponentFeature cf, int c) {
		double local_lambda, lambda_bg, first_lambda_ip, second_lambda_ip, third_lambda_ip, lambda_peak_ctrl, first_lambda_ctrl, second_lambda_ctrl, third_lambda_ctrl;

		Point pos = cf.getPosition();
		String chrom      = pos.getChrom();
		int chromLen      = gen.getChromLength(chrom);
		Region peakRegion = pos.expand(modelRange);

		int left_peak,  left_third_region,  left_second_region,  left_first_region;
		int right_peak, right_third_region, right_second_region, right_first_region;

		int num_first_lambda_ip = 0, num_second_lambda_ip = 0, num_third_lambda_ip = 0;
		int num_peak_ctrl = 0, num_first_lambda_ctrl = 0, num_second_lambda_ctrl = 0, num_third_lambda_ctrl = 0;

		double num_peak_ip;

		double ipCounts = condHitCounts.get(chrom).get(0).get(c);
		double ctrlCounts;
		if(controlDataExist)
			ctrlCounts = condHitCounts.get(chrom).get(1).get(c);
		else
			ctrlCounts = ipCounts;

		left_peak           = (int)Math.min(chromLen-1, peakRegion.getStart());
		left_third_region   = (int)Math.max(         0, pos.getLocation() - config.third_lambda_region_width/2);
		left_second_region  = (int)Math.max(         0, pos.getLocation() - config.second_lambda_region_width/2);
		left_first_region   = (int)Math.max(         0, pos.getLocation() - config.first_lambda_region_width/2);

		right_peak          = (int)Math.max(         0, peakRegion.getEnd());
		right_third_region  = (int)Math.min(chromLen-1, pos.getLocation() + config.third_lambda_region_width/2);
		right_second_region = (int)Math.min(chromLen-1, pos.getLocation() + config.second_lambda_region_width/2);
		right_first_region  = (int)Math.min(chromLen-1, pos.getLocation() + config.first_lambda_region_width/2);

		int smallestLeft = (int) Math.min(left_peak, Math.min(left_third_region, Math.min(left_second_region, left_first_region)));
		int largestRight = (int) Math.max(right_peak, Math.max(right_third_region, Math.max(right_second_region, right_first_region)));

		// Get the number of reads assigned to the peak - soft assignment
		num_peak_ip = cf.getEventReadCounts(c);

		// k = 0: '+' strand, k = 1: '-' strand
		for(int k = 0; k <= 1; k++) {
			
			int s_idx, e_idx;
			
			if(ipStrandFivePrimePos[k].length > 0) {
				// IP Channel
				s_idx = Arrays.binarySearch(ipStrandFivePrimePos[k], smallestLeft);
				e_idx = Arrays.binarySearch(ipStrandFivePrimePos[k], largestRight);
				if(s_idx < 0) { s_idx = -s_idx-1; }
				if(e_idx < 0) { e_idx = -e_idx-1; }
				s_idx = StatUtil.searchFrom(ipStrandFivePrimePos[k], ">=", smallestLeft, s_idx);
				e_idx = StatUtil.searchFrom(ipStrandFivePrimePos[k], "<=", largestRight, e_idx);

				for(int i = s_idx; i <= e_idx; i++) {
					if(left_first_region <= ipStrandFivePrimePos[k][i] && ipStrandFivePrimePos[k][i] <= right_first_region) {
						num_first_lambda_ip  += ipStrandFivePrimes[k].get(i).getCount();
						num_second_lambda_ip += ipStrandFivePrimes[k].get(i).getCount();
						num_third_lambda_ip  += ipStrandFivePrimes[k].get(i).getCount();
					}
					else if(left_second_region <= ipStrandFivePrimePos[k][i] && ipStrandFivePrimePos[k][i] <= right_second_region) {
						num_second_lambda_ip += ipStrandFivePrimes[k].get(i).getCount();
						num_third_lambda_ip  += ipStrandFivePrimes[k].get(i).getCount();
					}
					else if(left_third_region <= ipStrandFivePrimePos[k][i] && ipStrandFivePrimePos[k][i] <= right_third_region) {
						num_third_lambda_ip  += ipStrandFivePrimes[k].get(i).getCount();
					}
				}	
			}

			if(ctrlStrandFivePrimePos[k].length > 0) {
				// CTRL Channel
				s_idx = Arrays.binarySearch(ctrlStrandFivePrimePos[k], smallestLeft);
				e_idx = Arrays.binarySearch(ctrlStrandFivePrimePos[k], largestRight);
				if(s_idx < 0) { s_idx = -s_idx-1; }
				if(e_idx < 0) { e_idx = -e_idx-1; }
				s_idx = StatUtil.searchFrom(ctrlStrandFivePrimePos[k], ">=", smallestLeft, s_idx);
				e_idx = StatUtil.searchFrom(ctrlStrandFivePrimePos[k], "<=", largestRight, e_idx);
				for(int i = s_idx; i <= e_idx; i++) {
					if(left_peak <= ctrlStrandFivePrimePos[k][i] && ctrlStrandFivePrimePos[k][i] <= right_peak)
						num_peak_ctrl += ctrlStrandFivePrimes[k].get(i).getCount();

					if(left_first_region <= ctrlStrandFivePrimePos[k][i] && ctrlStrandFivePrimePos[k][i] <= right_first_region) {
						num_first_lambda_ctrl  += ctrlStrandFivePrimes[k].get(i).getCount();
						num_second_lambda_ctrl += ctrlStrandFivePrimes[k].get(i).getCount();
						num_third_lambda_ctrl  += ctrlStrandFivePrimes[k].get(i).getCount();
					}
					else if(left_second_region <= ctrlStrandFivePrimePos[k][i] && ctrlStrandFivePrimePos[k][i] <= right_second_region) {
						num_second_lambda_ctrl += ctrlStrandFivePrimes[k].get(i).getCount();
						num_third_lambda_ctrl  += ctrlStrandFivePrimes[k].get(i).getCount();
					}
					else if(left_third_region <= ctrlStrandFivePrimePos[k][i] && ctrlStrandFivePrimePos[k][i] <= right_third_region) {
						num_third_lambda_ctrl  += ctrlStrandFivePrimes[k].get(i).getCount();
					}
				}
			}
		}//end of for(int k = 0; k <= 1; k++) LOOP

		// Subtract the number of reads that are inside the peak region
		// We just want to compare it against the surrounding regions (excl. the peak specific signal)
		num_first_lambda_ip  = Math.max(0, num_first_lambda_ip  - (int)num_peak_ip);
		num_second_lambda_ip = Math.max(0, num_second_lambda_ip - (int)num_peak_ip);
		num_third_lambda_ip  = Math.max(0, num_third_lambda_ip  - (int)num_peak_ip);
		if(!controlDataExist) {
			num_peak_ctrl          = Math.max(0, num_peak_ctrl          - (int)num_peak_ip);
			num_first_lambda_ctrl  = Math.max(0, num_first_lambda_ctrl  - (int)num_peak_ip);
			num_second_lambda_ctrl = Math.max(0, num_second_lambda_ctrl - (int)num_peak_ip);
			num_third_lambda_ctrl  = Math.max(0, num_third_lambda_ctrl  - (int)num_peak_ip);
			ctrlCounts             = Math.max(0, ctrlCounts             - (int)num_peak_ip);
		}

		first_lambda_ip      = 1.0*num_first_lambda_ip*((double)peakRegion.getWidth()/(double)(config.first_lambda_region_width-peakRegion.getWidth()));
		second_lambda_ip     = 1.0*num_second_lambda_ip*((double)peakRegion.getWidth()/(double)(config.second_lambda_region_width-peakRegion.getWidth()));
		third_lambda_ip      = 1.0*num_third_lambda_ip*((double)peakRegion.getWidth()/(double)(config.third_lambda_region_width-peakRegion.getWidth()));

		double ip2ctrl_ratio   = ipCounts/ctrlCounts;
		int skippedRegionWidth = controlDataExist ? 0 : peakRegion.getWidth();
		lambda_peak_ctrl       = 1.0*num_peak_ctrl*ip2ctrl_ratio;
		first_lambda_ctrl      = 1.0*num_first_lambda_ctrl*ip2ctrl_ratio*((double)peakRegion.getWidth()/(double)(config.first_lambda_region_width-skippedRegionWidth));
		second_lambda_ctrl     = 1.0*num_second_lambda_ctrl*ip2ctrl_ratio*((double)peakRegion.getWidth()/(double)(config.second_lambda_region_width-skippedRegionWidth));
		third_lambda_ctrl      = 1.0*num_third_lambda_ctrl*ip2ctrl_ratio*((double)peakRegion.getWidth()/(double)(config.third_lambda_region_width-skippedRegionWidth));
		lambda_bg              = 1.0*ctrlCounts*ip2ctrl_ratio*((double)peakRegion.getWidth()/(double)(chromLen-peakRegion.getWidth()));


        //		Use exactly what they do in MACS statistical method
        //		double ip2ctrl_ratio = ipCounts/ctrlCounts;
        //		first_lambda_ip    = 1.0*num_first_lambda_ip*((double)peakRegion.getWidth()/(double)config.first_lambda_region_width);
        //		second_lambda_ip   = 1.0*num_second_lambda_ip*((double)peakRegion.getWidth()/(double)config.second_lambda_region_width);
        //		third_lambda_ip    = 1.0*num_third_lambda_ip*((double)peakRegion.getWidth()/(double)config.third_lambda_region_width);
        //		lambda_peak_ctrl   = 1.0*num_peak_ctrl*ip2ctrl_ratio;
        //		first_lambda_ctrl  = 1.0*num_first_lambda_ctrl*ip2ctrl_ratio*((double)peakRegion.getWidth()/(double)config.first_lambda_region_width);
        //		second_lambda_ctrl = 1.0*num_second_lambda_ctrl*ip2ctrl_ratio*((double)peakRegion.getWidth()/(double)config.second_lambda_region_width);
        //		third_lambda_ctrl  = 1.0*num_third_lambda_ctrl*ip2ctrl_ratio*((double)peakRegion.getWidth()/(double)config.third_lambda_region_width);
        //		lambda_bg          = 1.0*ctrlCounts*ip2ctrl_ratio*((double)peakRegion.getWidth()/(double)chromLen);

		if(controlDataExist)
			local_lambda = Math.max(lambda_bg, Math.max(first_lambda_ip, Math.max(second_lambda_ip, Math.max(third_lambda_ip, Math.max(lambda_peak_ctrl, Math.max(first_lambda_ctrl, Math.max(second_lambda_ctrl, third_lambda_ctrl)))))));
		else //skip first lambda regions (1K)
			local_lambda = Math.max(lambda_bg, Math.max(second_lambda_ip, Math.max(third_lambda_ip, Math.max(second_lambda_ctrl, third_lambda_ctrl))));

        return local_lambda;
	}//end of evalFeatureSignificance method

	/**
	 * Count non-specific reads  (all the reads not in the refinedRegions)
	 * for scaling experiment and control reads.
	 * Because the signal (specific peak region) can vary a lot in different conditions,
	 * the non-specific regions are a fair estimate for noise.
	 * We assume the noise in expt and control should be comparable, can be used to scale reads.
	 */
	public void countNonSpecificReads(List<ComponentFeature> compFeatures){
		if (!wholeGenomeDataLoaded){	
			ratio_non_specific_total=ratio_total;
			ComponentFeature.setNon_specific_ratio(ratio_total);
			return;
		}
        //		System.out.println("\nCounting non-specific read numbers ...\n");

		// construct the refined specific binding regions
		ArrayList<Region> refinedRegions = new ArrayList<Region>();
		for (ComponentFeature cf:compFeatures){
			refinedRegions.add(cf.getPosition().expand(0));
		}
		refinedRegions = mergeRegions(refinedRegions, true);
		
		//Count the specific read numbers
		int expt_test_region_total[]=new int[numConditions];
		int crtl_test_region_total[]=new int[numConditions];
		int expt_non_specific_total[]=new int[numConditions];
		int crtl_non_specific_total[]=new int[numConditions];
		int totalLength=0;	// total length of non-overlapping peak regions
		for(Region r : refinedRegions){
			totalLength += r.getWidth();
			for(int i=0; i<caches.size(); i++){
				Pair<ReadCache,ReadCache> e = caches.get(i);
				expt_test_region_total[i] += e.car().countHits(r);
				if(controlDataExist)
					crtl_test_region_total[i] += e.cdr().countHits(r);
			}
		}
		log(1, "\nTotal length of specific binding regions: "+totalLength);

		// non-specific = total - specific
		for(int i=0; i<numConditions; i++){
			Pair<ReadCache,ReadCache> e = caches.get(i);
			expt_non_specific_total[i]=(int)e.car().getHitCount()-expt_test_region_total[i];
			if(controlDataExist) {
				crtl_non_specific_total[i]=(int)e.cdr().getHitCount()-crtl_test_region_total[i];
				ratio_non_specific_total[i] = (double)expt_non_specific_total[i]/crtl_non_specific_total[i];
			}

	    	double noiseReadNum_per_kbp = (double) expt_non_specific_total[i] * 1000
	    		/ (config.mappable_genome_length - totalLength);

	    	StringBuilder sb = new StringBuilder();
	    	sb.append(conditionNames.get(i)+" data summary:\n");
	    	sb.append("\tIP total   \t\t" +(int)e.car().getHitCount());
	    	if (controlDataExist){
	    		sb.append("\n\tControl total\t\t" + (int)e.cdr().getHitCount() );
	    	    sb.append("\n\tIP/Control  \t\t" +String.format("%.3f", e.car().getHitCount()/e.cdr().getHitCount() ));
	    	    sb.append("\n\tIP non-specific\t\t" +expt_non_specific_total[i]);
	    	    sb.append("\n\tControl non-specific\t" +crtl_non_specific_total[i]);
                //	    	    sb.append("\nRatio non-specific\t" +String.format("%.3f",ratio_non_specific_total[i])+
	    	}
	    	sb.append("\n\tNoise reads per 1000bp\t" +String.format("%.2f",noiseReadNum_per_kbp)+"\n");
	    	log(1, sb.toString());
		}
		ComponentFeature.setNon_specific_ratio(ratio_non_specific_total);

        //		log(1, "countNonSpecificReads(): "+timeElapsed(tic));
	}

	private void addConfigString(String name, boolean value){
		configsb.append(name).append("\t").append(new Boolean(value).toString()).append("\n");
	}
	public String getConfigString(){
		return "BindingMixture Configurations\n"+configsb.toString();
	}

	public void addAnnotations(){
		try{
			//Add closest genes
			System.out.println("Adding gene annotations");
			addClosestGenes(signalFeatures);
	
			//Add other annotations
			System.out.println("Adding other annotations");
			addRegionAnnotations(signalFeatures);
		}
		catch(Exception e){
			System.err.println("Error in adding annotations.");
			e.printStackTrace(System.err);
		}
	}

	/* Default printing option
	 * for single condition, sort by p-value
	 * for multiple conditions, sort by location
	 */
	public void printFeatures(){
		// fs is sorted by location by default
		ArrayList<ComponentFeature> fs = new ArrayList<ComponentFeature>();
		for(Feature f : signalFeatures){
			ComponentFeature cf = (ComponentFeature)f;
			fs.add(cf);
		}
		//for single condition, sort by p-value
		if (!config.sort_by_location && this.numConditions==1){
			Collections.sort(fs, new Comparator<ComponentFeature>(){
                    public int compare(ComponentFeature o1, ComponentFeature o2) {
                        if(controlDataExist)
                            return o1.compareByPValue(o2);
                        else
                            return o1.compareByPValue_wo_ctrl(o2);
                    }
                });
		}	
		else
			Collections.sort(fs);
		String fname = outName+"_GPS_significant.txt";
		printFeatures(fname, fs);
	}

	public void printInsignificantFeatures(){
		if (insignificantFeatures==null || insignificantFeatures.size()==0)
			return;
		
		ArrayList<ComponentFeature> fs = new ArrayList<ComponentFeature>();
		for(Feature f : insignificantFeatures){
			fs.add((ComponentFeature)f);
		}
		//for single condition, sort by p-value
		if (!config.sort_by_location && this.numConditions==1){
			Collections.sort(fs, new Comparator<ComponentFeature>(){
                    public int compare(ComponentFeature o1, ComponentFeature o2) {
                        if(controlDataExist)
                            return o1.compareByPValue(o2);
                        else
                            return o1.compareByPValue_wo_ctrl(o2);
                    }
                });
		}
		else
			Collections.sort(fs);
		String fname = outName+"_GPS_insignificant.txt";
		printFeatures(fname, fs);
	}
	
	public void printFilteredFeatures(){
		if (filteredFeatures!=null && !filteredFeatures.isEmpty()){
			ArrayList<ComponentFeature> fs = new ArrayList<ComponentFeature>();
			for(Feature f : filteredFeatures){
				fs.add((ComponentFeature)f);
			}
			//for single condition, sort by IP count
			if (numConditions==1){
				Collections.sort(fs, new Comparator<ComponentFeature>(){
                        public int compare(ComponentFeature o1, ComponentFeature o2) {
                            return o1.compareByTotalResponsibility(o2);
                        }
                    });
			}
			String fname = outName+"_GPS_filtered.txt";
			printFeatures(fname, fs);
		}
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

		String fname = outName+"_GPS_significant.txt";
		printFeatures(fname, fs);
	}

	private void printFeatures(String fname, ArrayList<ComponentFeature> fs){
		try{
			FileWriter fw = new FileWriter(fname);
			boolean first=true;
			for(ComponentFeature f : fs){
				if(first){
					fw.write(f.headString_v1());
					first=false;
				}
				fw.write(f.toString_v1());
			}
			fw.close();
			
			// BED format
			if(config.outputBED){
				fname=fname.replaceAll("txt", "bed");
				fw = new FileWriter(fname);
				first=true;
				for(ComponentFeature f : fs){
					if(first){
						fw.write("track name=GPS_"+outName+" description=\"GPS Event Call\"\n");
                        //					fw.write(f.headString_BED());
						first=false;
					}
					fw.write(f.toBED());
				}
				fw.close();
			}	
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
			List<StrandedBase> bases_p= e.cdr().getStrandedBases(r, '+');  
			List<StrandedBase> bases_m= e.cdr().getStrandedBases(r, '-');  
			count+=StrandedBase.countBaseHits(bases_p)+StrandedBase.countBaseHits(bases_m);
		}
		return count;
	}	
	private float countCtrlReads(Region r, int cond){
		List<StrandedBase> bases_p= caches.get(cond).cdr().getStrandedBases(r, '+');  
		List<StrandedBase> bases_m= caches.get(cond).cdr().getStrandedBases(r, '-');  
		return StrandedBase.countBaseHits(bases_p)+StrandedBase.countBaseHits(bases_m);
	}


	public void writeDebugFile(){
		try{
			FileWriter fw = new FileWriter("Binding Mixture Debug.txt", false); //new file
			fw.write(log_all_msg.toString());
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	// show the message in STDOUT and log file
	// depending on "verbose" setting in property file
	private void log(int mode, String msg){
		if (constants.LOG_ALL)
			log_all_msg.append(msg);
		if ( config.bmverbose>=mode){
	    	System.out.println(msg);
			try{
				logFileWriter.write(msg+"\n");
				logFileWriter.flush();
			} catch (IOException e) {
				e.printStackTrace();
			}
			if (constants.LOG_ALL)
				log_all_msg.append("\n");
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

	public void printError() {

	}
	private void printTowerRegions(){
		// save the list of regions to file
		try{
			StringBuilder fileName = new StringBuilder(outName).append("_");
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
	private void printExcludedRegions(){
		if (excludedRegions.isEmpty())
			return;
		
		// save the list of regions to file
		try{
			StringBuilder fileName = new StringBuilder(outName).append("_ExcludedRegions.txt");
			FileWriter fw = new FileWriter(fileName.toString());

			StringBuilder txt = new StringBuilder();
			for (int i=0;i<excludedRegions.size();i++){
				txt.append(excludedRegions.get(i).toString()).append("\n");
			}
			fw.write(txt.toString());
			fw.close();
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}
	/*
	 * Calculate average square distance for read profile 
	 * 
	 */
	private double calcSqDiff(double[]profile){
		double sqDiff = 0;
		double[] m = model.getProbabilities();
		int nzPosCount = 0;
		double totalCount = 0;
		double nzProbSum = 0;
		for (int i=0;i<profile.length;i++){
			if (profile[i]!=0){
				nzPosCount++;
				totalCount+=profile[i];
				nzProbSum += m[i];
			}
		}
		if (nzPosCount==0)
			return config.shapeDeviation*2*0.75;
		for (int i=0;i<profile.length;i++){
			if (profile[i]!=0){
				double expected = totalCount*m[i]/nzProbSum;
				sqDiff += (profile[i]-expected)*(profile[i]-expected);
			}
		}
		return sqDiff/totalCount;
	}

    class GPSConstants {

        /****************************
         * Constants
         ***************************/
	
        public final static boolean LOG_ALL=false;

        // width for smoothing a read (used as a stddev for creating the Gaussian kernel of probability)
        public final static int READ_KERNEL_ESTIMATOR_WIDTH = 5;

        //Maximum number of components that EM can handle efficiently
        public final int MAX_NUM_COMPONENTS=1000;
        public final int OPTIMAL_NUM_COMPONENTS=100;
        public final int INIT_SPACING=5;
        // Maximum number of EM iterations
        public final int MAX_EM_ITER=10000;
        public final int MAX_EM_ML_ITER=1;
        //Run EM up until <tt>ANNEALING_ITER</tt> with smaller alpha based on the current iteration
        public final int ANNEALING_ITER=200;

        //EM convergence between the likelihood of the current and the previous step
        public final double EM_CONVERGENCE = 1e-5;
        public final double EM_ML_CONVERGENCE = 1e-8;

        public final int num_top_mfold_feats = 1000;
        public final int maxUpdateIters = 7;

        //	public final double lambda = 0;

        // Max number of reads to load from DB
        public final static int MAXREAD = 1000000;
        // true:  eliminated component in batch, as long as matching criteria in EM derivation
        // false: eliminate only the worse case components, re-distribute the reads of eliminated component
        public final static boolean BATCH_ELIMINATION = false;
        public final static boolean SMART_SPACING = true;		// dynamically determine init comopent spacing
        public final static boolean MAKE_HARD_ASSIGNMENT=false;
    }

    class GPSConfig {
        public boolean do_model_selection=false;
        public boolean use_joint_event = false;
        public boolean TF_binding = true;
        public boolean outputBED = false;
        public boolean testPValues = false;
        public boolean post_artifact_filter=false;
        public boolean filterEvents=false;
        public boolean kl_count_adjusted = false;
        public boolean sort_by_location=false;
        public boolean subtract_for_segmentation=false;
        public boolean exclude_unenriched = false;
        public int KL_smooth_width = 0;
        public int max_hit_per_bp = -1;
        /** percentage of candidate (enriched) peaks to take into account
         *  during the evaluation of non-specific signal */
        public double pcr = 0.0;
        public double q_value_threshold = 2.0;
        public double joint_event_distance = 500;
        public double alpha_factor = 3.0;
        public int top_events = 2000;
        public int min_event_count = 500;	// minimum num of events to update read distribution
        public int smooth_step = 30;
        public int window_size_factor = 3;	//number of model width per window
        public int min_region_width = 50;	//minimum width for select enriched region
        public double mappable_genome_length = 2.08E9; // mouse genome
        public double sparseness=6.0;
        public double fold = 3.0;
        public double kl_ic = -1.0;
        public double shapeDeviation = 0;
        public int gentle_elimination_factor = 2;	// factor to reduce alpha to a gentler pace after eliminating some component
        public int resolution_extend = 2;
        public int first_lambda_region_width  =  1000;
        public int second_lambda_region_width =  5000;
        public int third_lambda_region_width  = 10000;
        public boolean use_dynamic_sparseness = true;
        public boolean use_betaEM  = true;
        public boolean use_scanPeak  = true;
        public boolean refine_regions = false;		// refine the enrichedRegions for next round using EM results
        public int bmverbose=1;		// BindingMixture verbose mode
        public int base_reset_threshold = 200;	// threshold to set a base read count to 1
        public int windowSize;			// size for EM sliding window for splitting long regions
        //Run EM up until <tt>ML_ITER</tt> without using sparse prior
        public int ML_ITER=10;
        // the range to scan a peak if we know position from EM result
        public int SCAN_RANGE = 20;
        public int gentle_elimination_iterations = 5;

        public void parseArgs(String args[]) {
            Set<String> flags = Args.parseFlags(args);
            // default as false, need the flag to turn it on
            sort_by_location = flags.contains("sl");
            use_joint_event = flags.contains("refine_using_joint_event");
            post_artifact_filter = flags.contains("post_artifact_filter");
            kl_count_adjusted = flags.contains("adjust_kl");
            refine_regions = flags.contains("refine_regions");
            subtract_for_segmentation = flags.contains("subtract_ctrl_for_segmentation");
            outputBED = flags.contains("outBED");
            testPValues = flags.contains("testP");
            if (testPValues)
            	System.err.println("testP is " + testPValues);
            exclude_unenriched = flags.contains("ex_unenriched");
            // default as true, need the opposite flag to turn it off
            use_dynamic_sparseness = ! flags.contains("fa"); // fix alpha parameter
            use_betaEM = ! flags.contains("poolEM");
            filterEvents = !flags.contains("nf");	// not filtering of predicted events
            TF_binding = ! flags.contains("br");	// broad region, not TF data, is histone or pol II
            if (!TF_binding){
                use_joint_event = true;
                sort_by_location = true;
            }
            use_scanPeak = ! flags.contains("no_scanPeak");
            do_model_selection = !flags.contains("no_model_selection");
            mappable_genome_length = Args.parseDouble(args, "s", 2.08E9);	// size of mappable genome
            // Optional input parameter
            maxThreads = Args.parseInteger(args,"threads",java.lang.Runtime.getRuntime().availableProcessors());	// default to the # processors
            q_value_threshold = Args.parseDouble(args, "q", 2.0);	// q-value
            sparseness = Args.parseDouble(args, "a", 6.0);	// minimum alpha parameter for sparse prior
            alpha_factor = Args.parseDouble(args, "af", alpha_factor); // denominator in calculating alpha value
            fold = Args.parseDouble(args, "fold", fold); // minimum fold enrichment IP/Control for filtering
            shapeDeviation =  TF_binding?-0.45:-0.3;		// set default according to filter type    		
            shapeDeviation = Args.parseDouble(args, "sd", shapeDeviation); // maximum shapeDeviation value for filtering
            max_hit_per_bp = Args.parseInteger(args, "mrc", 0); //max read count per bp, default -1, estimate from data
            pcr = Args.parseDouble(args, "pcr", 0.0); // percentage of candidate (enriched) peaks to be taken into account during the evaluation of the non-specific slope
            window_size_factor = Args.parseInteger(args, "wsf", 3);
            second_lambda_region_width = Args.parseInteger(args, "w2", 5000);
            third_lambda_region_width = Args.parseInteger(args, "w3", 10000);
            joint_event_distance = Args.parseInteger(args, "j", 10000);		// max distance of joint events
            top_events = Args.parseInteger(args, "top", top_events);
            min_event_count = Args.parseInteger(args, "min", min_event_count);
            base_reset_threshold = Args.parseInteger(args, "reset", base_reset_threshold);
            min_region_width = Args.parseInteger(args, "min_region_width", 50);
            bmverbose = Args.parseInteger(args, "bmverbose", bmverbose);
            smooth_step = Args.parseInteger(args, "smooth", smooth_step);
            KL_smooth_width = Args.parseInteger(args, "kl_s_w", KL_smooth_width);
            kl_ic = Args.parseDouble(args, "kl_ic", kl_ic);
            resolution_extend = Args.parseInteger(args, "resolution_extend", resolution_extend);
            gentle_elimination_factor = Args.parseInteger(args, "gentle_elimination_factor", gentle_elimination_factor);
            // These are options for EM performance tuning
            // should NOT expose to user
            // therefore, still use UPPER CASE to distinguish
            ML_ITER = Args.parseInteger(args, "ML_ITER", ML_ITER);
            SCAN_RANGE = Args.parseInteger(args, "SCAN_RANGE", SCAN_RANGE);
        }
    }


    class GPSThread implements Runnable {

        private GPSMixture mixture;
        private GPSConstants constants;
        private GPSConfig config;
        private Collection<Region> regions;
        private Collection<Integer> processedRegionCount;
        private Collection<ComponentFeature> compFeatures;
        /**
         * <tt>HashMap</tt> containing the single event regions. <br>
         * Each <tt>key</tt> contains the (single event) region of interest and the corresponding <tt>value</tt>
         * the smaller sub-region (within the <tt>key</tt> region) containing the binding location
         * The purpose is that if we re-evaluate with new empirical distribution, we do not re-run EM for unary event.
         * The assumption is that changing empirical distribution will not have changed the existence of events,
         * it will only modify the position of the event slightly, so we can just do the scanning.
         */
        private HashMap<Region, Point> singleEventRegions=new HashMap<Region, Point>();
        //Number of non zero components of specified region
        private int nonZeroComponentNum=0;
        //Components representing the binding events in the analyzed window
        private ArrayList<BindingComponent> components = new ArrayList<BindingComponent>();
        //Spacing between the components
        private int componentSpacing=1;
        // Maximum number of components determined by the data
        private int componentMax;

        public GPSThread (Collection<Region> regions,
        				  Collection<Integer> processedRegionCount,
                          Collection<ComponentFeature> compFeatures,
                          GPSMixture mixture,
                          GPSConstants constants,
                          GPSConfig config) {
            this.regions = regions;
            this.processedRegionCount = processedRegionCount;
            this.mixture = mixture;
            this.constants = constants;
            this.config = config;
            this.compFeatures = compFeatures;
        }
        public void simpleRun(List<StrandedBase> bases, Region r) {
            components = new ArrayList<BindingComponent>();
            ArrayList<List<StrandedBase>> signals = new ArrayList<List<StrandedBase>>();
            signals.add(bases);

            Pair<double[][][], int[][][]> result = null;
            if (!mixture.doScanning){
                //Run EM and increase resolution
                initializeComponents(r, mixture.numConditions );
                while(nonZeroComponentNum>0){
                    double alpha = Math.max(Math.sqrt(StrandedBase.countBaseHits(bases))/config.alpha_factor, config.sparseness);
                    result = EMTrain(signals, alpha);
                    if(componentSpacing==1)
                        break;
                    updateComponentResolution(r, 1, componentSpacing);
                }
                setComponentResponsibilities(signals, result.car(), result.cdr());
            } else{		// scan it
                BindingComponent peak = mixture.scanPeak(signals, r);
                components.add(peak);
            }
            compFeatures.addAll(mixture.callFeatures(components));
        }
        public void run() {
            for (Region rr : regions) {
                mixture.log(2, rr.toString());
                try{
                    ArrayList<BindingComponent> comps= new ArrayList<BindingComponent>();
                    // Cut long regions into windowSize(1.5kb) sliding window (500bp overlap) to analyze
                    ArrayList<Region> windows = new ArrayList<Region>();
                    if (rr.getWidth()<=config.windowSize)
                        windows.add(rr);
                    else{
                        windows = mixture.splitWindows(rr, config.windowSize, mixture.modelWidth/2);
                    }
		 
                    // run EM for each window 
                    for (Region w : windows){
                        ArrayList<BindingComponent> result = analyzeWindow(w);
                        if (result!=null){
                            comps.addAll(result);
                        }
                    }
		
                    /* ****************************************************************
                     * fix sliding window boundary effect (version 1)
                     */
                    if (windows.size()>1 && comps.size()>0){
                        Collections.sort(comps);
                        // The whole region can be divided into subRegions with gap >= modelWidth
                        ArrayList<ArrayList<Region>> subRegions = new ArrayList<ArrayList<Region>>();
                        ArrayList<Region> rs = new ArrayList<Region>();
                        rs.add(comps.get(0).getLocation().expand(0));
                        subRegions.add(rs);
                        if (comps.size()>1){
                            for (int i=1;i<comps.size();i++){
                                if (comps.get(i).getLocation().distance(comps.get(i-1).getLocation())>mixture.modelWidth){
                                    rs = new ArrayList<Region>();
                                    subRegions.add(rs);
                                    rs.add(comps.get(i).getLocation().expand(0));
                                } else	{ // in same subregions
                                    rs.add(comps.get(i).getLocation().expand(0));
                                }
                            }
                        }
                        for (ArrayList<Region> subr : subRegions){
                            // expand with modelRange padding, and merge overlapped regions
                            // ==> Refined enriched regions
                            rs=mixture.mergeRegions(subr, true);
                            for (Region r:rs){
                                if (r.getWidth()==2*mixture.modelRange+1){	// for unary event, one of the windows should contain it in full, no need to evaluate again
                                    if (subr.size()>1){	// if multiple windows predict the same unary event
                                        Point p = r.getMidpoint();
                                        ArrayList<BindingComponent> duplicates = new ArrayList<BindingComponent>(); 
                                        for (BindingComponent b: comps){
                                            if (b.getLocation().equals(p))
                                                duplicates.add(b);
                                        }
                                        if (duplicates.size()>1){
                                            for (int i=1;i<duplicates.size();i++){
                                                comps.remove(duplicates.get(i));
                                            }
                                        }
                                    }
                                } else { // for joint events 
                                    // the regions in rs includes the influence paddings, remove it here
                                    int start=r.getStart();
                                    int end=r.getEnd();
                                    if (start!=0)
                                        start += mixture.modelRange;
                                    if (end != mixture.gen.getChromID(r.getChrom()))
                                        end -= mixture.modelRange;
                                    if (start>end){
                                        int tmp = start;
                                        end = tmp;
                                        start = end;
                                    }
                                    Region tightRegion = new Region(r.getGenome(), r.getChrom(), start, end);
                                    for (int i=0;i<windows.size()-1;i++){
                                        Region boundary = windows.get(i).getOverlap(windows.get(i+1));
                                        // if the predicted component overlaps with boundary of sliding window
                                        if (boundary.overlaps(tightRegion)){
                                            // remove old
                                            ArrayList<BindingComponent> toRemove = new ArrayList<BindingComponent>();
                                            for (BindingComponent m:comps){
                                                if (tightRegion.overlaps(m.getLocation().expand(0)))
                                                    toRemove.add(m);
                                            }
                                            comps.removeAll(toRemove);
                                            // re-process the boundary region
                                            int winSize = mixture.modelWidth*10;
                                            if (r.getWidth()<winSize){ // if the region is small, directly process it
                                                ArrayList<BindingComponent> result = analyzeWindow(r);
                                                if (result!=null){
                                                    comps.addAll(result);
                                                }
                                            } else {	// if the region is too long, split into windows (5kb size, 1kb overlap)
                                                ArrayList<Region> wins = mixture.splitWindows(r, winSize, mixture.modelWidth);
                                                // process each window, then fix boundary 
                                                ArrayList<ArrayList<BindingComponent>> comps_all_wins= new ArrayList<ArrayList<BindingComponent>>();
                                                for (Region w : wins){
                                                    ArrayList<BindingComponent> comp_win = new ArrayList<BindingComponent>();
                                                    ArrayList<BindingComponent> result = analyzeWindow(w);
                                                    if (result!=null){
                                                        comp_win.addAll(result);
                                                    }
                                                    comps_all_wins.add(comp_win);
                                                }
                                                /* ****************************************************************
                                                 * fix sliding window boundary effect (version 2)
                                                 * take the events in the left half of boundary from left window
                                                 * take the events in the right half of boundary from right window
                                                 * so that the events we reported have at least 500bp data as margin
                                                 * Note: this may result in slight inaccuracy comparing to re-process whole region
                                                 * 		 but it is the trade-off against running EM on a huge large region
                                                 * 		 and running for whole region might lead to inaccuracy too, because of too many components
                                                 */
                                                if (wins.size()>1){
                                                    for (int ii=0;ii<wins.size()-1;ii++){
                                                        int midCoor = wins.get(ii).getOverlap(wins.get(ii+1)).getMidpoint().getLocation();
                                                        ArrayList<BindingComponent> leftComps = comps_all_wins.get(ii);
                                                        toRemove = new ArrayList<BindingComponent>();//reset
                                                        for (BindingComponent m:leftComps){
                                                            if (m.getLocation().getLocation()>midCoor)
                                                                toRemove.add(m);
                                                        }
                                                        leftComps.removeAll(toRemove);
                                                        ArrayList<BindingComponent> rightComps = comps_all_wins.get(ii+1);
                                                        toRemove = new ArrayList<BindingComponent>();	//reset
                                                        for (BindingComponent m:rightComps){
                                                            if (m.getLocation().getLocation()<=midCoor)
                                                                toRemove.add(m);
                                                        }
                                                        rightComps.removeAll(toRemove);
                                                    }
                                                }
                                                for (ArrayList<BindingComponent>comps_win:comps_all_wins){
                                                    comps.addAll(comps_win);
                                                }
                                            }
                                            break;	// break from testing more boundary
                                        }
                                    }
                                }
                            }
                        }
                    }// END: fix sliding window boundary effect

                    /* ****************************************************************
                     * refine unary events and collect all the events as features
                     * this is last step because fixing boundary may result in some new unary events
                     */
                    if (comps.size()>0){
                        singleEventRegions.clear();
                        Collections.sort(comps);
                        // The whole region can be divided into subRegions with gap >= mixture.modelWidth
                        ArrayList<ArrayList<Region>> subRegions = new ArrayList<ArrayList<Region>>();
                        ArrayList<Region> rs = new ArrayList<Region>();
                        rs.add(comps.get(0).getLocation().expand(0));
                        subRegions.add(rs);
                        if (comps.size()>1){
                            for (int i=1;i<comps.size();i++){
                                if (comps.get(i).getLocation().distance(comps.get(i-1).getLocation())>mixture.modelWidth){
                                    rs = new ArrayList<Region>();
                                    subRegions.add(rs);
                                    rs.add(comps.get(i).getLocation().expand(0));
                                }
                                else	// in same subregions
                                    rs.add(comps.get(i).getLocation().expand(0));
                            }
                        }
                        for (ArrayList<Region> subr : subRegions){
                            ArrayList<BindingComponent> bs = new ArrayList<BindingComponent>();
                            if (subr.size()==1){
                                BindingComponent b;
                                // Scan Peak unary region
                                if(config.use_scanPeak) {
                                    Point p = subr.get(0).getMidpoint();
                                    b = mixture.scanPeak(p);
                                    if (b==null){
                                        continue;
                                    }
                                    b.setEMPosition(p);
                                    for (BindingComponent m:comps){
                                        if (p.getLocation()==m.getLocation().getLocation()) {
                                            b.setAlpha(m.getAlpha());	// inherit alpha value from previous EM component
                                            if(!config.use_betaEM) {
                                                for(int c = 0; c < mixture.numConditions; c++) { 
                                                    b.setConditionBeta(c, m.getConditionBeta(c)); 
                                                }
                                            }
                                            break;   // Once you found a previous component on the same location as b, there is no need to search all the others
                                        }
                                    }
                                } else {
                                    // Do not want to scan peak
                                    //		Run EM again on the unary region and take the component with the maximum strength
                                    ArrayList<BindingComponent> bl = analyzeWindow(subr.get(0).expand(mixture.modelRange, mixture.modelRange));
                                    if(bl == null || bl.size() == 0) { 
                                        continue; 
                                    } else if(bl.size() == 1) { 
                                        b = bl.get(0); 
                                    } else {
                                        double[] compStrength = new double[bl.size()];
                                        for(int k = 0; k < compStrength.length; k++) { compStrength[k] = bl.get(k).getMixProb(); }
                                        Pair<Double, TreeSet<Integer>> max_maxCompIdx = StatUtil.findMax(compStrength);
                                        b = bl.get(max_maxCompIdx.cdr().first());
                                    }
                                }
                                bs.add(b);
                                singleEventRegions.put(subr.get(0), b.getLocation());
                            } else {	// for joint events
                                rs=mixture.mergeRegions(subr, true);
                                for (Region r:rs){
                                    for (BindingComponent m:comps){
                                        if (r.overlaps(m.getLocation().expand(0)))
                                            bs.add(m);
                                    }
                                }
                            }
                            Collections.sort(bs);
                            if (bs.size()>=2){
                                ArrayList<BindingComponent> toRemove = new ArrayList<BindingComponent>();
                                for (int i=1;i<bs.size();i++){
                                    int distance = bs.get(i).getLocation().distance(bs.get(i-1).getLocation());
                                    if (distance==0){
                                        if (bs.get(i).getTotalSumResponsibility()>bs.get(i-1).getTotalSumResponsibility())
                                            toRemove.add(bs.get(i-1));
                                        else
                                            toRemove.add(bs.get(i));
                                    }
                                }
                                bs.removeAll(toRemove);
                            }
                            compFeatures.addAll(mixture.callFeatures(bs));
                        }
                    }
                    processedRegionCount.add(1); 
                } catch(Exception e){
                    System.err.println("ERROR: Java Exception when analyzing region "+rr.toString());
                    e.printStackTrace(System.err);
                    System.exit(-1);
                }
            }
        }

        private ArrayList<BindingComponent> analyzeWindow(Region w){

            ArrayList<List<StrandedBase>> signals = mixture.loadData_checkEnrichment(w);
            if (signals==null)
                return null;

            // We want to run EM only for potential overlapping regions
            // If after first round, we are sure the region contains unary event, we will just scan for peak
            if (!singleEventRegions.containsKey(w)) {
                // dynamically determine an alpha value for this sliding window
                double alpha = config.sparseness;
                if (config.use_dynamic_sparseness)
                    alpha = Math.max(mixture.estimateAlpha(w, signals), config.sparseness);

                Pair<double[][][], int[][][]> result = null;

                // Choose either Yuchun's beta EM or Giorgos's pool EM methods
                if(config.use_betaEM) {
                    //Run EM and increase resolution
                    initializeComponents(w, mixture.numConditions);
                    int lastResolution;

                    int countIters = 0;
                    while(nonZeroComponentNum>0){
                        lastResolution = componentSpacing;
                        //					int numAllComps = components.size();
                        // EM learning, components list will only contains non-zero components
                        result = EMTrain(signals, alpha);
                        // mixture.log(4, componentSpacing+" bp\t"+(int)nonZeroComponents+" components.");

                        countIters++;
                        //					System.out.printf("%s\tupdateIter\t%d\tnumComps\t%d\tnumNonZeroComps\t%d%n", w.toString(), countIters, numAllComps, components.size());
                        // increase resolution
                        if (componentSpacing==1)
                            break;
                        updateComponentResolution(w, mixture.numConditions, lastResolution);
                        if(componentSpacing==lastResolution)
                            break;
                    } 	// end of while (resolution)
                    if (nonZeroComponentNum==0)	return null;
                    setComponentResponsibilities(signals, result.car(), result.cdr());
                } else {
                    // Run MultiIndependentMixture (Temporal Coupling)

                    // Convert data in current region in primitive format
                    int[][] pos     = new int[mixture.numConditions][];  // position is relative to region's start
                    int[][] count   = new int[mixture.numConditions][];
                    char[][] strand = new char[mixture.numConditions][];
                    for(int t = 0; t < mixture.numConditions; t++) {
                        pos[t]    = new int[signals.get(t).size()];
                        count[t]  = new int[signals.get(t).size()];
                        strand[t] = new char[signals.get(t).size()];

                        for(int k = 0; k < signals.get(t).size(); k++) {
                            pos[t][k]    = signals.get(t).get(k).getCoordinate() - w.getStart();  // position is relative to region's start
                            count[t][k]  = (int) signals.get(t).get(k).getCount();
                            strand[t][k] = signals.get(t).get(k).getStrand();
                        }
                    }

                    // Create and just initialize the list to pass the while loop for the first time
                    List<BindingComponent> nonZeroComponents = new ArrayList<BindingComponent>();
                    nonZeroComponents.add(new BindingComponent(mixture.model, new Point(w.getGenome(), w.getChrom(), -1), mixture.numConditions));
                    List<BindingComponent> prevNonZeroComponents = new ArrayList<BindingComponent>(nonZeroComponents);
                    componentSpacing = initSpacing(w);
                    int[] compPos = mixture.setComponentPositions(w.getWidth(), componentSpacing);
                    int countUpdates = 0;
                    boolean hasConverged = false;
                    boolean compsAreTheSame = false;
                    double prev_loglik = Double.NEGATIVE_INFINITY;
                    MultiIndependentMixtureCounts prev_mim = new MultiIndependentMixtureCounts();
                    MultiIndependentMixtureCounts mim = new MultiIndependentMixtureCounts();
                    while( nonZeroComponents.size() > 0 && countUpdates < constants.maxUpdateIters ) {
                        // Set up params for TC to be consistent with those defined in BingingMixture
                        mim = new MultiIndependentMixtureCounts();
                        setMultiIndepMixtureParams(mim, alpha);
                        mim.exec(pos, count, strand, compPos, mixture.model);

                        nonZeroComponents.clear(); nonZeroComponents = null;
                        nonZeroComponents = mixture.determineNonZeroComps(w, mim, compPos, alpha);
                    
                        if(nonZeroComponents.size() == prevNonZeroComponents.size()) {
                            hasConverged = mim.convergence_check(mim.get_glob_loglik(), prev_loglik,  1e-5, true, -1e-6);
                            compsAreTheSame = mixture.checkCompExactMatching(nonZeroComponents, prevNonZeroComponents);
                        }

                        if(hasConverged || compsAreTheSame)
                            break;

                        countUpdates++;
                        prevNonZeroComponents = new ArrayList<BindingComponent>(nonZeroComponents);
                        prev_loglik = mim.get_glob_loglik();
                        prev_mim = mim;
                        compPos = updateComps(w, nonZeroComponents);
                    }

                    HashMap<Integer, double[][]> responsibilities = (HashMap<Integer, double[][]>) prev_mim.get_resp();
                    components.clear(); components = null;
                    components = (ArrayList<BindingComponent>) prevNonZeroComponents;
                    nonZeroComponentNum  = components.size();
                    if (nonZeroComponentNum==0)	return null;
                    setComponentResponsibilities(signals, responsibilities);
                }// end of else condition (running it with TC)

            } else{// if single event region, just scan it
                BindingComponent peak = mixture.scanPeak(singleEventRegions.get(w));
                components = new ArrayList<BindingComponent>();
                components.add(peak);
            }		
            return components;
        }//end of analyzeWindow method

        private int initSpacing(Region w) {
            if (w.getWidth()/50 > constants.MAX_NUM_COMPONENTS)
                System.err.println("Very large region, " + w.toString());
            int spacing = 1;
            componentMax = Math.max(constants.OPTIMAL_NUM_COMPONENTS, w.getWidth()/50);
            if(componentMax > 0 && constants.SMART_SPACING && w.getWidth() > componentMax)
                spacing = Math.max(2, (int)Math.round(1.0*w.getWidth()/componentMax));
            spacing = Math.max(w.getWidth()/constants.MAX_NUM_COMPONENTS, Math.min(spacing, constants.INIT_SPACING));
            return spacing;
        }//end of initSpacing method

        /**
         * Initializes the components. Spaces them evenly along the region for now.
         *
         * @param currReg
         * @param signals
         * @param weighted (uniform or weighted initialization?)
         */
        private void initializeComponents(Region currReg,int numCond){
            components = new ArrayList<BindingComponent>();

            //decide how many components
            componentMax = Math.max(constants.OPTIMAL_NUM_COMPONENTS, currReg.getWidth()/50);
            if(componentMax>0){
                if (constants.SMART_SPACING){
                    int spacing = componentMax>=currReg.getWidth() ? 1 : Math.max(2, (currReg.getWidth()/componentMax)+1);
                    componentSpacing = Math.max(currReg.getWidth()/constants.MAX_NUM_COMPONENTS, Math.min(spacing, constants.INIT_SPACING));
                } else {
                    componentSpacing = Math.max(currReg.getWidth()/constants.MAX_NUM_COMPONENTS, constants.INIT_SPACING);
                }
                int numComponents = currReg.getWidth()/componentSpacing;

                //Set up the components
                double totalMP =0;
                for(int i=0; i<numComponents; i++){
                    Point pos = new Point(mixture.gen, currReg.getChrom(), currReg.getStart()+(i*componentSpacing));
                    BindingComponent currComp = new BindingComponent(mixture.model, pos, numCond);
                    currComp.setMixProb(1);
                    totalMP+=currComp.getMixProb();
                    components.add(currComp);
                }
                //Normalize mixProbs
                for(BindingComponent b : components)
                    b.setMixProb(b.getMixProb()/totalMP);
            }
            nonZeroComponentNum=components.size();
        }//end of initializeComponents method

        /**
         * EM training
         *
         * Returns the responsibility of EM training
         *
         * purely matrix/array operations
         * After EM training, components list will only contains non-zero components
         */
        private Pair<double[][][], int[][][]>  EMTrain(ArrayList<List<StrandedBase>> signals, double alpha){
            int numComp = components.size();

            // H function and responsibility will be stored using an indirect indexing method
            // Because only the components within modelRange matters, we only store those components around the reads
            // and use a mapping array to keep track of the component index
            // This will reduce the memory requirement from N*numComp, to NInRange*numComp
            // and should save some running time because we only iterate the effective base-comps
            double[][]   counts= new double[mixture.numConditions][];	// Hit Count
            double[][][] h= new double[mixture.numConditions][][]; 		// H function
            double[][][] r= new double[mixture.numConditions][][];		// Responsibility
            int[][][] b2c= new int[mixture.numConditions][][]; 			// mapping from base to component
            int[][][] c2b= new int[mixture.numConditions][][]; 			// mapping from component to base
            double[][] b= new double[mixture.numConditions][];			// Beta
            double[] pi = new double[numComp];					// Pi

            // init pi
            for(int j=0;j<numComp;j++){
                BindingComponent comp = components.get(j);
                pi[j]= comp.getMixProb();
            }

            boolean no_data_bin = false;
            for(int c=0; c<mixture.numConditions; c++){
                List<StrandedBase> bases_old = signals.get(c);
                List<StrandedBase> bases = new ArrayList<StrandedBase>();
                if (componentSpacing==1 || no_data_bin) {
                    bases = bases_old;
                } else {
                    if (!bases_old.isEmpty()){
                        char strand = '+';
                        int pos = bases_old.get(0).getCoordinate()+componentSpacing/2;
                        float count = 0;
                        for (StrandedBase bb:bases_old){
                            if (bb.getStrand()!=strand){
                                if (count!=0)
                                    bases.add(new StrandedBase(strand, pos, count));
                                strand = '-';
                                pos = bb.getCoordinate()+componentSpacing/2;
                                count = 0;
                            }
                            if( bb.getCoordinate()>=pos-componentSpacing/2 &&
                                bb.getCoordinate()<=pos+componentSpacing-componentSpacing/2-1){
                                count += bb.getCount();
                            }else{
                                if (count!=0)
                                    bases.add(new StrandedBase(strand, pos, count));
                                count=bb.getCount();
                                pos = bb.getCoordinate()+componentSpacing/2;
                            }
                        }
                        if (count!=0)
                            bases.add(new StrandedBase(strand, pos, count));
                    }
                }
				
                int numBases = bases.size();
                //			System.out.println(StrandedBase.countBaseHits(bases));
                //			System.out.println(StrandedBase.countBaseHits(bases_old));
                double[] bc= new double[numComp];
                for(int j=0;j<numComp;j++)
                    bc[j]=1.0/mixture.numConditions;
                b[c] = bc;

                double[] countc= new double[numBases];
                for(int i=0;i<numBases;i++)
                    countc[i]=bases.get(i).getCount();
                counts[c]=countc;
			
                int[][] c2b_c = new int[numComp][];
                double[][] hc= new double[numComp][];
                double [] prob_comp = new double[numBases];
                for(int j=0;j<numComp;j++){
                    BindingComponent comp = components.get(j);
                    ArrayList<Integer> nzBases = new ArrayList<Integer>();
                    for(int i=0;i<numBases;i++){
                        StrandedBase base = bases.get(i);
                        int dist = base.getStrand()=='+' ? base.getCoordinate()-comp.getLocation().getLocation(): comp.getLocation().getLocation()-base.getCoordinate();
                        prob_comp[i] = mixture.model.probability(dist);
                        if (prob_comp[i]>1e-10){
                            nzBases.add(i);
                        }
                    }
                    c2b_c[j] = new int[nzBases.size()];
                    hc[j] = new double[nzBases.size()];				
                    for(int i=0;i<nzBases.size();i++){
                        c2b_c[j][i] = nzBases.get(i);
                        hc[j][i] = prob_comp[nzBases.get(i)];
                    }
                }
                h[c] = hc;
                c2b[c]=c2b_c;

                int[][] b2c_c = new int[numBases][];
                for(int i=0;i<numBases;i++){
                    ArrayList<Integer> nearComps = new ArrayList<Integer>();
                    for(int j=0;j<numComp;j++){
                        for( int ii: c2b_c[j])
                            if (ii==i){
                                nearComps.add(j);
                                break;
                            }
                    }
                    b2c_c[i] = new int[nearComps.size()];
                    for (int j=0;j<nearComps.size();j++)
                        b2c_c[i][j] = nearComps.get(j);
                }
                b2c[c] = b2c_c;
			
                //Initial Semi-E-step: initialize unnormalized responsibilities
                double[][] rc= new double[numComp][];
                for(int j=0;j<numComp;j++){
                    int[] baseIdx = c2b_c[j];
                    rc[j] = new double[baseIdx.length];
                    for(int i=0;i<baseIdx.length;i++){
                        rc[j][i] = hc[j][i]*pi[j]*bc[j];
                    }
                }
                r[c] = rc;
            }
            mixture.log(5, "\n"+componentSpacing+" bp:\t");
            //////////
            // Run EM steps
            //////////
            Pair<double[][][], int[][][]> result = EM_MAP(counts, h, r, b, b2c, c2b, pi, alpha);
            r = result.car();
            c2b = result.cdr();
		
            //////////
            // re-assign EM result back to the component objects
            //////////
            if (nonZeroComponentNum==0){
                components.clear();
                return new Pair<double[][][], int[][][]>(r, c2b);
            }
            ArrayList<BindingComponent> nonZeroComponents = new ArrayList<BindingComponent>();
            for(int j=0;j<numComp;j++){
                if (pi[j]>0){		// non-zero components
                    BindingComponent comp = components.get(j);
                    comp.setMixProb(pi[j]);
                    comp.setAlpha(alpha);
                    for(int c=0; c<mixture.numConditions; c++){
                        double[] bc = b[c];
                        comp.setConditionBeta(c, bc[j]);
                        comp.setOld_index(j);
                    }
                    nonZeroComponents.add(comp);
                }
            }
            components = nonZeroComponents;

            if (componentSpacing==1){
                for(int c=0; c<mixture.numConditions; c++){
                    for(int j=0;j<components.size();j++){	
                        double sum_resp = 0.0;	
                        int oldIndex = components.get(j).getOld_index();
                        int[] baseIdx = c2b[c][oldIndex];
                        for(int i=0;i<baseIdx.length;i++){
                            sum_resp += counts[c][baseIdx[i]]*r[c][oldIndex][i];
                        }
                        // Assign the summed responsibilities only to non-zero components
                        components.get(j).setCondSumResponsibility(c, sum_resp);
                    }
                }
            }
	
            return new Pair<double[][][], int[][][]>(r, c2b);
        }//end of EMTrain method


        /**
         * core EM steps, sparse prior (component elimination), multi-condition
         */
        private Pair<double[][][], int[][][]> EM_MAP(  	double[][]   counts,
                                                        double[][][] h,
                                                        double[][][] r,
                                                        double[][]   b,
                                                        int[][][] b2c,
                                                        int[][][] c2b,
                                                        double[] pi,
                                                        double alpha) {
            int numComp = pi.length;
            ArrayList<EM_State> models = new  ArrayList<EM_State> ();

            double lastLAP=0, LAP=0; // log posterior prob
            int t=0;
            double currAlpha = alpha/config.gentle_elimination_factor;
            boolean minElimination = false;
            int gentleCounts = 0; 	// count the iterations that we run on gentle mode after last elimination
            // when reach the threshold, increase currAlpha value
            // index of non-zero components, used to iterate components
            ArrayList<Integer> nzComps = new ArrayList<Integer>();	
            for (int j=0;j<pi.length;j++){
                nzComps.add(j);
            }
            mixture.log(5, (int)nonZeroComponentNum+" ");
            //Run EM while not converged
            for(t=0; t<constants.MAX_EM_ITER ; t++){

                lastLAP=LAP;

                //////////
                //Semi-E-step (presumes unnormalized responsibilities have been calculated)
                //Simply normalize responsibilities here
                //////////
                for(int c=0; c<mixture.numConditions; c++){
                    double[][] rc = r[c];
                    int numBases = counts[c].length;
                    double totalResp[] = new double[numBases];
                    // init
                    for(int i=0;i<numBases;i++){
                        totalResp[i] = 0;
                    }
                    // sum
                    for(int j:nzComps){
                        int[] baseIdx = c2b[c][j];
                        for(int i=0;i<baseIdx.length;i++)
                            totalResp[baseIdx[i]] += rc[j][i];
                    }
                    // normalize
                    for(int j:nzComps){
                        int[] baseIdx = c2b[c][j];
                        for(int i=0;i<baseIdx.length;i++)
                            if (totalResp[baseIdx[i]]>0)
                                rc[j][i] = rc[j][i]/totalResp[baseIdx[i]];
                    }
                }


                //////////
                //M-step
                //////////

                //Pi
                // ML: standard EM
                if (t<=config.ML_ITER || currAlpha==0){
                    for(int j:nzComps){
                        double r_sum=0;
                        for(int c=0; c<mixture.numConditions; c++){
                            int[] baseIdx = c2b[c][j];
                            for(int i=0;i<baseIdx.length;i++)
                                r_sum += r[c][j][i]*counts[c][baseIdx[i]];
                        }
                        pi[j]=r_sum;    // standard EM ML
                    }
                }
                // MAP: EM with sparce prior
                else {
                    if (constants.BATCH_ELIMINATION){
                        if (t >config.ML_ITER && t <= constants.ANNEALING_ITER)
                            currAlpha = alpha * (t-config.ML_ITER)/(constants.ANNEALING_ITER-config.ML_ITER);
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
                        for(int j:nzComps){
                            double r_sum=0;
                            for(int c=0; c<mixture.numConditions; c++){
                                int[] baseIdx = c2b[c][j];
                                for(int i=0;i<baseIdx.length;i++)
                                    r_sum += r[c][j][i]*counts[c][baseIdx[i]];
                            }

                            // component elimination
                            pi[j]=Math.max(0,r_sum-currAlpha);

                            // if component prob becomes 0, clear responsibility
                            if (pi[j]==0){
                                for(int c=0; c<mixture.numConditions; c++){
                                    for(int i=0;i<c2b[c][j].length;i++)
                                        r[c][j][i] = 0;
                                }
                            }
	                	
                        }
                    }else{ 	//batch_elimination is false
                        // eliminate only the worst cases
                        // redistribute to boost neighbor component by next round of EM
                        // in this case, we do not need annealing schedule for alpha
            		
                        double r_sum[]=new double[nzComps.size()];
                        for(int jnz=0;jnz<r_sum.length;jnz++){
                            for(int c=0; c<mixture.numConditions; c++){
                                int[] baseIdx = c2b[c][nzComps.get(jnz)];
                                for(int i=0;i<baseIdx.length;i++)
                                    r_sum[jnz] += r[c][nzComps.get(jnz)][i]*counts[c][baseIdx[i]];
                            }
                        }
                        if (gentleCounts>=config.gentle_elimination_iterations &&
                            currAlpha<alpha )
                            currAlpha=Math.min(alpha, currAlpha*2);
            		
                        // find the worst components
                        Pair<Double, TreeSet<Integer>> worst=null;	// worst components
                        if ( (!minElimination) && (componentSpacing!=1))
                            worst = mixture.findSmallestCases(r_sum, currAlpha);
                        else
                            worst = StatUtil.findMin(r_sum);
                        if (worst.car() > currAlpha){
                            // no component to be eliminated, update pi(j)
                            for(int jnz=0;jnz<r_sum.length;jnz++)
                                pi[nzComps.get(jnz)]=r_sum[jnz]-currAlpha;	// not normailzed yet
                            // stop Smallest cases elimination, only eliminate min from now on
                            minElimination = true;
                            gentleCounts ++;
                            //                     	System.out.println(t+":\t"+currAlpha+"\t iterating");
                        }else{
                            // eliminate worst case components, could be 1 or multiple components
                            // redistribute responsibilities in next E step
                            for (int jnz: worst.cdr()){
                                pi[nzComps.get(jnz)]=0;                        	
                                // clear responsibility
                                for(int c=0; c<mixture.numConditions; c++){
                                    for(int i=0; i<c2b[c][nzComps.get(jnz)].length;i++)
                                        r[c][nzComps.get(jnz)][i] = 0;
                                }
                            }                    	
                            // keep iterating on this Alpha value, until converge, then we raise it up to eliminate next one
                            // give EM a time to stabilize before eliminating the next components
                            //                    	System.out.println(t+":\t"+currAlpha+"\t elimination");
                            currAlpha = Math.max(worst.car(), alpha/config.gentle_elimination_factor/2);
                            gentleCounts = 0;
                            //                    	currAlpha = 0;
                        }
                    }
                }

                // update component count, normalize pi
                double totalPi=0;
                nzComps.clear();
                for(int j=0;j<pi.length;j++){
                    if (pi[j]!=0){
                        totalPi+=pi[j];
                        nzComps.add(j);
                    }
                }
                nonZeroComponentNum = nzComps.size();
                if (nonZeroComponentNum==0)
                    return new Pair<double[][][], int[][][]>(r, c2b);
            	
                if (totalPi!=0){
                    for(int j:nzComps){
                        pi[j]=pi[j]/totalPi;
                    }
                }

                //Beta parameters
                if(mixture.numConditions>1){
                    for(int j:nzComps){
                        double b_sum=0;
                        for(int c=0; c<mixture.numConditions; c++){
                            double sum_i=0;		// sum over i
                            int[] baseIdx = c2b[c][j];
                            for(int i=0;i<baseIdx.length;i++)
                                sum_i += r[c][j][i]*counts[c][baseIdx[i]];	

                            b[c][j]=sum_i;
                            b_sum+=sum_i;
                        }

                        // normalize across conditions
                        if (b_sum!=0){
                            for(int c=0; c<mixture.numConditions; c++){
                                b[c][j] /= b_sum;
                            }
                        }
                    } // for
                    // Beta clustering would go here.
                }

                //Semi-E-step:Calculate next un-normalized responsibilities
                for(int j:nzComps){
                    for(int c=0; c<mixture.numConditions; c++){
                        int[] baseIdx = c2b[c][j];
                        for(int i=0;i<baseIdx.length;i++)
                            r[c][j][i] = pi[j]*b[c][j]*h[c][j][i];
                    }
                }

                //Log-likelihood calculation
                double LL =0;
                // get the background probability if a read is out of range of the event
                double bgProb = Math.min(mixture.model.probability(mixture.model.getMax()),
                                         mixture.model.probability(mixture.model.getMin()));
                for(int c=0; c<mixture.numConditions; c++){
                    for(int i=0;i<counts[c].length;i++){
                        // for each read, each event will give a conditional prob or bg prob
                        double j_sum=0;
                        for(int j:nzComps){
                            boolean found = false;
                            for (int ii=0;ii<c2b[c][j].length;ii++)
                                if (i==c2b[c][j][ii] ){
                                    if (r[c][j][ii]!=0){
                                        j_sum += r[c][j][ii];
                                        found = true;
                                        break;
                                    }
                                }
                            if (!found)
                                j_sum += pi[j]*b[c][j]*bgProb;
                        }
                        if (j_sum!=0)
                            LL += Math.log(j_sum)*counts[c][i];
                    }
                }
                // log prior
                double LP=0;
                for(int j:nzComps)
                    if (pi[j]!=0)
                        LP+=Math.log(pi[j]);

                LP= -currAlpha*LP;
                LAP = LL+LP;

                //			System.out.println("EM: "+t+"\t"+currAlpha+"\t"+LAP+"\t"+lastLAP+"\t("+nzComps.size()+" non-zero components).");
                if (t<2 || Math.abs(LAP-lastLAP)>constants.EM_CONVERGENCE){
                    continue;
                }
                else{
                    //				System.out.println("Converged\t"+currAlpha);
                    // if converge on smallest alpha, raise one level
                    if (currAlpha<alpha/config.gentle_elimination_factor){
                        currAlpha = alpha/config.gentle_elimination_factor;
                        continue;
                    }
                    // if converge on a smaller alpha value
                    if (currAlpha<alpha){
                        currAlpha = alpha;		// raise alpha, it may come down again after eliminating the next comp
                        //					System.out.println(t+"::\t"+currAlpha);
                        continue;
                    }
                    // else: converge on full alpha value
                    if (componentSpacing==1 && config.do_model_selection && (!constants.BATCH_ELIMINATION)){
                        // BIC model selection to decide which solution is better
                        // 1. save the state
                        EM_State state = new EM_State(this.mixture.numConditions);
                        state.LL = LL;
                        state.numComponent = nonZeroComponentNum;
                        for (int c=0;c<mixture.numConditions;c++){
                            double[][]rc = r[c];
                            double[][]copy = new double[rc.length][];
                            for (int j=0;j<rc.length;j++){
                                copy[j] = rc[j].clone();
                            }
                            state.resp[c]= copy;
                        }
                        for (int c=0;c<mixture.numConditions;c++){
                            state.beta[c]=b[c].clone();
                        }
                        state.pi = pi.clone();
                        state.alpha = currAlpha;

                        models.add(state);
                        //2. more than one component left, raise alpha, continue EM
                        if (nonZeroComponentNum>1){
                            currAlpha *= 2;
                            continue;
                        } else	// single component left, done
                            break;
                    } else // if at 5bp resolution step, done with EM
                        break;
                }
            } //LOOP: Run EM while not converged

            //		System.out.println("numIters\t" + t + "\t MLIters\t" + config.ML_ITER + "\t annealIters\t" + constants.ANNEALING_ITER + "\tmaxIters\t" + constants.MAX_EM_ITER);

            if (componentSpacing==1 && config.do_model_selection && models.size()>1 && (!constants.BATCH_ELIMINATION)){
                // BIC model selection to decide which solution is better
                // 3. find best solution
                double totalBaseCount = 0;
                for(int c=0; c<mixture.numConditions; c++){
                    totalBaseCount += counts[c].length;
                }

                int best = 0;
                double bestBIC = models.get(0).BIC(totalBaseCount);
                for (int i=1;i<models.size();i++){
                    double bic = models.get(i).BIC(totalBaseCount);
                    //				System.out.println(String.format("%.3f\t%.3f\t", bestBIC, bic)+models.get(i).toString());
                    if (bestBIC <= bic){
                        best = i;
                        bestBIC = bic;
                    }
                }
                // copy the best model state back to memory
                EM_State bestModel = models.get(best);
                r=new double[mixture.numConditions][][];
                for (int c=0;c<mixture.numConditions; c++){
                    r[c]=bestModel.resp[c];
                }
                b=new double[mixture.numConditions][];
                for (int c=0;c<mixture.numConditions; c++){
                    b[c] = bestModel.beta[c];
                }
                nzComps.clear();
                for (int j=0;j<pi.length;j++){
                    pi[j] = bestModel.pi[j];
                    if (pi[j]!=0)
                        nzComps.add(j);
                    //				System.out.print(String.format("%.2f ", pi[j]));
                }
                nonZeroComponentNum = nzComps.size();
			
			
                //			System.out.println();
            }
            // normalize responsibilities here
            for(int c=0; c<mixture.numConditions; c++){
                double[][] rc = r[c];
                int numBases = counts[c].length;
                double totalResp[] = new double[numBases];
                // init
                for(int i=0;i<numBases;i++){
                    totalResp[i] = 0;
                }
                // sum
                for(int j:nzComps){
                    int[] baseIdx = c2b[c][j];
                    for(int i=0;i<baseIdx.length;i++)
                        totalResp[baseIdx[i]] += rc[j][i];
                }
                // normalize
                for(int j:nzComps){
                    int[] baseIdx = c2b[c][j];
                    for(int i=0;i<baseIdx.length;i++)
                        if (totalResp[baseIdx[i]]>0)
                            rc[j][i] = rc[j][i]/totalResp[baseIdx[i]];
                }
            }

            // make hard assignment
            //		if (componentSpacing==1 && constants.MAKE_HARD_ASSIGNMENT){


            mixture.log(4, "EM_MAP(): "+"\tt="+t+"\t"+
                        String.format("%.6f",LAP)+"\t("+(int)nonZeroComponentNum+" events)");
            return new Pair<double[][][], int[][][]>(r, c2b);

        }//end of EM_MAP method

        //Update the resolution of the components
        //Add new components in around non-zero components
        private void updateComponentResolution(Region currReg, int numCond, int lastResolution){
            ArrayList<BindingComponent> newComponents = new ArrayList<BindingComponent>();

            //First map the valid bases around non-zero components
            int numValidBases=0;
            int[] valid = new int[currReg.getWidth()];
            //		for(int v=0; v<currReg.getWidth(); v++){valid[v]=0;}
            for(BindingComponent b : components){
                for(int v=(int)Math.max(0, (b.getLocation().getLocation()-lastResolution*config.resolution_extend)-currReg.getStart());
                    v<=Math.min(currReg.getWidth()-1, (b.getLocation().getLocation()+lastResolution*config.resolution_extend)-currReg.getStart());
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
            if(componentSpacing>=lastResolution)
                componentSpacing = lastResolution -1;
		
            //Set up the components
            double totalMP =0;
            for(BindingComponent b: components){
                for(int v=(int)Math.max(0, (b.getLocation().getLocation()-lastResolution*config.resolution_extend+1)-currReg.getStart());
                    v<=Math.min(currReg.getWidth()-1, (b.getLocation().getLocation()+lastResolution*config.resolution_extend-1)-currReg.getStart());
                    v+=componentSpacing){
                    Point pos = new Point(mixture.gen, currReg.getChrom(), currReg.getStart()+v);
                    //Reusing the valid array for the sake of it
                    if(valid[v]!=-1){
                        BindingComponent currComp = new BindingComponent(mixture.model, pos, numCond);
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
        }//end of updateComponentResolution method

        // This is for beta EM method
        private void setComponentResponsibilities(ArrayList<List<StrandedBase>> signals, 
                                                  double[][][] responsibilities, int[][][] c2b) {
            // Set responsibility profile for each component (for kernel density and KL calculation)
            for(int j=0;j<components.size();j++){
                BindingComponent comp = components.get(j);
                int jr = comp.getOld_index();
                for(int c=0; c<mixture.numConditions; c++){
                    List<StrandedBase> bases = signals.get(c);
                    double[][] rc = responsibilities[c];

                    // store binding profile (read responsibilities in c condition) of this component
                    double[] profile_plus = new double[mixture.modelWidth];
                    double[] profile_minus = new double[mixture.modelWidth];
                    for(int i=0;i<c2b[c][jr].length;i++){
                        int base_idx = c2b[c][jr][i];
                        StrandedBase base = bases.get(base_idx);
                        if (rc[jr][i]>0){
                            try{
                                if (base.getStrand()=='+')
                                    profile_plus[base.getCoordinate()-comp.getLocation().getLocation()-mixture.model.getMin()]=rc[jr][i]*base.getCount();
                                else
                                    profile_minus[comp.getLocation().getLocation()-base.getCoordinate()-mixture.model.getMin()]=rc[jr][i]*base.getCount();
                            }
                            catch (Exception e){
                            }
                        }
                    }
                    comp.setReadProfile(c, profile_plus,  '+');
                    comp.setReadProfile(c, profile_minus, '-');
                }
            }
        }//end of setComponentResponsibilities method
	
        // This is for pool EM method
        private void setComponentResponsibilities(ArrayList<List<StrandedBase>> signals, HashMap<Integer, double[][]> responsibilities) {
            // Set responsibility profile for each component (for kernel density and KL calculation)
            for(int j=0;j<components.size();j++){
                BindingComponent comp = components.get(j);
                int jr = comp.getOld_index();
                for(int c=0; c<mixture.numConditions; c++){
                    List<StrandedBase> bases = signals.get(c);
                    double[][] rc = responsibilities.get(c);

                    // store binding profile (read responsibilities in c condition) of this component
                    double[] profile_plus = new double[mixture.modelWidth];
                    double[] profile_minus = new double[mixture.modelWidth];
                    for(int i=0;i<bases.size();i++){
                        StrandedBase base = bases.get(i);
                        if (rc[i][jr]>0){
                            try{
                                if (base.getStrand()=='+')
                                    profile_plus[base.getCoordinate()-comp.getLocation().getLocation()-mixture.model.getMin()]=rc[i][jr]*base.getCount();
                                else
                                    profile_minus[comp.getLocation().getLocation()-base.getCoordinate()-mixture.model.getMin()]=rc[i][jr]*base.getCount();
                            }
                            catch (Exception e){
                                System.err.println(comp.toString()+"\t"+base.getStrand()+"\t"+base.getCoordinate());
                                e.printStackTrace(System.err);
                            }
                        }
                    }
                    comp.setReadProfile(c, profile_plus,  '+');
                    comp.setReadProfile(c, profile_minus, '-');
                }
            }
        }//end of setComponentResponsibilities method

        private int[] updateComps(Region w, List<BindingComponent> comps) {
            int[] newComps;
            int spacing;
            int window = (int)Math.round(1.0*mixture.modelRange/(constants.maxUpdateIters-1));
            Set<Integer> newCompsSet = new TreeSet<Integer>();
            List<int[]> compGroupBounds = new ArrayList<int[]>();

            BindingComponent[] comps_arr = comps.toArray(new BindingComponent[0]);
            int[] compAssigns = new int[comps_arr.length];
            for(int j = 1; j < comps_arr.length; j++) {
                if(j != comps_arr.length-1) {
                    if(Math.abs(comps_arr[j].getLocation().getLocation() - comps_arr[j-1].getLocation().getLocation()) <= Math.abs(comps_arr[j].getLocation().getLocation() - comps_arr[j+1].getLocation().getLocation()))
                        compAssigns[j] = compAssigns[j-1];
                    else
                        compAssigns[j] = compAssigns[j-1]+1;
                }
                if(Math.abs(comps_arr[j].getLocation().getLocation() - comps_arr[j-1].getLocation().getLocation()) > 2*window)
                    compAssigns[j] = compAssigns[j-1]+1;

                if(j == comps_arr.length-1 && Math.abs(comps_arr[j].getLocation().getLocation() - comps_arr[j-1].getLocation().getLocation()) <= 2*window)
                    compAssigns[j] = compAssigns[j-1];
            }

            if(compAssigns.length > 0)
                compGroupBounds.add(new int[]{comps_arr[0].getLocation().getLocation(), comps_arr[0].getLocation().getLocation()});

            for(int j = 1; j < compAssigns.length; j++)
                if(compAssigns[j] != compAssigns[j-1])
                    compGroupBounds.add(new int[]{comps_arr[j].getLocation().getLocation(), comps_arr[j].getLocation().getLocation()});

            for(int j = 1; j < compAssigns.length; j++) {
                int start = compGroupBounds.get(compAssigns[j])[0];
                int end   = comps_arr[j].getLocation().getLocation();
                compGroupBounds.set(compAssigns[j], new int[]{start, end});
            }

            // If there is overlap between the component groups, modify the bounds properly
            for(int k = 1; k < compGroupBounds.size(); k++) {
                int start    = compGroupBounds.get(k)[0];
                int prev_end = compGroupBounds.get(k-1)[1];
                if(Math.abs(start-prev_end) < 2*window) {
                    start    += window;
                    prev_end -= window;
                    int end   = compGroupBounds.get(k)[1];
                    int prev_start = compGroupBounds.get(k-1)[0];
                    start = Math.min(start, end);
                    prev_end = Math.max(prev_start, prev_end);
                    compGroupBounds.set(k, new int[]{start, end});
                    compGroupBounds.set(k-1, new int[]{prev_start, prev_end});
                }
            }

            int numCandidBases = 0;
            for(int[] bounds:compGroupBounds) {
                int start = bounds[0];
                int end   = bounds[1];
                int left_bound  = Math.max(0, start-window-w.getStart());
                int right_bound = Math.min(w.getWidth()-1, end+window-w.getStart());
                numCandidBases += right_bound-left_bound+1;
            }
            spacing = componentMax >= numCandidBases ? 3 : Math.max(3, (int)Math.round(1.0*numCandidBases/componentMax));

            for(int[] bounds:compGroupBounds) {
                int start = bounds[0];
                int end   = bounds[1];
                for(int v =  Math.max(0, start-window-w.getStart());
                    v <= Math.min(w.getWidth()-1, end+window-w.getStart());
                    v += spacing) {
                    newCompsSet.add(v);
                }
                newCompsSet.add(Math.min(w.getWidth()-1, end+window-w.getStart()));
            }

            newComps = Utils.ref2prim(newCompsSet.toArray(new Integer[0]));
            return newComps;
        }//end of updateComps method

        private void setMultiIndepMixtureParams(MultiIndependentMixtureCounts mim, double alpha) {
            mim.set_trainVersion(1);
            mim.set_alpha(alpha);
            mim.set_maxIters(constants.MAX_EM_ITER);
            mim.set_ML_maxIters(config.ML_ITER);
            mim.set_anneal_maxIters(200);
            mim.set_convergence_thres(1e-8);
            mim.set_decreasing_thres(-1e-6);
            mim.set_check_increased(true);
            mim.set_isBatchElimOn(true);
            //		mim.set_isBatchElimOn(constants.BATCH_ELIMINATION);
        }//end of setMultiIndepMixtureParams method

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
            // for multi-condition, # of beta variables is (mixture.numConditions-1)*numComponents
            // n: is the number of data point, i.e. the base positions.
            double BIC(double n){
                return LL - (numComponent*2-1 + (mixture.numConditions-1)*numComponent )/2*Math.log(n);
            }
            public String toString(){
                return String.format("%.3f\t%.0f", LL, numComponent);
            }
        }

    }
}
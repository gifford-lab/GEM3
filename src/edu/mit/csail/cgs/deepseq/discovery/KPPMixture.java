package edu.mit.csail.cgs.deepseq.discovery;

import java.io.*;
import java.util.*;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import javax.imageio.ImageIO;

import cern.jet.random.Poisson;
import cern.jet.random.Binomial;
import cern.jet.random.engine.DRand;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.*;
import edu.mit.csail.cgs.deepseq.discovery.kmer.Kmer;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KmerSet;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KMAC.KmerGroup;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KMAC.MotifThreshold;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KMAC;
import edu.mit.csail.cgs.deepseq.features.*;
import edu.mit.csail.cgs.deepseq.multicond.MultiIndependentMixtureCounts;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.deepseq.utilities.ReadCache;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.Utils;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.models.data.DataFrame;
import edu.mit.csail.cgs.utils.models.data.DataRegression;
import edu.mit.csail.cgs.utils.probability.NormalDistribution;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public class KPPMixture extends MultiConditionFeatureFinder {
    private Config config;
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
	private boolean controlDataExist = false;
	// Ratio of each pair of IP/Ctrl for all conditions
	private boolean hasIpCtrlRatio = false;
	private double[] ratio_total;
	private double[] ratio_non_specific_total;
	private double background_proportion = -1;			//background read proportion in the candidate regions
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
	//Regions to be excluded, supplied by user, appended with un-enriched regions
	private ArrayList<Region> excludedRegions = new ArrayList<Region>();

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

    private StringBuilder configsb = new StringBuilder();
	private StringBuilder log_all_msg = new StringBuilder();
	private FileWriter logFileWriter;
	private File outputFolder = null;

	/** Kmer motif engine */
	private KMAC kmac;
	private boolean kmerPreDefined = false;
	private double kcm_threshold = -1;
	
	public KPPMixture(Genome g, 
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
		
		File outFile = new File(outName);
		outputFolder = outFile.getParentFile();
		File gem_outputs_folder;
		if (outputFolder!=null){
			if (!outputFolder.exists()){
				System.err.println("\nThe output file path is not correct: "+outFile.getAbsolutePath());
				cleanUpDataLoader();
	    		System.exit(-1);
			}
			gem_outputs_folder = new File(outputFolder, outFile.getName()+"_outputs");
		}
		else
			gem_outputs_folder = new File(outFile.getName()+"_outputs");
		gem_outputs_folder.mkdir();
		outName = new File(gem_outputs_folder, outFile.getName()).getAbsolutePath();	// re-direct outName prefix to a folder
		
		model.printToFile(outName+"_0_Read_distribution.txt");
		allModels.put(outName+"_0", model);
		
    	
		/* ***************************************************
		 * Print out command line options
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
    	 * Load options
    	 ***********************************/    
		try{
			config.parseArgs(args);   
		}
		catch (Exception e){
			e.printStackTrace();
			cleanUpDataLoader();
    		System.exit(-1);
		}
     // if mappable_genome_length is not provided, compute as 0.8 of total genome size
        if (config.mappable_genome_length<0){		
	        config.mappable_genome_length = 0.8 * gen.getGenomeSize();
	        System.out.println(String.format("\nMappable Genome Length is %,d.", (long)config.mappable_genome_length));
        }
        
    	if(config.second_lambda_region_width < config.first_lambda_region_width) {
    		System.err.println("\nThe first control region width (w2) has to be more than " + config.first_lambda_region_width + " bp.");
			cleanUpDataLoader();
    		System.exit(-1);
    	}    	
    	if(config.third_lambda_region_width < config.second_lambda_region_width) {
    		System.err.println("\nThe second control region width (w3) has to be more than " + config.second_lambda_region_width + " bp.");
			cleanUpDataLoader();
			System.exit(-1);
    	}
    	

    	/* *********************************
    	 * Load Kmer list
    	 ***********************************/
    	String kmerFile = Args.parseString(args, "kf", null);
    	if (kmerFile!=null){
    		kmerPreDefined = true;
			File kFile = new File(kmerFile);
			if(kFile.isFile()){
				KmerSet kc = new KmerSet(kFile);
	        	kmac = new KMAC(kc.getKmers(0), outName);
	        	kmac.setTotalSeqCounts(kc.posSeqCount, kc.negSeqCount);
	        	kmac.setConfig(config, outName);
	        	kcm_threshold = kc.ksmThreshold;
			}
    	}
    	
    	/* *********************************
    	 * load ChIP-Seq read data
    	 ***********************************/
    	//Determine the subset of genome to run EM
    	ArrayList<Region> subsetRegions = getSubsetRegions(args);
    	if (!subsetRegions.isEmpty())
    		config.cache_genome = false;
     	
    	//load read data
		this.conditionNames = conditionNames;
    	loadChIPSeqData(subsetRegions, args);
		
    	// check if there is miss control data
		if (controlDataExist){
			for (Pair<ReadCache, ReadCache> p:caches){
				if (p.cdr()==null){
					System.err.println("\nMissing control data to match "+p.car().getName());
					cleanUpDataLoader();
					System.exit(-1);
				}
			}
		}
		
		// exclude some regions
     	String excludedName = Args.parseString(args, "ex", "yes");
     	int exludedLength = 0;
     	if (!excludedName.equals("yes")){
     		excludedRegions = mergeRegions(CommonUtils.loadRegionFile(excludedName, gen), false);
     		log(1, "\nExclude " + excludedRegions.size() + " regions.\n");
			for(int c = 0; c < numConditions; c++) {
				caches.get(c).car().excludeRegions(excludedRegions);
				if(controlDataExist) {
					caches.get(c).cdr().excludeRegions(excludedRegions);
				}
			}   
			for (Region r:excludedRegions)
				exludedLength += r.getWidth();
     	}
     	// subtract excluded regions length from mappable genome size
     	config.mappable_genome_length -= exludedLength;
     	
		// print initial dataset counts
		for(int c = 0; c < numConditions; c++) {
			caches.get(c).car().displayStats();
			if(controlDataExist) {
				caches.get(c).cdr().displayStats();
			}
		}
		
		// Filtering/reset bases
		if(config.filterDupReads){
			applyPoissonFilter(false);
			applyPoissonFilter(true);
		}
		
		// Normalize conditions
		normExpts(caches);
		
		// esitimate sparseness minimum value from the whole genome coverage
		if (config.sparseness==-1){
			double totalReadCount=0;
			for(int i=0; i<numConditions; i++)
				totalReadCount += caches.get(i).car().getHitCount();	
			int expectedHitCount = calcExpectedHitCount(totalReadCount, config.poisson_alpha, modelWidth);
			config.sparseness = expectedHitCount;
			log(1, String.format("\nAt Poisson p-value %.1e, in a %dbp window, expect %d reads.\n", 
					config.poisson_alpha, modelWidth, expectedHitCount));
		}
		
		// estimate ratio and segment genome into regions
    	ratio_total=new double[numConditions];
    	ratio_non_specific_total = new double[numConditions];
        sigHitCounts=new double[numConditions];
        seqwin=100;
     	
		if (config.ip_ctrl_ratio>0){		// if ratio is provided
			for(int i=0; i<numConditions; i++){
				ratio_total[i]=config.ip_ctrl_ratio;
				ratio_non_specific_total[i]=config.ip_ctrl_ratio;
			}
		}
		else{
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
	        } else{	// want to analyze only specified regions, set default
	        	for(int i=0; i<numConditions; i++){
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
        
     	String subsetFormat = Args.parseString(args, "subFormat", "");
    	// if not provided region list, directly segment genome into enrichedRegions
		if (wholeGenomeDataLoaded || !subsetFormat.equals("Regions")){
			setRegions(selectEnrichedRegions(subsetRegions, true));
     		// ip/ctrl ratio by regression on non-enriched regions
			if (config.ip_ctrl_ratio==-1){
     			ArrayList<Region> temp = (ArrayList<Region>)restrictRegions.clone();
     			temp.addAll(excludedRegions);
    			calcIpCtrlRatio(mergeRegions(temp, false));
    			if(controlDataExist) {
    				for(int t = 0; t < numConditions; t++)
    					System.out.println(String.format("For condition %s, IP/Control = %.2f", conditionNames.get(t), ratio_non_specific_total[t]));
    			}
    		}
		} else{
			setRegions(subsetRegions);
		}
		
		log(1, "\nThe genome is segmented into "+restrictRegions.size()+" regions for analysis.");

        if (restrictRegions.isEmpty())
        	System.exit(-1);
        
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
	public KPPMixture(String[] args, boolean doScanning){
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
	
	protected void finalize() throws Throwable {
	    try {
	    	cleanUpDataLoader();        
	    } finally {
	        super.finalize();
	    }
	}

	private void cleanUpDataLoader(){
		for (int i=0;i<experiments.size();i++){
			Pair<DeepSeqExpt, DeepSeqExpt> pair = experiments.get(i);
			DeepSeqExpt ip = pair.car();
			DeepSeqExpt ctrl = pair.cdr();
			if (ip!=null)
				ip.closeLoaders();
			if (ctrl!=null)
				ctrl.closeLoaders();
		}
	}
	
	private void commonInit(String modelFile){
        constants = new GPSConstants();
        config = new Config();

		//Load Binding Model
		File pFile = new File(modelFile);
		if(!pFile.isFile()){
			System.err.println("\nCannot find read distribution file!");
			cleanUpDataLoader();
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
		// Compute the total length and read counts of enriched regions
		if (config.background_proportion==-1){
			if (background_proportion == -1){	// first round, do not estimate from candidate regions
				background_proportion = config.pi_bg_r0;
				log(2,String.format("Default initial noise proportion in the candidate regions = %.3f%n", background_proportion));
			}
			else{
				int length=0;
				int selectedReadCount=0;
				int totalReadCount=0;
				for (Region r: restrictRegions){
					length+= r.getWidth();
					selectedReadCount+=this.countIpReads(r);
				}
				for(int i=0; i<numConditions; i++){
					totalReadCount += caches.get(i).car().getHitCount();
				}
				background_proportion = 1.0*(totalReadCount-selectedReadCount)/(config.mappable_genome_length-length)*length/selectedReadCount;
				log(2,String.format("Estimated noise proportion in the candidate regions = %.3f%n", background_proportion));
			}
		}
		else{
			background_proportion = config.background_proportion;
			log(2,String.format("Pre-specified noise proportion in the candidate regions = %.3f%n", background_proportion));
		}
		// refresh the total profile sums every round
		profile_plus_sum = new double[modelWidth];
		profile_minus_sum = new double[modelWidth];

		signalFeatures.clear();
        Vector<ComponentFeature> compFeatures = new Vector<ComponentFeature>();
		Vector<KmerPP> allKmerHits = new Vector<KmerPP>();
		
        // prepare for progress reporting
        Vector<Integer> processRegionCount = new Vector<Integer>();		// for counting how many regions are processed by all threads
		int displayStep = (int) Math.pow(10, (int) (Math.log10(totalRegionCount)));
		TreeSet<Integer> reportTriggers = new TreeSet<Integer>();
		for (int i=1;i<=totalRegionCount/displayStep; i++){
			reportTriggers.add(i*displayStep);
		}
		reportTriggers.add(100);
		reportTriggers.add(1000);
		reportTriggers.add(10000);
		
		if (kmac!=null && kmac.isInitialized()){
			log(1,"Running EM with motif positional prior ...");
		}
		
		// create threads and run EM algorithms
		// the results are put into compFeatures
        Thread[] threads = new Thread[config.maxThreads];
        log(1,String.format("Running with %d threads ...", config.maxThreads));
        Vector<Region> regionsRunning = new Vector<Region>();		// object to pass info of currently running regions
        // regionsToRun is shared by all threads. Each thread will access it exclusively, lock the obj, get first region, remove it, then unlock.
        TreeSet<Region> regionsToRun = new TreeSet<Region>();
        regionsToRun.addAll(restrictRegions);
        for (int i = 0 ; i < threads.length; i++) {
            Thread t = new Thread(new GPS2Thread(regionsToRun,
            									processRegionCount,
            									regionsRunning,
                                                compFeatures,
                                                allKmerHits,
                                                this,
                                                constants,
                                                config, 
                                                true));
            t.start();
            threads[i] = t;
        }
        boolean anyrunning = true;
        int count = 0;
        HashSet<Region> heavyRegions = new HashSet<Region>();
        while (anyrunning) {
            anyrunning = false;
            try {
                Thread.sleep(1000);
            } catch (InterruptedException e) { }
            for (int i = 0; i < threads.length; i++) {
                if (threads[i].isAlive()) {
                    anyrunning = true;
                    if (count == processRegionCount.size()){
                    	try{
	                    	heavyRegions.addAll(regionsRunning);
	                    	if (config.verbose>1 && !regionsRunning.isEmpty())
	                    		System.out.println("Analyzing "+regionsRunning.elementAt(0).toString());
                    	}
                    	catch (ConcurrentModificationException e){
                    		//ignore
                    	}
                    }
                    break;
                }
            }    
            count = processRegionCount.size();
//            System.out.println(count);
            int trigger = totalRegionCount;
            if (!reportTriggers.isEmpty())
            	trigger = reportTriggers.first();
            if (count>trigger){
				System.out.println(trigger+"\t/"+totalRegionCount+"\t"+CommonUtils.timeElapsed(tic));
				reportTriggers.remove(reportTriggers.first());
            }
        }
        System.out.println(totalRegionCount+"\t/"+totalRegionCount+"\t"+CommonUtils.timeElapsed(tic));
        if (compFeatures.isEmpty()){
        	log(1, "No valid binding event was found.");
        	return signalFeatures;
        }
        processRegionCount.clear();
        compFeatures.trimToSize();
        
        log(3,String.format("%d threads have finished running", config.maxThreads));
        
        // print out the heavy regions
        if (!heavyRegions.isEmpty()){
        	ArrayList<Region> rs = new ArrayList<Region>();
        	rs.addAll(heavyRegions);
        	Collections.sort(rs);
        	StringBuilder sb0 = new StringBuilder();
        	for (Region r:rs)	
	        	sb0.append(r.toString()).append("\n");
	        CommonUtils.writeFile(outName+"_zzzzz_heavyRegions.txt", sb0.toString());
        }
        
        // print out kmer hit list
        if (config.kmer_print_hits){
	        Collections.sort(allKmerHits);
	        StringBuilder sb = new StringBuilder();
	        sb.append("Position\tKmer\tCount\tWeight\tpp\n");
	        for (KmerPP h:allKmerHits)
	        	sb.append(h.toString()).append("\n");
	        CommonUtils.writeFile(outName+"_Kmer_Hits.txt", sb.toString());
	        sb = null;
        }
        
		/* ********************************************************
		 * refine the bound regions that contain binding events
		 * ********************************************************/
		if (config.refine_regions){
			Collections.sort(compFeatures);
			ArrayList<Region> refinedRegions = new ArrayList<Region>();
			for (ComponentFeature cf:compFeatures){
				refinedRegions.add(cf.getPosition().expand(0));
			}
			refinedRegions.trimToSize();
			// expand with modelRange, and merge overlapped regions ==> Refined enriched regions
			this.restrictRegions=mergeRegions(refinedRegions, true);
		}

		log(2, "Finish predicting events: "+CommonUtils.timeElapsed(tic)+"\n");

		/* ********************************************************
		 * post EM processing
		 * ********************************************************/
		postEMProcessing(compFeatures);
		
		log(1, "Finish binding event prediction: "+CommonUtils.timeElapsed(tic)+"\n");

		return signalFeatures;
	}// end of execute method

	/**
	 * It performs statistical tests for determining significant peaks         <br>
	 * If a control is used, the current method is used. Otherwise, MACS proposed method is used.
	 * @param compFeatures
	 */
	private void postEMProcessing(List<ComponentFeature> compFeatures) {
		// use the refined regions to count non-specific reads
        /* don't do this any more.  all it does is set ratio_non_specific_total, which
           we no longer want to update because we compute it at the beginning of the run
           using the whole genome rather than just the unenriched regions
        */
//			countNonSpecificReads(compFeatures);
		
		// collect enriched regions to exclude to define non-specific region
		Collections.sort(compFeatures, new Comparator<ComponentFeature>() {
                public int compare(ComponentFeature o1, ComponentFeature o2) {
                    return o1.compareByTotalResponsibility(o2);
                }
            });
		// get median event strength
		double quarterStrength = compFeatures.get(compFeatures.size()/4).getTotalEventStrength();
//		System.out.println(String.format("Median event strength = %.1f\n",medianStrength));
		
		// only do this calculation at the first round, then throw out read data in non-specific regions (when we have control data for stat test)
		if (!hasIpCtrlRatio){		
			ArrayList<Region> exRegions = new ArrayList<Region>();
			for(int i = 0; i < compFeatures.size(); i++) {
				exRegions.add(compFeatures.get(i).getPosition().expand(modelRange));
			}
			for (Region ex:excludedRegions)
				if (ex!=null)
					exRegions.add(ex);	// also excluding the excluded regions (user specified + un-enriched)
			
			calcIpCtrlRatio(mergeRegions(exRegions, false));
			if(controlDataExist) {
				for(int c = 0; c < numConditions; c++)
					System.out.println(String.format("\nScaling condition %s, IP/Control = %.2f", conditionNames.get(c), ratio_non_specific_total[c]));
				System.out.println();
			}
			hasIpCtrlRatio = true;
			
			// delete read data in un-enriched region, we don't need them for binomial test
			// but if no control data, we need whole genome data to estimate lambda for Poisson test
//			if(controlDataExist) {
//				for(int c = 0; c < numConditions; c++) {
//					caches.get(c).car().deleteUnenrichedReadData(restrictRegions);				
//					caches.get(c).cdr().deleteUnenrichedReadData(restrictRegions);
//				}
//			}
		}
		
		/**
		 * Classify the events into IP-like and control-like groups
		 */
		if (config.classify_events)
			falseDiscoveryTest(compFeatures);
		
		/**
		 *  calculate p-values with or without control
		 */
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
				if (config.TF_binding){	// single event IP/ctrl only applies to TF
					if (cf.getQValueLog10(cond)>=config.q_value_threshold){
						significant = true;
						break;
					}
				}
				else		//TODO: if histone data, need to test as a set of events
					significant = true;
			}

			boolean notFiltered = false;
			// only filter high read count events (more likely to be artifacts) 
			if (config.filterEvents && cf.getTotalEventStrength()>quarterStrength){
				for (int cond=0; cond<numConditions; cond++){
					// if one condition is good event, this position is GOOD
					// logKL of event <= 2.5, and IP/control >= 4 --> good (true)
//					if (cf.getShapeDeviation(cond)<=config.shapeDeviation){
//						if (!controlDataExist){
//							notFiltered = true;
//							break;
//						}
//						else{
//							double ratio = cf.getEventReadCounts(cond)/cf.getScaledControlCounts(cond);
//							if ((ratio>=config.fold && cf.getAverageIpCtrlLogKL()>config.kl_ic) || (cf.getAverageIpCtrlLogKL()<config.kl_ic && ratio>=config.fold*2)){
//								notFiltered = true;
//								break;
//							}
//						}						
//					}
					// relax fold change requirement for event with good shape
					if (!controlDataExist){
						if (cf.getShapeDeviation(cond)<=config.shapeDeviation){
							notFiltered = true;
							break;
						}
					}
					else{
						double ratio = cf.getEventReadCounts(cond)/cf.getScaledControlCounts(cond);
						if ((ratio>=config.fold && cf.getShapeDeviation(cond)<=config.shapeDeviation) || 
								(ratio>=2 && cf.getShapeDeviation(cond)<=-0.5)){
							notFiltered = true;
							break;
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
		((ArrayList<Feature>)signalFeatures).trimToSize();
		filteredFeatures.trimToSize();
		insignificantFeatures.trimToSize();
		
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
		
		// post filtering for multi-conditions
		for(int c = 0; c < numConditions; c++)
			condSignalFeats[c] = condPostFiltering(signalFeatures, c);

		log(1, "---------------------------\n"+
			"Events discovered \nSignificant:\t"+signalFeatures.size()+
            "\nInsignificant:\t"+insignificantFeatures.size()+
            "\nFiltered:\t"+filteredFeatures.size()+"\n");
	}//end of post EM Processing

	/**
	 *  evaluate significance of each called events, calculate p-value from binomial distribution
	 */
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
	                if (ipCount==0){			// if one of the condition does not have reads, set p-value=1
	                	cf.setPValue_w_ctrl(1, cond);
	                	continue;
	                }
	                if (config.testPValues){
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
		
		                    poisson.setMean(config.minFoldChange * Math.max(scaledControlCount, totalIPCount[cond] * modelWidth / config.mappable_genome_length  ));
		                    pValuePoisson = 1 - poisson.cdf(ipCount) + poisson.pdf(ipCount);
		                } catch(Exception err){
		                    err.printStackTrace();
		                    System.err.println(cf.toString());
		                    throw new RuntimeException(err.toString(), err);
		                }

	                	cf.setPValue_w_ctrl(Math.max(Math.max(pValuePoisson,pValueBalance),Math.max(pValueControl,pValueUniform)), cond);
	                }
	                else
	                	cf.setPValue_w_ctrl(StatUtil.binomialPValue(scaledControlCount, scaledControlCount+ipCount), cond);
				}
			}
		} 
		/** compute Poisson p-value from IP only (similar to MACS) */
		if( (!controlDataExist) || (controlDataExist && config.strigent_event_pvalue)) {
			Collections.sort(compFeatures);				// sort by location
			
			createChromStats(compFeatures);
			
			Map<String, ArrayList<Integer>> chrom_comp_pair = new HashMap<String, ArrayList<Integer>>();
			for(int i = 0; i < compFeatures.size(); i++) {
				String chrom = compFeatures.get(i).getPosition().getChrom();
				if(!chrom_comp_pair.containsKey(chrom))
					chrom_comp_pair.put(chrom, new ArrayList<Integer>());			
				chrom_comp_pair.get(chrom).add(i);
			}
		
			for(String chrom:chrom_comp_pair.keySet()) {
				int chromLen = gen.getChromLength(chrom);
				for(int c = 0; c < numConditions; c++) {
					double chrom_lambda = ((double)condHitCounts.get(chrom).get(0).get(c))/chromLen*(modelRange*2+1);				
					for(int i:chrom_comp_pair.get(chrom)) {
						double thirdLambda = computeLambda(compFeatures, i, c, config.third_lambda_region_width);
						double secondLambda = computeLambda(compFeatures, i, c, config.second_lambda_region_width);					
						double local_lambda = Math.max(secondLambda,  Math.max(thirdLambda, chrom_lambda));
						ComponentFeature cf = compFeatures.get(i); 
						cf.setExpectedCounts(local_lambda, c);                        
						poisson.setMean(local_lambda);
	                    int count = (int)Math.ceil(cf.getEventReadCounts(c));
	                    double pValue = 1 - poisson.cdf(count) + poisson.pdf(count);
						cf.setPValue_wo_ctrl(pValue, c);
					}//end of for(int i:chrom_comp_pair.get(chrom)) LOOP	
				}				
			}
		}
		
		// calculate q-values, correction for multiple testing
		benjaminiHochbergCorrection(compFeatures);
		
		// find q-value cutoff by k-mer occurence
//		setQValueCutoff(compFeatures);
	}//end of evaluateConfidence method

	/** compute local lambda around event i, excluding nearby events */
	private double computeLambda(List<ComponentFeature> compFeatures, int i, int c, int length){
		ComponentFeature cf = compFeatures.get(i); 
		Region expandedRegion = cf.getPosition().expand(length/2);
		ArrayList<Region> peakRegions = new ArrayList<Region>();		// peaks in third_lambda_region
		peakRegions.add(cf.getPosition().expand(modelRange));
		int j=0;
		while(true){
			j++;
			if (i+j==compFeatures.size())
				break;
			Region r = compFeatures.get(i+j).getPosition().expand(modelRange);
			if (expandedRegion.overlaps(r))
				peakRegions.add(r);
			else
				break;
		}
		j=0;
		while(true){
			j--;
			if (i+j<0)
				break;
			Region r = compFeatures.get(i+j).getPosition().expand(modelRange);
			if (expandedRegion.overlaps(r))
				peakRegions.add(r);
			else
				break;
		}
		peakRegions = mergeRegions(peakRegions, false);
		expandedRegion = expandedRegion.combine(peakRegions.get(0)).combine(peakRegions.get(peakRegions.size()-1));
		double expandedCount = countIpReads(expandedRegion, c);
		int expandedLength = expandedRegion.getWidth();
		for (Region r:peakRegions){
			expandedLength-=r.getWidth();
			expandedCount-=countIpReads(r, c);
		}
		if (expandedCount==0||expandedLength==0)		// set to 0, so other local lambda will be picked instead
			return 0;
		else
			return expandedCount/expandedLength*(modelRange*2+1);
	}
//	
//	private void setQValueCutoff(List<ComponentFeature>compFeatures) {
//		KmerCluster pCluster = null;
//		for (KmerCluster kc:clusters){
//			if (kc.wm!=null){
//				pCluster = kc;
//				break;		// only use the first PWM, this should be the primary motif
//			}
//		}
//		if (pCluster == null)
//			return;
//		int N = compFeatures.size();
//		int motifHitCount[] = new int[N];		// the number of Motif Hits up to event #i
//		String[] seqs = new String[N];	// DNA sequences around binding sites
//
//		for(int i=0;i<N;i++){
//			Region posRegion = compFeatures.get(i).getPeak().expand(config.k_win/2);
//			seqs[i] = kmf.getSequenceUppercase(posRegion).toUpperCase();
//		}
//		for (int i=0;i<seqs.length;i++){
//			double score = WeightMatrixScorer.getMaxSeqScore(pCluster.wm, seqs[i]);
//			compFeatures.get(i).setEnrichedKmerHGPLog10(score);
////			boolean isHit = score>=pCluster.pwmThreshold;
////			int inc = (isHit)?1:0;
////			if (i==0)
////				motifHitCount[0] = inc;
////			else
////				motifHitCount[i] = inc + motifHitCount[i-1];
//		}
////		StringBuilder sb = new StringBuilder();
////		for (int i=0;i<N-1;i++){
////			int negHit = motifHitCount[N-1]-motifHitCount[i];
////			double hgp_log10 = kEngine.computeHGP(i+1, N-i-1, motifHitCount[i], negHit);
////			compFeatures.get(i).setEnrichedKmerHGPLog10(hgp_log10);
////			sb.append(String.format("%s\t%d\t%d\t%d\t%d\t%.1f\n", compFeatures.get(i).getPeak().toString(), i+1, N-i-1, motifHitCount[i], negHit, hgp_log10));
////		}	
////		CommonUtils.writeFile(outName+"_q_cutoff.txt", sb.toString());
//	}

	private void falseDiscoveryTest(List<ComponentFeature> ipFeatures){
		Vector<ComponentFeature> ctrlFeatures = predictEventsInControlData();
		StringBuilder sb = new StringBuilder();
		for (ComponentFeature cf:ipFeatures){
			sb.append(1+"\t").append("1:"+cf.getTotalEventStrength()).append("\t2:"+cf.getAverageLogKL());
		}
		for (ComponentFeature cf:ipFeatures){
			sb.append(-1+"\t").append("1:"+cf.getTotalEventStrength()).append("\t2:"+cf.getAverageLogKL());
		}
		CommonUtils.writeFile(outName+"_svmTrain.txt", sb.toString());
	}
	// run the same EM-KPP procedure on control data
	private Vector<ComponentFeature> predictEventsInControlData(){
		ArrayList<Region> regions = this.selectEnrichedRegions(new ArrayList<Region>(), false);
        Vector<ComponentFeature> ctrlFeatures = new Vector<ComponentFeature>();
		long tic = System.currentTimeMillis();
		int totalRegionCount = regions.size();
		if (totalRegionCount==0)
			return null;
       
        // prepare for progress reporting
        Vector<Integer> processRegionCount = new Vector<Integer>();		// for counting how many regions are processed by all threads
		int displayStep = (int) Math.pow(10, (int) (Math.log10(totalRegionCount)));
		TreeSet<Integer> reportTriggers = new TreeSet<Integer>();
		for (int i=1;i<=totalRegionCount/displayStep; i++){
			reportTriggers.add(i*displayStep);
		}
		reportTriggers.add(100);
		reportTriggers.add(1000);
		reportTriggers.add(10000);

		// create threads and run EM algorithms
		// the results are put into compFeatures
		Collection<KmerPP> allKmerHits = new Vector<KmerPP>();
        Thread[] threads = new Thread[config.maxThreads];
        log(1,String.format("Running EM on control data: creating %d threads", config.maxThreads));
        int regionsPerThread = regions.size()/threads.length;
        TreeSet<Region> regionsToRun = new TreeSet<Region>();
        regionsToRun.addAll(regions);
        Vector<Region> regionsRunning = new Vector<Region>();		// object to pass info of currently running regions
        for (int i = 0 ; i < threads.length; i++) {
//            ArrayList<Region> threadRegions = new ArrayList<Region>();
//            int nextStartIndex = (i+1)*regionsPerThread;
//            if (i==threads.length-1)		// last thread
//            	nextStartIndex = regions.size();
//            // get the regions as same chrom as possible, to minimize chrom sequence that each thread need to cache
//            for (int j = i*regionsPerThread; j < nextStartIndex; j++) {
//                threadRegions.add(regions.get(j));
//            }
            Thread t = new Thread(new GPS2Thread(regionsToRun,
            									processRegionCount,
            									regionsRunning,
                                                ctrlFeatures,
                                                allKmerHits,
                                                this,
                                                constants,
                                                config,
                                                false));
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
        
        System.out.println(ctrlFeatures.size()+" features found in control data.");
        
        log(1,String.format("%d threads have finished running", config.maxThreads));
        
        return ctrlFeatures;
	}
	
	// split a region into smaller windows if the width is larger than windowSize
	private ArrayList<Region> splitWindows(Region r, int windowSize, int overlapSize, boolean isIP){
		ArrayList<Region> windows=new ArrayList<Region>();
		if (r.getWidth()<=windowSize)
			return windows;
		List<StrandedBase> allBases = new ArrayList<StrandedBase>();
		for (int c=0; c<caches.size(); c++){
			ReadCache cache=null;
			if (isIP)
				cache = caches.get(c).car();
			else
				cache = caches.get(c).cdr();
			List<StrandedBase> bases = cache.getUnstrandedBases(r);
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
		
		// if IP/Control enrichment ratios are lower than cutoff for all ~300bp sliding windows in all conditions, 
		// skip this region, and record it in excludedRegions
		// as long as one of the sub window is enriched, the whole region is enriched, so it is quite conservative, 
		// but it can be useful to block out long continuous regions		
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
			if (!enriched && w!=null){
				excludedRegions.add(w);
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

	private ArrayList<Region> getSubsetRegions(String[] args){
		/* **************************************************
		 * Determine the subset of regions to run EM.
 		 * It can be specified as a file from command line.
		 * If no file pre-specified, estimate the enrichedRegions after loading the data.
		 * We do not want to run EM on regions with little number of reads
		 * **************************************************/
    	String subset_str = Args.parseString(args, "subs", null);
     	String subsetFormat = Args.parseString(args, "subFormat", "");
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
     	return subsetRegions;
	}
	
	/* ***************************************************
	 * Load ChIP-Seq data
	 * ***************************************************/
	private void loadChIPSeqData(ArrayList<Region> subsetRegions, String[] args) {
		this.caches = new ArrayList<Pair<ReadCache, ReadCache>>();
		this.numConditions = conditionNames.size();
		ComponentFeature.setConditionNames(conditionNames);
		condSignalFeats = new ArrayList[numConditions];
		for(int c = 0; c < numConditions; c++) { condSignalFeats[c] = new ArrayList<Feature>(); }
		long tic = System.currentTimeMillis();
		System.out.println("\nGetting 5' positions of all reads...");

		if (experiments.isEmpty()){		// special case, loading RSC file
			for (int i=0;i<numConditions;i++){				
				try {
					ReadCache ipCache = new ReadCache(gen, conditionNames.get(i)+"_IP  ");
					ipCache.readRSC(Args.parseString(args, "--rfexpt"+conditionNames.get(i), ""));
					
					ReadCache ctrlCache = null;
					String ctrlFile = Args.parseString(args, "--rfctrl"+conditionNames.get(i), null);
					if (ctrlFile!=null){
						ctrlCache = new ReadCache(gen, conditionNames.get(i)+"_CTRL");
						ctrlCache.readRSC(ctrlFile);
			    		controlDataExist = true;
					}
					this.caches.add(new Pair<ReadCache, ReadCache>(ipCache, ctrlCache));
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			wholeGenomeDataLoaded = true;
			
			return;
		}
		
		boolean fromReadDB = experiments.get(0).car().isFromReadDB();
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
							int chunkNum = count/constants.MAXREAD*2+1;
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
				ipCache.populateArrays(true);
                //				ipCache.displayStats();
				if (controlDataExist){
					ctrlCache.populateArrays(true);
                    //					ctrlCache.displayStats();
				}
			}
			else if (!fromReadDB){		// load from File
				ipCache.addAllFivePrimes(ip.getAllStarts());
				ipCache.populateArrays(true);
				if (controlDataExist){
					ctrlCache.addAllFivePrimes(ctrl.getAllStarts());
					ctrlCache.populateArrays(true);
				}
				wholeGenomeDataLoaded = true;
			}
			// cleanup			// clean up connection if it is readDB, cleanup data object if it is file
			ip.closeLoaders();
			ctrl.closeLoaders();
			ip=null;
			ctrl=null;
			System.gc();
			if (config.write_RSC_file){
				ipCache.writeRSC();
				ctrlCache.writeRSC();
			}
		} // for each condition

		if (fromReadDB){
			System.out.println("Finish loading data from ReadDB, " + CommonUtils.timeElapsed(tic));
			System.out.println();
		}
	}
	private void doBaseFiltering(){
		//  set max read count for bases to filter PCR artifacts (optional)
		if (config.max_hit_per_bp==0){
			// Mandatory base reset for extremely high read counts
			for(int i=0; i<numConditions; i++){
				Pair<ReadCache,ReadCache> e = caches.get(i);
				ArrayList<Pair<Point, Float>> f = e.car().resetHugeBases(config.base_reset_threshold);
				for (Pair<Point, Float> p:f)
					System.err.printf("%s IP=%.0f-->1 ",p.car().toString(), p.cdr());
				if (controlDataExist){
					f = e.cdr().resetHugeBases(config.base_reset_threshold);
					for (Pair<Point, Float> p:f)
						System.err.printf("%s CTRL=%.0f-->1 ",p.car().toString(), p.cdr());
				}
			}
			System.err.println();
			return;
		}
		
        for(int i=0; i<numConditions; i++){
			Pair<ReadCache,ReadCache> e = caches.get(i);
			double ipCount = e.car().getHitCount();
			// estimate max hit count per BP
			// if user supply using config.max_hit_per_bp, use it
			// if not supplied, take the max of default (3) and possion expected count
			if (config.max_hit_per_bp==-1){
				int maxPerBP = calcExpectedHitCount(ipCount, 1e-9, 1);
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
	
	// compute expected hit count by threshold of Poisson p-value
	// set average read in given region as lambda parameter for Poisson distribution
	private int calcExpectedHitCount(double totalReads, double threshold, int regionWidth){
		int countThres=0;
		DRand re = new DRand();
		Poisson P = new Poisson(0, re);
		double lambda = totalReads*regionWidth /config.mappable_genome_length;
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
		signals.trimToSize();
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
        GPS2Thread t = new GPS2Thread(null, null, null, out, null, this, constants,config, true);
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
	 * Calc the ratio of IP vs Control channel, exlcluding the specific regions
	 */
	private void calcIpCtrlRatio(ArrayList<Region> excludedRegions) {
		// linear regression to get the IP/control ratio
		// now we do not require whole genome data, because partial data could be run on 1 chrom, still enough data to esitmates
		if(controlDataExist) {
			if (config.ip_ctrl_ratio==-1){		// regression using non-specific regions
				ratio_non_specific_total = new double[numConditions];
				for(int t = 0; t < numConditions; t++){
					ratio_non_specific_total[t] = getSlope(t, t, "IP/CTRL", excludedRegions);
				}
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

	/** merge the overlapped regions<br>
	 * if "toExpandRegion"=true, expand each region on both side to leave enough space,
	 * to include every potential reads, then merge
	 */
	private ArrayList<Region> mergeRegions(ArrayList<Region> regions, 
                                             boolean toExpandRegion){
		if (toExpandRegion){
			for (int i=0;i<regions.size();i++){
				regions.set(i, regions.get(i).expand(modelRange, modelRange));
			}
		}
		return Region.mergeRegions(regions);
	}//end of mergeRegions method

	/**
	 * Select the regions where the method will run on.    <br>
	 * It will be either whole genome, chromosome(s) or extended region(s).
	 * This method can be run on either IP or Ctrl reads
	 * @param focusRegion
	 * @return
	 */
	private ArrayList<Region> selectEnrichedRegions(List<Region> focusRegions, boolean isIP){
		long tic = System.currentTimeMillis();
		ArrayList<Region> regions = new ArrayList<Region>();
		Map<String, List<Region>> chr2regions = new HashMap<String, List<Region>>();
		
		// sort regions by chrom for efficiency
		if(focusRegions.size() > 0) {
			for(Region focusRegion:focusRegions) {
				if(!chr2regions.containsKey(focusRegion.getChrom()))
					chr2regions.put(focusRegion.getChrom(), new ArrayList<Region>());
				chr2regions.get(focusRegion.getChrom()).add(focusRegion);
			}			
		}
		else {		// if no given regions, add each chrom as a region
			for (String chrom:gen.getChromList()){
				Region wholeChrom = new Region(gen, chrom, 0, gen.getChromLength(chrom)-1);
				chr2regions.put(chrom, new ArrayList<Region>());
				chr2regions.get(chrom).add(wholeChrom);
			}
		}
		
        Poisson poisson = new Poisson(1, new DRand());
        double totalConditionCount[] = new double[caches.size()];
        for (int i = 0; i < caches.size(); i++) {
        	if (isIP)
        		totalConditionCount[i] = caches.get(i).car().getHitCount();
        	else
        		totalConditionCount[i] = caches.get(i).cdr().getHitCount();
        }
        
		for(String chrom : chr2regions.keySet()) {
			for(Region focusRegion : chr2regions.get(chrom)) {
				List<Region> rs = new ArrayList<Region>();
				HashMap<Region, ArrayList<StrandedBase>> reg2bases = new HashMap<Region, ArrayList<StrandedBase>>();
				ArrayList<StrandedBase> allBases = new ArrayList<StrandedBase>();
				for (int c=0; c<caches.size(); c++){
					List<StrandedBase> bases = null;
					if (isIP)
						bases = caches.get(c).car().getUnstrandedBases(focusRegion);
					else
						bases = caches.get(c).cdr().getUnstrandedBases(focusRegion);
					
					if (bases==null || bases.size()==0){
						continue;
					}
					allBases.addAll(bases); // pool all conditions
				}
				allBases.trimToSize();
				Collections.sort(allBases);	// sort by location
				
				// cut the pooled reads into independent regions
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
							rs.add(r);
							ArrayList<StrandedBase> bases = new ArrayList<StrandedBase>();
							for (int p=start;p<=breakPoint;p++)
								bases.add(allBases.get(p));
							reg2bases.put(r, bases);
						}
						start = breakPoint+1;
					}
				}
				// the last region
				float count = 0;
				for(int m=start;m<allBases.size();m++){
					count += allBases.get(m).getCount();
				}
				if (count>=config.sparseness){
					Region r = new Region(gen, chrom, allBases.get(start).getCoordinate(), allBases.get(allBases.size()-1).getCoordinate());
					rs.add(r);
					ArrayList<StrandedBase> bases = new ArrayList<StrandedBase>();
					for (int p=start;p<allBases.size();p++)
						bases.add(allBases.get(p));
					reg2bases.put(r, bases);
				}

				// check regions, exclude un-enriched regions based on control counts.  If a region is too big (bigger
                // than modelWidth), break it into smaller pieces, keep the region if any of those pieces are enriched
				ArrayList<Region> toRemove = new ArrayList<Region>();
                List<Region> toTest = new ArrayList<Region>();
				for (Region r: rs){
					if (r.getWidth()<=config.min_region_width){
						toRemove.add(r);
						continue;
					}
                    boolean enriched = false;
                    toTest.clear();
                    if (r.getWidth()<=modelWidth){
                        toTest.add(r);
                    } else {
                        for (start = r.getStart(); start < r.getEnd(); start += modelWidth / 3) {
                            Region testr = new Region(r.getGenome(), r.getChrom(), start, start + modelWidth);
                            toTest.add(testr);
                        }
                    }
                    for (Region testr : toTest) {
                        for (int c=0;c<numConditions;c++){
                            int readCount = 0;
                            if (isIP)
                            	readCount = (int)countIpReads(testr,c); 
                            else
                            	readCount = (int)countCtrlReads(testr,c); 
                            
                            poisson.setMean(config.minFoldChange * totalConditionCount[c] * testr.getWidth() / config.mappable_genome_length);
                            double pval = 1 - poisson.cdf(readCount) + poisson.pdf(readCount);
                            if (pval <= .01) {
                                enriched = true;
                                break;
                            }
                            
                            if (isIP & controlDataExist) {		
                                double ctrlreads = countCtrlReads(testr,c);
                                poisson.setMean(config.minFoldChange * ctrlreads * this.ratio_total[c]);
                                pval = 1 - poisson.cdf(readCount) + poisson.pdf(readCount);
                                if (pval <= .01) {
                                    enriched = true;
                                    break;
                                }
                            }                            
                        }
                        if (enriched) { break ;}
                    }
                    if (!enriched){	// remove this region if it is not enriched in all conditions
                        toRemove.add(r);
                        reg2bases.remove(r);
                    } 
				}
				rs.removeAll(toRemove);		
				
				List<Region> smallRegions = new ArrayList<Region>();
				for (Region r:rs){
					int maxSize = 5000;
					if (config.TF_binding)
						maxSize = config.windowSize;
					start = r.getStart();
					int end = r.getEnd();
					if (r.getWidth()>maxSize){ // if the region is too large, break it further at the lowest coverage point
						// base count profile
						float[] profile = new float[end-start+1];
						for (StrandedBase b: reg2bases.get(r))
							profile[b.getCoordinate()-start] = profile[b.getCoordinate()-start]+b.getCount();
						reg2bases.remove(r);
						
						// moving average
						float[] movingAvg = new float[profile.length];
						int halfBin = 100;
						for (int p=0;p<=halfBin*2;p++)
							movingAvg[halfBin]=movingAvg[halfBin]+profile[p];
						for (int p=halfBin+1;p<profile.length-halfBin;p++){
							movingAvg[p]=movingAvg[p-1]-profile[p-1-halfBin]+profile[p+halfBin];
						}
						
						// for every maxSize region, start from modelWidth, find the lowest movingAvg point to break
						int subStart = halfBin;
						while( subStart<profile.length-maxSize+halfBin){
							int subEnd=0;
							float lowest = Float.MAX_VALUE;
							for (int p=subStart+modelWidth;p<subStart+maxSize-halfBin;p++){
								if (movingAvg[p]<lowest){
									subEnd = p;
									lowest = movingAvg[p];
								}
							}
							// if there is a region with same  lowest value, take the middle position
							int p = subEnd;
							for (p=subEnd;p<subStart+maxSize-halfBin;p++){
								if (movingAvg[p]!=lowest)
									break;
							}
							subEnd = (subEnd+p-1)/2;
							smallRegions.add(new Region(gen, chrom, (subStart==halfBin?0:subStart)+start, subEnd+start));
							subStart = subEnd+1;
						}
						if (subStart+start<end)
							smallRegions.add(new Region(gen, chrom, subStart+start, end));
					}
					else
						smallRegions.add(r);
				}
				if (!smallRegions.isEmpty())
					regions.addAll(smallRegions);
			}
		}
		
		log(3, "selectEnrichedRegions(): " + CommonUtils.timeElapsed(tic));
		regions.trimToSize();
		
		// display stats of enriched regions
		int[] binMins = {0,500,1000,2000,3000,5000};
		int[] counts = new int[binMins.length];
		for (Region r:regions){
			for (int i=binMins.length-1;i>=0;i--){
				if (r.getWidth()>binMins[i]){
					counts[i]++;
					break;
				}
			}
		}
		for (int i=0;i<binMins.length;i++){
			log(2, "[" + binMins[i] + " - " + (i==binMins.length-1?"...":binMins[i+1]) + "]\t" + counts[i]);
		}
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

		Collections.sort(comps);
		
		//Make feature calls (all non-zero components)
		// DO NOT remove events with IP < alpha (this happens because EM exit before convergence)
		// because for KPP EM learning, some event may get pseudo-count from kpp
		for(BindingComponent b : comps){
			if(b.getMixProb()>0 ){
				features.add(callFeature(b));
			}
		}

		// assign control reads to IP events for each pair of expt/ctrl (condition)
		// if no control data, skip
		if (controlDataExist){
			// expand to the region covering all events
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
								int idx = base.getCoordinate()-pos-model.getMin();
								if (idx>=modelWidth||idx<0){
									System.err.println("Invalid profile index "+idx+",\tpos "+pos+"\tbase +"+base.getCoordinate()+ " in region "+r.toString());
									continue;
								}
								ctrl_profile_plus[base.getCoordinate()-pos-model.getMin()]=assignment[i][j]*base.getCount();
							}
							else{
								int idx = pos-base.getCoordinate()-model.getMin();
								if (idx>=modelWidth||idx<0){
									System.err.println("Invalid profile index "+idx+",\tpos "+pos+"\tbase -"+base.getCoordinate()+ " in region "+r.toString());
									continue;
								}
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
                shapeDeviation[c]<=config.shapeDeviation &&
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
	private BindingComponent scanPeak(Point p, boolean isIP){
		Region scanRegion = p.expand(config.SCAN_RANGE);
		Region plusRegion = scanRegion.expand(-model.getMin(), model.getMax());
		Region minusRegion = scanRegion.expand(model.getMax(), -model.getMin());
		ArrayList<List<StrandedBase>> signals = new ArrayList<List<StrandedBase>>();
		//Load each condition's read hits
		int totalHitCounts = 0;
		for(int c=0;c<numConditions;c++){
			ReadCache cache = null;
			if (isIP)
				cache = caches.get(c).car();
			else
				cache = caches.get(c).cdr();
			List<StrandedBase> bases_p= cache.getStrandedBases(plusRegion, '+');  // reads of the current region - IP channel
			List<StrandedBase> bases_m= cache.getStrandedBases(minusRegion, '-');  // reads of the current region - IP channel

			bases_p.addAll(bases_m);
			totalHitCounts += StrandedBase.countBaseHits(bases_p);
			signals.add(bases_p);
		} // for loop

		if (totalHitCounts==0) // check for empty read region
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
	 * Apply a Possion filter for duplicate reads
	 */
	private void applyPoissonFilter(boolean printFilterMsg){
		// Filtering/reset bases
		if (printFilterMsg)
			System.out.println("\nFilter duplicate reads.");
		for(int c = 0; c < numConditions; c++) {
			caches.get(c).car().applyPoissonGaussianFilter(10e-3, 20);
			if(controlDataExist) {
				caches.get(c).cdr().applyPoissonGaussianFilter(10e-3, 20);
			}
		}	 		
		
		// print resulting dataset counts
		if (printFilterMsg){
			for(int c = 0; c < numConditions; c++) {
				caches.get(c).car().displayStats();
				if(controlDataExist) {
					caches.get(c).cdr().displayStats();
				}
			}
		}
	}
	/**
	 * This method normalizes conditions as follows:							<br>
	 * It evaluates the slope between each pair of an IP condition and a reference one (condition 1). 
	 * It does the same for the Ctrl channel.    <br>
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
	 * using the non-specific region (all the genome regions, excluding the input regions.)
	 * @param condX_idx index of condition (channel) where it will be used a dependent variable 
	 * @param condY_idx index of condition (channel) where it will be used an independent variable
	 * @param flag <tt>IP</tt> for pairs of IP conditions, <tt>CTRL</tt> for pairs of control conditions, <tt>IP/CTRL</tt> for ip/control pairs 
	 * @param excludedRegions regions to exclude for data for regression, this can be defined for different purposes
	 * @return slope, as the ratio of ( 1st dataset / 2nd dataset )
	 */
	private double getSlope(int condX_idx, int condY_idx, String flag, ArrayList<Region> excludedRegions) {
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
		for(Region r:excludedRegions) {
			String chrom = r.getChrom();
			if(!chrom_regions_map.containsKey(chrom))
				chrom_regions_map.put(chrom, new ArrayList<Region>());
			chrom_regions_map.get(chrom).add(r);
		}
		// for each chrom, construct non-specific regions, get read count
		for(String chrom:gen.getChromList()) {
			List<Region> chrom_non_specific_regs = new ArrayList<Region>();
			int chromLen = gen.getChromLength(chrom);
			// Get the excluded regions of this chrom sorted by location
			List<Region> chr_enriched_regs = new ArrayList<Region>();
			if(chrom_regions_map.containsKey(chrom)) {
				chr_enriched_regs = chrom_regions_map.get(chrom);
				Collections.sort(chr_enriched_regs);
			}
			else{	// We only estimate using the chrom that contains enriched data regions
				continue;
			}
				
			// Construct the non-excluded regions and check for overlapping with the excluded regions
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
				// Estimate the (ctrlCount, ipCount) pairs
				for(Region r:chrom_non_specific_regs) {
					double ctrlCounts_X = countCtrlReads(r, condX_idx);
					double ipCounts_Y   = countIpReads(r, condY_idx);
					if (ctrlCounts_X!=0 && ipCounts_Y!=0)
						scalePairs.add(new PairedCountData(ctrlCounts_X, ipCounts_Y));  // we want to see how many time ipCounts are larger from ctrlCounts
				}
			}
		}//end of for(String chrom:gen.getChromList()) LOOP
			
		// Calculate the slope for this condition
		slope = calcSlope(scalePairs, config.excluded_fraction, config.dump_regression);
		scalePairs.clear();
		return slope;
	}//end of getSlope method
		
	private static double calcSlope(List<PairedCountData> scalePairs, double excludedFraction, boolean dumpRegression) {
		double slope;
        if(scalePairs==null || scalePairs.size()==0) { return 1.0; }
        List<PairedCountData> selectedPairs = new ArrayList<PairedCountData>();
        // exclude the top and bottom fraction of each data series
        if (excludedFraction!=0){
        	double[]x = new double[scalePairs.size()];
        	double[]y = new double[scalePairs.size()];
        	for (int i=0;i<x.length;i++){
        		PairedCountData p = scalePairs.get(i);
        		x[i]=p.x;
        		y[i]=p.y;
        	}
        	Arrays.sort(x);
        	Arrays.sort(y);
        	double xLow = x[(int)(x.length*excludedFraction)];
        	double xHigh = x[(int)(x.length*(1-excludedFraction))-1];
        	double yLow = y[(int)(y.length*excludedFraction)];
        	double yHigh = y[(int)(y.length*(1-excludedFraction))-1];
        	for (int i=0;i<x.length;i++){
        		if (x[i]<=xHigh && x[i]>=xLow && y[i]<=yHigh && y[i]>=yLow)
        			selectedPairs.add(new PairedCountData(x[i],y[i]));
        	}
        }
        else{
        	selectedPairs = scalePairs;
        }
        
        if (dumpRegression){
        	for (PairedCountData p:selectedPairs)
        		System.out.println(p.x+"\t"+p.y);
        }
        	
        DataFrame df = new DataFrame(edu.mit.csail.cgs.deepseq.PairedCountData.class, selectedPairs.iterator());
        DataRegression r = new DataRegression(df, "y~x - 1");
        r.calculate();
        Map<String, Double> map = r.collectCoefficients();
        slope = map.get("x");
        return slope;
    }//end of calcSlope method

	
	double updateBindingModel(int left, int right, String roundLable){
		if (signalFeatures.size()<config.min_event_count){
			System.err.println("\nWarning: The read distribution is not updated, too few ("+signalFeatures.size()+"<"+config.min_event_count+") significant events.");
			return -100;
		}
		int width = left+right+1;
		
		int skipTopPeaks = (int) (0.01 * signalFeatures.size());
		ArrayList<ComponentFeature> cfs = new ArrayList<ComponentFeature>();
		for (Feature f: signalFeatures){
			if (skipTopPeaks>0){
				skipTopPeaks--;
				continue;
			}
			ComponentFeature cf = (ComponentFeature)f;
			// the events that are used to refine read distribution should be
			// having strength and shape at the upper half ranking
			if (config.use_joint_event || !cf.isJointEvent() )
				cfs.add(cf);
		}
		Collections.sort(cfs, new Comparator<ComponentFeature>(){
                public int compare(ComponentFeature o1, ComponentFeature o2) {
                    return o1.compareByTotalResponsibilityWithKmerMatch(o2);
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
						strengthThreshold = cf.getTotalEventStrength();
						break;
					}
				}
			}
		}
		if (strengthThreshold==-1)		// if not set, then we are using all events
			strengthThreshold=cfs.get(cfs.size()-1).getTotalEventStrength();

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
		BindingModel.minKL_Shift(model_plus, model_minus);

		//use single model for now
		List<Pair<Integer, Double>> dist = new ArrayList<Pair<Integer, Double>>();
		for (int i=0;i<width;i++){
			double sum = model_plus[i]+model_minus[i];
			sum = sum>=0?sum:2.0E-300;
			Pair<Integer, Double> p = new Pair<Integer, Double>(i-left, sum);
			dist.add(p);
		}

		double[] oldModel = model.getProbabilities();
		model = new BindingModel(dist);
		model.printToFile(roundLable+"_Read_distribution.txt");
		modelRange = model.getRange();
		modelWidth = model.getWidth();
		allModels.put(roundLable, model);

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
		System.setProperty("java.awt.headless", "true");
	    BufferedImage im = new BufferedImage(w, h, BufferedImage.TYPE_INT_ARGB);
	    Graphics2D g2 = (Graphics2D)im.getGraphics();
	    g2.setRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
	    g2.setColor(Color.white);
	    g2.fillRect(0, 0, w, h);	
	    g2.setColor(Color.gray);
	    g2.drawLine(20, h-margin, w-20, h-margin);		// x-axis
	    g2.drawLine(w/2, margin, w/2, h-margin);	// y-axis    
	    g2.setFont(new Font("Arial",Font.PLAIN,16));
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
		    g2.setFont(new Font("Arial",Font.PLAIN,20));
		    g2.drawString(new File(rounds.get(i)).getName(), w/2+50, i*25+margin+25);
	    }

	    try{
	    	ImageIO.write(im, "png", f);
	    }
	    catch(IOException e){
	    	System.err.println("Error in printing file "+filename);
	    }
	}


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

			totalIPCounts.put(chrom, counts);
			condHitCounts.put(chrom, currChromCondCounts);
		}//end of for(String chrom:gen.getChromList()) loop

	}//end of createChromStats method

	/**
	 * Count non-specific reads  (all the reads not in the refinedRegions)
	 * for scaling experiment and control reads.
	 * Because the signal (specific peak region) can vary a lot in different conditions,
	 * the non-specific regions are a fair estimate for noise.
	 * We assume the noise in expt and control should be comparable, can be used to scale reads.
	 */
	private void countNonSpecificReads(List<ComponentFeature> compFeatures){
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
	public void printFeatures(int round){
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
		String fname = outName+"_GEM_events.txt";
		printFeatures(fname, fs);
	}

	public void printInsignificantFeatures(int round){
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
		String fname = outName+"_GEM_insignificant.txt";
		printFeatures(fname, fs);
	}
	
	public void printFilteredFeatures(int round){
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
			String fname = outName+"_GEM_filtered.txt";
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
				if (!config.print_bound_seqs)
					f.setBoundSequence("");
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
						fw.write("track name=GEM_"+outName+" description=\"GEM Event Call\"\n");
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
	public void releaseMemory(){
		for (Feature f:signalFeatures){
			ComponentFeature cf = (ComponentFeature)f;
			cf.releaseMemory();
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
            count += e.car().countStrandedBases(r,'+');
            count += e.car().countStrandedBases(r,'-');
		}
		return count;
	}
	private float countIpReads(Region r, int cond){
        return caches.get(cond).car().countStrandedBases(r,'+') + 
            caches.get(cond).car().countStrandedBases(r,'-');
	}
	private float countCtrlReads(Region r){
		float count=0;
		for(Pair<ReadCache,ReadCache> e : caches){
            count += e.cdr().countStrandedBases(r,'+');
            count += e.cdr().countStrandedBases(r,'-');
		}
		return count;
	}	
	private float countCtrlReads(Region r, int cond){
        return caches.get(cond).cdr().countStrandedBases(r,'+') + 
            caches.get(cond).cdr().countStrandedBases(r,'-');
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
	public void log(int mode, String msg){
		if (constants.LOG_ALL)
			log_all_msg.append(msg);
		if ( config.verbose>=mode){
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

    /**
	 * Refine restrict regions
	 */
	public void refineRegions(){
		ArrayList<Feature> events = new ArrayList<Feature>();
		events.addAll(signalFeatures);
		if (insignificantFeatures!=null)
			events.addAll(insignificantFeatures);
		if (filteredFeatures!=null)
			events.addAll(filteredFeatures);
		ArrayList<ComponentFeature> compFeatures = new ArrayList<ComponentFeature>();
		for (Feature f : events){
			ComponentFeature cf = (ComponentFeature) f;
			for (int c=0;c<this.numConditions;c++){
				if (cf.getQValueLog10(c)> config.q_refine)		// relax to include more potential regions
					compFeatures.add(cf);
			}
		}
		Collections.sort(compFeatures);
		ArrayList<Region> refinedRegions = new ArrayList<Region>();
		for (ComponentFeature cf:compFeatures){
			refinedRegions.add(cf.getPosition().expand(0));
		}
		this.restrictRegions=mergeRegions(refinedRegions, true);
	}
		
	/**
     * Initalize the kmer engine, then run KMAC<br>
     * This is called only once for initial setup.<br>
     * It compact the cached sequence data and build the KMF<br>
     * If the return value is -1, KMF is not successful, should exit the program.
     */
    public int initKMAC(){
    	if (config.k==-1 && config.k_min==-1)
    		return -1;
		if (signalFeatures.isEmpty())
    		return -2;
    	System.out.println("Loading genome sequences ...");
		kmac = new KMAC(gen, config.cache_genome, config.use_db_genome, config.genome_path);
		long tic = System.currentTimeMillis();

		// setup lightweight genome cache
		if (!kmerPreDefined){
			ArrayList<Region> posituveRegions = new ArrayList<Region>();
			for (Region r: restrictRegions){
				posituveRegions.add(r.expand(config.k_win+modelRange, config.k_shift+config.k_win+modelRange));
			}
			posituveRegions = Region.mergeRegions(posituveRegions);
			int totalLength=0;
			for (Region r: posituveRegions){
				totalLength+=r.getWidth();
			}
			// get negative regions
			ArrayList<Region> negativeRegions = new ArrayList<Region>();

			// In proximal regions, but excluding binding regions
			int winSize = Math.max(config.k_win, config.k_win2); 
			for (Feature f:signalFeatures){
				String chr = f.getPeak().getChrom();
				int length = gen.getChromLength(chr)-1;
				int basis = f.getPeak().getLocation()+config.k_neg_dist;
				for (int i=0;i<config.k_negSeq_ratio;i++){
					int start = basis + (winSize+1)*i;
					if ( start+winSize>=length)
						continue;
					Region r = new Region(gen, chr, start, start+winSize);
					negativeRegions.add(r);
				}
			}
			Collections.sort(negativeRegions);
			ArrayList<Region> toRemove = new ArrayList<Region>();
			for (int i=1;i<negativeRegions.size();i++){
				if (negativeRegions.get(i).overlaps(negativeRegions.get(i-1)))
					toRemove.add(negativeRegions.get(i));
			}
			negativeRegions.removeAll(toRemove);
			
			double gc = kmac.setupRegionCache(posituveRegions, negativeRegions, config.k_neg_dist);
			if (config.gc==-1){
				config.setGC(gc);
			}
			if (config.verbose>1){
				System.out.println("Compact genome sequence cache to " + totalLength + " bps, "+CommonUtils.timeElapsed(tic));
			}
			else
				System.out.println("Done, "+CommonUtils.timeElapsed(tic));
			System.out.println(String.format("GC content=%.2f\n", config.gc));
		}
		
		return runKMAC(config.k_win);		
    }
    
    private ArrayList<ComponentFeature> getEvents(){
		ArrayList<ComponentFeature> events = new ArrayList<ComponentFeature>();
		int count = 1;
		for(Feature f : signalFeatures){
			if(count++>config.k_seqs)
				break;
			ComponentFeature cf = (ComponentFeature)f;
			events.add(cf);
		}
		int significantCount = signalFeatures.size();
		if (config.kmer_use_insig || significantCount<2000 ){
			for(Feature f : insignificantFeatures){
				if(count++>config.k_seqs)
					break;
				ComponentFeature cf = (ComponentFeature)f;
				events.add(cf);
			}
		}
		if (config.kmer_use_filtered || significantCount<2000 ){			
			for(Feature f : filteredFeatures){
				if(count++>config.k_seqs)
					break;
				ComponentFeature cf = (ComponentFeature)f;
				events.add(cf);
			}
		}	
		events.trimToSize();
		
		Collections.sort(events, new Comparator<ComponentFeature>(){
            public int compare(ComponentFeature o1, ComponentFeature o2) {
                if(controlDataExist)
                    return o1.compareByPValue(o2);
                else
                    return o1.compareByPValue_wo_ctrl(o2);
            }
        });
		
		return events;
    }   
    
    /**
     * Run the K-mer motif discovery procedure <br>
     * If the return value is -1, KMAC is not successful, should exit the program.
     */
    public int runKMAC( int winSize){
    	// set the parameters
    	kmac.setConfig(config, outName);
		int[] eventCounts = new int[]{signalFeatures.size(), insignificantFeatures.size(), filteredFeatures.size()};
   	
    	// load sequence from binding event positions
    	ArrayList<ComponentFeature> events = getEvents();    	
    	kmac.loadTestSequences(events, winSize);
    	log(1, String.format("Running KMAC motif discovery with %d+/%d- sequences.\n", kmac.getPositiveSeqs().length, kmac.getNegSeqCount()));
    	if (config.print_input_seqs)
    		kmac.printInputSequences(outName);
    	
    	// select best k value
		if (config.k_min!=-1){
			// compare different values of k to select most enriched k value
			int bestK = config.selectK_byTopKmer?kmac.selectK_byTopKmer(config.k_min, config.k_max, eventCounts):kmac.selectK(config.k_min, config.k_max, eventCounts);
			if (bestK!=0){
				config.k = bestK;
				if (config.allow_seed_inheritance){
					config.k_min = -1;		// prevent selecting k again
					config.k_max = -1;
				}
			}
			else{
				log(1, "The value of k can not be determined automatically.");
				return -1;
			}
		}
		else{
			config.k_win = winSize;
		}	
		
		// select enriched k-mers, cluster and align
		ArrayList<Kmer> kmers = kmac.selectEnrichedKmers(config.k);
		kmers = kmac.KmerMotifAlignmentClustering(kmers, -1, false, eventCounts);
		if (kmers ==null)		// No motif found, exit here!
			return -1;
		if (kmers.isEmpty()){
			System.err.print("Not able to find KSM motif");
			if (kmac.getPrimaryCluster().wm!=null){
				config.pp_use_kmer = false;
				System.err.println(" , use PWM as prior!");
			}else{
				System.err.println(" and PWM motif, exit here!");
				return -1;
			}
		}
		
		// use only primary cluster k-mers for search
		ArrayList<Kmer> primaryKmers = new ArrayList<Kmer>();
		for (Kmer km:kmers)
			if (km.getClusterId()==0)
				primaryKmers.add(km);
		kmac.updateEngine(primaryKmers);		
		kcm_threshold = kmac.getPrimaryCluster().ksmThreshold==null?0:kmac.getPrimaryCluster().ksmThreshold.score;
		return 0;
    }
    	
	public MotifThreshold estimateClusterKgsThreshold(ArrayList<Kmer> clusterKmers){
		kmac.updateEngine(clusterKmers);
		return kmac.optimizeKsmThreshold(outName, false);

	}
	
	public void estimateOverallKgsThreshold(){
		if (! kmac.isInitialized())
			return;
		
		MotifThreshold t = kmac.optimizeKsmThreshold(outName, true);
		if (t!=null)
			System.out.println(String.format("%.2f\t%d\t%d\t%.1f\n", t.score, t.posHit, t.negHit, t.hgp ));
	}
    
	/**
	 * Print the overlapping kmer counts with increasing number of k
	 * in a specified window around the binding events
	 */
	public void printOverlappingKmers(){
	    	ArrayList<ComponentFeature> events = getEvents();
	    	for (int k=config.k;k<config.k+5;k++){
			String name = outName+"_";//OK_win"+ (config.k*2);
		}
	}

    class GPSConstants {

        /****************************
         * Constants
         ***************************/
	
        public final boolean LOG_ALL=false;

        // width for smoothing a read (used as a stddev for creating the Gaussian kernel of probability)
        public final int READ_KERNEL_ESTIMATOR_WIDTH = 5;

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
        public final int MAXREAD = 1000000;
        // true:  eliminated component in batch, as long as matching criteria in EM derivation
        // false: eliminate only the worse case components, re-distribute the reads of eliminated component
        public final boolean BATCH_ELIMINATION = false;
        public final boolean SMART_SPACING = true;		// dynamically determine init comopent spacing
        public final boolean MAKE_HARD_ASSIGNMENT=false;
    }


    class GPS2Thread implements Runnable {
    	
        private KPPMixture mixture;
        private GPSConstants constants;
        private Config config;
        private TreeSet<Region> regions;
        private Collection<Integer> processedRegionCount;
        private Collection<Region> regionsRunning;
        private Collection<ComponentFeature> compFeatures;
        private Collection<KmerPP> allKmerHits;
        private boolean isIP;
        
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

        public GPS2Thread (TreeSet<Region> regions,
        				  Collection<Integer> processedRegionCount,
        				  Collection<Region> regionsRunning,
                          Collection<ComponentFeature> compFeatures,
                          Collection<KmerPP> allKmerHits,
                          KPPMixture mixture,
                          GPSConstants constants,
                          Config config,
                          boolean isIP) {
            this.regions = regions;
            this.processedRegionCount = processedRegionCount;
            this.regionsRunning = regionsRunning;
            this.mixture = mixture;
            this.constants = constants;
            this.config = config;
            this.compFeatures = compFeatures;
            this.allKmerHits = allKmerHits;
            this.isIP = isIP;
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
                    result = EMTrain(signals, null, alpha, new double[r.getWidth()]);
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
        	
        	SequenceGenerator<Region> seqgen = new SequenceGenerator<Region>();
        	if (config.use_db_genome)
        		seqgen.useLocalFiles(false);
        	if (config.cache_genome)
        		seqgen.useCache(true);
            while (!regions.isEmpty()) {
            	Region rr = null;
            	synchronized (regions){
	            	if (!regions.isEmpty()){
	            		rr = regions.first();
	            		regions.remove(rr);
            		}
	            	else
	            		break;
            	}
                mixture.log(3, rr.toString());
                regionsRunning.add(rr);
                try{
                    ArrayList<BindingComponent> comps= new ArrayList<BindingComponent>();
                    // Cut long regions into windowSize(1.5kb) sliding window (500bp overlap) to analyze
                    ArrayList<Region> windows = new ArrayList<Region>();
                    if (rr.getWidth()<=config.windowSize)
                        windows.add(rr);
                    else{
                        windows = mixture.splitWindows(rr, config.windowSize, mixture.modelWidth/2, isIP);
                    }
		 
                    // run EM for each window 
                    for (Region w : windows){
                        ArrayList<BindingComponent> result = analyzeWindow(w, seqgen);
                        if (result!=null){
                            comps.addAll(result);
                        }
                    }
		
                    /* ****************************************************************
                     * fix sliding window boundary effect (version 1)
                     * ****************************************************************/
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
                                            int winSize = mixture.modelWidth*5;
                                            if (r.getWidth()<winSize){ // if the region is small, directly process it
                                                ArrayList<BindingComponent> result = analyzeWindow(r, seqgen);
                                                if (result!=null){
                                                    comps.addAll(result);
                                                }
                                            } else {	// if the region is too long, split into windows (5kb size, 1kb overlap)
                                                ArrayList<Region> wins = mixture.splitWindows(r, winSize, mixture.modelWidth, isIP);
                                                // process each window, then fix boundary 
                                                ArrayList<ArrayList<BindingComponent>> comps_all_wins= new ArrayList<ArrayList<BindingComponent>>();
                                                for (Region w : wins){
                                                    ArrayList<BindingComponent> comp_win = new ArrayList<BindingComponent>();
                                                    ArrayList<BindingComponent> result = analyzeWindow(w, seqgen);
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
                     * if we have positional prior, use the EM result directly
                     * ****************************************************************/
                    if (kmac!=null){
                    	compFeatures.addAll(mixture.callFeatures(comps));
                    }
                    /* ****************************************************************
                     * if not positional prior, refine unary events by scanEvent()
                     * and collect all the resulting events as features
                     * this is last step because fixing boundary may result in some new unary events
                     * ****************************************************************/
                    else if (comps.size()>0){
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
                        for (ArrayList<Region> subr : subRegions){	// for each independent regions (may have multiple events)
                            ArrayList<BindingComponent> bs = new ArrayList<BindingComponent>();
                            if (subr.size()==1){
                                BindingComponent b;
                                // Scan Peak unary region
                                if(config.use_scanPeak) {
                                    Point p = subr.get(0).getMidpoint();
                                    b = mixture.scanPeak(p, isIP);
                                    if (b==null){
                                        continue;
                                    }
                                    b.setEMPosition(p);
                                    for (BindingComponent m:comps){
                                        if (p.getLocation()==m.getLocation().getLocation()) {
                                            b.setAlpha(m.getAlpha());	// inherit alpha value from previous EM component
                                            b.setNoiseFraction(m.getNoiseFraction());
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
                                    ArrayList<BindingComponent> bl = analyzeWindow(subr.get(0).expand(mixture.modelRange, mixture.modelRange), seqgen);
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
        			cleanUpDataLoader();
                    System.exit(-1);
                } finally{
                	regionsRunning.remove(rr);
                }
            }
        }

        private ArrayList<BindingComponent> analyzeWindow(Region w, SequenceGenerator<Region> seqgen){

            ArrayList<List<StrandedBase>> signals = mixture.loadData_checkEnrichment(w);
            if (signals==null)
                return null;

            ArrayList<List<StrandedBase>> bg_signals = null;
            if (mixture.controlDataExist && config.noise_distribution==2){		// if want to use control to set bg dist. 
            	bg_signals = loadBasesInWindow(w, "CTRL");
            }
            
            // We want to run EM only for potential overlapping regions
            // After first round, if we are sure the region contains unary event, we will just scan for peak
            if (!singleEventRegions.containsKey(w)) {
                // dynamically determine an alpha value for this sliding window
            	// when having a noies component, reduce alpha value further
                double alpha = config.sparseness;
                if (config.use_dynamic_sparseness){
                    alpha = Math.max(mixture.estimateAlpha(w, signals), config.sparseness);
                    if (config.noise_distribution!=0)		// have a noise model
                    	alpha /= config.alpha_fine_factor;
                }

                Pair<double[][][], int[][][]> result = null;

                // Choose either Yuchun's beta EM or Giorgos's pool EM methods
                if(config.use_betaEM) {

                	//construct the positional prior for each position in this region
                	double[] pp = new double[w.getWidth()+1];
                	KmerGroup[] pp_kmer = new KmerGroup[pp.length];
                	String seq = null;
                	if (kmac!=null && kmac.isInitialized()){
                		if (kmerPreDefined)		// if kmer is loaded from other sources, get fresh sequence 
                			seq = seqgen.execute(w).toUpperCase();
                		else					// otherwise, we have run KMF, get the cached sequences
                			seq = kmac.getSequenceUppercase(w);
                		
                		HashMap<Integer, KmerPP> hits = new HashMap<Integer, KmerPP>();
                		if (config.pp_use_kmer){		// use k-mer match to set KPP
	                		KmerGroup[] matchPositions = kmac.query(seq);
		                	if (config.print_PI)	
		                		System.out.println(seq);
		                	// Effectively, the top kmers will dominate, because we normalize the pp value
		                	for (KmerGroup g: matchPositions){
		                		// the posBS is the expected binding position
		                		int bindingPos = g.getPosBS();
		                		if (bindingPos>=pp.length || bindingPos<0)
		                			continue;
		                		if (g.getHgp()> -kcm_threshold)	// must past ksm threshold
		                			continue;
	
		                		double kmerCountSum = 0;
		                		if (config.use_weighted_kmer)
		                			kmerCountSum = g.getWeightedKmerStrength();
		                		else
		                			kmerCountSum = g.getGroupHitCount();
		                		
		                		// select the approach to generate pp from kmer count
		                		if (config.kpp_mode==0)
		                			pp[bindingPos] = kmerCountSum;
		                		else if (config.kpp_mode==1)
		                			pp[bindingPos] = kmerCountSum==0?0:Math.log(kmerCountSum);
		                		else if (config.kpp_mode==10)
		                			pp[bindingPos] = kmerCountSum==0?0:Math.log10(kmerCountSum);
		                		pp_kmer[bindingPos] = g;
		                		if (config.print_PI)	
			                		System.out.println(bindingPos+"\t"+g.getBestKmer().getKmerString());
		                		hits.put(bindingPos, new KmerPP(new Point(gen, w.getChrom(), w.getStart()+bindingPos), g, pp[bindingPos]));
		                	}
                		}
                		else {			// use PWM to set KPP
                			KMAC.KmerCluster cluster = kmac.getPrimaryCluster();
                			if (cluster!=null){
                				WeightMatrixScorer scorer = new WeightMatrixScorer(cluster.wm);
                				WeightMatrixScoreProfile profiler = scorer.execute(seq);
                				for (int i=0;i<profiler.length();i++){
                					double max = profiler.getMaxScore(i);
                					if (max<0)
                						continue;
                					char strand = profiler.getMaxStrand(i);
                					int pos = cluster.pos_BS_seed-cluster.pos_pwm_seed;
                					if (strand=='+')
                						pos = pos + i;
                					else
                						pos = i + cluster.wm.length()-1 + (-pos);
                					if (pos<0||pos>=seq.length())
                						continue;
                					if (max>pp[pos]){
                						pp[pos] = max;
                					}
                				}
                			}
                		}
                		
	                	// if kmer hits are 100bp apart, consider them as independent motif hit
	                	// therefore, the pp value are normalized to alpha in local region of 100bp apart
	                	// find all the dividing boundaries
	                	int GAP = 100;
	                	int prevBoundary = 0;
	                	ArrayList<Integer> boundaries=new ArrayList<Integer>() ;
	                	boundaries.add(0);
	                	for (int i=1;i<pp.length;i++){
	                		if (pp[i]>0){					// for non-zero pp
	                			if (i-prevBoundary>GAP){	// a new local
	                				boundaries.add(i);		// this is right exclusive
	                			}
	                			prevBoundary = i;
	                		}
	                	}
	                	if (prevBoundary!=pp.length-1)
	                		boundaries.add(pp.length-1);
	                	
	                	for (int i=0;i<boundaries.size()-1;i++){
	                		double total = 0;
	                		for (int j=boundaries.get(i);j<boundaries.get(i+1); j++){
	                			if (config.kpp_normalize_max)
			                		total = Math.max(total, pp[j]);	// scale so that (largest position prior == sparse prior)
	                			else
			                		total += pp[j];	// normalize the local region so that total positional prior pseudo-count equal to alpha (sparse prior)
	                		}
	                		for (int j=boundaries.get(i);j<boundaries.get(i+1); j++){
		                		if (pp[j]>0){
		                			pp[j] = pp[j]/total*(alpha*config.kpp_factor);
		                			if (config.pp_use_kmer)
		                				hits.get(j).pp = pp[j];
	                			}
	                		}
	                	}
	                	if (config.kmer_print_hits){
	                		allKmerHits.addAll(hits.values());
	                		hits.clear();
	                	}
                	}
                	
                    //Run EM and increase resolution
                    initializeComponents(w, mixture.numConditions);
                    int countIters = 0;
                    while(nonZeroComponentNum>0){                        
                        int numComp = components.size();
                        double[] pos_alpha = new double[numComp];						// positional alpha
                        if (kmac!=null){
	                        for (int i=0;i<numComp;i++){
	                        	BindingComponent b = components.get(i);
	                        	int bIdx = b.getLocation().getLocation()-w.getStart();
	                        	double maxPP = 0;
	                        	for (int j=0;j<componentSpacing;j++){
	                        		int idx = bIdx+j;
	                        		if (idx>=pp.length)
	                        			idx = pp.length-1;
	                        		maxPP = Math.max(maxPP,pp[idx]);
	                        	}
	                        	pos_alpha[i]=maxPP;
	                        }
                        }
                        else{
                        	for (int i=0;i<numComp;i++){
                        		pos_alpha[i]=0;
                        	}
                        }
                        // EM learning, components list will only contains non-zero components
                        result = EMTrain(signals, bg_signals, alpha, pos_alpha);
                        // mixture.log(4, componentSpacing+" bp\t"+(int)nonZeroComponents+" components.");

                        countIters++;
                        //					System.out.printf("%s\tupdateIter\t%d\tnumComps\t%d\tnumNonZeroComps\t%d%n", w.toString(), countIters, numAllComps, components.size());
                        // increase resolution
                        if (componentSpacing==1)
                            break;
                        int lastResolution = componentSpacing;
                        updateComponentResolution(w, mixture.numConditions, componentSpacing);
                        if(componentSpacing==lastResolution)
                            break;
                    } 	// end of while (resolution)
                    if (nonZeroComponentNum==0)	
                    	return null;
                    
                    setComponentResponsibilities(signals, result.car(), result.cdr());
                    if (kmac!=null && config.pp_use_kmer)
                    	setEventKmerGroup(pp_kmer, w.getStart(), seq);
                } else {
                    // Run MultiIndependentMixture (Temporal Coupling) -- PoolEM
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
                }// end of else condition (running it with MultiIndependentMixture TC)

            } else{// if single event region, just scan it
                BindingComponent peak = mixture.scanPeak(singleEventRegions.get(w), isIP);
                components = new ArrayList<BindingComponent>();
                components.add(peak);
            }		
            
            // discard components with less than alpha reads
            ArrayList<BindingComponent> toRemove = new ArrayList<BindingComponent>();
            for (BindingComponent c: components){
            	double discard_threshold = config.discard_subAlpha_components?c.getAlpha():5;
            	if (c.getTotalSumResponsibility()<discard_threshold)
            		toRemove.add(c);
            }
            components.removeAll(toRemove);

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
         * Initializes the components. Set event spacing evenly.
         *
         * @param currReg
         * @param signals
         * @param weighted (uniform or weighted initialization?)
         */
        private void initializeComponents(Region currReg, int numCond){
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
        private Pair<double[][][], int[][][]>  EMTrain(ArrayList<List<StrandedBase>> signals, ArrayList<List<StrandedBase>> bg_signals, double alpha, double[] pos_alpha){
            int numComp = components.size();
            // H function and responsibility will be stored using an indirect indexing method
            // Because only the components within modelRange matters, we only store those components around the reads
            // and use a mapping array to keep track of the component index
            // This will reduce the memory requirement from N*numComp, to NInRange*numComp
            // and should save some running time because we only iterate the effective base-comps
            double[][] counts= new double[mixture.numConditions][];	// Hit Count	[cond][base]
            double[][][] h= new double[mixture.numConditions][][]; 		// H function [cond][comp][base]
            double[][][] r= new double[mixture.numConditions][][];		// Responsibility
            int[][][] c2b= new int[mixture.numConditions][][]; 			// mapping from component to base
            double[][] b= new double[mixture.numConditions][];			// Beta
            double[] pi = new double[numComp];							// Pi
            double[][] prob_bg = new double[mixture.numConditions][];	// Background probability per base for each condition
            
            // init pi
            for(int j=0;j<numComp;j++){
                BindingComponent comp = components.get(j);
                pi[j]= comp.getMixProb();
            }

            boolean no_data_bin = false;
            int minPlus=Integer.MAX_VALUE, maxPlus=Integer.MIN_VALUE,minMinus=Integer.MAX_VALUE, maxMinus=Integer.MIN_VALUE;

            for(int c=0; c<mixture.numConditions; c++){
                List<StrandedBase> bases_old = signals.get(c);
                List<StrandedBase> bases = new ArrayList<StrandedBase>();
                if (componentSpacing==1 || no_data_bin) {
                    bases = bases_old;
                } else {		// merge read counts into data bin
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

                // compute the range of the reads, to estimate the init prob. for bg component
                for(int i=0;i<numBases;i++){
                    int coord = bases.get(i).getCoordinate();
                    if (bases.get(i).getStrand()=='+'){
                    	if (coord<minPlus)
                    		minPlus = coord;
                    	if (coord>maxPlus)
                    		maxPlus = coord;
                    }
                    else{
                    	if (coord<minMinus)
                    		minMinus = coord;
                    	if (coord>maxMinus)
                    		maxMinus = coord;
                    }
                }
                
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
                
                // background probability
                if (config.noise_distribution==2 && mixture.controlDataExist){ //config.noise_distribution==2 smoothed shape using control data
                	List<StrandedBase> ctrls = bg_signals.get(c);
                	double[] plusProfile = new double[Math.max(model.getRange(), maxPlus-minPlus+1)];
                	for (StrandedBase base_ctrl: ctrls){
                		if (base_ctrl.getStrand()=='+'){
	                		int offset = base_ctrl.getCoordinate() - minPlus;
	                		if (offset>=0 && offset<plusProfile.length)
	                			plusProfile[offset] = base_ctrl.getCount();
                		}
                	}
                	for (int i=0;i<plusProfile.length;i++)			// pseudo-count to deal with 0 ctrl reads
                		plusProfile[i] += StrandedBase.countBaseHits(ctrls)/plusProfile.length;
//                	System.out.println();
//                	System.out.println(CommonUtils.arrayToString(plusProfile, "%.4f"));
                	plusProfile = StatUtil.symmetricKernelSmoother(plusProfile, mixture.gaussian);
//                	System.out.println(CommonUtils.arrayToString(plusProfile, "%.4f"));
                	double[] minusProfile = new double[Math.max(model.getRange(), maxMinus-minMinus+1)];
                	for (StrandedBase base_ctrl: ctrls){
                		if (base_ctrl.getStrand()=='-'){
	                		int offset = base_ctrl.getCoordinate() - minMinus;
	                		if (offset>=0 && offset<minusProfile.length)
	                			minusProfile[offset] = base_ctrl.getCount();
                		}
                	}
                	for (int i=0;i<minusProfile.length;i++)			// pseudo-count to deal with 0 ctrl reads
                		minusProfile[i] += StrandedBase.countBaseHits(ctrls)/minusProfile.length;
                	minusProfile = StatUtil.symmetricKernelSmoother(minusProfile, mixture.gaussian);
                	
                	double [] prob_bg_c = new double[numBases];
                	for(int i=0;i<numBases;i++){
                		StrandedBase base = bases.get(i);
                		if (base.getStrand()=='+')
                			prob_bg_c[i] = plusProfile[base.getCoordinate()-minPlus];
                		else
                			prob_bg_c[i] = minusProfile[base.getCoordinate()-minMinus];
                	}
	                prob_bg[c] = prob_bg_c;
                }
                else if (config.noise_distribution!=0){	// uniform shape 
                    int range = Math.max(Math.max(maxPlus-minPlus+1, maxMinus-minMinus+1), model.getRange());
                    double prob = 1.0/range; 
		            double [] prob_bg_c = new double[numBases];
	                for(int i=0;i<numBases;i++)
	                	prob_bg_c[i] = prob;
	                prob_bg[c] = prob_bg_c;
                }
                 
            }
            mixture.log(5, "\n"+componentSpacing+" bp:\t");   
            
            //////////
            // Run EM steps
            //////////
            Pair<double[][][], Double> result = EM_MAP(counts, h, r, b, c2b, pi, alpha, pos_alpha, prob_bg);
            r = result.car();
            double noiseProb = result.cdr();
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
                    comp.setNoiseFraction(noiseProb);
                    for(int c=0; c<mixture.numConditions; c++){
                        double[] bc = b[c];
                        comp.setConditionBeta(c, bc[j]);
                        comp.setOld_index(j);
                    }
                    nonZeroComponents.add(comp);
                }
            }
            components = nonZeroComponents;

            // Assign the summed responsibilities only to non-zero components at 1bp resolution
            if (componentSpacing==1){
                for(int c=0; c<mixture.numConditions; c++){
                    for(int j=0;j<components.size();j++){	
                        double sum_resp = 0.0;	
                        int oldIndex = components.get(j).getOld_index();
                        int[] baseIdx = c2b[c][oldIndex];
                        for(int i=0;i<baseIdx.length;i++){
                            sum_resp += counts[c][baseIdx[i]]*r[c][oldIndex][i];
                        }
                        components.get(j).setCondSumResponsibility(c, sum_resp);
                    }
                }
            }
	
            return new Pair<double[][][], int[][][]>(r, c2b);
        }//end of EMTrain method


        /** 
         * core EM steps, sparse prior (component elimination), multi-condition
         * 
         * @param counts Read counts per base for each condition
         * @param h Prob of base given event location for each condition
         * @param r Responsibility of event to base for each condition
         * @param b Beta (proportion of pi) of event for each condition
         * @param c2b Mapping from component to base
         * @param pi Mixing probability of event
         * @param alpha Sparse prior (uniform negative Dirichlet prior)
         * @param pos_alpha Positional prior (positive count per event position)
         * @param prob_bg Background probability per base for each condition
         * @return
         */
        private Pair<double[][][], Double> EM_MAP(  	double[][]   counts,
                                                        double[][][] h,
                                                        double[][][] r,
                                                        double[][]   b,
                                                        int[][][] c2b,
                                                        double[] pi,
                                                        double alpha,
                                                        double[] pos_alpha, 
                                                        double[][] prob_bg) {
        	
        	boolean model_noise = config.noise_distribution !=0;
        	
        	double totalCounts = 0;
        	for (int i=0;i<counts.length;i++){
        		for (int j=0;j<counts[i].length;j++)
        			totalCounts += counts[i][j];
        	}
        	
            long tic=System.currentTimeMillis();
            ArrayList<EM_State> models = new  ArrayList<EM_State> ();

            // variable for bg (noise) component
            double[] b_bg = new double[mixture.numConditions];
            double pi_bg = 0;
            if (model_noise)
            	pi_bg = background_proportion;
            double pi_signal=1-pi_bg;
            double [][] r_bg = new double[mixture.numConditions][];	// [cond][base]
            for(int c=0; c<mixture.numConditions; c++){
            	r_bg[c] = new double[counts[c].length];
            	b_bg[c] = 1.0/mixture.numConditions;
            }
            if (model_noise){
            	for (int j=0;j<pi.length;j++){
            		pi[j]=pi[j] * pi_signal;
            	}
            	for(int c=0; c<mixture.numConditions; c++){
            		for (int j=0;j<pi.length;j++){
                        int[] baseIdx = c2b[c][j];
                        for(int i=0;i<baseIdx.length;i++)
                        	r[c][j][i] *= pi_signal;
                    }
            	}
            	for(int c=0; c<mixture.numConditions; c++){
            		for(int i=0;i<r_bg[c].length;i++)
            			r_bg[c][i] = prob_bg[c][i] * pi_bg * b_bg[c];
            	}
            }
            
            boolean hasPP = false;
            for (double p_a:pos_alpha){
            	if (p_a!=0){
            		hasPP = true;
            		break;
            	}
            }
            if (config.print_PI){
            	System.out.println(alpha+"\t"+CommonUtils.arrayToString(pos_alpha, "%.4f"));
            	System.out.println(-1+"\t"+CommonUtils.arrayToString(pi, "%.4f")+
            			(model_noise?String.format("\t%.2f",pi_bg):""));
            }
            
            double lastLAP=0, LAP=0; // log posterior prob
            int t=0;
            double currAlpha = alpha/config.gentle_elimination_factor;
            double effectiveAlpha = 0;	// the real Alpha applied to compute pi (may be different from currAlpha because of the dynamic elimination schedule)
            
            double minProb = 1.0/(pi.length*2);		// in ML speedup, threshold to eliminate ML components, length*2 to make sure it is smaller than 1/m
            double maxMinProb = currAlpha / totalCounts;		// the max minProb is bound by alpha
            
            boolean minElimination = false;
            int gentleCounts = 0; 	// count the iterations that we run on gentle mode after last elimination
            // when reach the threshold, increase currAlpha value
            // index of non-zero components, used to iterate components
            ArrayList<Integer> nzComps = new ArrayList<Integer>();	
            for (int j=0;j<pi.length;j++){
                nzComps.add(j);
            }
            mixture.log(5, (int)nonZeroComponentNum+" ");            
            
            //Run EM
//            System.out.println("maxMinProb="+maxMinProb);
            for(t=0; t<constants.MAX_EM_ITER ; t++){
//            	long toc = System.currentTimeMillis();
//            	System.out.println("t="+t+"\t"+minProb+"\t"+(toc-tic)+"\t"+nzComps.size());
//            	tic = toc;
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
                    if (model_noise && pi_bg!=0){
	                    for(int i=0;i<numBases;i++){
	                    	totalResp[i] += r_bg[c][i];
	                    }
                	}
                    
                    // normalize
                    for(int j:nzComps){
                        int[] baseIdx = c2b[c][j];
                        for(int i=0;i<baseIdx.length;i++)
                            if (totalResp[baseIdx[i]]>0)
                                rc[j][i] = rc[j][i]/totalResp[baseIdx[i]];
                    }
                    if (model_noise && pi_bg!=0){
	                    for(int i=0;i<numBases;i++){
	                    	r_bg[c][i] = r_bg[c][i]/totalResp[i];
	                    }
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
                        effectiveAlpha = currAlpha;
                        for(int j:nzComps){
                            double r_sum=0;
                            for(int c=0; c<mixture.numConditions; c++){
                                int[] baseIdx = c2b[c][j];
                                for(int i=0;i<baseIdx.length;i++)
                                    r_sum += r[c][j][i]*counts[c][baseIdx[i]];
                            }

                            // component elimination
                            pi[j]=Math.max(0,r_sum-currAlpha+pos_alpha[j]);

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
                            r_sum[jnz] += pos_alpha[nzComps.get(jnz)];	// adding positional prior as pseudo-count
                        }
                        if (gentleCounts>=config.gentle_elimination_iterations && currAlpha<alpha )
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
                            effectiveAlpha = currAlpha;
                            //                     	System.out.println(t+":\t"+currAlpha+"\t iterating");
                        }else{
                            // eliminate worst case components, could be 1 or multiple components
                            // not apply alpha here, redistribute responsibilities in next E step
                        	for(int jnz=0;jnz<r_sum.length;jnz++){
	                            if (worst.cdr().contains(jnz)){
	                                pi[nzComps.get(jnz)]=0;                        	
	                                // clear responsibility
	                                for(int c=0; c<mixture.numConditions; c++){
	                                    for(int i=0; i<c2b[c][nzComps.get(jnz)].length;i++)
	                                        r[c][nzComps.get(jnz)][i] = 0;
	                                }
	                            }
	                            else
	                            	pi[nzComps.get(jnz)]=r_sum[jnz];	// not normailzed, not apply alpha
                        	}
                        	effectiveAlpha = worst.car();
                            // keep iterating on this Alpha value, until converge, then we raise it up to eliminate next one
                            // give EM some time to stabilize before eliminating the next components
                            //                    	System.out.println(t+":\t"+currAlpha+"\t elimination");
                            currAlpha = Math.max(worst.car(), alpha/config.gentle_elimination_factor/2);
                            gentleCounts = 0;
                            //                    	currAlpha = 0;
                        }
                    }
                }
                // BG component computation is the same for either ML or MAP
                if (model_noise && pi_bg!=0){
                	double r_bg_sum=0;
                	for(int c=0; c<mixture.numConditions; c++){
                		for(int i=0;i<r_bg[c].length;i++)
                			r_bg_sum += r_bg[c][i]*counts[c][i];
                	}
//                	pi_bg = r_bg_sum;
                	pi_bg = r_bg_sum-effectiveAlpha;
                	if (pi_bg<0)
                		pi_bg = 0;
                }
//                if (config.print_PI)
//                	System.out.println(t+"Resp\t"+CommonUtils.arrayToString(pi, "%.4f")+
//                			(model_noise?String.format("\t%.4f",pi_bg):""));                
                // normalize pi
                double totalPi=0;
                for(int j=0;j<pi.length;j++){
                    if (pi[j]!=0){
                        totalPi+=pi[j];
                    }
                }      
                if (model_noise && pi_bg!=0){
                	totalPi+=pi_bg;
                }
                if (totalPi!=0){
                    for(int j:nzComps){
                        pi[j]=pi[j]/totalPi;
                    }                
                    if (model_noise && pi_bg!=0){
                    	pi_bg/=totalPi;
                    }

                }
                // EM speed up, eliminate components with probability less than initial avg prob.
                // Only do this for coase spacing, which has more components
                
                if (t<=config.ML_ITER && componentSpacing!=1 && config.ML_speedup){
                	boolean eliminated = false;
                	for (int j:nzComps){
                		if (pi[j]<minProb && pos_alpha[j]==0){
                			// eliminate component, clear responsibility
                			pi[j]=0;
                			eliminated = true;
                			for(int c=0; c<mixture.numConditions; c++){
                                for(int i=0; i<c2b[c][j].length;i++)
                                    r[c][j][i] = 0;
                            }
                		}
                	}
                	if (eliminated){
	                	 // normalize pi here again
	                    totalPi=0;
	                    for(int j:nzComps){
	                        if (pi[j]!=0){
	                            totalPi+=pi[j];
	                        }
	                    }      
	                    if (model_noise){
	                    	totalPi+=pi_bg;
	                    }       	
	                    if (totalPi!=0){
	                        for(int j:nzComps){
	                            pi[j]=pi[j]/totalPi;
	                        }
		                    if (model_noise){
		                    	pi_bg/=totalPi;
		                    }
	                    }
                	}
                	else{	// not eliminated, double minProb
                		if (minProb*2<maxMinProb)
                			minProb *= 2;
                	}
                }
                
                // update component count
                nzComps.clear();
                for(int j=0;j<pi.length;j++){
                    if (pi[j]!=0){
                        nzComps.add(j);
                    }
                }
                nonZeroComponentNum = nzComps.size();
                if (nonZeroComponentNum==0)
                    return new Pair<double[][][], Double>(r, pi_bg);
                
                if (config.print_PI)
                	System.out.println(t+"\t"+CommonUtils.arrayToString(pi, "%.4f")+
                			(model_noise?String.format("\t%.4f",pi_bg):""));

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
                    } 
                    if (model_noise && pi_bg!=0){
                    	double r_bg_sum=0;
                    	for(int c=0; c<mixture.numConditions; c++){
                    		double b_bg_c = 0;
                    		for(int i=0;i<r_bg[c].length;i++)
                    			b_bg_c += r_bg[c][i]*counts[c][i];
                    		b_bg[c] = b_bg_c;
                    		r_bg_sum += b_bg_c;
                    	}
                    	// normalize across conditions
                        if (r_bg_sum!=0){
                            for(int c=0; c<mixture.numConditions; c++){
                            	b_bg[c] /= r_bg_sum;
                            }
                        }
                    }
                    
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
                if (model_noise && pi_bg!=0){
	            	for(int c=0; c<mixture.numConditions; c++){
	            		for(int i=0;i<r_bg[c].length;i++)
	            			r_bg[c][i] = prob_bg[c][i] * pi_bg * b_bg[c];
	            	}
                }

                //Log-likelihood calculation
                double LL =0;
                // get the background probability if a read is out of range of the event
                double baselineProb = Math.min(mixture.model.probability(mixture.model.getMax()),
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
                                j_sum += pi[j]*b[c][j]*baselineProb;
                        }
                        if (model_noise && pi_bg!=0){
                        	j_sum += r_bg[c][i];
                        }
                        if (j_sum!=0)
                            LL += Math.log(j_sum)*counts[c][i];
                    }
                }
//                if (t==31)
//                	t=t;
                // log prior
                double LP=0;
                for(int j:nzComps)
                    if (pi[j]!=0)
                        LP+=(pos_alpha[j]-currAlpha)*Math.log(pi[j]);		// positional prior and sparse prior
                //TODO: sparse prior for noise component??
                if (model_noise && pi_bg!=0){
                	LP += -effectiveAlpha * Math.log(pi_bg);
                }
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
                        state.LAP = LAP;		// LAP for penalized likelihood
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
                        if (model_noise){
                        	for (int c=0;c<mixture.numConditions;c++){
                        		state.resp_bg[c]=r_bg[c].clone();
                        	}
                            state.beta_bg = b_bg.clone();
                            state.pi_bg = pi_bg;
                        }

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
                double bestBIC = models.get(0).BIC(totalBaseCount, hasPP);
                for (int i=1;i<models.size();i++){
                    double bic = models.get(i).BIC(totalBaseCount, hasPP);
//                   System.out.println(String.format("%.3f\t%.3f\t", bestBIC, bic)+models.get(i).toString());
                    if (bestBIC <= bic){
                        best = i;
                        bestBIC = bic;
                    }
                }
                // copy the best model state back to memory
                EM_State bestModel = models.get(best);
//                r=new double[mixture.numConditions][][];
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
                if (model_noise){
                	for (int c=0;c<mixture.numConditions;c++){
                		r_bg[c]=bestModel.resp_bg[c];
                	}
                    b_bg = bestModel.beta_bg;
                    pi_bg = bestModel.pi_bg;
                }
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
//                if (model_noise){			// Only normalize with all real components
//                    for(int i=0;i<numBases;i++){
//                    	totalResp[i] += r_bg[c][i];
//                    }
//            	}
                
                // normalize
                for(int j:nzComps){
                    int[] baseIdx = c2b[c][j];
                    for(int i=0;i<baseIdx.length;i++)
                        if (totalResp[baseIdx[i]]>0)
                            rc[j][i] = rc[j][i]/totalResp[baseIdx[i]];
                }
//                if (model_noise){
//                    for(int i=0;i<numBases;i++){
//                    	r_bg[c][i] = r_bg[c][i]/totalResp[i];
//                    }
//            	}
            }

            mixture.log(4, "EM_MAP(): "+"\tt="+t+"\t"+
                        String.format("%.6f",LAP)+"\t("+(int)nonZeroComponentNum+" events)");

            return new Pair<double[][][], Double>(r, pi_bg);
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

        // matched EM resulted binding components with the kmer prior
        //TODO: maybe we should search closest (<k/4) pp_kmer because some position may end up out-competed by nearby positions
        private void setEventKmerGroup(KmerGroup[] pp_kmer, int startPos, String seq){
        	for (BindingComponent b:components){
        		int kIdx = b.getLocation().getLocation()-startPos;
        		b.setKmerGroup(pp_kmer[kIdx]);
        		if (seq!=null){
	        		int left = kIdx - config.k_win/2;
	        		int right = kIdx + config.k_win/2;
	        		if (left<0){
	        			left = 0;
	        		}
	        		if (right+1>seq.length()){
	        			right = seq.length()-1;
	        		}
	        		// if left and right are adjusted, the seq length is not config.k_win+1 anymore
	        		// Then we can not assume the middle is binding position
	        		String bs = seq.substring(left,right+1);
	        		if (b.getKmerGroup()!=null){
	        			if (!bs.contains(b.getKmerGroup().getBestKmer().getKmerString())){
	        				bs = SequenceUtils.reverseComplement(bs);
	        				b.setKmerStrand('-');
	        			}
	        			if (!bs.contains(b.getKmerGroup().getBestKmer().getKmerString())){
	        				b.setKmerGroup(null);
	        				b.setKmerStrand('0');
	        			}
	        		}
	    			b.setBoundSequence(bs);
        		}
        	}
        }
        
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

        // this method is used in Giorgos's PoolEM method
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
            double[][] resp_bg;
            double[]beta_bg;
            double pi_bg;

            double LAP;
            double numComponent;
            EM_State(int numCond){
                resp=new double[numCond][][];
                beta=new double[numCond][];
                resp_bg=new double[numCond][];
                beta_bg=new double[numCond];
            }
            // BIC=LAP-#param/2*ln(n)
            // # param: GPS: Each component has 2 parameters, mixing prob and position, thus "*2";
            // "-1" comes from the fact that total mix prob sum to 1.
            // for multi-condition, # of beta variables is (mixture.numConditions-1)*numComponents
            // n: is the number of data point, i.e. the base positions of reads.

            // BIC_GEM
            // Each component has 3 parameters, position prior, mixing prob and position, thus "*3";
            double BIC(double n, boolean hasPP){
            	double num_parameters = numComponent*(hasPP?3:2)-1 + (mixture.numConditions-1)*numComponent;
            	if (config.noise_distribution!=0)
            		num_parameters += mixture.numConditions; //# of beta_bg is (mixture.numConditions-1), #pi_bg is 1
            	if (config.bic)
            		return LAP - num_parameters/2*Math.log(n);
            	else
            		return LAP - num_parameters; // AIC
            }
            public String toString(){
                return String.format("%.3f\t%.0f", LAP, numComponent);
            }
        }
    }
    
    class KmerPP implements Comparable<KmerPP>{
    	Point coor;
    	KmerGroup kmerMatches;
    	double pp;
    	public KmerPP(Point coor, KmerGroup matches, double pp) {
			this.coor = coor;
			this.kmerMatches = matches;
			this.pp = pp;
		}
		public String toString(){
    		return String.format("%s\t%d\t%.0f\t%.3f\t%.3f", coor.getLocationString(), kmerMatches.getKmers().size(), kmerMatches.getGroupHitCount(),kmerMatches.getWeightedKmerStrength(), pp);
    	}
		public int compareTo(KmerPP h) {
			return coor.compareTo(h.coor);
		}
    }

 }
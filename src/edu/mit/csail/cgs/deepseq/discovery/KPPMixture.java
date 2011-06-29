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

import net.sf.samtools.util.SequenceUtil;

import cern.jet.random.Poisson;
import cern.jet.random.Binomial;
import cern.jet.random.engine.DRand;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.*;
import edu.mit.csail.cgs.deepseq.features.*;
import edu.mit.csail.cgs.deepseq.multicond.MultiIndependentMixtureCounts;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.deepseq.utilities.ReadCache;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.ewok.verbs.motifs.Kmer;
import edu.mit.csail.cgs.ewok.verbs.motifs.KmerEngine;
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
import edu.mit.csail.cgs.utils.strings.multipattern.AhoCorasick;


class KPPMixture extends MultiConditionFeatureFinder {
	private final char[] LETTERS = {'A','C','G','T'};
	private final int MAXLETTERVAL = Math.max(Math.max(Math.max('A','C'),Math.max('T','G')),
            Math.max(Math.max('a','c'),Math.max('t','g'))) + 1;
	
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
	private boolean controlDataExist = false;
	// Ratio of each pair of IP/Ctrl for all conditions
	private boolean hasIpCtrlRatio = false;
	private double[] ratio_total;
	private double[] ratio_non_specific_total;
	
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

	/** Kmer motif engine
	 **/
	private KmerEngine kEngine;
	private boolean kmerPreDefined = false;
	
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
		model.printToFile(outName+"_0_Read_distribution.txt");
		allModels.put(outName+"_0", model);

    	/* *********************************
    	 * Load Kmer list
    	 ***********************************/
    	String kmerFile = Args.parseString(args, "kf", null);
    	if (kmerFile!=null){
    		kmerPreDefined = true;
			File kFile = new File(kmerFile);
			if(kFile.isFile()){
				try {
					ArrayList<Kmer> kmers = new ArrayList<Kmer>(); 
					BufferedReader reader = new BufferedReader(new FileReader(kFile.getName()));
			        String line;
			        while ((line = reader.readLine()) != null) {
			            line = line.trim();
			            String[] words = line.split("\\s+");
			            try {
				            Kmer kmer = new Kmer(words[0], Integer.parseInt(words[1]));
				            kmers.add(kmer);
		            	}
		            	catch (NumberFormatException nfe) {	// ignore if not a number, such as header
		            		continue;
		            	}
			        }
			        if (!kmers.isEmpty()){
			        	kEngine = new KmerEngine(kmers, outName);
			        }
				} catch (FileNotFoundException e) {
					e.printStackTrace();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
    	}
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
        config.parseArgs(args);    
     // if mappable_genome_length is not provided, compute as 0.8 of total genome size
        if (config.mappable_genome_length<0){		
	        config.mappable_genome_length = 0.8 * gen.getGenomeSize();
	        System.out.println(String.format("\nMappable Genome Length is %,d.", (long)config.mappable_genome_length));
        }
        
    	if(config.second_lambda_region_width < config.first_lambda_region_width) {
    		System.err.println("\nThe first control region width (w2) has to be more than " + config.first_lambda_region_width + " bp.");
    		System.exit(-1);
    	}    	
    	if(config.third_lambda_region_width < config.second_lambda_region_width) {
    		System.err.println("\nThe second control region width (w3) has to be more than " + config.second_lambda_region_width + " bp.");
    		System.exit(-1);
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
    	loadChIPSeqData(subsetRegions);
		
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
     		// ip/ctrl ratio by regression on non-enriched regions
			if (config.ip_ctrl_ratio==-1){
     			setRegions(selectEnrichedRegions(subsetRegions, true));
     			ArrayList<Region> temp = (ArrayList<Region>)restrictRegions.clone();
     			temp.addAll(excludedRegions);
    			calcIpCtrlRatio(mergeRegions(temp, false));
    			if(controlDataExist) {
    				for(int t = 0; t < numConditions; t++)
    					System.out.println(String.format("For condition %s, IP/Control = %.2f", conditionNames.get(t), ratio_non_specific_total[t]));
    			}
    		}
			setRegions(selectEnrichedRegions(subsetRegions, true));
		} else{
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
		
		if (kEngine!=null && kEngine.isInitialized()){
			System.out.println("\nRunning EM with k-mer positional prior ...\n");
		}
		
		// create threads and run EM algorithms
		// the results are put into compFeatures
        Thread[] threads = new Thread[maxThreads];
        log(1,String.format("Creating %d threads ...", maxThreads));
        int regionsPerThread = restrictRegions.size()/threads.length;
        for (int i = 0 ; i < threads.length; i++) {
            ArrayList<Region> threadRegions = new ArrayList<Region>();
            int nextStartIndex = (i+1)*regionsPerThread;
            if (i==threads.length-1)		// last thread
            	nextStartIndex = restrictRegions.size();
            // get the regions as same chrom as possible, to minimize chrom sequence that each thread need to cache
            for (int j = i*regionsPerThread; j < nextStartIndex; j++) {
                threadRegions.add(restrictRegions.get(j));
            }
            Thread t = new Thread(new GPS2Thread(threadRegions,
            									processRegionCount,
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
        
        log(3,String.format("%d threads have finished running", maxThreads));

        // print out kmer hit list
        if (config.kmer_print_hits){
	        Collections.sort(allKmerHits);
	        StringBuilder sb = new StringBuilder();
	        sb.append("Position\tKmer\tCount\tWeight\tpp\n");
	        for (KmerPP h:allKmerHits)
	        	sb.append(h.toString()).append("\n");
	        CommonUtils.writeFile(outName+"_Kmer_Hits.txt", sb.toString());
        }
		/* ********************************************************
		 * refine the specific regions that contain binding events
		 * ********************************************************/
		if (config.refine_regions){
			Collections.sort(compFeatures);
			ArrayList<Region> refinedRegions = new ArrayList<Region>();
			for (ComponentFeature cf:compFeatures){
				refinedRegions.add(cf.getPosition().expand(0));
			}
			// expand with modelRange, and merge overlapped regions ==> Refined enriched regions
			this.restrictRegions=mergeRegions(refinedRegions, true);
		}

		log(2, "Finish predicting events: "+CommonUtils.timeElapsed(tic)+"\n");

		/* ********************************************************
		 * post EM processing
		 * ********************************************************/
		postEMProcessing(compFeatures);
		
		log(1, "Finish prediction: "+CommonUtils.timeElapsed(tic)+"\n");

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
		double medianStrength = compFeatures.get(compFeatures.size()/2).getTotalEventStrength();
//		System.out.println(String.format("Median event strength = %.1f\n",medianStrength));
		
		// only do this calculation at the first round, then throw out read data in non-specific regions (when we have control data for stat test)
		if (!hasIpCtrlRatio){		
			ArrayList<Region> exRegions = new ArrayList<Region>();
			for(int i = 0; i < compFeatures.size(); i++) {
				exRegions.add(compFeatures.get(i).getPosition().expand(modelRange));
			}
			exRegions.addAll(excludedRegions);	// also excluding the excluded regions (user specified + un-enriched)
			
			calcIpCtrlRatio(mergeRegions(exRegions, false));
			if(controlDataExist) {
				for(int c = 0; c < numConditions; c++)
					System.out.println(String.format("\nScaling condition %s, IP/Control = %.2f", conditionNames.get(c), ratio_non_specific_total[c]));
				System.out.println();
			}
			hasIpCtrlRatio = true;
			
			// delete read data in un-enriched region, we don't need them for binomial test
			// but if no control data, we need whole genome data to estimate lambda for Poisson test
			if(controlDataExist) {
				for(int c = 0; c < numConditions; c++) {
					caches.get(c).car().deleteUnenrichedReadData(restrictRegions);				
					caches.get(c).cdr().deleteUnenrichedReadData(restrictRegions);
				}
			}
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
			if (config.filterEvents && cf.getTotalEventStrength()>medianStrength){
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
	
	                    poisson.setMean(config.minFoldChange * Math.max(scaledControlCount, totalIPCount[cond] * modelWidth / config.mappable_genome_length  ));
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
		
		// find q-value cutoff by k-mer occurence
		setQValueCutoff(compFeatures);
	}//end of evaluateConfidence method

	private void setQValueCutoff(List<ComponentFeature>compFeatures) {
		int N = compFeatures.size();
		int kmerHitCount[] = new int[N];		// the number of kmerHits up to event #i
		for (int i=0;i<N;i++){
			ComponentFeature cf = compFeatures.get(i);
			Kmer kmer = cf.getKmer();
			int inc = (kmer==null||kmer.getSeqHitCount()<=1)?0:1;
			if (i==0)
				kmerHitCount[0] = inc;
			else
				kmerHitCount[i] = inc + kmerHitCount[i-1];
		}
		for (int i=0;i<N-1;i++){
//			System.out.print(String.format("%d\t%d\t%d\t%d\t", kmerHitCount[i], N, kmerHitCount[N-1], i+1));
			double hgp = 1-StatUtil.hyperGeometricCDF_cache(kmerHitCount[i], N, kmerHitCount[N-1], i+1);
			compFeatures.get(i).setEnrichedKmerHGPLog10(-Math.log10(hgp));
//			System.out.println(String.format("%.1f\t%.3f", compFeatures.get(i).getEnrichedKmerHGPLog10(), 1-(kmerHitCount[i]/(i+1.0))));
		}		
	}

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
        Thread[] threads = new Thread[maxThreads];
        log(1,String.format("Running EM on control data: creating %d threads", maxThreads));
        int regionsPerThread = regions.size()/threads.length;
        for (int i = 0 ; i < threads.length; i++) {
            ArrayList<Region> threadRegions = new ArrayList<Region>();
            int nextStartIndex = (i+1)*regionsPerThread;
            if (i==threads.length-1)		// last thread
            	nextStartIndex = regions.size();
            // get the regions as same chrom as possible, to minimize chrom sequence that each thread need to cache
            for (int j = i*regionsPerThread; j < nextStartIndex; j++) {
                threadRegions.add(regions.get(j));
            }
            Thread t = new Thread(new GPS2Thread(threadRegions,
            									processRegionCount,
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
        
        log(1,String.format("%d threads have finished running", maxThreads));
        
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
		
		// if IP/Control enrichment ratios are lower than cutoff for all 500bp sliding windows in all conditions, 
		// skip this region, and record it in excludedRegions
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
	private void loadChIPSeqData(ArrayList<Region> subsetRegions){
		this.caches = new ArrayList<Pair<ReadCache, ReadCache>>();
		this.numConditions = experiments.size();
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
        GPS2Thread t = new GPS2Thread(null, null, out, null, this, constants,config, true);
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

	/* 
	 * Calc the ratio of IP vs Control channel, exlcluding the specific regions
	 */
	private void calcIpCtrlRatio(ArrayList<Region> specificRegions) {
		// linear regression to get the IP/control ratio
		// now we do not require whole genome data, because partial data could be run on 1 chrom, still enough data to esitmates
		if(controlDataExist) {
			if (config.ip_ctrl_ratio==-1){		// regression using non-specific regions
				ratio_non_specific_total = new double[numConditions];
				for(int t = 0; t < numConditions; t++){
					ratio_non_specific_total[t] = getSlope(t, t, "IP/CTRL", specificRegions);
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

	// merge the overlapped regions
	// if "toExpandRegion"=true, expand each region on both side to leave enough space,
	// to include every potential reads, then merge
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
				List<StrandedBase> allBases = new ArrayList<StrandedBase>();
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
							// if the average read count per modelWidth is less than config.sparseness/2, find sparse point to further split
							rs.add(r);
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
                    if (!enriched){	// remove this region if it is not enriched in any condition
                        toRemove.add(r);
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
	 * Apply a Possion filter for duplicate reads
	 */
	private void applyPoissonFilter(boolean printFilterMsg){
		// Filtering/reset bases
		if (printFilterMsg)
			System.out.println("\nApply Poisson filter for duplicate reads.");
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

	
	double updateBindingModel(int left, int right){
		if (signalFeatures.size()<config.min_event_count){
			System.err.println("Warning: The read distribution is not updated, too few ("+signalFeatures.size()+"<"+config.min_event_count+") significant events.");
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
            count += e.car().countStrandedBases(r,'+');
            count += e.car().countStrandedBases(r,'-');
		}
		return count;
	}
	private float countIpReads(Region r, int cond){
        return caches.get(cond).car().countStrandedBases(r,'+') + 
            caches.get(cond).car().countStrandedBases(r,'+');
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
            caches.get(cond).cdr().countStrandedBases(r,'+');
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
	
	/**
	 * Some old methods, not used any more
	 */
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
     * Initalize the kmer engine
     * This is called only once for initial setup.
     * It compact the cached sequence data and build the kmerEngine
     */
    public void initKmerEngine(){
    	if (config.k==-1 && config.k_min==-1)
    		return;
		
    	System.out.println("Loadning genome sequences ...");
		kEngine = new KmerEngine(gen, config.cache_genome);
		long tic = System.currentTimeMillis();
		
		// refine restrictRegions, to reduce memory and run time
		ArrayList<Feature> events = new ArrayList<Feature>();
		events.addAll(signalFeatures);
		events.addAll(insignificantFeatures);
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
		// setup lightweight genome cache
		if (!kmerPreDefined){
			ArrayList<Region> expandedRegions = new ArrayList<Region>();
			for (Region r: restrictRegions){
				expandedRegions.add(r.expand(config.k_win+modelRange, config.k_shift+config.k_win+modelRange));
			}
			expandedRegions = this.mergeRegions(expandedRegions, false);
			int totalLength=0;
			for (Region r: expandedRegions){
				totalLength+=r.getWidth();
			}
			// get negative regions
			ArrayList<Region> negativeRegions = new ArrayList<Region>();
// Random genome-wide
//			cern.jet.random.engine.RandomEngine randomEngine = new cern.jet.random.engine.MersenneTwister();
//			for (Feature f:signalFeatures){
//				String chr = f.getPeak().getChrom();
//				int length = gen.getChromLength(chr)-config.k_win-1;
//				// for each event, get random negative regions from the same chromosome
//				for (int i=0;i<config.negative_ratio;i++){
//					double rand = randomEngine.nextDouble();
//					int start = (int)(rand*length);
//					Region r = new Region(gen, chr, start, start+config.k_win);
//					negativeRegions.add(r);
//				}
//			}
//			negativeRegions = Region.filterOverlapRegions(negativeRegions, expandedRegions);
			// In proximal regions, but excluding binding regions
			for (Feature f:signalFeatures){
				String chr = f.getPeak().getChrom();
				int length = gen.getChromLength(chr)-config.k_win-1;
				// for each event, get random negative regions from the same chromosome
				for (int i=1;i<=config.negative_ratio;i++){
					int start = f.getPeak().getLocation()+config.k_neg_dist*i;
					if ( start+config.k_win>=length)
						continue;
					Region r = new Region(gen, chr, start, start+config.k_win);
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
			
			kEngine.setupRegionCache(expandedRegions, negativeRegions);
			if (config.bmverbose>1)
				System.out.println("Compact genome sequence cache to " + totalLength + " bps, "+CommonUtils.timeElapsed(tic));
		}
		
		buildEngine();
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
		if (config.kmer_use_insig){
			for(Feature f : insignificantFeatures){
				if(count++>config.k_seqs)
					break;
				ComponentFeature cf = (ComponentFeature)f;
				events.add(cf);
			}
		}
		if (config.kmer_use_filtered){			
			for(Feature f : filteredFeatures){
				if(count++>config.k_seqs)
					break;
				ComponentFeature cf = (ComponentFeature)f;
				events.add(cf);
			}
		}
		return events;
    }   
    
    private ArrayList<Point> getEventPoints(){
    	ArrayList<Point> points = new ArrayList<Point>();
    	ArrayList<ComponentFeature> events = getEvents();
		for(ComponentFeature cf : events){
			points.add(cf.getPeak());
		}
		return points;
    }
    
    /**
     * Build k-mer engine.<br>
     * Select enriched k-mers, align k-mers, setup k-mer engine
     */
    public void buildEngine(){
    	ArrayList<Point> points = getEventPoints();
		// compare different values of k to select most enriched k value
		int k=0;
		double max_value=0;
		if (config.k_min!=-1){
			int eventCounts[] = new int[config.k_max-config.k_min+1];
			for (int i=0;i<eventCounts.length;i++){
				ArrayList<Kmer> kms = kEngine.selectEnrichedKmers(i+config.k_min, points, config.k_win, config.hgp, config.k_fold, outName+"_OK_win"+ (config.k_win));
				if (kms.isEmpty())
					eventCounts[i] = 0;
				else{
					Collections.sort(kms, new Comparator<Kmer>(){
					    public int compare(Kmer o1, Kmer o2) {
					    		return o1.compareByHGP(o2);
					    }
					});
					eventCounts[i] = kms.get(0).getSeqHitCount();
				}
			}
			for (int i=0;i<eventCounts.length-1;i++){
				if (eventCounts[i]==0)
					continue;
				if (eventCounts[0]/eventCounts[i]>=2)
					break;
				
				double percent = (eventCounts[i]-eventCounts[i+1])*100.0/(double)eventCounts[i];
				if (percent>max_value){
					k=i+config.k_min;
					max_value = percent;
				}
				System.out.println(String.format("k=%d, decrease=%.2f%%\teventCount=%d", i+config.k_min, percent, eventCounts[i]));
			}
			System.out.println(String.format("selected k=%d, max decrease=%.2f%%", k, max_value));
			config.k = k;
		}
//		if (config.k_shift==-1)				// default k_shift value to be k
//			config.k_shift = config.k;
		
		ArrayList<Kmer> kmers = kEngine.selectEnrichedKmers(config.k, points, config.k_win, config.hgp, config.k_fold, outName+"_OK_win"+ (config.k_win));
		if (config.align_overlap_kmer)
			kmers = alignOverlappedKmers(kmers, getEvents());
		kEngine.updateEngine(kmers, outName+"_OK_win"+ (config.k_win), false);		
    }
    
    /** 
     * Align overlapped k-mers using the positive sequences containing these k-mers<br>
     * Use the seed kmer position as ancher, specify the position of kmers and seqs as the position relative to the seed kmer start<br>
     * seq_seed = - seed_seq = - (index of start of seed kmer in seq)<br>
     * kmer_seed is kmer position relative to seed k-mer, also defined as k-mer shift . For the seed kmer or its mismatch, kmer_seed = 0<br>
     * kmer_seed = kmer_seq + seq_seed<br>
     * seq_seed = kmer_seed - kmer_seq<br>
     * Thus, we get a consensus (most frequent) of "kmer_seed" using the kmer's occurrences in aligned sequences<br>
     * Then, we align the unaligned sequences (also containing this k-mer) using the consensus value of kmer_seed.
     * Aligning k-mers is done through such alternating aligning kmer and aligning sequences.
     * @param kmers
     * @param events
     * @return
     */
	private ArrayList<Kmer> alignOverlappedKmers(ArrayList<Kmer> kmers, ArrayList<ComponentFeature> events){
		if (kmers.size()==0)
			return kmers;
		String[] seqs_old = kEngine.getPositiveSeqs();
		ArrayList<Kmer> allAlignedKmers = new ArrayList<Kmer>();
		// build the kmer search tree
		AhoCorasick oks = new AhoCorasick();
		HashMap<String, Kmer> str2kmer = new HashMap<String, Kmer>();
		for (Kmer km: kmers){
			str2kmer.put(km.getKmerString(), km);
			oks.add(km.getKmerString().getBytes(), km);
	    }
		oks.prepare();
		
		// index kmer->seq, seq->kmer
		ArrayList<HashSet<Kmer>> seq2kmer = new ArrayList<HashSet<Kmer>>();
		HashMap <Kmer, HashSet<Integer>> kmer2seq = new HashMap <Kmer, HashSet<Integer>>();
		for (int i=0;i<seqs_old.length;i++){
			String seq = seqs_old[i];
			HashSet<Kmer> results = KmerEngine.queryTree (seq, oks);
			if (results.isEmpty()){
				seq2kmer.add(null);
			}
			else{
				for (Kmer km: results){		
					if (!kmer2seq.containsKey(km)){
						kmer2seq.put(km, new HashSet<Integer>());
					}
					kmer2seq.get(km).add(i);
				}
				seq2kmer.add(results);
			}
		}
		
		// specify the position of kmers and seqs as the position relative to the seed kmer start
		Collections.sort(kmers, new Comparator<Kmer>(){
		    public int compare(Kmer o1, Kmer o2) {
		    		return o1.compareByHGP(o2);
		    }
		});

		// cluster and align
		final int STRAND = 1000;		// extra bp add to indicate negative strand match of kmer
		final int UNALIGNED = 999;
		int clusterID = 0;
		StringBuilder alignedKmer_sb = new StringBuilder();
		ArrayList<KmerCluster> clusters = new ArrayList<KmerCluster>();
		
		while(!kmers.isEmpty()){
			KmerCluster cluster = new KmerCluster();
			cluster.clusterId = clusterID;
			clusters.add(cluster);
			
			String[] seqs = seqs_old.clone();
			int posSeqs[] = new int[seqs.length]; 		// the position of sequences
			String seqAlignRefs[] = new String[seqs.length];
			boolean isPlusStrands[] = new boolean[seqs.length];
			// init posSeqs, so each new kmer cluter align with all the sequences 
			for (int i=0;i<posSeqs.length;i++){
				posSeqs[i] = UNALIGNED;
				isPlusStrands[i] = true;
			}		
			ArrayList<Kmer> alignedKmers = new ArrayList<Kmer>();
			
			/** get seed kmer and its mismatch k-mers, align sequences */
			Kmer seed = kmers.get(0);
			cluster.seedKmer = seed;
			ArrayList<Kmer> seedFamily = getSeedKmerFamily(kmers, seed);
			// align the containing sequences
			for (Kmer km:seedFamily){
				alignSequences(km, kmer2seq.get(km), seqs, posSeqs, isPlusStrands, seqAlignRefs);
			}
			alignedKmers.addAll(seedFamily);
			kmers.removeAll(seedFamily);
			boolean seedFamilyMismatchUsed = false;
			
			/** greedly align kmers until no kmer can be aligned */
			WeightMatrix wm = null;
			String pfmStr = "";

			// newly aligned kmers
			ArrayList<Kmer> aligned_new = new ArrayList<Kmer>();
			aligned_new.addAll(seedFamily);
			
			while(!kmers.isEmpty()){
				for (Kmer km:kmers){
					int kmer_seed = 0;			// the shift of this kmer w.r.t. seed kmer
					HashSet<Integer> hits = kmer2seq.get(km);
					
					/** use k-mers already in the aligned sequences, find the consensus k-mer position */
					ArrayList<Integer> posKmer = new ArrayList<Integer>(); // the occurences of this kmer w.r.t. seed kmer
					for (int seqId:hits){
						if (posSeqs[seqId] != UNALIGNED){		// aligned seqs
							int pos = seqs[seqId].indexOf(km.getKmerString());
							if (pos==-1){
								pos = seqs[seqId].indexOf(km.getKmerRC());
								if (pos!=-1){
									posKmer.add(posSeqs[seqId]+pos+STRAND);		// STRAND --> match of KmerRC
								}
							}
							else{
								posKmer.add(posSeqs[seqId]+pos);
							}
						}
					}
					
					// find the most frequent kmerPos
					if (posKmer.size()>km.getSeqHitCount()/2){		// the kmer positions aligned must be at least half of total kmer count
						Pair<int[], int[]> sorted = StatUtil.sortByOccurences(posKmer);
						int counts[] = sorted.cdr();
						int posSorted[] = sorted.car();
						int maxCount = counts[counts.length-1];
						if (maxCount<Math.round(posKmer.size()*config.kmer_freq_pos_ratio))// the most freq position count must be large enough
							continue;
						ArrayList<Integer> maxPos = new ArrayList<Integer>();
						ArrayList<Boolean> isPositive = new ArrayList<Boolean>();
						for (int i=counts.length-1;i>=0;i--){
							if (counts[i]==maxCount){
								int p = posSorted[i];
								if (p>STRAND/2){			// if match of KmerRC
									maxPos.add(p-STRAND);
									isPositive.add(false);
								}
								else{
									maxPos.add(p);
									isPositive.add(true);
								}
							}	
							else		// do not need to count for non-max elements
								break;
						}
						if (maxPos.size()>1){		// if tie with 1+ positions, get the one closest to seed kmer
							int min = Integer.MAX_VALUE;
							int minIdx = 0;
							for (int i=0;i<maxPos.size();i++){
								int distance = Math.abs(maxPos.get(i));
								if (distance<min){
									minIdx = i;
									min = distance;
								}
							}
							kmer_seed = maxPos.get(minIdx);
							if (Math.abs(kmer_seed)>config.k_shift)		// if the kmer is too far away
								continue;
							if (!isPositive.get(minIdx)){	// if match of KmerRC
								km.RC();
							}
						}
						else{
							kmer_seed = maxPos.get(0);
							if (Math.abs(kmer_seed)>config.k_shift)		// if the kmer is too far away
								continue;
							if (!isPositive.get(0))	// if match of KmerRC
								km.RC();
						}
						km.setAlignString(maxCount+"/"+posKmer.size());
					}
					else {
						continue;	// continue for next kmer
					}
					km.setShift(kmer_seed);
					aligned_new.add(km);
					
					/** use kmer shift to align the containing sequences */
					alignSequences(km, hits, seqs, posSeqs, isPlusStrands, seqAlignRefs);
					
				} //for (Kmer km:kmers)
				
				kmers.removeAll(aligned_new);
				alignedKmers.addAll(aligned_new);
				Collections.sort(aligned_new);
				
				/** use aligned k-mer to align mismatch k-mers */
				if (config.use_kmer_mismatch){
					ArrayList<Kmer> mmaligned = new ArrayList<Kmer>();			// kmers aligned by using mismatch
					for (Kmer km: kmers){
						String seq = km.getKmerString();
						for (Kmer akm: aligned_new){
							if (this.mismatch(akm.getKmerString(), seq)==1){
								km.setShift(akm.getShift());
								km.setKmerStartOffset(akm.getKmerStartOffset());
								km.setAlignString("MM:"+akm.getKmerString());
								mmaligned.add(km);
								break;
							}
							if (this.mismatch(akm.getKmerRC(), seq)==1){
								km.RC();
								km.setShift(akm.getShift());
								km.setKmerStartOffset(akm.getKmerStartOffset());
								km.setAlignString("MM:"+akm.getKmerString());
								mmaligned.add(km);
								break;
							}
						}
					}
					// use kmer shift to align the containing sequences 
					for (Kmer km: mmaligned){
						alignSequences(km, kmer2seq.get(km), seqs, posSeqs, isPlusStrands, seqAlignRefs);
					}
					kmers.removeAll(mmaligned);
					aligned_new.clear();				// clear aligned kmers that are used for mismatch search
					aligned_new.addAll(mmaligned);
					alignedKmers.addAll(mmaligned);
				} //if (config.use_kmer_mismatch)
				
				/** build PWM to continue grow cluster */
				int alignedSeqCount=0;
				for (int i=0;i<posSeqs.length;i++){
					if (posSeqs[i] != UNALIGNED)
						alignedSeqCount++;
				}
				if (alignedSeqCount>=config.kmer_cluster_seq_count){
					
	//				boolean pwm_aligned=false;
					double[][] pfm = new double[config.k_win+1][MAXLETTERVAL];
					ArrayList<String> alignedSeqs = new ArrayList<String>();
					for (int i=0;i<posSeqs.length;i++){
						// get aligned sequences
						int pos = posSeqs[i];
						if (pos == UNALIGNED)
							continue;
						int kmMid_seq = -pos+config.k/2;
						if (!isPlusStrands[i])			// minus strand
							kmMid_seq = config.k_win - kmMid_seq;// - (1-config.k%2);				// adjust for odd/even value of k
						Point center = new Point(gen, events.get(i).getPeak().getChrom(), 
								events.get(i).getPeak().getLocation()-(config.k_win/2)+kmMid_seq);
						String seq = kEngine.getSequence(center.expand(config.k_win/2));
						if (!isPlusStrands[i])
							seq=SequenceUtils.reverseComplement(seq);
//						alignedSeqs.add(seq);		// for debugging
						// count base frequencies
			 			double strength = config.use_strength?events.get(i).getTotalEventStrength():1;
			    		for (int p=0;p<config.k_win+1;p++){
			    			char base = seq.charAt(p);
			    			pfm[p][base] +=strength;
			    		}	    	
			    	}
//					for (String s:alignedSeqs)
//						System.out.println(s);
					
					double[][] pwm = pfm.clone();
					for (int i=0;i<pfm.length;i++)
						pwm[i]=pfm[i].clone();
					
			    	// make the PWM
					Pair<Integer, Integer> ends = constructPWM(pwm);
					int leftIdx = ends.car();
					int rightIdx = ends.cdr();
					if (rightIdx-leftIdx+1>config.k/2){		// pwm is long enough
				    	float[][] matrix = new float[rightIdx-leftIdx+1][MAXLETTERVAL];   
				    	for(int p=leftIdx;p<=rightIdx;p++){
				    		for (int b=0;b<LETTERS.length;b++){
				    			matrix[p-leftIdx][LETTERS[b]]=(float) pwm[p][LETTERS[b]];
				    		}
				    	}
				    	float[][] pfm_trim = new float[rightIdx-leftIdx+1][MAXLETTERVAL];   
				    	for(int p=leftIdx;p<=rightIdx;p++){
				    		for (int b=0;b<LETTERS.length;b++){
				    			pfm_trim[p-leftIdx][LETTERS[b]]=(float) pfm[p][LETTERS[b]];
				    		}
				    	}
				    	cluster.pfmString = makeTRANSFAC (pfm_trim, String.format("DE %s_%d_c%d\n", outName, clusterID, alignedSeqCount));
				    	cluster.sequenceCount = alignedSeqCount;
//				    	CommonUtils.writeFile(outName+"_OK_"+clusterID+"_PFM.txt", pfmStr);
				    	
				    	wm = new WeightMatrix(matrix);
				    	cluster.wm = wm;
				    	// Check the quality of new PWM: hyper-geometric p-value test using the positive and negative sequences
				    	cluster.pwmThreshold = kEngine.estimatePwmThreshold(wm, config.wm_factor, outName);	
				    	if (cluster.pwmThreshold<0){
				    		System.out.println("alignOverlapKmers: PWM "+WeightMatrix.getMaxLetters(wm)+" is non-specific.");
				    	}	
				    	else{
				    		int count_pwm_aligned = 0;
					    	if (cluster.pwmThreshold<wm.getMaxScore()*config.wm_factor)
					    		cluster.pwmThreshold = wm.getMaxScore()*config.wm_factor;
					        WeightMatrixScorer scorer = new WeightMatrixScorer(wm);
					    	HashMap<String, Integer> pwmAlignedKmerStr = new HashMap<String, Integer>();	// kmerString -> count
					    	int pos_pwm_seed = leftIdx-(config.k_win/2-config.k/2);
					    	cluster.pos_pwm_seed = pos_pwm_seed;
					    	HashMap<String, ArrayList<Integer>> pwmKmerStr2seq = new HashMap<String, ArrayList<Integer>>();	// kmerString -> list of sequence ids
					    	for (int i=0;i<posSeqs.length;i++){
							  int pos = posSeqs[i];
							  if (pos != UNALIGNED)					// align the unaligned seqs with pwm
								continue;
					    	  String seq = seqs[i];
					    	  if (seq.length()<wm.length()-1)
					    		  continue;
					    	      	  
					          WeightMatrixScoreProfile profiler = scorer.execute(seq);
					          double maxSeqScore = Double.NEGATIVE_INFINITY;
					          int maxScoringShift = 0;
					          char maxScoringStrand = '+';
					          for (int j=0;j<profiler.length();j++){
					        	  double score = profiler.getMaxScore(j);
					        	  if (maxSeqScore<score){
					        		  maxSeqScore = score;
					        		  maxScoringShift = j;
					        		  maxScoringStrand = profiler.getMaxStrand(j);
					        	  }
					          }
					          // if a sequence pass the motif score, reset the kmer to the binding site
					          if (maxSeqScore >= cluster.pwmThreshold){
								if (maxScoringStrand =='-'){
									seqs[i] = SequenceUtils.reverseComplement(seqs[i]);
									isPlusStrands[i] = false;
									maxScoringShift = seqs[i].length()-maxScoringShift-wm.length();
									// i.e.  (seq.length()-1)-maxScoringShift-(wm.length()-1);
								}
								int start = -(pos_pwm_seed-maxScoringShift);
								int end = start+config.k;
								if (start<0 || end > seqs[i].length()){		// if the seed kmer alignment exceed the sequence, skip
									if (maxScoringStrand =='-'){	// reverse the change
										seqs[i] = SequenceUtils.reverseComplement(seqs[i]);
										isPlusStrands[i] = true;
									}
									continue;
								}
								posSeqs[i] = pos_pwm_seed-maxScoringShift;
								String kmerStr = seqs[i].substring(-posSeqs[i], end);
								if (pwmAlignedKmerStr.containsKey(kmerStr))
									pwmAlignedKmerStr.put(kmerStr, pwmAlignedKmerStr.get(kmerStr)+1);
								else
									pwmAlignedKmerStr.put(kmerStr, 1);
								if (pwmKmerStr2seq.containsKey(kmerStr))
									pwmKmerStr2seq.get(kmerStr).add(i);
								else{
									ArrayList<Integer> list = new ArrayList<Integer>();
									list.add(i);
									pwmKmerStr2seq.put(kmerStr, list);
								}
								seqAlignRefs[i] = "PWM:"+WeightMatrix.getMaxLetters(wm);
								count_pwm_aligned ++;
					          }
					        }	// each unaligned sequence
					    	for (String kmStr: pwmAlignedKmerStr.keySet()){
					    		Kmer km = null;
					    		String kmCR = SequenceUtil.reverseComplement(kmStr);
					    		if (str2kmer.containsKey(kmStr)){	// if existing k-mers
					    			km = str2kmer.get(kmStr);
					    			kmers.remove(km);
					    		}
					    		else if (str2kmer.containsKey(kmCR)){	
					    			km = str2kmer.get(kmCR);
					    			km.RC();
					    			kmers.remove(km);
					    		}
					    		else {								// new found k-mers
					    			km = new Kmer(kmStr, pwmAlignedKmerStr.get(kmStr));
					    		}
				    			km.setShift(0);
				    			km.setAlignString("PWM:"+WeightMatrix.getMaxLetters(wm));
				    			aligned_new.add(km);
				    			alignedKmers.add(km);
				    			// connect kmer back to the sequence
				    			for (int i:pwmKmerStr2seq.get(kmStr)){
				    				if (seq2kmer.get(i)==null)
				    					seq2kmer.set(i, new HashSet<Kmer>());
				    				seq2kmer.get(i).clear();
									seq2kmer.get(i).add(km);
				    			}
					    	}
					    	if (config.bmverbose>1)
					    		System.out.println("PWM "+WeightMatrix.getMaxLetters(wm)+" align "+count_pwm_aligned+" sequences and "+pwmAlignedKmerStr.size()+" k-mers.");
				    	}
					}
				}
		    	
				if (aligned_new.isEmpty()){		// if no new mismatch or pwm aligned kmer, then no new aligned sequence
					break;			// stop this cluster, start a new one
				}
			} // greedily growing cluster
	        
			/** use all aligned sequences to find expected binding sites */
	    	// average all the binding positions to decide the expected binding position
	    	// weighted by strength if "use_strength" is true
			StringBuilder sb = new StringBuilder();
			double sum_offsetXstrength = 0;
	    	double sum_strength = 0;
	    	int leftmost = Integer.MAX_VALUE;
			for (int i=0;i<posSeqs.length;i++){
				int pos = posSeqs[i];
				if (pos < leftmost )
					leftmost = pos;					
			}
			for (int i=0;i<posSeqs.length;i++){
				int pos = posSeqs[i];
				if (pos == UNALIGNED)
					continue;
				if (config.bmverbose>1)
					sb.append(CommonUtils.padding(-leftmost+pos, '-')+seqs[i]+"\t\t"+seqAlignRefs[i]+"\n");
	 			double strength = config.use_strength?events.get(i).getTotalEventStrength():1;
    			sum_offsetXstrength += strength*(config.k_win/2+pos);
        		sum_strength += strength;
	    	}
			if (config.bmverbose>1)
				CommonUtils.writeFile(outName+"_seqs_aligned_"+seed.getKmerString()+".txt", sb.toString());
	    	int bPos=StatUtil.round(sum_offsetXstrength/sum_strength);		// mean BS position relative to seed k-mer start
	    	cluster.bindingPosition = bPos - cluster.pos_pwm_seed;
	    	
	    	alignedKmer_sb.append("Cluster #"+clusterID+"\n");
	    	int leftmost_km = Integer.MAX_VALUE;
			for (Kmer km: alignedKmers){
				km.setKmerStartOffset(km.getShift()-bPos);
				if (km.getKmerStartOffset()<leftmost_km)
					leftmost_km = km.getKmerStartOffset();
			}	    	
			for (Kmer km: alignedKmers){
				alignedKmer_sb.append(km.getKmerStartOffset()+"\t"+CommonUtils.padding(-leftmost_km+km.getKmerStartOffset(), '-')+km.toOverlapString()+"\t"+km.getAlignString()+"\n");
			}
			cluster.kmerCount = alignedKmers.size();
			allAlignedKmers.addAll(alignedKmers);
			clusterID++;
		} // each cluster
		
		CommonUtils.writeFile(outName+"_OK_aligned.txt", alignedKmer_sb.toString());
		
		// take the best kmer from each sequence
		TreeSet<Kmer> bestKmers = new TreeSet<Kmer>();
		StringBuilder best_sb = new StringBuilder();
		int leftmost_km = Integer.MAX_VALUE;
		for (int i=0;i<seq2kmer.size();i++){
			ArrayList<Kmer> kms = new ArrayList<Kmer>();
			if (seq2kmer.get(i)==null)
				continue;
			kms.addAll(seq2kmer.get(i));
			// specify the position of kmers and seqs as the position relative to the seed kmer start
			Collections.sort(kms);
			Kmer km = kms.get(0);
			bestKmers.add(km);
			if (km.getKmerStartOffset()<leftmost_km)
				leftmost_km = km.getKmerStartOffset();
		}
		for (Kmer km:bestKmers){
			best_sb.append(km.getKmerStartOffset()+"\t"+CommonUtils.padding(-leftmost_km+km.getKmerStartOffset(), '-')+km.toOverlapString()+"\t"+km.getAlignString()+"\n");
		}
		CommonUtils.writeFile(outName+"_OK_best_aligned.txt", best_sb.toString());
		
		// output cluster information, PFM, and PWM
		StringBuilder pfm_sb = new StringBuilder();		
		for (KmerCluster c:clusters){
    		WeightMatrix wm = c.wm;
    		if (wm==null || c.sequenceCount<config.kmer_cluster_seq_count)
    			continue;
    		System.out.println(String.format("-------------------------------\n%s k-mer cluster #%d, from %d k-mers, %d binding events.", outName, c.clusterId, c.kmerCount, c.sequenceCount));
    		int pos = c.bindingPosition;
    		if (pos>0)
    			System.out.println(CommonUtils.padding(pos, ' ')+"|\n"+ WeightMatrix.printMatrixLetters(wm));
    		else
    			System.out.println(WeightMatrix.printMatrixLetters(wm));
    		System.out.println(String.format("PWM threshold: %.2f/%.2f", c.pwmThreshold, c.wm.getMaxScore()));
			pfm_sb.append(c.pfmString);
		}
		CommonUtils.writeFile(outName+"_OK_PFM.txt", pfm_sb.toString());
//		return allAlignedKmers;
		ArrayList<Kmer> result = new ArrayList<Kmer>();
		result.addAll(bestKmers);
		return result;
	}

    private ArrayList<Kmer> getSeedKmerFamily(ArrayList<Kmer> kmers, Kmer seed) {
    	ArrayList<Kmer> family = new ArrayList<Kmer>();
    	String seedKmerStr = seed.getKmerString();
    	String seedKmerRC = seed.getKmerRC();
    	for (Kmer kmer: kmers){
	    	if (kmer.hasString(seedKmerStr) || mismatch(seedKmerStr, kmer.getKmerString())<=1+config.k*0.1){
	    		kmer.setShift(0);
	    		family.add(kmer);
	    		kmer.setAlignString(seedKmerStr);
	    	}
	    	else if (kmer.hasString(seedKmerRC)||mismatch(seedKmerRC, kmer.getKmerString())<=1+config.k*0.1){
	    		kmer.setShift(0);
	    		kmer.RC();
	    		family.add(kmer);
	    		kmer.setAlignString(seedKmerStr);
	    	}
    	}
		return family;
	}
    
	/** use kmer shift to align the containing sequences */
    private void alignSequences(Kmer km, HashSet<Integer> seqIds, String[]seqs, int[] posSeqs, boolean[] isPlusStrands, String[] seqAlignRefs){
		final int UNALIGNED = 999;
    	for (int seqId:seqIds){
			if (posSeqs[seqId] == UNALIGNED){		// not aligned yet
				// align using kmer positions
				int kmer_seq = seqs[seqId].indexOf(km.getKmerString());
				if (kmer_seq<0){
					seqs[seqId] = SequenceUtils.reverseComplement(seqs[seqId]);
					isPlusStrands[seqId] = false;
					kmer_seq = seqs[seqId].indexOf(km.getKmerString());
				}
				posSeqs[seqId] = -kmer_seq + km.getShift();
				seqAlignRefs[seqId] = km.getKmerString();
			}
		}
	}

    private Pair<Integer, Integer> constructPWM(double[][] pwm){
    	// normalize, compare to background, and log2
    	double[] ic = new double[pwm.length];						// information content
    	for (int p=0;p<pwm.length;p++){
    		int sum=0;
    		for (int b=0;b<LETTERS.length;b++){
    			char base = LETTERS[b];
    			if (pwm[p][base]==0)
    				pwm[p][base]=0.001; 
    			sum += pwm[p][base];
    		}
    		for (int b=0;b<LETTERS.length;b++){
    			char base = LETTERS[b];
    			double f = pwm[p][base]/sum;						// normalize freq
    			pwm[p][base] = Math.log(f/config.bg[b])/Math.log(2.0);		//log base 2
    			ic[p] += f*pwm[p][base];
    		}
    	}
    	// make a WeightMatrix object
    	int leftIdx, rightIdx;
    	if (config.trim_simple){	// trim low ic ends (simple method)
	    	leftIdx=ic.length-1;
	    	for (int p=0;p<ic.length;p++){
	    		if (ic[p]>=config.ic_trim){
	    			leftIdx = p;
	    			break;
	    		}
	    	}
	    	rightIdx=0;
	    	for (int p=ic.length-1;p>=0;p--){
	    		if (ic[p]>=config.ic_trim){    			
	    			rightIdx=p;
	    			break;
	    		}
	    	}
    	}
	    else{	// trim low ic ends (more sophisticated method)
	    	// To avoid situations where a remote position happens to pass the ic threshold
	    	leftIdx=-1;
	    	double score = 0;
	    	for (int p=0;p<ic.length;p++){
	    		if (ic[p]>config.ic_trim){
	    			score ++;
	    		}
	    		else{
	    			score -= 0.3;
	    		}
	    		if (score<0 && p-leftIdx<config.k/2){
	    			score=0;
	    			leftIdx=p;
	    		}
	    	}
	    	leftIdx++;
	    	
	    	rightIdx=ic.length;
	    	score = 0;
	    	for (int p=ic.length-1;p>=0;p--){
	    		if (ic[p]>config.ic_trim){
	    			score ++;
	    		}
	    		else{
	    			score -= 0.3;
	    		}
	    		if (score<0 && rightIdx-p<config.k/2){
	    			score=0;
	    			rightIdx=p;
	    		}
	    	}
	    	rightIdx--;
	    }
    	
    	if (rightIdx-leftIdx+1<=config.k/2){
    		System.out.println("makePWM: PWM is too short, stop here.");
    		if (config.bmverbose>1){
		    	StringBuilder sb = new StringBuilder("Information contents of aligned positions\n");
		    	for (int p=0;p<ic.length;p++){
		    		sb.append(String.format("%d\t%.1f\t%s\n", p, ic[p], (p==leftIdx||p==rightIdx)?"<--":""));
		    	}
		    	System.out.println(sb.toString());
    		}
    	}
    	return new Pair<Integer, Integer>(leftIdx, rightIdx);
    }
    
	// update kmerEngine with the predicted kmer-events
	public void updateKmerEngine(boolean makePFM){
		long tic = System.currentTimeMillis();
		if (config.k==-1)
			return;
		ArrayList<ComponentFeature> compFeatures = new ArrayList<ComponentFeature>();
		for(Feature f : signalFeatures)
			compFeatures.add((ComponentFeature)f);
		
		// reload the test sequences, and purify kmers
		ArrayList<Point> events = getEventPoints();
		kEngine.loadTestSequences(events, config.k_win);
		purifyKmers(compFeatures);
		
		if (makePFM){
			updateKmersWithPWM(compFeatures);
			purifyKmers(compFeatures);
		}
		
		consolidateKmers(compFeatures);
		
		ArrayList<Kmer> kmers = countKmers(compFeatures);	
		
		log(1, "Kmers ("+kmers.size()+") updated, "+CommonUtils.timeElapsed(tic));
		kEngine.updateEngine(kmers, outName, false);
	}

	private void updateKmersWithPWM(ArrayList<ComponentFeature> compFeatures){
	    	
	    	// cluster binding sequences
	    	ArrayList<ComponentFeature> unalignedFeatures = new ArrayList<ComponentFeature>();
	    	ArrayList<ComponentFeature> nullKmerFeatures = new ArrayList<ComponentFeature>();
	    	for(ComponentFeature cf : compFeatures){
	    		if (cf.getKmer()!=null )
	    			unalignedFeatures.add(cf);
	    		else
	    			nullKmerFeatures.add(cf);
	    	}
	    	ArrayList<MotifCluster> clusters = clusterEventSequences(unalignedFeatures, nullKmerFeatures);
	    	
	    	// build WeightMatrix from each aligned group of sequences
	    	StringBuilder sb_pfm = new StringBuilder();
	    	//print out the kmers
	    	StringBuilder sb_kmer = new StringBuilder();
	    	sb_kmer.append("Kmer").append("\t").append("KmerRC").append("\t").append("seqCt").append("\t")
	    	  .append("gShift").append("\t").append("Alignment").append("\t").append("Reference").append("\n");
	    	//print out the kmer strings in fasta format
//	    	StringBuilder sb_fa = new StringBuilder();
	    	int goodClusterCount = 0;
	    	for (MotifCluster cluster : clusters){
	    		ArrayList<ComponentFeature> alignedFeatures = cluster.alignedFeatures;
	    		ArrayList<Integer> motifPos = cluster.motifStartInSeq;
	    		if ( alignedFeatures.size()<=config.kmer_cluster_seq_count)
	    			continue;
	    		
	    		goodClusterCount++;
	    		cluster.clusterId = goodClusterCount;
	    		System.out.println(String.format("-------------------------------\n%s motif cluster #%d, from %d binding events.", outName, goodClusterCount, alignedFeatures.size()));
	    		WeightMatrix wm = cluster.matrix;
	    		if (wm==null)
	    			continue;
	    		double maxScore = wm.getMaxScore();
	    		int pos = cluster.bindingPosition;
	    		if (pos>0)
	    			System.out.println(CommonUtils.padding(pos, ' ')+"|\n"+ WeightMatrix.printMatrixLetters(wm));
	    		else
	    			System.out.println(WeightMatrix.printMatrixLetters(wm));
	    		System.out.println(String.format("PWM threshold: %.2f", cluster.pwmThreshold));
	    		// use PWM to scan null-kmer events, to discover kmers that was not included in the inital set
	    		// TODO: maybe should favor PWM hit that is close to the event location
	    		WeightMatrixScorer scorer = new WeightMatrixScorer(wm);
	        	ArrayList<ComponentFeature> temp = new ArrayList<ComponentFeature>();
	        	ArrayList<Kmer> newKmers = new ArrayList<Kmer>();
	        	
	        	HashMap<String, Kmer> alignedKmerMap = new HashMap<String, Kmer> ();
	        	ArrayList<Kmer> dups = new ArrayList<Kmer>();
	        	for (Kmer km : cluster.alignedKmers){
	        		String kmStr = km.getKmerString();
	        		if (alignedKmerMap.containsKey(kmStr)){
	        			System.err.println("nullKmerFeatures: "+kmStr+" duplicated!");
	        			alignedKmerMap.get(kmStr).mergeKmer(km);
	        			dups.add(km);
	        		}
	        		else
	        			alignedKmerMap.put(kmStr, km);
	        	}
	        	cluster.alignedKmers.removeAll(dups);
	        	
	    		for (ComponentFeature nf : nullKmerFeatures){
					if (nf.getBoundSequence().length()<config.k)
						continue;
					
	    			Pair<Integer, Double> hit = scanPWM(nf.getBoundSequence(), wm, scorer);	// PWM hit in the bound sequence
	    			int hitPos = hit.car();									// motif hit start pos
	    			double score = hit.cdr();
//	    			if (score < maxScore*config.wm_factor)
	    			if (score < cluster.pwmThreshold)
	    				continue;
	    			  
	    			if (hitPos<0){											// if match on '-' strand
	    				nf.flipBoundSequence();
	    				hitPos = scanPWM(nf.getBoundSequence(), wm, scorer).car();		
	    			}
	    			String seq = nf.getBoundSequence();
					if (seq.length()==config.k){
						// check if the k-mer is from negative set
//						if (kEngine.isNegativeKmer(seq))
//							continue;
						Kmer kmer = new Kmer(seq, 1);
						//pos_kmer - pos_motif_binding
						kmer.setKmerStartOffset(0-(hitPos+cluster.bindingPosition));
						kmer.incrStrength(nf.getTotalEventStrength());
						kmer.setAlignString("nullPWM:"+WeightMatrix.getMaxLetters(wm));
	    				newKmers.add(kmer);
	    				nf.setKmer(kmer);
	    				alignedFeatures.add(nf);
	    				assert(hitPos>=0);
	    				motifPos.add(hitPos);
	    				temp.add(nf);
					}		    			
					else{
						// offset between PWM hit and Kmer start, this should work no matter which one is wider		    			
						int offset = cluster.bindingPosition-config.k/2;	
		    			int left = hitPos+offset;
		    			int right = left+config.k;
		    			if (right>seq.length() || left<0){
		    				// this sometimes happens when the PWM is shorter than kmer, just ignore it
	//	    				System.out.println("Warning: Get Kmer for nullKmerFeatures: seqLen="+seq.length()+" <"+left+", "+right+">");
		    				continue;
		    			}
	    				String kmerStr = seq.substring(left, right);
	    				// check if the k-mer is from negative set
//	    				if (kEngine.isNegativeKmer(kmerStr))
//	    					continue;
	    				Kmer kmer = new Kmer(kmerStr, 1);
						//pos_kmer - pos_motif_binding
//	    				kmer.setKmerStartOffset(left-(hitPos+cluster.bindingPosition));
	    				kmer.setKmerStartOffset(-config.k/2);
	    				kmer.incrStrength(nf.getTotalEventStrength());
						kmer.setAlignString("nullPWM:"+WeightMatrix.getMaxLetters(wm));
	    				newKmers.add(kmer);
	    				nf.setKmer(kmer);
	    				alignedFeatures.add(nf);
	    				motifPos.add(hitPos);
	    				temp.add(nf);
	    			}
	    		}
	    		if (!temp.isEmpty()){
	    			ArrayList<Kmer> realNewKmers = new ArrayList<Kmer>();
		    		nullKmerFeatures.removeAll(temp);
		    		for (Kmer kmer: newKmers){
		    			String kmStr = kmer.getKmerString();
		    			if (alignedKmerMap.containsKey(kmStr)){
		    				alignedKmerMap.get(kmStr).mergeKmer(kmer);
		    			}
		    			else{
		    				cluster.alignedKmers.add(kmer);
		    				alignedKmerMap.put(kmStr, kmer);
		    				realNewKmers.add(kmer);
		    			}
		    		}
		    		
		    		clusterByNewKmers(cluster, nullKmerFeatures, realNewKmers);
		    		System.out.println("Rescued kmers, now cluster contains "+cluster.alignedFeatures.size()+" events.");  		    		
	    		}
	    	}
	    	System.out.println("-------------------------------\nTotal clusters of motifs: " + goodClusterCount);
	    	
	    	goodClusterCount=0;
	    	for (MotifCluster cluster : clusters){
	    		ArrayList<ComponentFeature> alignedFeatures = cluster.alignedFeatures;
	    		ArrayList<Integer> motifPos = cluster.motifStartInSeq;
	    		if (alignedFeatures.size()<=config.kmer_cluster_seq_count)
	    			continue;
	    		else
	    			sb_kmer.append(outName+"_motif cluster "+cluster.clusterId+", from "+alignedFeatures.size()+" binding events.\n");
	    		
	    		// update the kmer offset information (kmer start relative from the binding position of PWM)
	    		// First get offset for the seed kmer
	    		ArrayList<Integer> seedOffsets = new ArrayList<Integer>();
	    		for (int f=0;f<alignedFeatures.size();f++){
	            	ComponentFeature cf = alignedFeatures.get(f);
	            	if (cf.getKmer()==cluster.seedKmer){
	            		String seq = cf.getBoundSequence();
		            	int start = seq.indexOf(cluster.seedKmer.getKmerString());
		            	if (start==-1)
		            		continue;
		            	seedOffsets.add(start - motifPos.get(f) - cluster.bindingPosition);
	            	}
	    		}
	    		int seedOffset =0;
	    		for (int offset: seedOffsets){
	    			seedOffset+=offset;	    			
	    		}
	    		seedOffset = StatUtil.round(seedOffset/seedOffsets.size());
	    		// then update the Kmer aligned by SeedKmer or Overlap using seed offset
	    		for (Kmer kmer:cluster.alignedKmers){
	    			if (kmer.getAlignString()==null)
	    				continue;
	    			if (!kmer.getAlignString().contains("PWM")){
	    				kmer.setKmerStartOffset(kmer.getShift()+seedOffset);
	    			}
	    		}
	    		// for kmers aligned by PWM scanning	    		
	            TreeMap<Kmer, ArrayList<Integer>> kmerOffsets = new TreeMap<Kmer, ArrayList<Integer>> ();
	    		for (int f=0;f<alignedFeatures.size();f++){	    			
	            	ComponentFeature cf = alignedFeatures.get(f);
	            	Kmer kmer = cf.getKmer();
	    			if (kmer.getAlignString()==null || !kmer.getAlignString().contains("PWM"))		// skip other kmers
	    				continue;
	    			
            		String seq = cf.getBoundSequence();
	            	int start = seq.indexOf(kmer.getKmerString());
	            	if (start==-1)
	            		continue;
	            	if (!kmerOffsets.containsKey(kmer)){
	            		kmerOffsets.put(kmer, new ArrayList<Integer>());
	            	}
	            	kmerOffsets.get(kmer).add(start - motifPos.get(f) - cluster.bindingPosition);
	    		}	
	    		for (Kmer kmer:kmerOffsets.keySet()){
	    			double sum = 0;
	    			for (int i:kmerOffsets.get(kmer)){
	    				sum += i;
	    			}
	    			kmer.setKmerStartOffset(StatUtil.round(sum/kmerOffsets.get(kmer).size()));
	    			kmer.setGroup(cluster.clusterId);
	    		}
	            
	            // make PFM
	    		if (cluster.matrix==null)
	    			continue;
	    		sb_pfm.append(getPFMString(alignedFeatures, motifPos, cluster.matrix.length(), cluster.clusterId));
	    		
	    		// print aligned kmers with their shift information
	    		// print fasta file with kmer string (maybe useful to do multiple alignment)
	            ArrayList<Kmer> akmers = cluster.alignedKmers;
	            Collections.sort(akmers);

				int kk=0;
	    		for (Kmer km: akmers){
					String shiftedKmer = CommonUtils.padding(Math.max(0, km.getKmerStartOffset()+config.k), '-').concat(km.getKmerString());
					sb_kmer.append(km.getKmerString()).append("\t").append(km.getKmerRC()).append("\t").append(km.getSeqHitCount()).append("\t")
					  .append(km.getKmerStartOffset()).append("\t").append(shiftedKmer).append("\t").append(km.getAlignString()).append("\n");		
					// fasta
//					sb_fa.append(">"+outName+"_"+cluster.clusterId+"_"+kk).append("\n");
//					sb_fa.append(km.getKmerString()).append("\n");		
		    		kk++;				
		    	}
				sb_kmer.append("\n");
	    	}
	    	CommonUtils.writeFile(outName+"_PFM.txt", sb_pfm.toString());
	    	sb_kmer.append("Total kmer groups: " + clusters.size());
	    	CommonUtils.writeFile(outName+"_Kmers_Aligned.txt", sb_kmer.toString());
//	    	CommonUtils.writeFile(outName+"_Kmers.fa", sb_fa.toString());
		}

	/**
	 * Print the overlapping kmer counts with increasing number of k
	 * in a specified window around the binding events
	 */
	public void printOverlappingKmers(){
		for (int k=config.k;k<config.k+5;k++){
			String name = outName+"_OK_win"+ (config.k*2);
	    	ArrayList<Point> events = getEventPoints();
			kEngine.buildEngine(k, events, config.k*2, config.hgp, config.k_fold, name);
		}
	}
	
	/** Test kmers against positive/negative sequence set<br>
     * 	if kmer is not enriched, componentFeature is modified to have a null kmer
     */
    private void purifyKmers(ArrayList<ComponentFeature> compFeatures){
    	HashMap<Kmer, ArrayList<ComponentFeature>> kmer2cfs = new HashMap<Kmer, ArrayList<ComponentFeature>>();
    	for(ComponentFeature cf : compFeatures){
    		Kmer kmer = cf.getKmer();
    		if (kmer==null)
    			continue;
    		if (!kmer2cfs.containsKey(kmer))
    			kmer2cfs.put(kmer, new ArrayList<ComponentFeature>());
    		kmer2cfs.get(kmer).add(cf);
    	}
    	HashMap<String, ArrayList<Kmer>> str2kmers = new HashMap<String, ArrayList<Kmer>>();
    	for (Kmer kmer:kmer2cfs.keySet()){
    		String kmerStr = kmer.getKmerString();
    		if (!str2kmers.containsKey(kmerStr))
    			str2kmers.put(kmerStr, new ArrayList<Kmer>());
    		str2kmers.get(kmerStr).add(kmer);
    	}
    	ArrayList<String> kmerStrings = new ArrayList<String>();
    	kmerStrings.addAll(str2kmers.keySet());
    	
    	double[] hgps = kEngine.computeHGPs(kmerStrings);
    	for (int i=0;i<hgps.length;i++){
    		for (Kmer kmer:str2kmers.get(kmerStrings.get(i))){
//    			kmer.setHgp(hgps[i]);			//TODO: verify if it is good to update hgp here
    			if (hgps[i]>config.hgp){
    			// update the Kmer reference in the binding event object
        			for (ComponentFeature cf : kmer2cfs.get(kmer))
        				cf.setKmer(null);
        		}
    		}
    	}
    }
	
	/** Consolidate kmers that have same String, sum the count and strength, average the shift
     * There are maybe kmers created as different object at different time, but they have same String
     * Test kmers against negative kmer set
     * This method will modify compFeatures to update the kmer information
     */
    private void consolidateKmers(ArrayList<ComponentFeature> compFeatures){
    	HashMap<Kmer, ArrayList<ComponentFeature>> kmer2cfs = new HashMap<Kmer, ArrayList<ComponentFeature>>();
    	for(ComponentFeature cf : compFeatures){
    		Kmer kmer = cf.getKmer();
    		if (kmer==null)
    			continue;
    		if (!kmer2cfs.containsKey(kmer))
    			kmer2cfs.put(kmer, new ArrayList<ComponentFeature>());
    		kmer2cfs.get(kmer).add(cf);
    	}
    	HashMap<String, ArrayList<Kmer>> str2kmers = new HashMap<String, ArrayList<Kmer>>();
    	for (Kmer kmer:kmer2cfs.keySet()){
    		kmer.setAlignString(null);		// reset, so it does not carry over to next round
    		kmer.setShift(-1000);
    		kmer.setGroup(-1);
    		String kmerStr = kmer.getKmerString();
    		if (!str2kmers.containsKey(kmerStr))
    			str2kmers.put(kmerStr, new ArrayList<Kmer>());
    		str2kmers.get(kmerStr).add(kmer);
    	}
    	for (String ks: str2kmers.keySet()){
    		ArrayList<Kmer> kmers = str2kmers.get(ks);
    		if (kmers.size()==1)
    			continue;
    		int count=0;
    		double strength = 0;
    		int offset = 0;
    		ArrayList<Integer> allOffsets = new ArrayList<Integer>();
    		for (Kmer kmer:kmers){
    			count += kmer.getSeqHitCount();
    			strength += kmer.getStrength();
    			allOffsets.add(kmer.getKmerStartOffset());
    			offset += kmer.getKmerStartOffset(); 
    		}

    		offset = StatUtil.round(offset/allOffsets.size());			// mean shift
//    		Collections.sort(allOffsets);
//    		shift = allOffsets.get(allOffsets.size()/2);			// median shift
    		Kmer cKmer = new Kmer(ks, count);						// consolidated kmer
    		cKmer.incrStrength(strength);
    		cKmer.setKmerStartOffset(offset);
    		// update the Kmer reference in the binding event object
    		for (Kmer kmer:kmers){
    			for (ComponentFeature cf : kmer2cfs.get(kmer))
    				cf.setKmer(cKmer);
    		}
    	}
    }
    
    /**
	 * Count kmers in all binding events
	 * Assuming the Kmer object reference in the componentFeature is unique for each kmer string
	 * Need to run consolidateKmers to ensure this.
	 */
	private ArrayList<Kmer> countKmers(ArrayList<ComponentFeature> compFeatures){
		HashMap<Kmer, Integer> kmer2count = new HashMap<Kmer, Integer>();
		HashMap<Kmer, Double> kmer2strength = new HashMap<Kmer, Double>();
		for(ComponentFeature cf : compFeatures){
			Kmer kmer = cf.getKmer();
			if (kmer==null)
				continue;
			if (kmer2strength.containsKey(kmer))	
				kmer2strength.put(kmer, kmer2strength.get(kmer)+cf.getTotalEventStrength());
			else
				kmer2strength.put(kmer, cf.getTotalEventStrength());
	
			if (kmer2count.containsKey(kmer))				
				kmer2count.put(kmer, kmer2count.get(kmer)+1);
			else
				kmer2count.put(kmer, 1);
		}
		ArrayList<Kmer> kmers = new ArrayList<Kmer>();
		for (Kmer kmer:kmer2count.keySet()){
			kmer.setSeqHitCount(kmer2count.get(kmer));
			kmer.setNegCount(-1);
			kmer.setStrength(kmer2strength.get(kmer));
			kmer.setGroup(-1);								// clear group ID 
			kmers.add(kmer);
		}	
		return kmers;
	}

	/**
     * Cluster the bound sequences at binding events based on Kmer and PWM learning
     */
    private ArrayList<MotifCluster> clusterEventSequences(ArrayList<ComponentFeature> unalignedFeatures, 
    		ArrayList<ComponentFeature> nullKmerFeatures){

    	ArrayList<MotifCluster> clusters = new ArrayList<MotifCluster>();
    	if (unalignedFeatures.isEmpty())
    		return clusters;
    	
    	fixInvalidKmerString(unalignedFeatures, nullKmerFeatures);
    	
    	// greedy grouping of bound sequences until all are grouped
    	while(!unalignedFeatures.isEmpty()){
    		
    		MotifCluster cluster = growSeqCluster(unalignedFeatures);
    		
    		if (cluster!=null){
    			if (cluster.alignedFeatures.size()<config.kmer_cluster_seq_count)	{	// if the cluster is too small, set as nullKmer to process later
    				for (ComponentFeature cf : cluster.alignedFeatures){
    					cf.setKmer(null);
    					nullKmerFeatures.add(cf);
    				}
    			}
    			else	
    				clusters.add(cluster);
    			
    			fixInvalidKmerString(unalignedFeatures, nullKmerFeatures);
    		}
    		else
    			break;
    	}
    	return clusters;
      }
    
    /**
     * fix Kmers that are not contained in boundSequence of the events
     */     
    private void fixInvalidKmerString(ArrayList<ComponentFeature> unalignedFeatures, 
    		ArrayList<ComponentFeature> nullKmerFeatures){
    	// error checking
    	ArrayList<ComponentFeature> toNullList = new ArrayList<ComponentFeature>();
    	for (ComponentFeature cf:unalignedFeatures){   
			if (!cf.getBoundSequence().contains(cf.getKmer().getKmerString())){
				// if kmer is not in sequence, flip and try again
				cf.flipBoundSequence();
				if (!cf.getBoundSequence().contains(cf.getKmer().getKmerString())){
					// if still not working, set kmer to null, put it in nullCluster.
					cf.setKmer(null);
					toNullList.add(cf);
				}
			}
    	}
    	unalignedFeatures.removeAll(toNullList);
    	nullKmerFeatures.addAll(toNullList);
    }
    
      /** 
       * greedily grow a cluster from the top count kmer
       * @param unalignedFeatures
       * @return
       */
      private MotifCluster growSeqCluster(ArrayList<ComponentFeature> unalignedFeatures){
    	
    	if (unalignedFeatures.size()==0)
    		return null;
    	
    	MotifCluster motifCluster = new MotifCluster();
    	ArrayList<ComponentFeature> alignedFeatures = new ArrayList<ComponentFeature>();
    	// PWM hit start position in the bound sequence
    	ArrayList<Integer> motifStartInSeq = new ArrayList<Integer>();	
    	motifCluster.alignedFeatures = alignedFeatures;
    	motifCluster.motifStartInSeq = motifStartInSeq;
    	
    	// aligne kmers that are similar to seedKmer
    	clusterBySeedKmer(motifCluster, unalignedFeatures);
    	if (alignedFeatures.size()<config.kmer_cluster_seq_count)		// do not have enough data to build PWM for further analysis
    		return motifCluster;
    	
    	boolean noMore=false;
    	WeightMatrix wm = makePWM( motifCluster, true);
    	if (wm==null){
    		return motifCluster;
    	}
//    	System.out.println("After clusterByKmer()\n" +
//    			CommonUtils.padding(motifCluster.bindingPosition, ' ')+"|\n"+
//    			WeightMatrix.printMatrixLetters(wm));
    	while(!noMore){
    		// Save the current state
    		MotifCluster old = motifCluster.clone();
    		
    		ArrayList<ComponentFeature> unaligned_old = new ArrayList<ComponentFeature>();
    		unaligned_old.addAll(unalignedFeatures);
    		
    		noMore = clusterByPWM(motifCluster, unalignedFeatures, config.bg);
    		wm = makePWM( motifCluster, true);
    		if (wm==null){		
    			// if the pwm is not good, return the previous result		
    			System.out.println("growSeqCluster: pwm is not good, take previous one.");
    			String bp = (motifCluster.bindingPosition>=0)? 
    					(CommonUtils.padding(motifCluster.bindingPosition, ' ')+"|\n"):"";
            	System.out.println( bp + WeightMatrix.printMatrixLetters(old.matrix));
            	unalignedFeatures.clear();
    			unalignedFeatures.addAll(unaligned_old);
    	    	for (ComponentFeature cf:unalignedFeatures){
    	    		// reverse some kmer-seq mismatch effect from those unsuccessful clustering manipulations
    				if (cf.getBoundSequence().indexOf(cf.getKmer().getKmerString())==-1)		
    					cf.flipBoundSequence();
    	    	}
    			return old;
    		}
    		else{
//            	System.out.println("After clusterByPWM()\n" +
//            			CommonUtils.padding(motifCluster.bindingPosition, ' ')+"|\n"+
//            			WeightMatrix.printMatrixLetters(wm));    			
    		}
        	if (wm.length()<=config.k*3/4)		// PWM is too short for further analysis
        		return motifCluster;
    	}
    	
    	return motifCluster;
      }
      //greedily grow a cluster from the top count kmer only by matching kmers
	    private void clusterBySeedKmer(MotifCluster motifCluster, 
	    	ArrayList<ComponentFeature> unalignedFeatures){
	    	ArrayList<ComponentFeature> alignedFeatures = motifCluster.alignedFeatures;
	    	// initially motifStartInSeq is used to store the position that the seedKmer may 
	    	// align to on this sequence, so that we can align the sequences to learn initial PWM
	    	// only after learning the PWM, we update these to be the PWM motif start position
	    	ArrayList<Integer> motifStartInSeq = motifCluster.motifStartInSeq;
	    	
	    	ArrayList<Kmer> kmers = countKmers(unalignedFeatures);
	    	if (kmers.size()<1)
	    		return;
//	    	Collections.sort(kmers);
			Collections.sort(kmers, new Comparator<Kmer>(){
			    public int compare(Kmer o1, Kmer o2) {
			    		return o1.compareByHGP(o2);
			    }
			});
	    	for (Kmer kmer:kmers)
	    		kmer.setAlignString(null);
	    	
	    	motifCluster.seedKmer = kmers.get(0);
	    	ArrayList<Kmer> alignedKmers = new ArrayList<Kmer>();
	    	String seedKmerStr = motifCluster.seedKmer.getKmerString();
	    	String seedKmerRC = motifCluster.seedKmer.getKmerRC();
	    	
	    	// map kmers to the events
	    	HashMap<Kmer, ArrayList<ComponentFeature>> kmer2cf = new HashMap<Kmer, ArrayList<ComponentFeature>>();
	    	int invalidKmerCount = 0;
	    	for(ComponentFeature cf : unalignedFeatures){
	    		Kmer kmer = cf.getKmer();
	    		assert(kmer!=null);
	    		if (!(cf.getBoundSequence().contains(kmer.getKmerString())||
	    				cf.getBoundSequence().contains(kmer.getKmerRC())) ){
	    			invalidKmerCount++;
	    			continue;
	    		}
	    		if (!kmer2cf.containsKey(kmer))
	    			kmer2cf.put(kmer, new ArrayList<ComponentFeature>());
	    		kmer2cf.get(kmer).add(cf);
	    	}
	    	if (invalidKmerCount>0){
	    		System.out.println("clusterBySeedKmer: "+invalidKmerCount+" out of "+unalignedFeatures.size()+" kmers can not match boundSequence.");
	    	}
	    	
	    	for (Kmer kmer:kmer2cf.keySet()){
	    		// perfect match or 1 mismatch (2 mismatch if k>10), no shift
	    		if (kmer.hasString(seedKmerStr) || mismatch(seedKmerStr, kmer.getKmerString())<=1+config.k*0.1){
	    			for (ComponentFeature cf:kmer2cf.get(kmer)){
	    				int idx = cf.getBoundSequence().indexOf(kmer.getKmerString());
	    				if (idx==-1)
	    					System.err.println("growByKmer: kmer "+kmer.getKmerString()+" is not in seq "+cf.getBoundSequence());
	    				else{
	    					alignedFeatures.add(cf);
			    			motifStartInSeq.add(idx);
		    			}
	    			}
		    		alignedKmers.add(kmer);
		    		kmer.setShift(0);
		    		kmer.setKmerStartOffset(motifCluster.seedKmer.getKmerStartOffset()+kmer.getShift());
		    		continue;
	    		}
	    		else if (kmer.hasString(seedKmerRC)||mismatch(seedKmerRC, kmer.getKmerString())<=1+config.k*0.1){
	    			for (ComponentFeature cf:kmer2cf.get(kmer)){
	    				cf.flipBoundSequence();
	    				int idx = cf.getBoundSequence().indexOf(kmer.getKmerRC());
	    				if (idx==-1)
	    					System.err.println("growByKmer: kmerRC "+kmer.getKmerRC()+" is not in seq "+cf.getBoundSequence());
	    				else{
	    					alignedFeatures.add(cf);
			    			motifStartInSeq.add(idx);
		    			}
	    			}
	    			alignedKmers.add(kmer);
	    			kmer.RC();
		    		kmer.setShift(0);
		    		kmer.setKmerStartOffset(motifCluster.seedKmer.getKmerStartOffset()+kmer.getShift());
	    			continue;
	    		}
	    		else{
	    		// match or 1 mismatch after 1 shift
		    		String kmerStr = kmer.getKmerString().substring(1);
		    		String kmerRC = kmer.getKmerRC().substring(1);
		    		String ref = seedKmerStr.substring(0, seedKmerStr.length()-1);
		    		if (mismatch(ref, kmerStr)<=1){
		    			for (ComponentFeature cf:kmer2cf.get(kmer)){
		    				int idx = cf.getBoundSequence().indexOf(kmer.getKmerString());
		    				if (idx==-1)
		    					System.err.println("growByKmer: kmer "+kmer.getKmerString()+" is not in seq "+cf.getBoundSequence());
		    				else{
		    					alignedFeatures.add(cf);
				    			motifStartInSeq.add(idx+1);
			    			}
		    			}
		    			alignedKmers.add(kmer);
			    		kmer.setShift(-1);
			    		kmer.setKmerStartOffset(motifCluster.seedKmer.getKmerStartOffset()+kmer.getShift());
		    			continue;
		    		}
		    		else if(mismatch(ref, kmerRC)<=1){	// if match RC, flip kmer and seq
		    			for (ComponentFeature cf:kmer2cf.get(kmer)){
		    				cf.flipBoundSequence();
		    				int idx = cf.getBoundSequence().indexOf(kmer.getKmerRC());
		    				if (idx==-1)
		    					System.err.println("growByKmer: kmerRC "+kmer.getKmerRC()+" is not in seq "+cf.getBoundSequence());
		    				else{
		    					alignedFeatures.add(cf);
				    			motifStartInSeq.add(idx+1);
			    			}
		    			}
		    			alignedKmers.add(kmer);
		    			kmer.RC();
		    			kmer.setShift(-1);
			    		kmer.setKmerStartOffset(motifCluster.seedKmer.getKmerStartOffset()+kmer.getShift());
		    			continue;
		    		}
		    		// shift to the other direction
		    		kmerStr = kmer.getKmerString().substring(0, seedKmerStr.length()-1);
		    		kmerRC = kmer.getKmerRC().substring(0, seedKmerStr.length()-1);
		    		ref = seedKmerStr.substring(1);
		    		if (mismatch(ref, kmerStr)<=1){
		    			for (ComponentFeature cf:kmer2cf.get(kmer)){
		    				int idx = cf.getBoundSequence().indexOf(kmer.getKmerString());
		    				if (idx==-1)
		    					System.err.println("growByKmer: kmer "+kmer.getKmerString()+" is not in seq "+cf.getBoundSequence());
		    				else{
		    					alignedFeatures.add(cf);
				    			motifStartInSeq.add(idx-1);
			    			}
		    			}
		    			alignedKmers.add(kmer);
		    			kmer.setShift(1);
			    		kmer.setKmerStartOffset(motifCluster.seedKmer.getKmerStartOffset()+kmer.getShift());
		    			continue;
		    		}
		    		else if(mismatch(ref, kmerRC)<=1){	// if match RC, flip kmer and seq
		    			for (ComponentFeature cf:kmer2cf.get(kmer)){
		    				cf.flipBoundSequence();
		    				int idx = cf.getBoundSequence().indexOf(kmer.getKmerRC());
		    				if (idx==-1)
		    					System.err.println("growByKmer: kmerRC "+kmer.getKmerRC()+" is not in seq "+cf.getBoundSequence());
		    				else{
		    					alignedFeatures.add(cf);
				    			motifStartInSeq.add(idx-1);
			    			}
		    			}
		    			alignedKmers.add(kmer);
		    			kmer.RC();
		    			kmer.setShift(1);
			    		kmer.setKmerStartOffset(motifCluster.seedKmer.getKmerStartOffset()+kmer.getShift());
		    			continue;
		    		}
	    		}
	    	} // for each kmer
	    	
	    	// update sequences
	    	unalignedFeatures.removeAll(alignedFeatures);
	    	// update motifCluster
	    	if (motifCluster.alignedKmers==null)
	    		motifCluster.alignedKmers = new ArrayList<Kmer>();    	
	    	motifCluster.alignedKmers.addAll(alignedKmers);
	    	
	    	for (Kmer kmer:alignedKmers){
	    		kmer2cf.remove(kmer);
	    		kmer.setAlignString("Seed:"+seedKmerStr);
	    	}
	    	
	    	// extend to overlapped kmers from all aligned kmers, iteratively
	    	long tic = System.currentTimeMillis();
	    	clusterByOverlapKmer(motifCluster, unalignedFeatures, kmer2cf);
//	    	System.out.println("ClusterByOverlapKmer, "+CommonUtils.timeElapsed(tic));
	    }

	//greedily extend aligned kmers to overlapped kmers
	    private void clusterByOverlapKmer(MotifCluster motifCluster, ArrayList<ComponentFeature> unalignedFeatures,
	    	HashMap<Kmer, ArrayList<ComponentFeature>> kmer2cf){
	    	ArrayList<Kmer> alignedKmers = new ArrayList<Kmer>();
	    	alignedKmers.addAll(motifCluster.alignedKmers);
	    	ArrayList<Kmer> newAlignedKmers = new ArrayList<Kmer>();
	    	while(!alignedKmers.isEmpty()){
		    	for (Kmer km : alignedKmers){
		    		String kmStr = km.getKmerString();
		    		for (int i=0;i<kmStr.length()-config.k_overlap;i++){
		    			String overlap = kmStr.substring(i, i+config.k_overlap);
		    			ArrayList<Kmer> overlapped = new ArrayList<Kmer>();
		    			for (Kmer kmer:kmer2cf.keySet()){
		    				if (km.getSeqHitCount() < kmer.getSeqHitCount())				// only extend overlap from more confident kmers
		    					continue;
		    				int idx = kmer.getKmerString().indexOf(overlap);
		    				if (idx!=-1){	// the kmer overlaps
		    					int newShift = km.getShift()+i-idx;
		    					if (newShift<=config.k/2){						// the overlap should not be from too many overlap extension
			    					kmer.setShift(newShift);
			    					overlapped.add(kmer);
			    					String rcString = kmer.getKmerRC();
			    					String kmerString = kmer.getKmerString();
			    					for (ComponentFeature cf:kmer2cf.get(kmer)){
			    						if (cf.getBoundSequence().contains(rcString))
			    							cf.flipBoundSequence();
			    						else if (!cf.getBoundSequence().contains(kmerString)){
					    					System.err.println("clusterByOverlapKmer: kmer "+kmerString+" is not in seq "+cf.getBoundSequence());
					    					continue;
			    						}
				    					motifCluster.alignedFeatures.add(cf);
				    					unalignedFeatures.remove(cf);
				    					int idx_seq = cf.getBoundSequence().indexOf(kmerString);
				    					motifCluster.motifStartInSeq.add(idx_seq-kmer.getShift());				    					
					    			}
			    					newAlignedKmers.add(kmer);
			    		    		kmer.setAlignString("Overlap:"+kmStr);
		    					}
		    				}
		    				else{
			    				int idx2 = kmer.getKmerRC().indexOf(overlap);
			    				if (idx2!=-1){	// the kmer overlaps to RC
			    					String rcString = kmer.getKmerString();
			    					kmer.RC();
			    					String kmerString = kmer.getKmerString();
			    					int newShift = km.getShift()+i-idx2;
			    					if (newShift<=config.k/2){						// the overlap should not be from too many overlap extension
				    					kmer.setShift(newShift);
				    					overlapped.add(kmer);
				    					for (ComponentFeature cf:kmer2cf.get(kmer)){
				    						if (cf.getBoundSequence().contains(rcString))
				    							cf.flipBoundSequence();
				    						else if (!cf.getBoundSequence().contains(kmerString)){
						    					System.err.println("clusterByOverlapKmer: kmerRC "+kmerString+" is not in seq "+cf.getBoundSequence());
						    					continue;
				    						}
					    					motifCluster.alignedFeatures.add(cf);
					    					unalignedFeatures.remove(cf);
					    					int idx_seq = cf.getBoundSequence().indexOf(kmerString);
					    					motifCluster.motifStartInSeq.add(idx_seq-kmer.getShift());
						    			}
				    					newAlignedKmers.add(kmer);
				    		    		kmer.setAlignString("OverlapRC:"+kmStr);
			    					}
			    				}
		    				}
		    			}
				    	for (Kmer kmer:overlapped)
				    		kmer2cf.remove(kmer);
		    		}
		    	}
		    	alignedKmers.clear();
		    	motifCluster.alignedKmers.addAll(newAlignedKmers);
		    	alignedKmers.addAll(newAlignedKmers);

		    	newAlignedKmers.clear();
	    	}
	    }
	    
      //continue to greedily grow the sequence cluster by building a PWM
      private boolean clusterByPWM(MotifCluster motifCluster, ArrayList<ComponentFeature> unalignedFeatures, double[] bg){
      	ArrayList<ComponentFeature> alignedFeatures = motifCluster.alignedFeatures;
    	ArrayList<Integer> motifStartInSeq = motifCluster.motifStartInSeq;
    	WeightMatrix wm = motifCluster.matrix;
    	
    	boolean noMore = true;
    	ArrayList<ComponentFeature> temp = new ArrayList<ComponentFeature>();
    	// merge kmers with duplicate string
    	HashMap<String, Kmer> alignedKmers = new HashMap<String, Kmer> ();
    	ArrayList<Kmer> dups = new ArrayList<Kmer>();
    	for (Kmer km: motifCluster.alignedKmers){
    		String kmStr = km.getKmerString();
    		if (alignedKmers.containsKey(kmStr)){
//    			System.err.println("ClusterByPWM: "+kmStr+" duplicated!");
    			alignedKmers.get(kmStr).mergeKmer(km);
    			dups.add(km);
    		}
    		else
    			alignedKmers.put(kmStr, km);
    	}
    	motifCluster.alignedKmers.removeAll(dups);
    	
    	double threshold = motifCluster.pwmThreshold;	
        WeightMatrixScorer scorer = new WeightMatrixScorer(wm);
    	ArrayList<Kmer> newAlignedKmers = new ArrayList<Kmer>();

        for (ComponentFeature cf: unalignedFeatures){
    	  String seq = cf.getBoundSequence();
    	  if (seq==null||seq.length()<wm.length()-1){
    		  temp.add(cf);
    		  continue;
    	  }
    	      	  
          WeightMatrixScoreProfile profiler = scorer.execute(seq);
          double maxSeqScore = Double.NEGATIVE_INFINITY;
          int maxScoringShift = 0;
          char maxScoringStrand = '+';
          for (int i=0;i<profiler.length();i++){
        	  double score = profiler.getMaxScore(i);
        	  if (maxSeqScore<score){
        		  maxSeqScore = score;
        		  maxScoringShift = i;
        		  maxScoringStrand = profiler.getMaxStrand(i);
        	  }
          }
          // if a sequence pass the motif score, reset the kmer to the binding site
          if (maxSeqScore >= threshold){
			if (maxScoringStrand =='-'){
			  cf.flipBoundSequence();
			  maxScoringShift = seq.length()-maxScoringShift-wm.length();	
			  // i.e.  (seq.length()-1)-maxScoringShift-(wm.length()-1);	
			}
			
			// get the k-mer
			seq = cf.getBoundSequence();
			int left = maxScoringShift+motifCluster.bindingPosition-config.k/2;
			int right = left+config.k;
			if (left<0||right>seq.length())		// if kmer matched to the end of bound sequence, skip
				continue;
			String kmerStr = cf.getBoundSequence().substring(left, right);
			
			// check if the k-mer is from negative set
//			if (kEngine.isNegativeKmer(kmerStr))
//				continue;
			
			alignedFeatures.add(cf);
			motifStartInSeq.add(maxScoringShift);
			temp.add(cf);
			
			if (alignedKmers.containsKey(kmerStr)){
				Kmer kmer = alignedKmers.get(kmerStr);
				kmer.incrSeqHitCount();
				kmer.incrStrength(cf.getTotalEventStrength());
				cf.setKmer(kmer);
			}
			else{
				Kmer kmer =  new Kmer(kmerStr, 1);
				kmer.setKmerStartOffset(left-maxScoringShift-motifCluster.bindingPosition);
				kmer.incrStrength(cf.getTotalEventStrength());
				alignedKmers.put(kmerStr, kmer);
				newAlignedKmers.add(kmer);
				kmer.setAlignString("PWM:"+WeightMatrix.getMaxLetters(wm));
				cf.setKmer(kmer);
			}
			noMore = false;
          }
        }

        unalignedFeatures.removeAll(temp);
    	temp.clear();	
    	
    	clusterByNewKmers(motifCluster, unalignedFeatures, newAlignedKmers);
    	
    	motifCluster.alignedKmers.addAll(newAlignedKmers);
    	
//    	System.err.println("Grow by PWM is done!");
    	return noMore;
      }
    	
    /** Used the newly PWM aligned kmers to cluster more sequences that contain these kmers
     *  This is to ensure all the sequences containing this kmer will be included in same cluster
     *  so that the kmers will be be split into different clusters.
    	For Kmers found by seed kmer and overlapped kmers, do not need to run this
    	This is only for Kmers found by PWM, because they are found through the bound sequence, not kmer
    */
    private boolean clusterByNewKmers(MotifCluster motifCluster, ArrayList<ComponentFeature> unalignedFeatures, Collection<Kmer> newKmers){
    	if (newKmers.isEmpty())
    		return true;
    	
//    	for (Kmer kmer:newKmers)
//    		System.out.println(kmer.getKmerString());
	
    	boolean noMore = true;
    	ArrayList<Kmer> kmers = new ArrayList<Kmer>();
    	kmers.addAll(newKmers);
    	KmerEngine engine = new KmerEngine(kmers, null);
    	ArrayList<ComponentFeature> temp = new ArrayList<ComponentFeature>();
		
    	for (ComponentFeature cf: unalignedFeatures){
      	  String seq = cf.getBoundSequence();
      	  if (seq==null)
      		  continue;
      	  
      	  HashMap<Integer, ArrayList<Kmer>> hits = engine.query(seq);
      	  if (hits.isEmpty())
      		  continue;
      	  
      	  noMore = false;
      	  
      	  // if multiple kmer match, take the hit nearest the middle of sequence
      	  int distToMiddle = seq.length();
      	  int nearstHit = seq.length();
      	  for (int pos:hits.keySet()){
      		  if (distToMiddle>Math.abs(Math.abs(pos)-seq.length()/2)){
      			distToMiddle = Math.abs(Math.abs(pos)-seq.length()/2);
      			nearstHit = pos;
      		  }
      	  }
      	  if (nearstHit==seq.length())
      		  continue;
      	  
      	  // if multiple kmers on that hit position, select the kmer with largest count (just to break the tie)
      	  Kmer selected = null;
      	  int count = 0;
      	  for (Kmer kmer: hits.get(nearstHit)){
      		  if (count<kmer.getSeqHitCount()){
      			  count = kmer.getSeqHitCount();
      			  selected = kmer;
      		  }
      	  }
      	  
      	  if (nearstHit<0){
      		  cf.flipBoundSequence();
      		  nearstHit=cf.getBoundSequence().indexOf(selected.getKmerString())-selected.getKmerStartOffset();
      	  }
      	  int motifStartPos = nearstHit-motifCluster.bindingPosition;
      	  if (motifStartPos<0){
      		  // this may happen if the binding position in motif is quite off-center
//      		  System.err.println("clusterByNewKmers: motifStartInSeq<0, centerHit="+centerHit+
//      				  ", motifCluster.bindingPosition="+motifCluster.bindingPosition+", kmer="+
//      				  selected.getKmerString()+ "\t" + cf.toString_v1());
      		  continue;
      	  }
      	  motifCluster.alignedFeatures.add(cf);
      	  motifCluster.motifStartInSeq.add(motifStartPos);
//      	  System.out.println(cf.getBoundSequence().substring(motifStartPos));
      	  cf.setKmer(selected);
      	  selected.incrSeqHitCount();
      	  selected.incrStrength(cf.getTotalEventStrength());
      	  temp.add(cf);
    	}
    	
    	unalignedFeatures.removeAll(temp);
    	
    	return noMore;
    }
    /**
     *  Scan the bound sequence of componentFeature using weight matrix
     *  @return  Pair of values, the start position of highest scoring PWM hit and the score
     *  The position will be negative if the match is on '-' strand    
     */
    private Pair<Integer, Double> scanPWM(String sequence, WeightMatrix wm, WeightMatrixScorer scorer){
		if (sequence==null||sequence.length()<wm.length()-1){
			return new Pair<Integer, Double>(-1,-1.0);
		}
		WeightMatrixScoreProfile profiler = scorer.execute(sequence);
		double maxSeqScore = Double.NEGATIVE_INFINITY;
		int maxScoringShift = 0;
		char maxScoringStrand = '+';
		for (int i=0;i<profiler.length();i++){
			double score = profiler.getMaxScore(i);
			if (maxSeqScore<score){
				maxSeqScore = score;
				maxScoringShift = i;
				maxScoringStrand = profiler.getMaxStrand(i);
			}
		}
  
		if (maxScoringStrand =='-'){
			maxScoringShift = -maxScoringShift;		
		}
		return new Pair<Integer, Double>(maxScoringShift, maxSeqScore);
    }
    
    /**
     *  Scan the bound sequence of componentFeature using weight matrix
     *  It scans outwards from the middle point, until a match pass the scoreTarget
     *  @return  Pair of values, the start position of nearest PWM hit and the score
     *  The position will be negative if the match is on '-' strand    
     *  If no match pass the scoreTarget, return -999 as position. The called need to check for this.
     */
    private Pair<Integer, Double> scanPWMoutwards(String sequence, WeightMatrix wm, WeightMatrixScorer scorer, int middle, double scoreTarget){
		if (sequence==null||sequence.length()<wm.length()-1){
			return new Pair<Integer, Double>(-999,-1.0);
		}
		WeightMatrixScoreProfile profiler = scorer.execute(sequence);
		double goodScore = 0;
		int goodScoreShift = -999;
		char goodScoreStrand = '+';
		int maxRange = Math.max(middle, profiler.length()-middle);
		for (int i=0;i<maxRange;i++){
			int idx = middle+i;
			if (idx<0 || idx>=sequence.length()-wm.length())
				continue;
			double score = profiler.getMaxScore(idx);
			if (score>=scoreTarget){
				goodScore = score;
				goodScoreShift = idx;
				goodScoreStrand = profiler.getMaxStrand(idx);
				break;
			}
			idx = middle-i;
			if (idx<0 || idx>=sequence.length()-wm.length())
				continue;
			score = profiler.getMaxScore(idx);
			if (score>=scoreTarget){
				goodScore = score;
				goodScoreShift = idx;
				goodScoreStrand = profiler.getMaxStrand(idx);
				break;
			}
		}
  
		if (goodScoreStrand =='-'){
			goodScoreShift = -goodScoreShift;		
		}
		return new Pair<Integer, Double>(goodScoreShift, goodScore);
    }

	// build PWM from Kmer cores, then update the motifStartInSeq using new PWM
    // The binding site is averaged from the GPS binding positions (relative to the PWM)
    private WeightMatrix makePWM(MotifCluster motifCluster, boolean isFromKmers){   	
    	ArrayList<ComponentFeature> alignedFeatures = motifCluster.alignedFeatures;
    	ArrayList<Integer> motifStartInSeq = motifCluster.motifStartInSeq;
    	// count
    	int seqLen = config.k_win+1;
    	int leftMost = seqLen;
    	int length = seqLen;	// the shortest length starting from the shift position
    	ArrayList<Integer> starts = new ArrayList<Integer>();
    	ArrayList<Integer> lengths = new ArrayList<Integer>();
    	for (int i=0;i<alignedFeatures.size();i++){
    		ComponentFeature cf = alignedFeatures.get(i);
    		String seq = cf.getBoundSequence();
    		int pos_motif = motifStartInSeq.get(i);			
    		if (pos_motif<0){
    			//TODO: still a few miss here
//    			System.err.println("Warning: makePWM(), pos_motif "+pos_motif+"<0,"+cf.toString_v1());
    			continue;
    		}

    		if (seq.length()-pos_motif<config.k){	
    			continue;
    		}
    		starts.add(pos_motif);
    		lengths.add(seq.length()-pos_motif);
    	}
    	Collections.sort(starts);
    	leftMost = starts.get(starts.size()/100);
    	Collections.sort(lengths);
    	length = lengths.get(lengths.size()/10);
    	
    	double sum_offsetXstrength = 0;
    	double sum_strength = 0;
    	double[][] pwm = new double[length][MAXLETTERVAL];
    	ArrayList<Integer> allOffsets = new ArrayList<Integer>();
    	StringBuilder pwm_seqs = new StringBuilder();
    	for (int i=0;i<alignedFeatures.size();i++){
 			double strength = config.use_strength?alignedFeatures.get(i).getTotalEventStrength():1;
 	 		String seq = alignedFeatures.get(i).getBoundSequence();
     		int pos_motif = motifStartInSeq.get(i);	
    		if (pos_motif<0 ||seq.length()-pos_motif<config.k)
    			continue;
			int start =pos_motif-leftMost;
    		if(start<0 || start+length>seq.length())
    			continue;
			String subSeq = seq.substring(start, start+length);
			if (seqLen==seq.length()){		// if length is not k_win+1, binding position may not be in middle
    			sum_offsetXstrength += strength*(config.k_win/2-start);
        		sum_strength += strength;
        		allOffsets.add(config.k_win/2-start);
			}
    		for (int p=0;p<length;p++){
    			char base = subSeq.charAt(p);
    			pwm[p][base] +=strength;
    		}
    		pwm_seqs.append(subSeq).append("\n");
    	}
    	
    	// average all the GPS binding positions to decide the binding position in the PWM
    	// weighted by strength if "use_strength" is true
    	int bPos=StatUtil.round(sum_offsetXstrength/sum_strength);		// mean
    	Collections.sort(allOffsets);
    	if (config.print_kmer_bPos){
	    	StringBuilder sb = new StringBuilder();
	    	for (int offset:allOffsets){
	    		sb.append(offset).append("\n");
	    	}
	    	ArrayList<Kmer> kmers = motifCluster.alignedKmers;
//	    	Collections.sort(kmers);
			Collections.sort(kmers, new Comparator<Kmer>(){
			    public int compare(Kmer o1, Kmer o2) {
			    		return o1.compareByHGP(o2);
			    }
			});
	    	String name = outName+"_PWM_"+kmers.get(0).getKmerString()+"_kmerShifts.txt";
	    	CommonUtils.writeFile(name, sb.toString());
    	}
//    	int bPos=allOffsets.get(allOffsets.size()/2);				// median

    	// make the PWM
		Pair<Integer, Integer> ends = constructPWM(pwm);
		int leftIdx = ends.car();
		int rightIdx = ends.cdr();
		
		bPos -= leftIdx;	// adjust PWM binding positin with the left_end trim
    	
    	// if the pwm is not good, return null. The operations in this method so far 
    	// does not change the state of componentFeatures or motifCluster, so we can discard this pwm and take previous result
    	if (rightIdx-leftIdx+1<=config.k/2){
    		motifCluster.isGood = false;
    		return null;
    	}
    	
    	float[][] matrix = new float[rightIdx-leftIdx+1][MAXLETTERVAL];   
    	for(int p=leftIdx;p<=rightIdx;p++){
    		for (int b=0;b<LETTERS.length;b++){
    			matrix[p-leftIdx][LETTERS[b]]=(float) pwm[p][LETTERS[b]];
    		}
    	}
    	WeightMatrix wm = new WeightMatrix(matrix);
		if (bPos>0)
			wm.consensus = CommonUtils.padding(bPos, ' ')+"|\n"+ WeightMatrix.printMatrixLetters(wm);
		else
			wm.consensus = WeightMatrix.printMatrixLetters(wm);    	
		
		// Check the quality of new PWM: hyper-geometric p-value test using the positive and negative sequences
    	double threshold = kEngine.estimatePwmThreshold(wm, config.wm_factor, outName);
    	threshold = Math.max(threshold, wm.getMaxScore()*config.wm_factor);
    	// if the pwm is not good, return null. The operations in this method so far 
    	// does not change the state of componentFeatures or motifCluster, so we can discard this pwm and take previous result
//    	if (threshold<0){
//    		motifCluster.isGood = false;
//    		System.out.println("makePWM: PWM is non-specific, stop here.");
//    		return null;
//    	}
  
    	/************************************************************************ 
    	 * update the motifStartInSeq for aligned sequences, w.r.t. PWM
    	 ************************************************************************/
        WeightMatrixScorer scorer = new WeightMatrixScorer(wm);

        // first get the offset of seed kmer relative to start of pwm
        ArrayList<Integer> seedOffsets = new ArrayList<Integer>();
		for (int p=0;p<alignedFeatures.size();p++){
        	ComponentFeature cf = alignedFeatures.get(p);
        	if (cf.getKmer()==motifCluster.seedKmer){
        		String seq = cf.getBoundSequence();
        		int pos_kmer = seq.indexOf(cf.getKmer().getKmerString());
            	int pos_pwm = scanPWMoutwards(seq, wm, scorer, pos_kmer, threshold).car();
//            	int pos_pwm = scanPWMoutwards(seq, wm, scorer, pos_kmer, wm.getMaxScore()*config.wm_factor).car();
            	if (pos_pwm==-999)
            		continue;
            	seedOffsets.add(pos_kmer-Math.abs(pos_pwm));
        	}
		}
		if (seedOffsets.isEmpty()){
    		System.out.println(String.format("makePWM: seed kmer do not pass PWM threshold %.2f, stop here.", threshold));
			return null;
		}
		
		int seedOffset = 0;
		for (int offset: seedOffsets)
			seedOffset += offset;
		seedOffset = StatUtil.round(seedOffset/seedOffsets.size());
		
        for (int p=0;p<alignedFeatures.size();p++){
        	ComponentFeature cf = alignedFeatures.get(p);
        	Kmer kmer = cf.getKmer();
        	if (kmer.getAlignString()==null)
        		continue;
        	if (!kmer.getAlignString().contains("PWM")){
        		// set motifStartInSeq relative to seed kmer
        		String seq = cf.getBoundSequence();
        		int pos_kmer = seq.indexOf(kmer.getKmerString());
        		motifStartInSeq.set(p, pos_kmer-seedOffset-kmer.getShift());		// update
        		continue;
        	}
        	String seq = cf.getBoundSequence();
//        	int pos = scanPWMoutwards(seq, wm, scorer, seq.indexOf(kmer.getKmerString()), threshold).car();
        	int pos = scanPWMoutwards(seq, wm, scorer, seq.indexOf(kmer.getKmerString()), wm.getMaxScore()*config.wm_factor).car();
        	// TODO: maybe should scan the kmer, need to consider the length of kmer and PWM
        	if (pos==-999){	// no match pass the target score, just get the best match
        		pos = scanPWM(cf.getBoundSequence(), wm, scorer).car();
        	}
        	if(pos<0){
        		cf.flipBoundSequence();
  			  	pos = seq.length()-wm.length()+pos;	
  			  	cf.getKmer().RC();
        	}
    		motifStartInSeq.set(p, pos);		// update
        }
        motifCluster.matrix = wm;
        motifCluster.pwmThreshold = threshold;
        motifCluster.bindingPosition = bPos;
        motifCluster.isGood = true;
        
    	return wm;
    }

    // make PFM from aligned sequences
	  private String getPFMString(ArrayList<ComponentFeature> alignedFeatures, ArrayList<Integer> motifStartInSeq, int length, int num){
	
		float[][] pfm = new float[length][MAXLETTERVAL];
		for (int i=0;i<alignedFeatures.size();i++){
			ComponentFeature cf = alignedFeatures.get(i);
			int pos_motif = motifStartInSeq.get(i);	
			if (pos_motif<0||pos_motif+length>cf.getBoundSequence().length())
				continue;
			String seq = cf.getBoundSequence().substring(pos_motif, pos_motif+length);		
			for (int p=0;p<length;p++){
				char base = seq.charAt(p);
				float strength = (float) (config.use_strength?alignedFeatures.get(i).getTotalEventStrength():1);
				pfm[p][base] +=strength;
			}
		}
		
		return makeTRANSFAC (pfm, String.format("DE %s_%d_c%d\n", outName, num, alignedFeatures.size()));		
	  }
	  
	  private String makeTRANSFAC (float[][] pfm, String header){
		// make string in TRANSFAC format
		StringBuilder msb = new StringBuilder();
		msb.append(header);
		for (int p=0;p<pfm.length;p++){
			msb.append(p+1).append(" ");
			int maxBase = 0;
			float maxCount=0;
			for (int b=0;b<LETTERS.length;b++){
				msb.append(String.format("%d ", (int)pfm[p][LETTERS[b]]));
				if (maxCount<pfm[p][LETTERS[b]]){
					maxCount=pfm[p][LETTERS[b]];
					maxBase = b;
				}
			}
			msb.append(LETTERS[maxBase]).append("\n");
		}
		msb.append("XX\n\n");
		return msb.toString();
	  }

	private int mismatch(String ref, String seq){
	  if (ref.length()!=seq.length())
		  return -1;
	  int mismatch = 0;
	  for (int i=0;i<ref.length();i++){
		  if (ref.charAt(i)!=seq.charAt(i))
			  mismatch ++;
	  }
	  return mismatch;
	}

	private class MotifCluster{
		int clusterId;
		Kmer seedKmer;
		ArrayList<ComponentFeature> alignedFeatures=new ArrayList<ComponentFeature>();
		// PWM hit start position in the bound sequence
		ArrayList<Integer> motifStartInSeq = new ArrayList<Integer>();	
		ArrayList<Kmer> alignedKmers = new ArrayList<Kmer>();;	
		WeightMatrix matrix;
		double pwmThreshold;
		int bindingPosition;
		boolean isGood=false;
		
		MotifCluster(){}	// empty constructor;

		protected MotifCluster clone(){
			MotifCluster cluster = new MotifCluster();
			cluster.seedKmer = this.seedKmer;
			cluster.alignedFeatures.addAll(this.alignedFeatures);
			cluster.motifStartInSeq.addAll(this.motifStartInSeq);
			cluster.alignedKmers.addAll(this.alignedKmers);
			cluster.matrix = this.matrix;
			cluster.bindingPosition = this.bindingPosition;
			cluster.pwmThreshold = this.pwmThreshold;
			return cluster;
		}
	}

	private class KmerCluster{
		int clusterId;
		Kmer seedKmer;
		String pfmString;
		WeightMatrix wm;
		double pwmThreshold;
		int pos_pwm_seed;
		int bindingPosition;
		int sequenceCount;
		int kmerCount;
		
		KmerCluster(){}	// empty constructor;

		protected KmerCluster clone(){
			KmerCluster cluster = new KmerCluster();
			cluster.clusterId = this.clusterId;
			cluster.seedKmer = this.seedKmer;
			cluster.pfmString = this.pfmString;
			cluster.wm = this.wm;
			cluster.pwmThreshold = this.pwmThreshold;
			cluster.pos_pwm_seed = this.pos_pwm_seed;
			cluster.bindingPosition = this.bindingPosition;
			cluster.sequenceCount = this.sequenceCount;
			cluster.kmerCount = this.kmerCount;
			
			return cluster;
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

    class GPSConfig {
		public boolean trim_simple=false;
		public boolean do_model_selection=false;
		public boolean classify_events = false;
        public boolean use_joint_event = false;
        public boolean TF_binding = true;
        public boolean outputBED = false;
        public boolean kmer_print_hits = false;
        public boolean testPValues = false;
        public boolean post_artifact_filter=false;
        public boolean filterEvents=true;
        public boolean filterDupReads=true;
        public boolean kl_count_adjusted = false;
        public boolean sort_by_location=false;
        public boolean exclude_unenriched = false;
        public boolean dump_regression = false;
        public boolean use_strength = false;
        public boolean print_kmer_bPos = false;
      	public int KL_smooth_width = 0;
        public int max_hit_per_bp = -1;
        
        // k-mer related
        public int k = -1;			// the width of kmer
        public int k_min = -1;		// the minimum value of k
        public int k_max= -1;		// the maximum value of k        
        public int k_seqs = 50000;	// the top number of event to get underlying sequences for initial Kmer learning 
        public int k_win = 60;		// the window around binding event to search for kmers
        public int k_neg_dist = 200;// the distance of the window for negative sequences from binding sites 
        public int k_shift = 99;	// the max shift from seed kmer when aligning the kmers     
        public int k_overlap = 7;	// the number of overlapped bases to assemble kmers into PWM    
        public int kpp_mode = 0;	// different mode to convert kmer count to positional prior alpha value
        public double hgp = 1e-3; 	// p-value threshold of hyper-geometric test for enriched kmer 
        public double k_fold = 2;	// the minimum fold of kmer count in positive seqs vs negative seqs
        public double gc = 0.42;	// GC content in the genome
        public double[] bg;			// background frequency based on GC content
        public double wm_factor = 0.5;		// The threshold relative to the maximum PWM score, for including a sequence into the cluster 
        public double ic_trim = 0.4;		// The information content threshold to trim the ends of PWM
        public double kmer_freq_pos_ratio = 0.8;	// The fraction of most frequent k-mer position in aligned sequences
        public int kmer_cluster_seq_count = 100;	// minimum number of sequences to be reported as a cluster, to build a PWM (for overlapping kmer)
        public int negative_ratio = 1; 		// The ratio of negative sequences to positive sequences
        public boolean kmer_use_insig = false;
        public boolean kmer_use_filtered = false;
        public boolean use_kmer_mismatch = true;
      	public boolean align_overlap_kmer=true;
        
        public double ip_ctrl_ratio = -1;	// -1: using non-specific region for scaling, -2: total read count for scaling, positive: user provided ratio
        public double q_value_threshold = 2.0;	// -log10 value of q-value
        public double q_refine = 1.5;
        public double joint_event_distance = 500;
        public double alpha_factor = 3.0;
        public double excluded_fraction = 0.05;	// top and bottom fraction of region read count to exclude for regression
        public int top_events = 2000;
        public int min_event_count = 500;	// minimum num of events to update read distribution
        public int smooth_step = 30;
        public int window_size_factor = 3;	//number of model width per window
        public int min_region_width = 50;	//minimum width for select enriched region
        public double mappable_genome_length = -1; // defalut is to compute
        public double sparseness=6.0;
        public double fold = 3.0;
        public double kl_ic = 0.0;
        public double shapeDeviation = 0;
        public int gentle_elimination_factor = 2;	// factor to reduce alpha to a gentler pace after eliminating some component
        public int resolution_extend = 2;
        public int first_lambda_region_width  =  1000;
        public int second_lambda_region_width =  5000;
        public int third_lambda_region_width  = 10000;
        public boolean use_dynamic_sparseness = true;
        public boolean use_betaEM = true;
        public boolean use_scanPeak  = true;
        public boolean refine_regions = false;		// refine the enrichedRegions for next round using EM results
        public boolean cache_genome = true;			// cache the genome sequence
        public int bmverbose=1;		// BindingMixture verbose mode
        public int base_reset_threshold = 200;	// threshold to set a base read count to 1
        public int windowSize;			// size for EM sliding window for splitting long regions
        //Run EM up until <tt>ML_ITER</tt> without using sparse prior
        public int ML_ITER=10;
        // the range to scan a peak if we know position from EM result
        public int SCAN_RANGE = 20;
        public int gentle_elimination_iterations = 5;
        public double minFoldChange = 1.5; // minimum fold change for a significant event.  applied after event discovery during the p-value filtering stage

        public void parseArgs(String args[]) {
            Set<String> flags = Args.parseFlags(args);
            // default as false, need the flag to turn it on
            classify_events = flags.contains("classify");
            sort_by_location = flags.contains("sl");
            use_joint_event = flags.contains("refine_using_joint_event");
            post_artifact_filter = flags.contains("post_artifact_filter");
            kl_count_adjusted = flags.contains("adjust_kl");
            refine_regions = flags.contains("refine_regions");
            outputBED = flags.contains("outBED");
            testPValues = flags.contains("testP");
            if (testPValues)
            	System.err.println("testP is " + testPValues);
            exclude_unenriched = flags.contains("ex_unenriched");
            dump_regression = flags.contains("dump_regression");
            use_strength = flags.contains("use_strength");
            kmer_print_hits = flags.contains("kmer_print_hits");
            kmer_use_insig = flags.contains("kmer_use_insig");
            kmer_use_filtered = flags.contains("kmer_use_filtered");
          
                // default as true, need the opposite flag to turn it off
            use_dynamic_sparseness = ! flags.contains("fa"); // fix alpha parameter
            use_betaEM = ! flags.contains("poolEM");
            filterEvents = !flags.contains("nf");	// not filtering of predicted events
            filterDupReads = !flags.contains("nrf");	// no read filtering of duplicate reads
            TF_binding = ! flags.contains("br");	// broad region, not TF data, is histone or pol II
            if (!TF_binding){
                use_joint_event = true;
                sort_by_location = true;
            }
            use_scanPeak = ! flags.contains("no_scanPeak");
            do_model_selection = !flags.contains("no_model_selection");
            use_kmer_mismatch = !flags.contains("no_kmm");
            align_overlap_kmer = !flags.contains("no_aok");

            mappable_genome_length = Args.parseDouble(args, "s", mappable_genome_length);	// size of mappable genome
           
            // Optional input parameter
            k = Args.parseInteger(args, "k", k);
            k_min = Args.parseInteger(args, "k_min", k_min);
            k_max = Args.parseInteger(args, "k_max", k_max);
            k_seqs = Args.parseInteger(args, "k_seqs", k_seqs);
            k_win = Args.parseInteger(args, "k_win", k_win);
            k_neg_dist = Args.parseInteger(args, "k_neg_dist", k_neg_dist);
            k_shift = Args.parseInteger(args, "k_shift", k_shift);
            k_overlap = Args.parseInteger(args, "k_overlap", Math.max(k_overlap, StatUtil.round(k*0.75)));
            kpp_mode = Args.parseInteger(args, "kpp_mode", kpp_mode);
            k_fold = Args.parseDouble(args, "k_fold", k_fold);
            gc = Args.parseDouble(args, "gc", gc);
        	bg = new double[4]; bg[0]=0.5-gc/2; bg[1]=gc/2; bg[2]=bg[1]; bg[3]=bg[0];
            wm_factor = Args.parseDouble(args, "wmf", wm_factor);
            ic_trim = Args.parseDouble(args, "ic", ic_trim);
            hgp = Args.parseDouble(args, "hgp", hgp);
            kmer_freq_pos_ratio = Args.parseDouble(args, "kmer_freq_pos_ratio", kmer_freq_pos_ratio);
            kmer_cluster_seq_count = Args.parseInteger(args, "cluster_seq_count", kmer_cluster_seq_count);

            ip_ctrl_ratio = Args.parseDouble(args, "icr", ip_ctrl_ratio);
            maxThreads = Args.parseInteger(args,"t",java.lang.Runtime.getRuntime().availableProcessors());	// default to the # processors
            q_value_threshold = Args.parseDouble(args, "q", q_value_threshold);	// q-value
            q_refine = Args.parseDouble(args, "q2", q_refine);	// q-value for refine regions
            sparseness = Args.parseDouble(args, "a", 6.0);	// minimum alpha parameter for sparse prior
            alpha_factor = Args.parseDouble(args, "af", alpha_factor); // denominator in calculating alpha value
            fold = Args.parseDouble(args, "fold", fold); // minimum fold enrichment IP/Control for filtering
            shapeDeviation =  TF_binding?-0.4:-0.3;		// set default according to filter type    		
            shapeDeviation = Args.parseDouble(args, "sd", shapeDeviation); // maximum shapeDeviation value for filtering
            max_hit_per_bp = Args.parseInteger(args, "mrc", 0); //max read count per bp, default -1, estimate from data
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
            excluded_fraction = Args.parseDouble(args, "excluded_fraction", excluded_fraction);
            kl_ic = Args.parseDouble(args, "kl_ic", kl_ic);
            resolution_extend = Args.parseInteger(args, "resolution_extend", resolution_extend);
            gentle_elimination_factor = Args.parseInteger(args, "gentle_elimination_factor", gentle_elimination_factor);
            minFoldChange = Args.parseDouble(args,"min_fold_change",minFoldChange);
            // These are options for EM performance tuning
            // should NOT expose to user
            // therefore, still use UPPER CASE to distinguish
            ML_ITER = Args.parseInteger(args, "ML_ITER", ML_ITER);
            SCAN_RANGE = Args.parseInteger(args, "SCAN_RANGE", SCAN_RANGE);
        }
    }
    class GPS2Thread implements Runnable {

        private KPPMixture mixture;
        private GPSConstants constants;
        private GPSConfig config;
        private Collection<Region> regions;
        private Collection<Integer> processedRegionCount;
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

        public GPS2Thread (Collection<Region> regions,
        				  Collection<Integer> processedRegionCount,
                          Collection<ComponentFeature> compFeatures,
                          Collection<KmerPP> allKmerHits,
                          KPPMixture mixture,
                          GPSConstants constants,
                          GPSConfig config,
                          boolean isIP) {
            this.regions = regions;
            this.processedRegionCount = processedRegionCount;
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
                    result = EMTrain(signals, alpha, new double[r.getWidth()]);
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
        	if (config.cache_genome)
        		seqgen.useCache(true);
            for (Region rr : regions) {
                mixture.log(2, rr.toString());
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
                                            int winSize = mixture.modelWidth*10;
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
                    if (kEngine!=null){
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
                    System.exit(-1);
                }
            }
        }

        private ArrayList<BindingComponent> analyzeWindow(Region w, SequenceGenerator<Region> seqgen){

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

                	//construct the positional prior for each position in this region
                	double[] pp = new double[w.getWidth()+1];
                	Kmer[] pp_kmer = new Kmer[pp.length];
                	String seq = null;
                	if (kEngine!=null && kEngine.isInitialized()){
                		if (kmerPreDefined)
                			seq = seqgen.execute(w).toUpperCase();
                		else
                			seq = kEngine.getSequence(w);
	            	    HashMap<Integer, ArrayList<Kmer>> kmerHits = kEngine.query(seq);
// TODO: allowing more motif in the region	                	
//	                	double kmerCount_max = kEngine.getMaxCount();
//	                	double kmerCount_min = kEngine.getMinCount();
//	                	for (int pos: kmerHits.keySet()){
//	                		// the pos is the start position, hence +k/2
//	                		//if pos<0, then the reverse compliment of kmer is matched
//	                		int bindingPos = Math.abs(pos)+this.config.k/2;
//	                		int kmerCount = kmerHits.get(pos).getSeqHitCount();
//	                		// select the approach to generate pp from kmer count
//		                	// scale positional prior pseudo-count to be in the range of [kpp_min, kpp_max]*sparseness
//	                		if (config.kc2pp==0)
//	                			pp[bindingPos] = ((kmerCount-kmerCount_min)*(config.kpp_max-config.kpp_min)/(kmerCount_max-kmerCount_min)+config.kpp_min)*config.sparseness;
//	                		else if (config.kc2pp==1)
//	                			pp[bindingPos] = ((kmerCount==0?0:Math.log(kmerCount)-Math.log(kmerCount_min))*(config.kpp_max-config.kpp_min)/(Math.log(kmerCount_max)-Math.log(kmerCount_min))+config.kpp_min)*config.sparseness;
//	                		
//	                		pp_kmer[bindingPos] = kmerHits.get(pos);
	                		
	                	// Effectively, the top kmers will dominate, because we normalize the pp value
                		HashMap<Integer, KmerPP> hits = new HashMap<Integer, KmerPP>();
	                	for (int pos: kmerHits.keySet()){
	                		// the pos is the motif middle position
	                		//if pos<0, then the reverse compliment of kmer is matched
	                		int bindingPos = Math.abs(pos);
	                		if (bindingPos>=pp.length){
//	                			System.err.println("KPP: bindingPos " + bindingPos + " out of bound ("+pp.length+") at "+w.toString());
	                			continue;
	                		}
//	                		// The nearest kmer will take counts from all others (on the same position)
	                		// This may resulted in splitting some kmer into overlapping kmers, so we use largest kmer as below
//	                		double kmerCountSum = 0;
//	                		Kmer nearestKmer = null;
//	                		int distance = 1000;
//	                		for (Kmer kmer:kmerHits.get(pos)){		
//	                			double count=0;
//		                		if (config.use_strength)
//		                			count = kmer.getStrength();
//		                		else
//		                			count = kmer.getSeqHitCount();
//		                		kmerCountSum+=count;
//		                		if (distance>Math.abs(kmer.getKmerShift()+config.k/2)){
//		                			distance = Math.abs(kmer.getKmerShift()+config.k/2);
//		                			nearestKmer = kmer;
//		                		}		                			
//	                		}
	                		// The largest kmer will take counts from all others (on the same position)	
	                		double kmerCountSum = 0;
	                		Kmer largetKmer = null;
	                		double largetstCount = 0;
	                		for (Kmer kmer:kmerHits.get(pos)){		
	                			double count=0;
		                		if (config.use_strength){
		                			count = kmer.getStrength();
		                			if (count<1)				// on first kpp round, kmers do not have strength value, use count here
		                				count = kmer.getSeqHitCount();
		                		}
		                		else
		                			count = kmer.getSeqHitCount();
		                		kmerCountSum+=count;
		                		if (largetstCount<count){
		                			largetstCount = count;
		                			largetKmer = kmer;
		                		}		                			
	                		}
	                		// select the approach to generate pp from kmer count
	                		if (config.kpp_mode==0)
	                			pp[bindingPos] = kmerCountSum;
	                		else if (config.kpp_mode==1)
	                			pp[bindingPos] = kmerCountSum==0?0:Math.log(kmerCountSum);
	                		else if (config.kpp_mode==10)
	                			pp[bindingPos] = kmerCountSum==0?0:Math.log10(kmerCountSum);
	                		pp_kmer[bindingPos] = largetKmer;
	                		hits.put(bindingPos, new KmerPP(new Point(gen, w.getChrom(), w.getStart()+bindingPos), largetKmer, pp[bindingPos]));
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
	                	
	                	// normalize the local region so that total positional prior pseudo-count equal to alpha (sparse prior)
	                	// alternatively, we can scale so that (largest position prior == sparse prior)
	                	for (int i=0;i<boundaries.size()-1;i++){
	                		int total = 0;
	                		for (int j=boundaries.get(i);j<boundaries.get(i+1); j++){
		                		total += pp[j];
	                		}
	                		for (int j=boundaries.get(i);j<boundaries.get(i+1); j++){
		                		if (pp[j]>0){
		                			pp[j] = pp[j]/total*(alpha-1);
		                			hits.get(j).pp = pp[j];
	                			}
	                		}
	                	}
	                	if (config.kmer_print_hits)
	                		allKmerHits.addAll(hits.values());
                	}
                	
                    //Run EM and increase resolution
                    initializeComponents(w, mixture.numConditions);
                    int countIters = 0;
                    while(nonZeroComponentNum>0){                        
                        int numComp = components.size();
                        double[] p_alpha = new double[numComp];						// positional alpha
                        if (kEngine!=null){
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
	                        	p_alpha[i]=maxPP;
	                        }
                        }
                        // EM learning, components list will only contains non-zero components
                        result = EMTrain(signals, alpha, p_alpha);
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
                    if (nonZeroComponentNum==0)	return null;
                    setComponentResponsibilities(signals, result.car(), result.cdr());
                    if (kEngine!=null)
                    	setComponentKmers(pp_kmer, w.getStart(), seq);
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
        private Pair<double[][][], int[][][]>  EMTrain(ArrayList<List<StrandedBase>> signals, double alpha, double[] p_alpha){
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
            double[] pi = new double[numComp];							// Pi

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
            Pair<double[][][], int[][][]> result = EM_MAP(counts, h, r, b, b2c, c2b, pi, alpha, p_alpha);
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
                                                        double alpha,
                                                        double[] p_alpha) {
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
                            pi[j]=Math.max(0,r_sum-currAlpha+p_alpha[j]);

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
                            r_sum[jnz] += p_alpha[nzComps.get(jnz)];	// adding positional prior as pseudo-count
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
                        LP+=(p_alpha[j]-currAlpha)*Math.log(pi[j]);		// positional prior and sparse prior

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
                        state.LL = LAP;		//TODO: LL or LAP?
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

        // matched EM resulted binding components with the kmer prior
        //TODO: maybe we should search closest (<k/4) pp_kmer because some position may end up out-competed by nearby positions
        private void setComponentKmers(Kmer[] pp_kmer, int startPos, String seq){
        	for (BindingComponent b:components){
        		int kIdx = b.getLocation().getLocation()-startPos;
        		b.setKmer(pp_kmer[kIdx]);
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
	        		if (b.getKmer()!=null){
	        			if (!bs.contains(b.getKmer().getKmerString()))
	        				bs = SequenceUtils.reverseComplement(bs);
	        			if (!bs.contains(b.getKmer().getKmerString()))
	        				b.setKmer(null);
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
    
    class KmerPP implements Comparable<KmerPP>{
    	Point coor;
    	Kmer kmer;
    	double pp;
    	public KmerPP(Point coor, Kmer kmer, double pp) {
			this.coor = coor;
			this.kmer = kmer;
			this.pp = pp;
		}
		public String toString(){
    		return String.format("%s\t%s\t%d\t%.3f\t%.3f", coor.getLocationString(), kmer.getKmerString(), kmer.getSeqHitCount(),kmer.getStrength(), pp);
    	}
		public int compareTo(KmerPP h) {
			return coor.compareTo(h.coor);
		}
    }
 }
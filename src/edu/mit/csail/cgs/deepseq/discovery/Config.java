package edu.mit.csail.cgs.deepseq.discovery;

import java.util.Set;

import edu.mit.csail.cgs.tools.utils.Args;

public class Config {
	public boolean multl_thread_hgp = false;
	public boolean trim_simple=false;
	public boolean print_PI = false;
	public boolean classify_events = false;
    public boolean use_joint_event = false;
    public boolean outputBED = false;
    public boolean outputNarrowPeak = false;
    public boolean outputMEME = false;
    public boolean outputHOMER = false;
    public boolean outputJASPAR = false;
    public boolean write_RSC_file = false;
    public boolean write_genetrack_file = false;
    public boolean kmer_print_hits = false;
    public boolean testPValues = false;
    public boolean post_artifact_filter=false;
    public boolean kl_count_adjusted = false;
    public boolean sort_by_location=false;
    public boolean dump_regression = false;
    public boolean use_kmer_strength = false;
    public boolean print_kmer_bPos = false;    
    public boolean discard_subAlpha_components=false;			// discard the component whose strength is less than alpha    
    public boolean process_all_regions = false;
    public boolean refine_window_boundary = false;
    public boolean is_branch_point_data = false;

    public boolean TF_binding = true;
    public boolean exclude_unenriched = true;
    public boolean filterEvents=true;
    public boolean filterDupReads=true;
	public boolean do_model_selection=true;
	/** strand_type <br>
	 * run event finding and motif discovery as <br>
	 * 0) unstranded, for typical ChIP-seq data <br>
	 * 1) only single strand, for BranchPoint data (Burge Lab), run motif discovery only on event strand, 
	 * 2) ChIP-exo (if wanting to find peak boundary), run motif discovery on both strands 
	 */
    public int strand_type = 0;		
  	public int KL_smooth_width = 0;
    public int max_hit_per_bp = -1;
    public int maxThreads;		// default to #CPU
    // k-mer related
    public int k = -1;			// the width of kmer
    public int k_min = -1;		// the minimum value of k
    public int k_max= -1;		// the maximum value of k        
    public String seed = null;
    public int seq_weight_type = 3;	// 0: no weighting, 1: strength weighting, 2: sqrt(strength) weighting 3: ln(strength) weighting
    
    /** number of top k-mers selected from density clustering to run KMAC */
    public int k_top = 3;
    /** kmer distance cutoff, kmers with smaller or equal distance are consider neighbors when computing local density, in density clustering */
    public int dc = 3;
    /** delta value cutoff, k-mers with equal or higher delta values are used for selecting cluster centers */
    public int delta = 2;
    
    public int k_seqs = 5000;	// the top number of event to get underlying sequences for initial Kmer learning 
    public int k_win = 61;		// the window around binding event to search for kmers
    public int k_win2 = 101;	// the window around binding event to search for maybe secondary motifs (in later rounds)
    public int k_win_f = 4;		// k_win = k_win_f * k
    public int k_neg_dist = 300;// the distance of the nearest edge of negative region from binding sites 
    public int k_negSeq_ratio = 2; 		// The ratio of cache negative sequences to positive sequences
    public int k_shift = 99;	// the max shift from seed kmer when aligning the kmers     
    public int max_cluster = 20;
    public int topKmer_trials = 5;	// the number of initial k-mers to try
    public float k_mask_f = 1;	// the fraction of PWM to mask
    public int kpp_mode = 0;	// different mode to convert kmer count to positional prior alpha value
    public double hgp = -3; 	// p-value threshold of hyper-geometric test for enriched motif 
    public double kmer_hgp = -3; 	// p-value threshold of hyper-geometric test for enriched kmer 
    public double k_fold = 2;	// the minimum fold of kmer count in positive seqs vs negative seqs
    public double gc = -1;	// GC content in the genome			//0.41 for human, 0.42 for mouse
    public double[] bg= new double[4];	// background frequency based on GC content
    public double wm_factor = 0.6;		// The threshold relative to the maximum PWM score, for including a sequence into the cluster 
    public double fpr = 0.15;		// The false positive rate for partial ROC
    public double motif_relax_factor = 1;		// A factor to multiply the motif PWM threshold, used for GEM motif positional prior
    public double ic_trim = 0.4;		// The information content threshold to trim the ends of PWM
//    public double kmer_freq_pos_ratio = 0.8;	// The fraction of most frequent k-mer position in aligned sequences
    public double motif_hit_factor = 0.005;
    public double motif_hit_factor_report = 0.05;
//    public double kmer_set_overlap_ratio = 0.5;
    public double motif_hit_overlap_ratio = 0.5;
    public double repeat_fraction=1;		// ignore lower case letter and N in motif discovery if less than _fraction_ of sequence
    public int kmer_remove_mode = 0;
    public int seed_range = 3;
    public double kmer_inRange_fraction = 0.3;		// the fraction of kmer in the seed_range out of all k-mer hit count
    public double kmer_consistent_fraction = 0.5;		// the fraction of consistently aligned kmers in the seed_range
    public boolean use_grid_search = true;
    public boolean optimize_pwm_threshold = true;
    public boolean optimize_kmer_set = true;
    public boolean kmer_use_insig = false;
    public boolean kmer_use_filtered = false;
    public boolean use_weighted_kmer = true;		// strength weighted k-mer count
    public boolean use_pos_kmer = true;				// position weighted k-mer count
    public boolean k_neg_shuffle = false;
    public boolean k_neg_dinu_shuffle = false;		// di-nuleotide shuffle
    public int rand_seed = 0;
    public double neg_pos_ratio = 1;					// number of shuffled negative / positive seqs
   	public boolean re_align_kmer = false;
   	public boolean use_kmer_mismatch = true;
   	public boolean use_seed_family = true;		// start the k-mer alignment with seed family (kmers with 1 or 2 mismatch)
   	/** Align and cluster motif using KSM */
   	public boolean use_ksm = true;	
   	public boolean use_gapped = true;
   	public boolean use_sub_kmers = true;
 	public boolean estimate_ksm_threshold = true;
  	public boolean kpp_normalize_max = true;
  	public boolean pp_use_kmer = true;			// position prior using k-mer(true) or PWM(false)
  	public boolean progressive_PWM_trim = true;
  	public double kpp_factor = 0.8;
  	public int kpp_nmotifs = 1;	// number of motifs to use for kpp setup
  	public double noise = 0.0;
    public boolean print_aligned_seqs = false;
    public boolean print_input_seqs = false;
    public boolean print_all_kmers = false;
    public boolean print_bound_seqs = false;
    public boolean re_train = false;
	public boolean refine_pwm = false;
	public boolean refine_ksm = false;	// refine the KSM at the end of KMAC using un-masked sequences
	public boolean refine_final_motifs = false;	// refine the final motifs
	public boolean optimize_KG_kmers = false;	
    public boolean print_pwm_fdr = false;
    /** whether to use K-mer Set Model to evaluate improvement of new cluster, default to use PWM */
    public boolean evaluate_by_ksm = false;	
    public boolean evaluate_by_both = false;	
    
    public boolean use_strength_weight = true;	// use binding event strength to weight 
    public boolean use_pos_weight = false;		// use binding position profile to weight motif site
    public boolean allow_seed_reset=true;		// reset primary motif if secondary motif is more enriched
    public boolean allow_seed_inheritance=true;	// allow primary seed k-mer to pass on to the next round of GEM
    public boolean filter_pwm_seq = true;
//    public boolean k_select_seed = false;
    public boolean pwm_align_new = true;		// use PWM to align only un-aligned seqs (vs. all sequences)
    public boolean strigent_event_pvalue = true;// stringent: binomial and poisson, relax: binomial only
    public boolean use_db_genome = false;// get the sequence from database, not from file
    public boolean k_mask_1base = false;
    public boolean selectK_byTopKmer = false;
    
    public double ip_ctrl_ratio = -1;	// -1: using non-specific region for scaling, -2: total read count for scaling, positive: user provided ratio
    public double q_value_threshold = 2.0;	// -log10 value of q-value
    public double q_refine = -1;
    public double alpha_factor = 3.0;
    public double alpha_fine_factor = 2.0;
    public double excluded_fraction = 0.05;	// top and bottom fraction of region read count to exclude for regression
    public int top_events = 2000;
    public int min_event_count = 500;	// minimum num of events to update read distribution
    public double top_fract_to_skip = 0.01;		// fraction of top events to skip for updating read distribution
    public int smooth_step = 30;
    public int window_size_factor = 3;	//number of model width per window
    public int min_region_width = 50;	//minimum width for select enriched region
    public int min_event_distance = 1;		//minimum distance for nearby events (if less, events are merged)
    public int noise_distribution = 1;	// the read distribution for noise component, 0 NO, 1 UNIFORM, 2 SMOOTHED CTRL
    public String out_name="out";
    
    public double mappable_genome_length = -1; // default is to compute
    public double background_proportion = -1;	// default is to compute
    public double pi_bg_r0 = 0.02;
    public double sparseness=-1;
    public double poisson_alpha=1e-3; 				// the Poisson p-value for estimating alpha
    public double fold = 2.5;
    public double kl_ic = -2.0;
    public double shapeDeviation;
    public int gentle_elimination_factor = 2;	// factor to reduce alpha to a gentler pace after eliminating some component
    public int resolution_extend = 2;
    public int second_lambda_region_width =  5000;
    public int third_lambda_region_width  = 10000;
    public boolean bic = false;				// use BIC or AIC for model selection
//    public boolean model_noise = false;		// have a noise component for background reads
    public boolean ML_speedup = false;		
    public boolean use_dynamic_sparseness = true;
    public boolean use_scanPeak  = true;
    public boolean refine_regions = false;		// refine the enrichedRegions for next round using EM results
    public boolean print_stranded_read_distribution = false;	
    public boolean cache_genome = true;			// cache the genome sequence
    public String genome_path = null;
    
    public int verbose=1;		// BindingMixture verbose mode
    public int base_reset_threshold = 200;	// threshold to set a base read count to 1
    public int windowSize;	// size for EM window, for splitting long regions, set modelWidth * config.window_size_factor in KPPMixture
    //Run EM up until <tt>ML_ITER</tt> without using sparse prior
    public int ML_ITER=10;
    // the range to scan a peak if we know position from EM result
    public int SCAN_RANGE = 20;
    public int gentle_elimination_iterations = 5;
    public double minFoldChange = 1.5; // minimum fold change for a significant event.  applied after event discovery during the p-value filtering stage

    public void parseArgs(String args[]) throws Exception {
        Set<String> flags = Args.parseFlags(args);
        is_branch_point_data = flags.contains("bp");
        if (is_branch_point_data){
        	strand_type = 1;	
        	min_region_width = 1;
        	min_event_distance = 3;
        	smooth_step = 0;
//        	kmer_inRange_fraction = 0.3;
        	noise_distribution = 0;
        	window_size_factor = 10;			//TODO: why 10??
        	alpha_factor = 0.7;
        	motif_relax_factor = 0.5;		// relax to half of the threshold
        }

        // default as false, need the flag to turn it on
        classify_events = flags.contains("classify");
        print_PI = flags.contains("print_PI");
        sort_by_location = flags.contains("sl");
        use_joint_event = flags.contains("use_joint_event");
        post_artifact_filter = flags.contains("post_artifact_filter");
        kl_count_adjusted = flags.contains("adjust_kl");
        outputBED = flags.contains("outBED");
        outputNarrowPeak = flags.contains("outNP");
        outputMEME = flags.contains("outMEME");
        outputHOMER = flags.contains("outHOMER");
        outputJASPAR = flags.contains("outJASPAR");
        write_RSC_file = flags.contains("writeRSC");
        write_genetrack_file = flags.contains("write_genetrack_file");
        testPValues = flags.contains("testP");
        if (testPValues)
        	System.err.println("testP is " + testPValues);
        dump_regression = flags.contains("dump_regression");
        progressive_PWM_trim = !flags.contains("npt");		// no progressive trimming of PWM
        use_kmer_strength = flags.contains("use_kmer_strength");
        kmer_print_hits = flags.contains("kmer_print_hits");
        kmer_use_insig = flags.contains("kmer_use_insig");
        k_neg_shuffle = flags.contains("k_neg_shuffle");
        k_neg_dinu_shuffle = flags.contains("k_neg_dinu_shuffle");
        rand_seed = Args.parseInteger(args, "rand_seed", rand_seed);
        neg_pos_ratio = Args.parseDouble(args, "neg_pos_ratio", neg_pos_ratio);
        re_align_kmer = flags.contains("rak");
        print_aligned_seqs = flags.contains("print_aligned_seqs");
        print_input_seqs = flags.contains("print_input_seqs");
        print_all_kmers = flags.contains("print_all_kmers");
        print_bound_seqs = flags.contains("print_bound_seqs");
        re_train = flags.contains("re_train");
        refine_pwm = flags.contains("refine_pwm");
        refine_ksm = flags.contains("refine_ksm");
        refine_final_motifs = flags.contains("refine_final_motifs");
        optimize_KG_kmers = flags.contains("optimize_KG_kmers");
        print_pwm_fdr = flags.contains("print_pwm_fdr");      
        use_db_genome = flags.contains("use_db_genome");
        evaluate_by_ksm = flags.contains("evaluate_by_ksm");
        evaluate_by_both = flags.contains("evaluate_by_both");
        k_mask_1base = flags.contains("k_mask_1base");
        bic = flags.contains("bic"); 					// BIC or AIC
        discard_subAlpha_components = flags.contains("no_sub_alpha");
        refine_regions = flags.contains("refine_regions");
        print_stranded_read_distribution = flags.contains("print_stranded_read_distribution");
        process_all_regions = flags.contains("process_all_regions");
        refine_window_boundary = flags.contains("refine_window_boundary");
        use_pos_weight = flags.contains("use_pos_weight");
        
        // default as true, need the opposite flag to turn it off
        exclude_unenriched = !flags.contains("not_ex_unenriched");
        use_dynamic_sparseness = ! flags.contains("fa"); // fix alpha parameter
        filterEvents = !flags.contains("nf");	// no filtering for predicted events
        filterDupReads = !flags.contains("nrf");	// no read filtering of duplicate reads
        TF_binding = ! flags.contains("br");	// broad region, not TF data, is histone or pol II
        if (!TF_binding){
            use_joint_event = true;
            sort_by_location = true;
        }
        ML_speedup = !flags.contains("no_fast_ML");
        use_scanPeak = ! flags.contains("no_scanPeak");
        do_model_selection = !flags.contains("no_model_selection");
        use_kmer_mismatch = !flags.contains("no_kmm");
        use_seed_family = !flags.contains("no_seed_family");
        use_ksm = !flags.contains("no_ksm");
        use_gapped = !flags.contains("ng");
        use_sub_kmers = !flags.contains("ns");
        pp_use_kmer = !flags.contains("pp_pwm");
        estimate_ksm_threshold = !flags.contains("no_ksm_threshold");
        use_weighted_kmer = !flags.contains("no_weighted_kmer");
        use_pos_kmer = !flags.contains("no_pos_kmer");
        optimize_pwm_threshold = !flags.contains("not_optimize_pwm_threshold");
        optimize_kmer_set = !flags.contains("not_optimize_kmer_set");
        use_grid_search = !flags.contains("no_grid_search");
        allow_seed_reset = !flags.contains("no_seed_reset");
        selectK_byTopKmer = flags.contains("selectK_byTopKmer");	
        if (selectK_byTopKmer)														// overwrite allow_seed_reset
        	allow_seed_reset = false;
        allow_seed_inheritance = !flags.contains("no_seed_inheritance");
        pwm_align_new = !flags.contains("pwm_align_all");
        filter_pwm_seq = !flags.contains("pwm_seq_asIs");
        strigent_event_pvalue = !flags.contains("relax");

        mappable_genome_length = Args.parseDouble(args, "s", mappable_genome_length);	// size of mappable genome
        background_proportion = Args.parseDouble(args, "pi_bg", background_proportion);	// proportion of background read signal
        pi_bg_r0 = Args.parseDouble(args, "pi_bg_r0", pi_bg_r0);	// proportion of background read signal for round 0
        
        // Optional input parameter
        genome_path = Args.parseString(args, "genome", genome_path);
        out_name = Args.parseString(args, "out_name", out_name);
        k = Args.parseInteger(args, "k", k);
        if (k==-1){
	        k_min = Args.parseInteger(args, "k_min", k_min);
	        k_max = Args.parseInteger(args, "k_max", k_max);
	        k_min = Args.parseInteger(args, "kmin", k_min);				// fail-safe
	        k_max = Args.parseInteger(args, "kmax", k_max);
	        if (k_max<k_min)
	        	System.err.println("\n\nWARNING: k_max value is smaller than k_min !!!\n\n");
        }
        else{
        	k_min = k;
        	k_max = k;	
        }
        seed = Args.parseString(args, "seed", null);
        if (seed!=null){
        	k = seed.length();
        	k_min = k;
        	k_max = k;
        	allow_seed_reset = false;
        }
        k_top = Args.parseInteger(args, "k_top", k_top);
        dc = Args.parseInteger(args, "dc", dc);
        delta = Args.parseInteger(args, "delta", delta);
        k_seqs = Args.parseInteger(args, "k_seqs", k_seqs);
        seq_weight_type = Args.parseInteger(args, "swt", seq_weight_type);
        k_win = Args.parseInteger(args, "k_win", k_win);
        k_win_f = Args.parseInteger(args, "k_win_f", k_win_f);
        k_neg_dist = Args.parseInteger(args, "k_neg_dist", k_neg_dist);
        k_negSeq_ratio = Args.parseInteger(args, "k_neg_ratio", k_negSeq_ratio);
        k_shift = Args.parseInteger(args, "k_shift", k_shift);
        max_cluster = Args.parseInteger(args, "max_cluster", max_cluster);
        topKmer_trials = Args.parseInteger(args, "topKmer_trials", topKmer_trials);
        k_mask_f = Args.parseFloat(args, "k_mask_f", k_mask_f);
        kpp_mode = Args.parseInteger(args, "kpp_mode", kpp_mode);
        k_fold = Args.parseDouble(args, "k_fold", k_fold);
        gc = Args.parseDouble(args, "gc", gc);
        if (gc>0)
        	setGC(gc);
        wm_factor = Args.parseDouble(args, "wmf", wm_factor);
        fpr = Args.parseDouble(args, "fpr", fpr);
        motif_relax_factor = Args.parseDouble(args, "mrf", motif_relax_factor);
        ic_trim = Args.parseDouble(args, "ic", ic_trim);
        hgp = Args.parseDouble(args, "hgp", hgp);
        kmer_hgp = Args.parseDouble(args, "kmer_hgp", kmer_hgp);
//        kmer_freq_pos_ratio = Args.parseDouble(args, "kmer_freq_pos_ratio", kmer_freq_pos_ratio);
//        kmer_cluster_seq_count = Args.parseInteger(args, "cluster_seq_count", kmer_cluster_seq_count);
        kpp_factor = Args.parseDouble(args, "kpp_factor", kpp_factor);
        kpp_nmotifs = Args.parseInteger(args, "kpp_nmotifs", kpp_nmotifs);
        noise = Args.parseDouble(args, "noise", noise);
        motif_hit_factor = Args.parseDouble(args, "pwm_hit_factor", motif_hit_factor);
        kmer_inRange_fraction = Args.parseDouble(args, "kmer_aligned_fraction", kmer_inRange_fraction);
        kmer_consistent_fraction = Args.parseDouble(args, "kmer_consistent_fraction", kmer_consistent_fraction);
//        kmer_set_overlap_ratio = Args.parseDouble(args, "kmer_set_overlap_ratio", kmer_set_overlap_ratio);
        repeat_fraction = Args.parseDouble(args, "repeat_fraction", repeat_fraction);
        seed_range = Args.parseInteger(args, "seed_range", seed_range);
        kmer_remove_mode = Args.parseInteger(args, "kmer_shift_remove", kmer_remove_mode);
               
        ip_ctrl_ratio = Args.parseDouble(args, "icr", ip_ctrl_ratio);
        maxThreads = Args.parseInteger(args,"t",java.lang.Runtime.getRuntime().availableProcessors());	// default to the # processors
        q_value_threshold = Args.parseDouble(args, "q", q_value_threshold);	// q-value
        q_refine = Args.parseDouble(args, "q2", q_refine);	// q-value for refine regions
        if (q_refine==-1)
        	q_refine = q_value_threshold*0.5;
        else{
        	if (q_refine>q_value_threshold){
        		System.err.println("q2>q");
        		throw new Exception("Invalide command line option: q2>q");
        	}
        }
        	
        sparseness = Args.parseDouble(args, "a", sparseness);	// minimum alpha parameter for sparse prior
        poisson_alpha = Args.parseDouble(args, "pa", poisson_alpha);	
        alpha_factor = Args.parseDouble(args, "af", alpha_factor); // denominator in calculating alpha value from read count in the region
        alpha_fine_factor = Args.parseDouble(args, "aff", alpha_fine_factor); // Divider in calculating alpha value when having a noise component
        noise_distribution = Args.parseInteger(args, "nd", noise_distribution);
        
        fold = Args.parseDouble(args, "fold", fold); // minimum fold enrichment IP/Control for filtering
        shapeDeviation =  TF_binding?-0.3:-0.2;		// set default according to filter type    		
        shapeDeviation = Args.parseDouble(args, "sd", shapeDeviation); // maximum shapeDeviation value for filtering
        max_hit_per_bp = Args.parseInteger(args, "mrc", 0); //max read count per bp, default -1, estimate from data
        window_size_factor = Args.parseInteger(args, "wsf", window_size_factor);
        second_lambda_region_width = Args.parseInteger(args, "w2", second_lambda_region_width);
        third_lambda_region_width = Args.parseInteger(args, "w3", third_lambda_region_width);
//        joint_event_distance = Args.parseInteger(args, "j", joint_event_distance);		// max distance of joint events
        top_events = Args.parseInteger(args, "top", top_events);
        min_event_count = Args.parseInteger(args, "min", min_event_count);
        top_fract_to_skip = Args.parseDouble(args, "skip_frac", top_fract_to_skip);
        base_reset_threshold = Args.parseInteger(args, "reset", base_reset_threshold);
        min_region_width = Args.parseInteger(args, "min_region_width", min_region_width);
        verbose = Args.parseInteger(args, "v", verbose);
        smooth_step = Args.parseInteger(args, "smooth", smooth_step);
        strand_type = Args.parseInteger(args, "strand_type", strand_type);
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
    public void setGC(double gc){
    	this.gc = gc;
    	bg[0]=0.5-gc/2; 
    	bg[1]=gc/2; 
    	bg[2]=bg[1]; 
    	bg[3]=bg[0];
    }
}

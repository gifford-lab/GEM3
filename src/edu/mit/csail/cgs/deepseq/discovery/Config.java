package edu.mit.csail.cgs.deepseq.discovery;

import java.util.Set;

import edu.mit.csail.cgs.tools.utils.Args;

public class Config {
	public boolean trim_simple=false;
	public boolean print_PI = false;
	public boolean do_model_selection=false;
	public boolean classify_events = false;
    public boolean use_joint_event = false;
    public boolean TF_binding = true;
    public boolean outputBED = false;
    public boolean write_RSC_file = false;
    public boolean kmer_print_hits = false;
    public boolean testPValues = false;
    public boolean post_artifact_filter=false;
    public boolean filterEvents=true;
    public boolean filterDupReads=true;
    public boolean kl_count_adjusted = false;
    public boolean sort_by_location=false;
    public boolean exclude_unenriched = true;
    public boolean dump_regression = false;
    public boolean use_event_strength = false;
    public boolean use_event_rank = false;
    public boolean use_kmer_strength = false;
    public boolean print_kmer_bPos = false;
  	public int KL_smooth_width = 0;
    public int max_hit_per_bp = -1;
    public int maxThreads;		// default to #CPU
    // k-mer related
    public int k = -1;			// the width of kmer
    public int k_min = -1;		// the minimum value of k
    public int k_max= -1;		// the maximum value of k        
    public String seed = null;
    public int k_seqs = 5000;	// the top number of event to get underlying sequences for initial Kmer learning 
    public int k_win = 61;		// the window around binding event to search for kmers
    public int k_win2 = 101;	// the window around binding event to search for motifs (in later rounds)
    public int k_win_f = 4;		// k_win = k_win_f * k
    public int k_neg_dist = 300;// the distance of the nearest edge of negative region from binding sites 
    public int k_negSeq_ratio = 2; 		// The ratio of cache negative sequences to positive sequences
    public int k_shift = 99;	// the max shift from seed kmer when aligning the kmers     
    public int max_cluster = 50;
//    public int k_overlap = 7;	// the number of overlapped bases to assemble kmers into PWM    
    public float k_mask_f = 1;	// the fraction of PWM to mask
    public int kpp_mode = 0;	// different mode to convert kmer count to positional prior alpha value
    public double hgp = -3; 	// p-value threshold of hyper-geometric test for enriched kmer 
    public double k_fold = 3;	// the minimum fold of kmer count in positive seqs vs negative seqs
    public double gc = -1;	// GC content in the genome			//0.41 for human, 0.42 for mouse
    public double[] bg= new double[4];	// background frequency based on GC content
    public double wm_factor = 0.6;		// The threshold relative to the maximum PWM score, for including a sequence into the cluster 
    public double ic_trim = 0.4;		// The information content threshold to trim the ends of PWM
    public double kmer_freq_pos_ratio = 0.8;	// The fraction of most frequent k-mer position in aligned sequences
    public double motif_hit_factor = 0.005;
    public double motif_hit_factor_report = 0.05;
    public double kmer_set_overlap_ratio = 0.5;
    public double repeat_fraction=1;		// ignore lower case letter and N in motif discovery if less than _fraction_ of sequence
    public int kmer_remove_mode = 0;
    public int seed_range = 3;
    public double kmer_aligned_fraction = 0.5;		// the fraction of kmer in the seed_range
    public boolean select_seed = false;
    public boolean use_grid_search = true;
    public boolean kmer_use_insig = false;
    public boolean kmer_use_filtered = false;
    public boolean use_weighted_kmer = true;		// strength weighted k-mer count
    public boolean use_pos_kmer = true;				// position weighted k-mer count
    public boolean k_neg_shuffle = false;
    public boolean k_neg_dinu_shuffle = false;		// di-nuleotide shuffle
   	public boolean re_align_kmer = false;
   	public boolean use_kmer_mismatch = true;
   	public boolean use_seed_family = true;		// start the k-mer alignment with seed family (kmers with 1 or 2 mismatch)
   	public boolean use_ksm = true;				// align with KSM (together with PWM)
 	public boolean estimate_ksm_threshold = true;
  	public boolean kpp_normalize_max = true;
  	public boolean pp_use_kmer = true;			// position prior using k-mer(true) or PWM(false)
  	public double kpp_factor = 0.8;
  	public double noise = 0.0;
    public boolean print_aligned_seqs = false;
    public boolean print_input_seqs = false;
    public boolean print_all_kmers = false;
    public boolean print_bound_seqs = false;
    public boolean re_train = false;
	public boolean refine_pwm = false;
    public boolean print_pwm_fdr = false;
    public boolean evaluate_by_kcm = false;		// whether to use K-mer class model to evaluate improvement of new cluster, default to use PWM
    public boolean use_strength_weight = true;	// use binding event strength to weight 
    public boolean use_pos_weight = true;		// use binding position profile to weight motif site
    public boolean allow_single_family =true;	// allow the kmer family only contains seed, i.e. no mismatch kmers
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
    public double joint_event_distance = 500;
    public double alpha_factor = 3.0;
    public double excluded_fraction = 0.05;	// top and bottom fraction of region read count to exclude for regression
    public int top_events = 2000;
    public int min_event_count = 500;	// minimum num of events to update read distribution
    public int smooth_step = 30;
    public int window_size_factor = 4;	//number of model width per window
    public int min_region_width = 50;	//minimum width for select enriched region
    public double mappable_genome_length = -1; // default is to compute
    public double sparseness=-1;
    public double poisson_alpha=1e-4; 				// the Poisson p-value for estimating alpha
    public double fold = 2.5;
    public double kl_ic = -2.0;
    public double shapeDeviation;
    public int gentle_elimination_factor = 2;	// factor to reduce alpha to a gentler pace after eliminating some component
    public int resolution_extend = 2;
    public int first_lambda_region_width  =  1000;
    public int second_lambda_region_width =  5000;
    public int third_lambda_region_width  = 10000;
    public boolean bic = false;				// use BIC or AIC for model selection
    public boolean use_dynamic_sparseness = true;
    public boolean use_betaEM = true;
    public boolean use_scanPeak  = true;
    public boolean refine_regions = false;		// refine the enrichedRegions for next round using EM results
    public boolean cache_genome = true;			// cache the genome sequence
    public String genome_path = null;
    
    public int verbose=1;		// BindingMixture verbose mode
    public int base_reset_threshold = 200;	// threshold to set a base read count to 1
    public int windowSize;			// size for EM sliding window for splitting long regions, set modelWidth * config.window_size_factor in KPPMixture
    //Run EM up until <tt>ML_ITER</tt> without using sparse prior
    public int ML_ITER=10;
    // the range to scan a peak if we know position from EM result
    public int SCAN_RANGE = 20;
    public int gentle_elimination_iterations = 5;
    public double minFoldChange = 1.5; // minimum fold change for a significant event.  applied after event discovery during the p-value filtering stage

    public void parseArgs(String args[]) throws Exception {
        Set<String> flags = Args.parseFlags(args);
        // default as false, need the flag to turn it on
        classify_events = flags.contains("classify");
        print_PI = flags.contains("print_PI");
        sort_by_location = flags.contains("sl");
        use_joint_event = flags.contains("refine_using_joint_event");
        post_artifact_filter = flags.contains("post_artifact_filter");
        kl_count_adjusted = flags.contains("adjust_kl");
        refine_regions = flags.contains("refine_regions");
        outputBED = flags.contains("outBED");
        write_RSC_file = flags.contains("writeRSC");
        testPValues = flags.contains("testP");
        if (testPValues)
        	System.err.println("testP is " + testPValues);
        dump_regression = flags.contains("dump_regression");
        use_event_strength = flags.contains("use_event_strength");
        use_event_rank = flags.contains("use_event_rank");
        use_kmer_strength = flags.contains("use_kmer_strength");
        kmer_print_hits = flags.contains("kmer_print_hits");
        select_seed = flags.contains("select_seed");
        kmer_use_insig = flags.contains("kmer_use_insig");
        kmer_use_filtered = flags.contains("kmer_use_filtered");

        k_neg_shuffle = flags.contains("k_neg_shuffle");
        k_neg_dinu_shuffle = flags.contains("k_neg_dinu_shuffle");
        re_align_kmer = flags.contains("rak");
        print_aligned_seqs = flags.contains("print_aligned_seqs");
        print_input_seqs = flags.contains("print_input_seqs");
        print_all_kmers = flags.contains("print_all_kmers");
        print_bound_seqs = flags.contains("print_bound_seqs");
        re_train = flags.contains("re_train");
        refine_pwm = flags.contains("refine_pwm");
        print_pwm_fdr = flags.contains("print_pwm_fdr");      
        use_db_genome = flags.contains("use_db_genome");
        evaluate_by_kcm = flags.contains("evaluate_by_kcm");
        k_mask_1base = flags.contains("k_mask_1base");
        bic = flags.contains("bic"); 					// BIC or AIC
        
        // default as true, need the opposite flag to turn it off
        exclude_unenriched = !flags.contains("not_ex_unenriched");
        use_dynamic_sparseness = ! flags.contains("fa"); // fix alpha parameter
        use_betaEM = ! flags.contains("poolEM");
        filterEvents = !flags.contains("nf");	// no filtering for predicted events
        filterDupReads = !flags.contains("nrf");	// no read filtering of duplicate reads
        TF_binding = ! flags.contains("br");	// broad region, not TF data, is histone or pol II
        if (!TF_binding){
            use_joint_event = true;
            sort_by_location = true;
        }
        use_scanPeak = ! flags.contains("no_scanPeak");
        do_model_selection = !flags.contains("no_model_selection");
        use_kmer_mismatch = !flags.contains("no_kmm");
        use_seed_family = !flags.contains("no_seed_family");
        use_ksm = !flags.contains("no_ksm");
        pp_use_kmer = !flags.contains("pp_pwm");
        estimate_ksm_threshold = !flags.contains("no_ksm_threshold");
        use_strength_weight = !flags.contains("no_weight");
        use_pos_weight = !flags.contains("no_pos_weight");
        use_weighted_kmer = !flags.contains("no_weighted_kmer");
        use_pos_kmer = !flags.contains("no_pos_kmer");
        use_grid_search = !flags.contains("no_grid_search");
        allow_single_family = !flags.contains("no_single_family");
        allow_seed_reset = !flags.contains("no_seed_reset");
        selectK_byTopKmer = !flags.contains("selectK_byPWM");	
        if (selectK_byTopKmer)														// overwrite allow_seed_reset
        	allow_seed_reset = false;
        allow_seed_inheritance = !flags.contains("no_seed_inheritance");
        pwm_align_new = !flags.contains("pwm_align_all");
        filter_pwm_seq = !flags.contains("pwm_seq_asIs");
        strigent_event_pvalue = !flags.contains("relax");

        mappable_genome_length = Args.parseDouble(args, "s", mappable_genome_length);	// size of mappable genome
       
        // Optional input parameter
        genome_path = Args.parseString(args, "genome", genome_path);
        k = Args.parseInteger(args, "k", k);
        k_min = Args.parseInteger(args, "k_min", k_min);
        k_max = Args.parseInteger(args, "k_max", k_max);
        k_min = Args.parseInteger(args, "kmin", k_min);				// fail-safe
        k_max = Args.parseInteger(args, "kmax", k_max);
        seed = Args.parseString(args, "seed", null);
        if (seed!=null){
        	k = seed.length();
        	k_min = -1;
        	k_max = -1;
        }
        k_seqs = Args.parseInteger(args, "k_seqs", k_seqs);
        k_win = Args.parseInteger(args, "k_win", k_win);
        k_win_f = Args.parseInteger(args, "k_win_f", k_win_f);
        k_neg_dist = Args.parseInteger(args, "k_neg_dist", k_neg_dist);
        k_negSeq_ratio = Args.parseInteger(args, "k_neg_ratio", k_negSeq_ratio);
        k_shift = Args.parseInteger(args, "k_shift", k_shift);
        max_cluster = Args.parseInteger(args, "max_cluster", max_cluster);
        k_mask_f = Args.parseFloat(args, "k_mask_f", k_mask_f);
        kpp_mode = Args.parseInteger(args, "kpp_mode", kpp_mode);
        k_fold = Args.parseDouble(args, "k_fold", k_fold);
        gc = Args.parseDouble(args, "gc", gc);
        if (gc>0)
        	setGC(gc);
        wm_factor = Args.parseDouble(args, "wmf", wm_factor);
        ic_trim = Args.parseDouble(args, "ic", ic_trim);
        hgp = Args.parseDouble(args, "hgp", hgp);
        kmer_freq_pos_ratio = Args.parseDouble(args, "kmer_freq_pos_ratio", kmer_freq_pos_ratio);
//        kmer_cluster_seq_count = Args.parseInteger(args, "cluster_seq_count", kmer_cluster_seq_count);
        kpp_factor = Args.parseDouble(args, "kpp_factor", kpp_factor);
        noise = Args.parseDouble(args, "noise", noise);
        motif_hit_factor = Args.parseDouble(args, "pwm_hit_factor", motif_hit_factor);
        kmer_aligned_fraction = Args.parseDouble(args, "kmer_aligned_fraction", kmer_aligned_fraction);
        kmer_set_overlap_ratio = Args.parseDouble(args, "kmer_set_overlap_ratio", kmer_set_overlap_ratio);
        repeat_fraction = Args.parseDouble(args, "repeat_fraction", repeat_fraction);
        seed_range = Args.parseInteger(args, "seed_range", seed_range);
        kmer_remove_mode = Args.parseInteger(args, "kmer_shift_remove", kmer_remove_mode);
               
        ip_ctrl_ratio = Args.parseDouble(args, "icr", ip_ctrl_ratio);
        maxThreads = Args.parseInteger(args,"t",java.lang.Runtime.getRuntime().availableProcessors());	// default to the # processors
        q_value_threshold = Args.parseDouble(args, "q", q_value_threshold);	// q-value
        q_refine = Args.parseDouble(args, "q2", q_refine);	// q-value for refine regions
        if (q_refine==-1)
        	q_refine = q_value_threshold*0.75;
        else{
        	if (q_refine>q_value_threshold){
        		System.err.println("q2>q");
        		throw new Exception("Invalide command line option: q2>q");
        	}
        }
        	
        sparseness = Args.parseDouble(args, "a", sparseness);	// minimum alpha parameter for sparse prior
        poisson_alpha = Args.parseDouble(args, "pa", poisson_alpha);	
        alpha_factor = Args.parseDouble(args, "af", alpha_factor); // denominator in calculating alpha value
        fold = Args.parseDouble(args, "fold", fold); // minimum fold enrichment IP/Control for filtering
        shapeDeviation =  TF_binding?-0.3:-0.2;		// set default according to filter type    		
        shapeDeviation = Args.parseDouble(args, "sd", shapeDeviation); // maximum shapeDeviation value for filtering
        max_hit_per_bp = Args.parseInteger(args, "mrc", 0); //max read count per bp, default -1, estimate from data
        window_size_factor = Args.parseInteger(args, "wsf", 3);
        second_lambda_region_width = Args.parseInteger(args, "w2", second_lambda_region_width);
        third_lambda_region_width = Args.parseInteger(args, "w3", third_lambda_region_width);
        joint_event_distance = Args.parseInteger(args, "j", 10000);		// max distance of joint events
        top_events = Args.parseInteger(args, "top", top_events);
        min_event_count = Args.parseInteger(args, "min", min_event_count);
        base_reset_threshold = Args.parseInteger(args, "reset", base_reset_threshold);
        min_region_width = Args.parseInteger(args, "min_region_width", 50);
        verbose = Args.parseInteger(args, "v", verbose);
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
    public void setGC(double gc){
    	this.gc = gc;
    	bg[0]=0.5-gc/2; 
    	bg[1]=gc/2; 
    	bg[2]=bg[1]; 
    	bg[3]=bg[0];
    }
}

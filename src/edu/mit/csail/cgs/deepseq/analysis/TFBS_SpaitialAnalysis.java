package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;

import edu.mit.csail.cgs.clustering.affinitypropagation.CorrelationSimilarity;
import edu.mit.csail.cgs.datasets.chipseq.ChipSeqLocator;
import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedPoint;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrixImport;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.DeepSeqExpt;
import edu.mit.csail.cgs.deepseq.StrandedBase;
import edu.mit.csail.cgs.deepseq.features.ComponentFeature;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils.NarrowPeak;
import edu.mit.csail.cgs.deepseq.utilities.ReadCache;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSPeak;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.utils.stats.StatUtil;
/**
 * Compute the spatial relationship among multiple TFs' binding sites<br>
 * Find the clusters of multiTF binding and the relative positions of each TF sites
 * @author Yuchun
 *
 */
public class TFBS_SpaitialAnalysis {
	private int signal_radius = 500;
	String[] args = null;
	Set<String> flags =null;
	Genome genome=null;
	ArrayList<String> expts = new ArrayList<String>();
	ArrayList<String> tf_names = new ArrayList<String>();
	ArrayList<String> motif_names = new ArrayList<String>();
	ArrayList<String> readdb_names = new ArrayList<String>();
	ArrayList<String> indirect_tf_expts = new ArrayList<String>();
	/** Mapping of TFSS(sequence specific) GEM expt ids (from direct binding id to indirect binding id */
	HashMap<Integer, Integer> directid2indirectid = null;
	
	ArrayList<WeightMatrix> pwms = new ArrayList<WeightMatrix>();
	ArrayList<String> kmers = new ArrayList<String>();
	ArrayList<ArrayList<Region>> annoRegions = new ArrayList<ArrayList<Region>>();		// annotation regions that do not use to merge regions, but count overlaps, such as histone mark
	ArrayList<String> annoLabels = new ArrayList<String>();		// annotation labels, corresponding to the regions
	ArrayList<ArrayList<Site>> all_sites = new ArrayList<ArrayList<Site>>();
	ArrayList<ArrayList<Site>> all_indirect_sites = new ArrayList<ArrayList<Site>>();
	ArrayList<Point> all_TSS;
	ArrayList<Cluster> all_clusters;
	ArrayList<Region> rs = null; 		// the user supplied regions
	int[][]tss_signals = null;
	
	double gc = 0.42;//mouse		gc=0.41 for human
	int profile_range = 100;		// the range (-x, +x) around 0 spacing position, the cluster distance parameter must be >= this range to have correct result
	int distance = 50;		// distance between TFBS within a cluster
	int top = -1;			// use only number of top ranking events for analysis ( default -1, use all events), this applies to all factors
	double q = 2;			// the q-value (-log10) cutoff for top ranking events for analysis ( default 2, use all events), this applies to all factors
	int anno_expand_distance = 500;		// the distance to expand a peak cluster when overlapping with annotations
	int range = 1000;		// the range around anchor site to search for targets
	int cluster_motif_padding = 100;  // the padding distance added to the cluster range for motif searching
	int exclude_range = 500;		// the range around anchor site to exclude same site
	int min_site = 1;				// the minimum number of sites for the cluster to be printed out
	int width = 3;
	int height = 3;
	double wm_factor = 0.6;	// PWM threshold, as fraction of max score
	double cutoff = 0.3;	// corr score cutoff
	double spacing_cutoff = 4;	// cutoff for spacing count Poisson p-value (corrected) 
	File dir;
	boolean no_gem_pwm = false;
	boolean oldFormat =  false;
	boolean useDirectBindingOnly = false;
	boolean print_uci_matlab_format = false;
	boolean print_hdp_format = false;
	boolean print_matrix = false;
	boolean print_hdp_rc_format = false;
	boolean print_matrix_rc = false;
	boolean print_full_format = false;
	boolean print_TMT_format = false;	
	boolean print_spacing_profile = false;
	
	private SequenceGenerator<Region> seqgen;
	boolean dev = false;
	boolean zero_or_one = false;	// for each TF, zero or one site per cluster, no multiple sites
	boolean out_subset = false;		// output the call coords in selective subset of regions
	String outPrefix = "out";
	String tfss_file;				// Sequence-specific TFs
	String tss_file;
	String piq_file;				// for peak files other than GEM calls (in PIQ format)
	String point_file;				// for peak/motif files (in CGS format)
	String pwm_file;
	String kmer_file;
	String anno_region_file;
	String query_region_file;
	String tss_signal_file;
	String cluster_file;
	String cluster_key_file;
	String exclude_sites_file;
	String anchor_string;	// the id of TF to anchor the sites/regions/sequences
	String target_string;	// the id of target TF, and its spacing range
	String sort_string;		// id:left:right, the id of TF to sort the sites/regions/sequences according to its offset from anchor
	
	// command line option:  (the folder contains GEM result folders) 
	// --dir C:\Data\workspace\gse\TFBS_clusters --species "Mus musculus;mm9" --r 2 --pwm_factor 0.6 --info expt.info.txt [--no_cache --old_format] 
	public static void main(String[] args) {
		TFBS_SpaitialAnalysis analysis = new TFBS_SpaitialAnalysis(args);
		int round = Args.parseInteger(args, "r", 2); 	//GEM output round (1 for GPS, 2 for GEM)
		int type = Args.parseInteger(args, "type", 999);
		ArrayList<ArrayList<Site>> clusters=null;
		switch(type){
		case 999:	// default: simplified file loading for RPD public code
			analysis.loadBindingEventsSimple();
			clusters = analysis.mergeTfbsClusters();
			analysis.outputTFBSclusters(clusters);
			break;
		case 0:
			analysis.loadBindingEvents_old();
			clusters = analysis.mergeTfbsClusters();
			analysis.outputTFBSclusters(clusters);
			break;
		case 1:
			analysis.loadClusterAndTssSignals();
			analysis.computeCorrelations();
			break;
		case 2:
			analysis.loadBindingEvents_old();
			analysis.computeTfbsSpacingDistribution();
			break;
		case 3:		// to print all the binding sites and motif positions in the clusters for downstream spacing/grammar analysis
			// 1. GEM: java edu.mit.csail.cgs.deepseq.analysis.TFBS_SpaitialAnalysis --species "Mus musculus;mm10"  --type 3 --dir /cluster/yuchun/www/guo/mES  --info mES.info.txt  --r $round --pwm_factor 0.6 --distance ${distance} --min_site ${min} --out $analysis
			// 2. NarrowPeak: --g /Users/yguo/data/projects/shared_data/GEM_support/hg19.info --type 3 --tf_peak_file /Users/yguo/proj/rpd/161118_HDP_ENCODE_TF_K562/k562.peak.test.txt --format bed --out k562.test
			System.out.println("Type 3: output binding and/or motif sites ...");
			analysis.loadEventsAndMotifs(round);
			clusters = analysis.mergeTfbsClusters();
			analysis.outputBindingAndMotifSites(clusters);
			break;
		case 31:		// to print all the binding sites and motif positions in the SPECIFIED regions for co-binding analysis
			// java edu.mit.csail.cgs.deepseq.analysis.TFBS_SpaitialAnalysis --species "Mus musculus;mm10"  --type 31 --dir /cluster/yuchun/www/guo/mES  --info mES.info.txt  --r $round --pwm_factor 0.6 --distance ${distance} --min_site ${min} --out $analysis
			analysis.loadEventsAndMotifs(round);
			clusters = analysis.addTfbs2Clusters();
			analysis.outputTFBSclusters(clusters);
			break;
		case 32:		// to print all the binding sites and motif positions in the SPECIFIED regions for spacing/grammar analysis
			// java edu.mit.csail.cgs.deepseq.analysis.TFBS_SpaitialAnalysis --species "Mus musculus;mm10"  --type 32 --dir /cluster/yuchun/www/guo/mES  --info mES.info.txt  --r $round --pwm_factor 0.6 --distance ${distance} --min_site ${min} --out $analysis
			analysis.loadEventsAndMotifs(round);
			clusters = analysis.addTfbs2Clusters();
			analysis.outputBindingAndMotifSites(clusters);
			break;
		case 33:		// to print all the binding sites and motif positions in the SPECIFIED regions for co-binding/spacing/grammar analysis
			// java edu.mit.csail.cgs.deepseq.analysis.TFBS_SpaitialAnalysis --species "Mus musculus;mm10"  --type 33 --dir /cluster/yuchun/www/guo/mES  --info mES.info.txt  --r $round --pwm_factor 0.6 --distance ${distance} --min_site ${min} --out $analysis
			analysis.loadEventsAndMotifs(round);
			clusters = analysis.addTfbs2Clusters();
			analysis.outputTFBSclusters(clusters);
			analysis.outputBindingAndMotifSites(clusters);
			break;
		case 4:		// spacing histogram
			// java edu.mit.csail.cgs.deepseq.analysis.TFBS_SpaitialAnalysis --species "Mus musculus;mm10" --type 4 --out $analysis --cluster 0_BS_Motif_clusters.${analysis}.d${distance}.min${min}.txt --key Key_flie --profile_range 200
			analysis.printSpacingHistrograms();
			break;
		case 5:		// anchor/sorted sites for matlab plot and sequences
			// java edu.mit.csail.cgs.deepseq.analysis.TFBS_SpaitialAnalysis --g C:\Data\workspace\GEM_support\mm10.info  --type 5 --out mES --cluster "C:\Data\workspace\gse\mES_spacing\0_BS_Motif_clusters.mES.d200.min1.txt" --anchor B16:-90:10 --target B04:-26:-26 --sort B03 --height 3
			analysis.plotAlignedSites();
			break;
		case 11:	// old code
			analysis.loadClusterAndTSS();
			analysis.computeCorrelations_db();
			break;		
		case -1:
			analysis.mergedTSS();
			break;
		case -2:
			analysis.printTssSignal();
			break;
		}
		System.out.println("\nDone!");
	}
	
	public TFBS_SpaitialAnalysis(String[] args){
		genome = CommonUtils.parseGenome(args);

		this.args = args;
		this.flags = Args.parseFlags(args);
		outPrefix = Args.parseString(args, "out", outPrefix);
		zero_or_one = flags.contains("zoo");
		dev = flags.contains("dev");
		out_subset = flags.contains("out_subset");
		oldFormat = flags.contains("old_format");
		no_gem_pwm = !flags.contains("gem_pwm");
//		indirect_binding = flags.contains("indirect_binding");
		useDirectBindingOnly = flags.contains("direct");
		print_uci_matlab_format = flags.contains("uci_matlab");
		print_hdp_format = flags.contains("hdp");
		print_matrix = flags.contains("matrix");
		print_hdp_rc_format = flags.contains("hdprc");
		print_matrix_rc = flags.contains("matrixrc");
		print_full_format = flags.contains("full");
		print_TMT_format = flags.contains("TMT");
		print_spacing_profile = flags.contains("spacing");
		dir = new File(Args.parseString(args, "dir", "."));
		expts = new ArrayList<String>();
		tf_names = new ArrayList<String>();
		motif_names = new ArrayList<String>();
		String tf_info_file = Args.parseString(args, "gem", null);
		int type = Args.parseInteger(args, "type", 999);
		
		if (type==999){		// default for public version (RMD paper)
			tf_info_file = Args.parseString(args, "tf_peak_file", null);
			print_hdp_format = true;
			print_matrix = true;
			print_full_format = true;
		}
		if (type==3){		// default for public version 
			tf_info_file = Args.parseString(args, "tf_peak_file", null);
		}
		if (tf_info_file!=null){
			ArrayList<String> info = CommonUtils.readTextFile(tf_info_file);
			// expt name | short name | factor type | display name | motif | readdb name
			for (String txt: info){
				if ( ! (txt.equals("") || txt.startsWith("#")) ){
					String[] f = txt.trim().split("\t");
					if (f.length == 6){
						expts.add(f[0]);
						tf_names.add(f[3]);
						motif_names.add(f[4]);
						readdb_names.add(f[5]);
					}
					else if (f.length == 2){		// 2 columns, TF.name<TAB>GEM.path
						tf_names.add(f[0]);
						expts.add(f[1]);
						motif_names.add("N.A.");
					}
				}
			}
		}
		System.out.println("Loading datasets ...");
		top = Args.parseInteger(args, "top", top);
		q = Args.parseDouble(args, "q", q);
		anchor_string = Args.parseString(args, "anchor", anchor_string);	// the id of TF/PWM/Kmer to anchor the sites/regions/sequences
		target_string = Args.parseString(args, "target", target_string);	// the id of TF/PWM/Kmer to anchor the sites/regions/sequences
		sort_string = Args.parseString(args, "sort", sort_string);			// the id of TF/PWM/Kmer to sort the sites/regions/sequences

		width = Args.parseInteger(args, "width", width);
		height = Args.parseInteger(args, "height", height);
		profile_range = Args.parseInteger(args, "profile_range", profile_range);
		
		tfss_file = Args.parseString(args, "tfss_file", null);
		if (tfss_file!=null){
			directid2indirectid = new HashMap<Integer, Integer>();
			indirect_tf_expts = CommonUtils.readTextFile(tfss_file);
			for (int i=0;i<indirect_tf_expts.size();i++){
				String e = indirect_tf_expts.get(i);
				if (type==999){	//public version
					int idx = tf_names.indexOf(e);
					if (idx>=0){		// add fake expt and iNames for the indirect TFSS sites
						tf_names.add("i_"+tf_names.get(idx));
						expts.add(expts.get(idx));
						directid2indirectid.put(idx, tf_names.size()-1);
					}
				}
				else{
					int idx = expts.indexOf(e);
					if (idx>=0){		// add fake expt and iNames for the indirect TFSS sites
						expts.add(e);
						tf_names.add("i_"+tf_names.get(idx));
						readdb_names.add(readdb_names.get(idx));
						directid2indirectid.put(idx, expts.size()-1);
					}
				}
			}
		}
		
		tss_file = Args.parseString(args, "tss", null);
		signal_radius = Args.parseInteger(args, "signal_radius", signal_radius);
		tss_signal_file = Args.parseString(args, "tss_signal", null);
		cluster_file = Args.parseString(args, "cluster", null);
		cluster_key_file = Args.parseString(args, "key", null);
		piq_file = Args.parseString(args, "piq", null);				// additional peaks, e.g. PIQ calls
		point_file = Args.parseString(args, "point", null);				// additional peaks, e.g. PIQ calls
		pwm_file = Args.parseString(args, "pwms", null);				// additional pwms
		kmer_file = Args.parseString(args, "kmers", null);
		anno_region_file = Args.parseString(args, "anno_regions", null);
		query_region_file = Args.parseString(args, "regions", null);
		exclude_sites_file = Args.parseString(args, "ex", null);
		distance = Args.parseInteger(args, "distance", distance);
		anno_expand_distance = Args.parseInteger(args, "anno_expand_distance", anno_expand_distance);
		range = Args.parseInteger(args, "range", range);
		cluster_motif_padding = Args.parseInteger(args, "cluster_motif_padding", cluster_motif_padding);
		exclude_range = Args.parseInteger(args, "exclude", exclude_range);
		min_site = Args.parseInteger(args, "min_site", min_site);
		wm_factor = Args.parseDouble(args, "pwm_factor", wm_factor);
		cutoff = Args.parseDouble(args, "cutoff", cutoff);
		spacing_cutoff = Args.parseDouble(args, "spacing_cutoff", spacing_cutoff);
		gc = Args.parseDouble(args, "gc", gc);
		String genome_dir = Args.parseString(args, "genome", null);
		seqgen = new SequenceGenerator<Region>();
		if (genome_dir!=null)
			seqgen.setGenomePath(genome_dir);
		seqgen.useCache(!flags.contains("no_cache"));
	}

	/**
	 * This method loads events/region for topic model analysis.<br>
	 */
	private void loadBindingEventsSimple(){
		boolean isBED = Args.parseString(args, "format", "GEM").equalsIgnoreCase("BED");
		ArrayList<Region> ex_regions = new ArrayList<Region>();
		if(exclude_sites_file!=null){
			ex_regions = CommonUtils.loadCgsRegionFile(exclude_sites_file, genome);
		}
		for (int tf=0;tf<tf_names.size();tf++){
			if (tf_names.get(tf).startsWith("i_"))		// names start with i_ are artificially created id for TFSS indirect binding
				continue;
			boolean tfss = false;		// this is false for non-tfss factors (e.g. Pol2), or when not considering direct/indirect
			if (directid2indirectid!=null && directid2indirectid.containsKey(tf))
				tfss = true;
			String expt = expts.get(tf);

			System.out.print(String.format("TF#%d: loading %s", tf, expt));
			
			if (isBED){
				ArrayList<NarrowPeak> ps = CommonUtils.load_narrowPeak(genome, expt, true);
				ArrayList<Site> sites = new ArrayList<Site>();
			eachpeak:	for (int i=0;i<ps.size();i++){
					NarrowPeak p = ps.get(i);
					Site site = new Site();
					site.tf_id = tf;
					site.event_id = i;
					site.signal = p.signal;
					site.motifStrand = '*';
					site.bs = p.summit;
					
					// skip site in the ex_regions
					for (Region rr: ex_regions){
						if (rr.contains(site.bs))
							continue eachpeak;
					}
					sites.add(site);
				}
				System.out.println(",\t n="+sites.size());
				Collections.sort(sites);
				all_sites.add(sites);
				continue;
			}
			
			// GEM format
			try{
				List<GPSPeak> gpsPeaks = GPSParser.parseGPSOutput(expt, genome);
				ArrayList<Site> sites = new ArrayList<Site>();
				ArrayList<Site> indirectSites = new ArrayList<Site>();
			eachpeak:	for (int i=0;i<gpsPeaks.size();i++){
					GPSPeak p = gpsPeaks.get(i);
					Site site = new Site();
					site.tf_id = tf;
					site.event_id = i;
					site.signal = p.getStrength();
					site.motifStrand = p.getKmerStrand();
					site.bs = (Point)p;
					
					// skip site in the ex_regions
					for (Region r: ex_regions){
						if (r.contains(site.bs))
							continue eachpeak;
					}
					if (tfss && site.motifStrand=='*'){
						// for factor that are TFSS, thus considering indirect
						site.tf_id = directid2indirectid.get(site.tf_id);
						indirectSites.add(site);
					}
					else
						sites.add(site);
				}
				if (tfss)	
					System.err.println(",\t d="+sites.size()+",\t i="+indirectSites.size());
				else
					System.err.println(",\t n="+sites.size());
				
				Collections.sort(sites);
				all_sites.add(sites);
				if (tfss_file!=null && !indirectSites.isEmpty()){
					Collections.sort(indirectSites);
					all_indirect_sites.add(indirectSites);
				}
			}
			catch (IOException e){
				System.out.println(expt+" does not have valid GPS/GEM event call file.");
				System.exit(1);
			}
		}
		all_sites.addAll(all_indirect_sites);

	}
	/**
	 * This method loads events/region for topic model analysis.<br>
	 */
	private void loadBindingEvents_old(){
		ArrayList<Region> ex_regions = new ArrayList<Region>();
		if(exclude_sites_file!=null){
			ex_regions = CommonUtils.loadCgsRegionFile(exclude_sites_file, genome);
		}
		for (int tf=0;tf<tf_names.size();tf++){
			if (tf_names.get(tf).startsWith("i_"))		// names start with i_ are artificially created id for TFSS indirect binding
				continue;
			boolean tfss = false;		// this is false for non-tfss factors (e.g. Pol2), or when not considering direct/indirect
			if (directid2indirectid!=null && directid2indirectid.containsKey(tf))
				tfss = true;
			String expt = expts.get(tf);

			System.err.print(String.format("TF#%d: loading %s", tf, expt));
			
			File dir2= new File(dir, expt);
			if (!oldFormat)
				dir2= new File(dir2, expt+"_outputs");
			
			// load binding event files, TFSS use GEM calls, nonTFSS (or don't care) use GPS calls
			File gpsFile = new File(dir2, expt+"_"+ (tfss?2:1) +
					(oldFormat?"_GPS_significant.txt":"_GEM_events.txt"));
			String filePath = gpsFile.getAbsolutePath();

			try{
				List<GPSPeak> gpsPeaks = GPSParser.parseGPSOutput(filePath, genome);
				ArrayList<Site> sites = new ArrayList<Site>();
				ArrayList<Site> indirectSites = new ArrayList<Site>();
			eachpeak:	for (int i=0;i<gpsPeaks.size();i++){
					GPSPeak p = gpsPeaks.get(i);
					Site site = new Site();
					site.tf_id = tf;
					site.event_id = i;
					site.signal = p.getStrength();
					site.motifStrand = p.getKmerStrand();
					site.bs = (Point)p;
					
					// skip site in the ex_regions
					for (Region r: ex_regions){
						if (r.contains(site.bs))
							continue eachpeak;
					}
					if (tfss && site.motifStrand=='*'){
						// for factor that are TFSS, thus considering indirect
						site.tf_id = directid2indirectid.get(site.tf_id);
						indirectSites.add(site);
					}
					else
						sites.add(site);
				}
				if (tfss)	
					System.err.println(",\t d="+sites.size()+",\t i="+indirectSites.size());
				else
					System.err.println(",\t n="+sites.size());
				
				Collections.sort(sites);
				all_sites.add(sites);
				if (tfss_file!=null && !indirectSites.isEmpty()){
					Collections.sort(indirectSites);
					all_indirect_sites.add(indirectSites);
				}
			}
			catch (IOException e){
				System.out.println(expt+" does not have a valid GPS/GEM event call file.");
				e.printStackTrace(System.err);
				System.exit(1);
			}
		}
		all_sites.addAll(all_indirect_sites);

	}
	
	/**
	 * This method loads events/motifs for spacing analysis.<br>
	 * it loads the event positions, also loads GEM motifs by option, <br>
	 * loads PIQ calls/motifs by --piq, <br>
	 * it also loads additional motifs specified by --pwms, or --kmers<br>
	 * @param round GEM output round (1 for GPS, 2 for GEM)
	 */
	private void loadEventsAndMotifs(int round){
		boolean isBED = Args.parseString(args, "format", "GEM").equalsIgnoreCase("BED");
		ArrayList<Region> exclude_regions = new ArrayList<Region>();
		if(exclude_sites_file!=null){
			exclude_regions = CommonUtils.loadCgsRegionFile(exclude_sites_file, genome);
		}
		ArrayList<Region> queryRegions = new ArrayList<Region>();
		if(query_region_file!=null){
			queryRegions = CommonUtils.load_BED_regions(genome, query_region_file).car();
		}

		for (int tf=0;tf<tf_names.size();tf++){
			String expt = expts.get(tf);
			File dir2= new File(dir, expt);
			dir2= new File(dir2, expt+"_outputs");
			System.err.print(String.format("GEM #%d: loading %s", tf, expt));
			
			if (isBED){
				ArrayList<NarrowPeak> ps = CommonUtils.load_narrowPeak(genome, expt, true);
				ArrayList<Site> sites = new ArrayList<Site>();
			eachpeak:	for (int i=0;i<ps.size();i++){
					NarrowPeak p = ps.get(i);
					// skip site in the ex_regions
					for (Region rr: exclude_regions){
						if (rr.contains(p.summit))
							continue eachpeak;
					}
					Site site = new Site();
					site.tf_id = tf;
					site.event_id = i;
					site.signal = p.signal;
					site.motifStrand = '*';
					site.bs = p.summit;
					sites.add(site);
				}
				System.out.println(",\t n="+sites.size());
				Collections.sort(sites);
				all_sites.add(sites);
				continue;
			}
			
			// GEM format
			// load binding event files 
			File gpsFile = new File(dir2, expt+"_"+ (round>=2?round:1) + "_GEM_events.txt");
			String filePath = gpsFile.getAbsolutePath();
			try{
				List<GPSPeak> gpsPeaks = GPSParser.parseGPSOutput(filePath, genome);
				ArrayList<Site> sites = new ArrayList<Site>();
				eachpeak:	for (int i=0;i<gpsPeaks.size();i++){
					if (top!=-1 && sites.size()>=top){		// only use top ranking events for analysis
						break eachpeak;
					}
					GPSPeak p = gpsPeaks.get(i);
					if (q!=2 && p.getQV_lg10()<q){		// only use events with high q value
						continue;
					}					
					
					// skip site in the ex_regions and those not in queryRegion
					// TODO: this loop can be more efficient
					if(query_region_file!=null){
						boolean isPeakInRegion = false;
						for (Region r:queryRegions){
							if (r.contains(p)){
								isPeakInRegion = true;
								break;
							}
						}
						if (!isPeakInRegion)
							continue eachpeak;
					}
					for (Region r: exclude_regions){
						if (r.contains(p))
							continue eachpeak;
					}
						
					Site site = new Site();
					site.tf_id = tf;
					site.event_id = i;
					site.signal = p.getStrength();
					site.motifStrand = p.getKmerStrand();
					site.bs = (Point)p;					
					sites.add(site);
				}				
					
				System.err.println(", n="+sites.size());
				Collections.sort(sites);
				all_sites.add(sites);
			}
			catch (IOException e){
				System.out.println(expt+" does not have valid GPS/GEM event call file.");
				System.exit(1);
			}
			// load motif files
			if (no_gem_pwm)
				continue;
			WeightMatrix wm = null;
			final String suffix = expt+"_"+ (round>=2?round:1) +"_PFM";
			File[] files = dir2.listFiles(new FilenameFilter(){
				public boolean accept(File arg0, String arg1) {
					if (arg1.startsWith(suffix))
						return true;
					else
						return false;
				}
			});
			if (files.length==0){
				System.out.println(expt+" does not have a motif PFM file.");
				pwms.add(null);
			}
			else{				// if we have valid PFM file
				wm = CommonUtils.loadPWM_PFM_file(files[0].getAbsolutePath(), gc);
				pwms.add( wm );
			}
		} // for each expt
		
		// load PIQ calls
		if (piq_file!=null){
			ArrayList<String> lines = CommonUtils.readTextFile(piq_file); 
			int newTFID = expts.size();
			for (int l=0;l<lines.size();l++){
				if (lines.get(l).startsWith("#"))
					continue;
				String[] f = lines.get(l).split("\t");
				expts.add(f[0]);
				tf_names.add(f[1]);
				motif_names.add(f[2]);

				System.err.print(String.format("Peak#%d: loading %s", l, f[0]));
				ArrayList<Site> sites = new ArrayList<Site>();
				ArrayList<String> piqs = CommonUtils.readTextFile(f[3]);
				eachpiq:	for (int i=0;i<piqs.size();i++){
					String line = piqs.get(i);
					if (!line.contains("score")){
						String[] q = line.split(",");
						Point p = new Point(genome, q[1].replaceAll("\"", "").replaceAll("chr", ""), Integer.parseInt(q[2]));
					
						// skip site in the ex_regions and those not in queryRegion
						// TODO: this loop can be more efficient
						if(query_region_file!=null){
							boolean isPeakInRegion = false;
							for (Region r:queryRegions){
								if (r.contains(p)){
									isPeakInRegion = true;
									break;
								}
							}
							if (!isPeakInRegion)
								continue eachpiq;
						}
						for (Region r: exclude_regions){
							if (r.contains(p))
								continue eachpiq;
						}

						Site site = new Site();
						site.tf_id = newTFID+l;
						site.event_id = i;
						site.signal = Double.parseDouble(q[5]);
						site.motifStrand = '+';
						if (q.length>7){
							site.shapeOnlyScore = Double.parseDouble(q[7]);
							site.readScore = Double.parseDouble(q[8]);
						}
						site.bs = p;					
						for (Region r: exclude_regions){
							if (r.contains(site.bs))
								continue eachpiq;
						}
						sites.add(site);					
					}
				}		
				piqs = CommonUtils.readTextFile(f[4]);			// load PIQ calls on minus strand (RC.call)
				eachpiq2:	for (int i=0;i<piqs.size();i++){
					String line = piqs.get(i);
					if (!line.contains("score")){
						String[] q = line.split(",");
						Point p = new Point(genome, q[1].replaceAll("\"", "").replaceAll("chr", ""), Integer.parseInt(q[2]));
						
						// skip site in the ex_regions and those not in queryRegion
						// TODO: this loop can be more efficient
						if(query_region_file!=null){
							boolean isPeakInRegion = false;
							for (Region r:queryRegions){
								if (r.contains(p)){
									isPeakInRegion = true;
									break;
								}
							}
							if (!isPeakInRegion)
								continue eachpiq2;
						}
						for (Region r: exclude_regions){
							if (r.contains(p))
								continue eachpiq2;
						}
						Site site = new Site();
						site.tf_id = newTFID+l;
						site.event_id = -i;
						site.signal = Double.parseDouble(q[5]);
						site.motifStrand = '-';
						if (q.length>7){
							site.shapeOnlyScore = Double.parseDouble(q[7]);
							site.readScore = Double.parseDouble(q[8]);
						}
						site.bs = p;					
						for (Region r: exclude_regions){
							if (r.contains(site.bs))
								continue eachpiq2;
						}
						sites.add(site);					
					}
				}	// for
				System.err.println(", n="+sites.size());
				Collections.sort(sites);
				all_sites.add(sites);	
			}
		}		
		
		// load additional pwms
		if (pwm_file!=null){
			pwms.addAll( CommonUtils.loadPWMs_PFM_file(pwm_file, gc) );
		}
		
		// load kmers
		if (kmer_file!=null){
			kmers.addAll( CommonUtils.readTextFile(kmer_file));
		}
		
		System.out.println(String.format("In total, loaded %d GEM/other peak files, %d kmers and %d PWMs.", all_sites.size(), kmers.size(), pwms.size()));
		
		StringBuilder sb = new StringBuilder();
		int digits = (int) Math.ceil(Math.log10(expts.size()));
		for (int i=0;i<expts.size();i++){
			sb.append(String.format("B%0"+digits+"d\t%s\t%s\t%s\n", i, expts.get(i), tf_names.get(i), motif_names.get(i)));
		}
		if (!pwms.isEmpty()){
			digits = (int) Math.ceil(Math.log10(pwms.size()));
			for (int i=0;i<pwms.size();i++){
				sb.append(String.format("M%0"+digits+"d\t%s\t%s\t%s\n", i, pwms.get(i).getName(),pwms.get(i).getName(),pwms.get(i).getName()));
			}
		}
		if (!kmers.isEmpty()){
			digits = (int) Math.ceil(Math.log10(kmers.size()));
			for (int i=0;i<kmers.size();i++){
				sb.append(String.format("K%0"+digits+"d\t%s\t%s\t%s\n", i, kmers.get(i),kmers.get(i),kmers.get(i)));
			}
		}
		CommonUtils.writeFile("0_BS_Motif_clusters."+outPrefix+"_keys.txt", sb.toString());
	}
	/**
	 * Data structure for a binding site (with absolute genome coordinate and strand)
	 *
	 */
	class Site implements Comparable<Site>{
		int tf_id;
		int event_id;		// original event id in the list of this TF binding data
		Point bs;
		int id;				// global id of all BS
		double signal;
		char motifStrand;								// motif match strand
		double shapeOnlyScore;
		double readScore;
		public int compareTo(Site s) {					// ascending coordinate
			return(bs.compareTo(s.bs));
		}
	}

	class Cluster{
		Region region;
		ArrayList<Integer> TFIDs = new ArrayList<Integer>();
		ArrayList<Double> TF_Signals = new ArrayList<Double>();
		ArrayList<Boolean> TF_hasMotifs = new ArrayList<Boolean>();
	}
	
	/**
	 * Merge all the sites in all_sites (binding calls or nearest motif match positions)
	 * into clusters of positions. Organized by chromosomes.
	 * @return
	 */
	private ArrayList<ArrayList<Site>> mergeTfbsClusters(){
		// classify sites by chrom
		System.out.println("Merging binding/motif sites into non-overlaping regions (clusters).");
		TreeMap<String, ArrayList<Site>> chrom2sites = new TreeMap<String, ArrayList<Site>>();
		for (ArrayList<Site> sites:all_sites){
			for (Site s:sites){
				String chr = s.bs.getChrom();
				if (!chrom2sites.containsKey(chr))
					chrom2sites.put(chr, new ArrayList<Site>());
				chrom2sites.get(chr).add(s);
			}
		}
		
		// sort sites and form clusters
		ArrayList<ArrayList<Site>> clusters = new ArrayList<ArrayList<Site>>();		
		ArrayList<Site> cluster = new ArrayList<Site>();
		for (String chr: chrom2sites.keySet()){
			ArrayList<Site> sites = chrom2sites.get(chr);
			Collections.sort(sites);

			cluster.add(sites.get(0));
			for (int i=1;i<sites.size();i++){
				Site s = sites.get(i);
				Site p = cluster.get(cluster.size()-1);		//previous
				if (s.bs.getLocation()-p.bs.getLocation()<distance)
					cluster.add(s);
				else{
					cluster.trimToSize();
					clusters.add(cluster);
					cluster = new ArrayList<Site>();
					cluster.add(s);
				}					
			}
			// finish the chromosome
			cluster.trimToSize();
			clusters.add(cluster);
			cluster = new ArrayList<Site>();
		}
		// finish all the sites
		cluster.trimToSize();
		clusters.add(cluster);

		// output
		
		// remove if less than min_site cutoff
		ArrayList<ArrayList<Site>> newClusters = new ArrayList<ArrayList<Site>>();	
		for (ArrayList<Site> c :clusters){
			if (c.size()>=min_site)				
				newClusters.add(c);
		}
		clusters.clear();
		clusters = newClusters;
		
		return clusters;
	}
		
	/**
	 * Add all the sites/motifs in all_sites that are in the query regions into clusters of positions. 
	 * @return
	 */
	private ArrayList<ArrayList<Site>> addTfbs2Clusters(){
		rs = CommonUtils.load_BED_regions(genome, query_region_file).car();
		System.out.println("Assign binding/motif sites to "+rs.size()+" regions (clusters).");
		
		// classify sites by chrom, so that the sorting space is smaller
		TreeMap<String, ArrayList<Site>> chrom2sites = new TreeMap<String, ArrayList<Site>>();
		for (ArrayList<Site> sites:all_sites){
			for (Site s:sites){
				String chr = s.bs.getChrom();
				if (!chrom2sites.containsKey(chr))
					chrom2sites.put(chr, new ArrayList<Site>());
				chrom2sites.get(chr).add(s);
			}
		}
		for (String chr: chrom2sites.keySet()){
			ArrayList<Site> sites = chrom2sites.get(chr);
			sites.trimToSize();
			Collections.sort(sites);
		}
		
		// add sites to regions
		// NOTE: this is not very efficient for large data set, because it just loops over everything in a chromosome
		ArrayList<ArrayList<Site>> clusters = new ArrayList<ArrayList<Site>>();		
		for (Region r: rs){
			ArrayList<Site> sites = chrom2sites.get(r.getChrom());
			ArrayList<Site> cluster = new ArrayList<Site>();
			if (sites==null || sites.isEmpty()){
				clusters.add(cluster);		// add a empty cluster
				continue;
			}
			for (Site s: sites){
				if (r.contains(s.bs))
					cluster.add(s);			
			}
			cluster.trimToSize();
			clusters.add(cluster);
		}
		clusters.trimToSize();
		
		if (out_subset){
			// separate sites by their tf_id
			HashMap<Integer,ArrayList<Site>> subsets = new HashMap<Integer,ArrayList<Site>>();
			for (ArrayList<Site> c: clusters){
				for (Site s: c){
					if (!subsets.containsKey(s.tf_id))
						subsets.put(s.tf_id, new ArrayList<Site>());
					subsets.get(s.tf_id).add(s);
				}
			}
			for (int id: subsets.keySet()){
				StringBuilder sb = new StringBuilder();
				for (Site s: subsets.get(id)){
					sb.append(String.format("%s\t%.4f\t%.4f\t%.4f\n", new StrandedPoint(s.bs, s.motifStrand).toString(),
							s.signal, s.shapeOnlyScore, s.readScore));
					
				}
				CommonUtils.writeFile(outPrefix+"."+id+".txt", sb.toString());
			}
		}
		return clusters;
	}
		
	/** 
	 * Output all the binding sites in the clusters, for topic modeling analysis or clustering analysis<br>
	 * This method also include pseudo sites from overlapping annotation regions (such as histone mark broad peaks, TSSs, DHS, etc.)<br>
	 * The site clusters may contain empty clusters, neet to check the length of the sites.
	 * @param clusters
	 */
	private void outputTFBSclusters(ArrayList<ArrayList<Site>> clusters){
		
		// load annotation regions (for example, histone marks, genome segmentation states), 
		// these are not used for region merging, but used to count overlaps of the merged region with the anno regions
		if (anno_region_file!=null){
			ArrayList<String> lines = CommonUtils.readTextFile(anno_region_file);
			for (String s:lines){
				String f[]=s.split("\t");
				annoLabels.add(f[0]);
				System.out.print("Loading "+f[0]+" data ...");
				ArrayList<Region> rs = null;
				if (f[1].equalsIgnoreCase("BED"))
					rs = CommonUtils.load_BED_regions(genome, f[2]).car();					
				if (f[1].equalsIgnoreCase("CGS"))
					rs = CommonUtils.loadCgsRegionFile(f[2], genome);
				System.out.println("    "+rs.size()+" regions have been loaded.");
				Collections.sort(rs);		// need to be sorted to use Region.computeOverlapLength();
				annoRegions.add(rs);
			}
		}
		else{
			anno_expand_distance = 0;
		}
		// update clusters to include annotation info as a pseudo site
		int tf_count = expts.size();
		for (String s:annoLabels)
			tf_names.add(s);
		
		for (int j=0;j<clusters.size();j++){
			ArrayList<Site> c = clusters.get(j);
			int numSite = c.size();
//			if (numSite==0)
//				continue;
			Region r = null;
			Point lastSiteCoord = null;
			if (rs!=null){
				r = rs.get(j);
				lastSiteCoord = new Point(genome, r.getChrom(), r.getEnd());
			}
			else{
				r = new Region(genome, c.get(0).bs.getChrom(), c.get(0).bs.getLocation(), c.get(numSite-1).bs.getLocation()).expand(anno_expand_distance, anno_expand_distance);
				lastSiteCoord = c.get(numSite-1).bs;
			}
			for (int i=0;i<annoRegions.size();i++){
				ArrayList<Region> rs = annoRegions.get(i);
				int length = Region.computeOverlapLength(r, rs);
				// If not TSS2kb, overlap with the expanded distance
				// If TSS2kb,  overlap WITHOUT the expanded distance
				if ((length>0 && !annoLabels.get(i).equalsIgnoreCase("tss2kb")) || length>anno_expand_distance){		
					// add a pseudo site
					Site site = new Site();
					site.tf_id = tf_count+i;
					site.event_id = 0;
					site.signal = 0;
					site.motifStrand = '*';
					site.bs = lastSiteCoord;
					c.add(site);
				}
			}
		}
		
		if (print_full_format){
			StringBuilder sb = new StringBuilder();
			sb.append("#Region\tLength\t#Sites\tTFs\tTFIDs\tSignals\tPos\tMotifs\t#Motif\n");
			for (int j=0;j<clusters.size();j++){
				ArrayList<Site> c = clusters.get(j);
				int numSite = c.size();
				if (numSite==0){		// no sites, must be from addTfbs2Clusters(), add a empty line
					sb.append(rs.get(j).toString()).append("\t").append(rs.get(j).getWidth()).append("\t").
						append(numSite).append("\t").append("NA\tNA\tNA\tNA\tNA\t0\n");
					continue;
				}
				Region r = null;
				if (rs!=null)
					r = rs.get(j);
				else
					r = new Region(genome, c.get(0).bs.getChrom(), c.get(0).bs.getLocation(), c.get(numSite-1).bs.getLocation()).expand(anno_expand_distance, anno_expand_distance);
				StringBuilder sb_tfs = new StringBuilder();
				StringBuilder sb_tfids = new StringBuilder();
				StringBuilder sb_tf_signals = new StringBuilder();
				StringBuilder sb_tf_positions = new StringBuilder();
				StringBuilder sb_tf_motifs = new StringBuilder();
				int totalMotifs = 0;
				for (Site s:c){
					sb_tfs.append(tf_names.get(s.tf_id)).append(",");
					sb_tfids.append(s.tf_id).append(",");
					sb_tf_signals.append(String.format("%d", Math.round(s.signal))).append(",");
					sb_tf_positions.append(s.bs.getLocation()-r.getStart()).append(",");
					sb_tf_motifs.append(s.motifStrand).append(",");
					totalMotifs += s.motifStrand=='*'?0:1;
				}
				if (sb_tfs.length()!=0){
					sb_tfs.deleteCharAt(sb_tfs.length()-1);
					sb_tfids.deleteCharAt(sb_tfids.length()-1);
					sb_tf_signals.deleteCharAt(sb_tf_signals.length()-1);
					sb_tf_positions.deleteCharAt(sb_tf_positions.length()-1);
					sb_tf_motifs.deleteCharAt(sb_tf_motifs.length()-1);
				}
				sb.append(r.toString()).append("\t").append(r.getWidth()).append("\t").append(numSite).append("\t").
				append(sb_tfs.toString()).append("\t").append(sb_tfids.toString()).append("\t").
				append(sb_tf_signals.toString()).append("\t").append(sb_tf_positions.toString()).append("\t").
				append(sb_tf_motifs.toString()).append("\t").append(totalMotifs)
				.append("\n");
			}
	
			CommonUtils.writeFile("0_BS_clusters."+outPrefix+".d"+distance+".min"+min_site+".full.txt", sb.toString());
		}
		if (print_TMT_format){	// Stanford TMT format
			StringBuilder sb = new StringBuilder();
			for (int j=0;j<clusters.size();j++){
				ArrayList<Site> c = clusters.get(j);
				int numSite = c.size();
				if (numSite==0)
					continue;
				Region r = null;
				if (rs!=null)
					r = rs.get(j);
				else
					r = new Region(genome, c.get(0).bs.getChrom(), c.get(0).bs.getLocation(), c.get(numSite-1).bs.getLocation()).expand(anno_expand_distance, anno_expand_distance);
				StringBuilder sb_tfs = new StringBuilder().append(r.toString()).append("\t").append(numSite).append("\t");
				for (Site s:c){
					sb_tfs.append(tf_names.get(s.tf_id)).append(" ");
				}
				if (sb_tfs.length()!=0){
					sb_tfs.deleteCharAt(sb_tfs.length()-1);
				}
				sb.append(sb_tfs.toString()).append("\n");
			}
	
			CommonUtils.writeFile("0_BS_clusters."+outPrefix+".d"+distance+".min"+min_site+".TMT.txt", sb.toString());
		}
		if (print_uci_matlab_format){	// UCI Matlab Topic Modeling Toolbox 1.4 format
			StringBuilder sb = new StringBuilder();
			int[] factorSiteCount = new int[tf_names.size()];

			int docID = 1;
			for (ArrayList<Site> c :clusters){
				for (int s=0;s<c.size();s++){
					Site site = c.get(s);
					factorSiteCount[site.tf_id]++;
				}
				for (int f=0;f<factorSiteCount.length;f++){
					if (factorSiteCount[f]>0){
						sb.append(docID).append(" ").append(f+1).append(" ").append(factorSiteCount[f]).append("\n");
						factorSiteCount[f]=0;// reset to 0 for next cluster
					}
				}
				docID++;
			}
			CommonUtils.writeFile("0_BS_clusters."+outPrefix+".d"+distance+".min"+min_site+".UCI.txt", sb.toString());
		}
		
		if (print_hdp_format){	// Blei HDP (lda-c) format, SITE count for each TF
			StringBuilder sb = new StringBuilder();
			int[] factorSiteCount = new int[tf_names.size()];

			for (ArrayList<Site> c :clusters){
				for (int s=0;s<c.size();s++){
					Site site = c.get(s);
					factorSiteCount[site.tf_id]++;
				}
				int uniqueTermCount=0;
				for (int count:factorSiteCount){
					if (count!=0)
						uniqueTermCount++;
				}
				sb.append(uniqueTermCount);
				for (int f=0;f<factorSiteCount.length;f++){
					if (factorSiteCount[f]>0){
						sb.append(" ").append(f).append(":").append(factorSiteCount[f]);
						factorSiteCount[f]=0;// reset to 0 for next cluster
					}
				}
				sb.append("\n");
			}
			CommonUtils.writeFile("0_BS_clusters."+outPrefix+".d"+distance+".min"+min_site+".HDP.txt", sb.toString());
		}
		
		if (print_hdp_rc_format){	// Blei HDP (lda-c) format, read count for each TF
			StringBuilder sb = new StringBuilder();
			double[] factorReadCount = new double[tf_names.size()];

			for (ArrayList<Site> c :clusters){
				for (int s=0;s<c.size();s++){
					Site site = c.get(s);
					factorReadCount[site.tf_id]+=site.signal;
				}
				int uniqueTermCount=0;
				for (double count:factorReadCount){
					if (count!=0)
						uniqueTermCount++;
				}
				sb.append(uniqueTermCount);
				for (int f=0;f<factorReadCount.length;f++){
					if (factorReadCount[f]>0){
						sb.append(" ").append(f).append(":").append((int)factorReadCount[f]);
						factorReadCount[f]=0;// reset to 0 for next cluster
					}
				}
				sb.append("\n");
			}
			CommonUtils.writeFile("0_BS_clusters."+outPrefix+".d"+distance+".min"+min_site+".HDP_RC.txt", sb.toString());
		}
		
		if (print_hdp_format||print_uci_matlab_format){
			StringBuilder sb = new StringBuilder();
			for (int i=0;i<tf_names.size();i++)
				sb.append(tf_names.get(i)).append("\n");
			CommonUtils.writeFile("0_BS_clusters."+outPrefix+".Dictioinary.txt", sb.toString());
		}
		
		if (print_matrix){	// Print region-tfSiteCount matrix, can be used for clustering analysis
			StringBuilder sb = new StringBuilder();
			int[] factorSiteCount = new int[tf_names.size()];
			sb.append("#Region    ").append("\t");
			for (int i=0;i<tf_names.size();i++)
				sb.append(tf_names.get(i)).append("\t");
			CommonUtils.replaceEnd(sb, '\n');
			
			for (int j=0;j<clusters.size();j++){
				ArrayList<Site> c = clusters.get(j);
				int numSite = c.size();

				Region r = null;
				if (rs!=null)
					r = rs.get(j);
				else
					r = new Region(genome, c.get(0).bs.getChrom(), c.get(0).bs.getLocation(), c.get(numSite-1).bs.getLocation()).expand(anno_expand_distance, anno_expand_distance);
				sb.append(r.toString()).append("\t");
				for (int s=0;s<c.size();s++){
					Site site = c.get(s);
					factorSiteCount[site.tf_id]++;
				}
				for (int f=0;f<factorSiteCount.length;f++){
					sb.append(factorSiteCount[f]).append("\t");
					factorSiteCount[f]=0;// reset to 0 for next cluster
				}
				CommonUtils.replaceEnd(sb, '\n');
			}
			CommonUtils.writeFile("0_BS_clusters."+outPrefix+".d"+distance+".min"+min_site+".factorCount_matrix.txt", sb.toString());
		}
		
		if (print_matrix_rc){	// Print region-tfSiteCount matrix, can be used for clustering analysis
			StringBuilder sb = new StringBuilder();
			double[] factorReadCount = new double[tf_names.size()];
			sb.append("#Region    ").append("\t");
			for (int i=0;i<tf_names.size();i++)
				sb.append(tf_names.get(i)).append("\t");
			CommonUtils.replaceEnd(sb, '\n');
			
			for (int j=0;j<clusters.size();j++){
				ArrayList<Site> c = clusters.get(j);
				int numSite = c.size();
				Region r = null;
				if (rs!=null)
					r = rs.get(j);
				else
					r = new Region(genome, c.get(0).bs.getChrom(), c.get(0).bs.getLocation(), c.get(numSite-1).bs.getLocation()).expand(anno_expand_distance, anno_expand_distance);
				sb.append(r.toString()).append("\t");
				for (int s=0;s<c.size();s++){
					Site site = c.get(s);
					factorReadCount[site.tf_id]+=site.signal;
				}
				for (int f=0;f<factorReadCount.length;f++){
					sb.append((int)factorReadCount[f]).append("\t");
					factorReadCount[f]=0;// reset to 0 for next cluster
				}
				CommonUtils.replaceEnd(sb, '\n');
			}
			CommonUtils.writeFile("0_BS_clusters."+outPrefix+".d"+distance+".min"+min_site+".factorCount_matrix_RC.txt", sb.toString());
		}
	}
	
	/**
	 * Output all the binding sites, and all the motif matches
	 * @param clusters
	 */
	private void outputBindingAndMotifSites (ArrayList<ArrayList<Site>> clusters){
		StringBuilder sb = new StringBuilder();
		sb.append("#Region+Padding\tRegion\tClusterId\tLength\t#Binding\t#PWMs\t#Kmers\tSite:Positions\tPadded_Sequence\n");
		ArrayList<WeightMatrixScorer> scorers = new ArrayList<WeightMatrixScorer>();
		ArrayList<Integer> wmLens = new ArrayList<Integer>();
		ArrayList<Double> wmThresholds = new ArrayList<Double>();
		for (WeightMatrix wm: pwms){
			scorers.add(new WeightMatrixScorer(wm));
			wmLens.add(wm.length());
			wmThresholds.add(wm.getMaxScore()*wm_factor);
		}

		int digits_binding = (int) Math.ceil(Math.log10(expts.size()));
		int digits_pwm = (int) Math.ceil(Math.log10(pwms.size()));
		int digits_kmer = (int) Math.ceil(Math.log10(kmers.size()));
		
		boolean motif_kmer_scan = pwms.size()+kmers.size() > 0;
		for (int id=0;id<clusters.size();id++){
			ArrayList<Site> bindingSites = clusters.get(id);
			int numSite = bindingSites.size();
			Region r = null;
			if (rs!=null)
				r = rs.get(id);
			else
				r = new Region(genome, bindingSites.get(0).bs.getChrom(), bindingSites.get(0).bs.getLocation(), bindingSites.get(numSite-1).bs.getLocation());	
			// Position 0: the start position of the padded cluster region
			Region region = r.expand(cluster_motif_padding,cluster_motif_padding);
			int start = region.getStart();
			
			// Binding list
			ArrayList<Integer> bindingIds = new ArrayList<Integer>(); 	// binding factor id
			ArrayList<Integer> bindingPos = new ArrayList<Integer>();	// position in the padded region
			ArrayList<Integer> eventIds = new ArrayList<Integer>();		// event id in the original site list
			ArrayList<Character> strands = new ArrayList<Character>();	// strand of site match
			ArrayList<Double> bindingStrength = new ArrayList<Double>();// binding strength, read count of the event
			for (Site s:bindingSites){
				bindingIds.add(s.tf_id);
				bindingPos.add(s.bs.getLocation()-start);
				eventIds.add(s.event_id);
				bindingStrength.add(s.signal);
				strands.add(s.motifStrand);
			}
			
			// PWM motif matches
			ArrayList<Integer> pwmMatchIds = new ArrayList<Integer>();
			ArrayList<Integer> pwmMatchPos = new ArrayList<Integer>();
			ArrayList<Integer> kmerMatchIds = new ArrayList<Integer>();
			ArrayList<Integer> kmerMatchPos = new ArrayList<Integer>();
			String seq = null;
			if (motif_kmer_scan && numSite>0){
				seq = seqgen.execute(region).toUpperCase();
			
				// scan motif matches in the cluster region sequence
				for (int i=0;i<pwms.size();i++){
					ArrayList<Integer> matchPos = CommonUtils.getAllPWMHit(seq, wmLens.get(i), scorers.get(i), wmThresholds.get(i));
					for (int p:matchPos){
						pwmMatchIds.add(i);
						pwmMatchPos.add(p);
					}
				}
				
				// K-mer matches
				for (int i=0;i<kmers.size();i++){
					String kmer = kmers.get(i);
					ArrayList<Integer> matchPos = CommonUtils.getAllKmerHit(seq, kmer);
					for (int p:matchPos){
						kmerMatchIds.add(i);
						kmerMatchPos.add(p);
					}
				}
			}
			sb.append(String.format("%s\t%s\t%d\t%d\t%d\t%d\t%d\t", region.toString(), r.toString(),
					id, r.getWidth(), numSite, pwmMatchIds.size(), kmerMatchIds.size()));
			if (bindingIds.size()+pwmMatchIds.size()+kmerMatchIds.size()==0){
				sb.append("\tNA");
			}
			else{
				for (int i=0;i<bindingIds.size();i++)
					sb.append(String.format("B%0"+digits_binding+"d:%d:%d:%s:%.1f ", bindingIds.get(i), bindingPos.get(i), eventIds.get(i), strands.get(i), bindingStrength.get(i)));
				for (int i=0;i<pwmMatchIds.size();i++)
					sb.append(String.format("M%0"+digits_pwm+"d:%d ", pwmMatchIds.get(i), pwmMatchPos.get(i)));
				for (int i=0;i<kmerMatchIds.size();i++)
					sb.append(String.format("K%0"+digits_kmer+"d:%d ", kmerMatchIds.get(i), kmerMatchPos.get(i)));
				sb.append("\t").append(seq);
			}
			sb.append("\n");
		}

		CommonUtils.writeFile("0_BS_Motif_clusters."+outPrefix+".d"+distance+".min"+min_site+".txt", sb.toString());
	}
	
	/**
	 * This method output the data for plotting pairwise spacing constraints (histogram and matrix)
	 */
	private void printSpacingHistrograms(){
		// The input file is the output from the type 3 method, outputBindingAndMotifSites()
		// Some filtering may be done (e.g. grep) to select desired regions/clusters.
		if (cluster_key_file==null){
			System.err.println("Please privide the experiment key file.");
			System.exit(-1);
		}
		ArrayList<String> lines0 = CommonUtils.readTextFile(cluster_key_file);
		HashMap<String, String> name_maps = new HashMap<String, String>();
		HashMap<String, String> motif_maps = new HashMap<String, String>();
		for (String l:lines0){
			String[] f = l.split("\t");
			name_maps.put(f[0], f[2]);
			motif_maps.put(f[0], f[3]);
		}
		
		ArrayList<String> lines = CommonUtils.readTextFile(cluster_file);
		ArrayList<BMCluster> clusters = new ArrayList<BMCluster>();
		TreeMap<String, TreeMap<String, SpacingProfile>> profiles = new TreeMap<String, TreeMap<String, SpacingProfile>>();
		TreeMap<String, TreeMap<String, Integer>> overlaps = new TreeMap<String, TreeMap<String, Integer>>();
		TreeMap<String,Integer> anchor_counts = new TreeMap<String,Integer>();		// total site count of this anchor
		for (String l:lines){	// each line is a cluster merged from nearby TFBS
			if (l.startsWith("#"))
				continue;
			
			BMCluster bmc = new BMCluster();			
			String[] f = l.split("\t");
			Region r = Region.fromString(genome, f[0]);
			bmc.cluster_id = Integer.parseInt(f[2]);
			String[] positions = f[7].trim().split(" ");
			bmc.anchor = new Point(r.getGenome(), r.getChrom(), r.getStart());
			
			// parse all sites
			ArrayList<RSite> sites = new ArrayList<RSite>();		// all the sites in this region (cluster/line)
			bmc.sites = sites;
			for (String ps: positions){
				if (ps.equals(""))
					continue;
				String[] bs = ps.split(":");
				RSite s = new RSite();
				s.id = Integer.parseInt(bs[0].substring(1));
				s.type = bs[0].charAt(0);
				s.label = bs[0];
				switch (s.type){
				case 'B':
					s.pos = Integer.parseInt(bs[1]);
					s.strand = bs[3].charAt(0);
					break;
				case 'M':
				case 'K':
					s.pos = Integer.parseInt(bs[1]);
					if (s.pos<0){				// for PWM and kmer, the strand is encoded as negative pos value
						s.pos =  -s.pos; 
						s.strand = '-';
					}
					else
						s.strand = '+';
					break;
				}
				sites.add(s);
			}
			
			// compute all pairwise spacings
			// anchor is at 0 position, then count the instances of target at each spacing, separate same/opposite strand
			HashSet<String> overlap_target_labels = new HashSet<String>();
			for (RSite anchor: sites){
				String ancStr = anchor.label;
				if (!profiles.containsKey(ancStr)){
					profiles.put(ancStr, new TreeMap<String,SpacingProfile>());
					overlaps.put(ancStr, new TreeMap<String,Integer>());
					anchor_counts.put(ancStr, 0);
				}
				TreeMap<String,SpacingProfile> profiles_anchor = profiles.get(ancStr);
				TreeMap<String,Integer> overlaps_anchor = overlaps.get(ancStr);
				anchor_counts.put(ancStr, anchor_counts.get(ancStr)+1);
				
				for (RSite target: sites){
					String tarStr = target.label;
					if (!profiles_anchor.containsKey(tarStr))
						profiles_anchor.put(tarStr, new SpacingProfile());
					SpacingProfile pf = profiles_anchor.get(tarStr);
					Pair<Character, Integer> p = getProfileIndex(anchor, target);
					if (p==null)
						continue;
					switch(p.car()){
					case 's':
						pf.profile_same[p.cdr()]+=1;
						break;
					case 'd':
						pf.profile_diff[p.cdr()]+=1;
						break;
					case 'u':
						pf.profile_unknown[p.cdr()]+=1;
					}
					int offset = target.pos - anchor.pos;
					if (offset<=profile_range && offset>=-profile_range)
						overlap_target_labels.add(tarStr);						
				}	// each site as target
				
				for (String s:overlap_target_labels){
					if (!overlaps_anchor.containsKey(s))
						overlaps_anchor.put(s, 1);
					else
						overlaps_anchor.put(s, overlaps_anchor.get(s)+1);
				}
			} // each site as anchor
			clusters.add(bmc);
		}// for each line
					
		StringBuilder sb_count = new StringBuilder(); 		// the count of tallest bar of spacings
		StringBuilder sb_offset = new StringBuilder("Anchor\t");		// the offset of tallest bar of spacings
		StringBuilder sb_overlap = new StringBuilder();
		StringBuilder sb_overlap_fraction = new StringBuilder("Overlap fraction\nAnchor\t");
		sb_count.append("\nThe following tables are similar to the pairwise spacing matrix in GEM paper.\n");
		sb_count.append("Anchor\t");
		sb_overlap.append("Number of events in rows that are covered by events in column.\n");
		sb_overlap.append("Anchor").append("\t");

		for (String ancStr: profiles.keySet()){
			sb_count.append(name_maps.get(ancStr)+"\t");
			sb_offset.append(name_maps.get(ancStr)+"\t");
			sb_overlap.append(name_maps.get(ancStr)+"\t");
			sb_overlap_fraction.append(name_maps.get(ancStr)+"\t");
		}
		CommonUtils.replaceEnd(sb_count, '\n');
		CommonUtils.replaceEnd(sb_offset, '\n');
		CommonUtils.replaceEnd(sb_overlap, '\n');
		CommonUtils.replaceEnd(sb_overlap_fraction, '\n');
		sb_overlap.append("Total").append("\t");
		for (String ancStr: profiles.keySet()){
			sb_overlap.append(anchor_counts.get(ancStr)+"\t");
		}
		CommonUtils.replaceEnd(sb_overlap, '\n');	
		
		StringBuilder sb_strong_spacings = new StringBuilder("A-Code\tT-Code\tAnchor\tTarget\tOrient\tPositn\tMaxCt\tBgCt\tP-value\tFold\n");	
		double bonferronni_factor = profiles.size()*profiles.size()*profile_range/2;
		DRand re = new DRand();
		Poisson poissonEngine = new Poisson(0, re);
		for (String ancStr: profiles.keySet()){
			
			StringBuilder sb_profiles = new StringBuilder(name_maps.get(ancStr)+"\t");			
			for (int i=-profile_range;i<=profile_range;i++)
				sb_profiles.append(i+"\t");
			CommonUtils.replaceEnd(sb_profiles, '\n');
			
			TreeMap<String,SpacingProfile> profiles_anchor = profiles.get(ancStr);
			sb_count.append(name_maps.get(ancStr)+"\t");
			sb_offset.append(name_maps.get(ancStr)+"\t");
			sb_overlap.append(name_maps.get(ancStr)+"\t");
			sb_overlap_fraction.append(name_maps.get(ancStr)+"\t");
			for (String tarStr: profiles.keySet()){
				SpacingProfile pf = profiles_anchor.get(tarStr);
				if (pf==null)
					pf = new SpacingProfile();
				
				double sum = 0;
				for (int i=0;i<profile_range/2;i++){
					sum += pf.profile_same[i];
				}
				for (int i=profile_range+profile_range/2+1;i<pf.profile_same.length;i++){
					sum += pf.profile_same[i];
				}
				for (int i=0;i<profile_range/2;i++){
					sum += pf.profile_diff[i];
				}
				for (int i=profile_range+profile_range/2+1;i<pf.profile_same.length;i++){
					sum += pf.profile_diff[i];
				}
				for (int i=0;i<profile_range/2;i++){
					sum += pf.profile_unknown[i];
				}
				for (int i=profile_range+profile_range/2+1;i<pf.profile_same.length;i++){
					sum += pf.profile_unknown[i];
				}
				double mean = sum/profile_range;	// the background per-base spacing count, average from 100-200 positions
				poissonEngine.setMean(mean+1);		// add 1 as pseudo-count when computing Poisson p-value
			
				boolean has_strong_spacing = false;
				ProfileStats stats_same = getProfileStats(pf.profile_same, ancStr.equals(tarStr));
				// p-value as the tail of Poisson, add the pdf term to gain numeric precision
				double qvalue_log10 = -Math.log10((1-poissonEngine.cdf(stats_same.max+1)+poissonEngine.pdf(stats_same.max+1))*bonferronni_factor);
				if (qvalue_log10>spacing_cutoff){
					sb_strong_spacings.append(String.format("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%.1f\t%.1f\t%.1f\n", ancStr,tarStr,name_maps.get(ancStr),name_maps.get(tarStr),1,stats_same.max_idx.first()-profile_range,stats_same.max,mean,qvalue_log10,stats_same.max/mean));
					has_strong_spacing = true;
				}
				ProfileStats stats_diff = getProfileStats(pf.profile_diff, ancStr.equals(tarStr));
				qvalue_log10 = -Math.log10((1-poissonEngine.cdf(stats_diff.max+1)+poissonEngine.pdf(stats_diff.max+1))*bonferronni_factor);
				if (qvalue_log10>spacing_cutoff){
					sb_strong_spacings.append(String.format("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%.1f\t%.1f\t%.1f\n", ancStr,tarStr,name_maps.get(ancStr),name_maps.get(tarStr),-1,stats_diff.max_idx.first()-profile_range,stats_diff.max,mean,qvalue_log10,stats_diff.max/mean));
					has_strong_spacing = true;
				}
				ProfileStats stats_unknown = getProfileStats(pf.profile_unknown, ancStr.equals(tarStr));
				qvalue_log10 = -Math.log10((1-poissonEngine.cdf(stats_unknown.max+1)+poissonEngine.pdf(stats_unknown.max+1))*bonferronni_factor);
				if (qvalue_log10>spacing_cutoff){
					sb_strong_spacings.append(String.format("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%.1f\t%.1f\t%.1f\n", ancStr,tarStr,name_maps.get(ancStr),name_maps.get(tarStr),0,stats_unknown.max_idx.first()-profile_range,stats_unknown.max,mean,qvalue_log10,stats_unknown.max/mean));
					has_strong_spacing = true;
				}
				
				ProfileStats max_all=null;
				if (stats_same.max>=stats_diff.max)
					max_all = stats_same;
				else 
					max_all = stats_diff;
				if (max_all.max<stats_diff.max)
					max_all = stats_unknown;
				sb_count.append(max_all.max+"\t");
				if (has_strong_spacing && max_all.max>10){
					for (int idx:max_all.max_idx){
						sb_offset.append((idx-profile_range)+",");
					}
					CommonUtils.replaceEnd(sb_offset, '\t');
				}
				else
					sb_offset.append("999\t");
				Integer ov = overlaps.get(ancStr).get(tarStr);
				sb_overlap.append((ov==null?0:ov.intValue())+"\t");
				sb_overlap_fraction.append(String.format("%.1f", (ov==null?0:ov.intValue()*100.0/anchor_counts.get(ancStr)))+"\t");
				
				sb_profiles.append(name_maps.get(tarStr)+"_s\t").append(CommonUtils.arrayToString(pf.profile_same)).append("\n");
				sb_profiles.append(name_maps.get(tarStr)+"_d\t").append(CommonUtils.arrayToString(pf.profile_diff)).append("\n");
				sb_profiles.append(name_maps.get(tarStr)+"_u\t").append(CommonUtils.arrayToString(pf.profile_unknown)).append("\n");
			}
			CommonUtils.replaceEnd(sb_count, '\n');
			CommonUtils.replaceEnd(sb_offset, '\n');
			CommonUtils.replaceEnd(sb_overlap, '\n');	
			CommonUtils.replaceEnd(sb_overlap_fraction, '\n');
			
			CommonUtils.writeFile(outPrefix+"_"+ancStr+"_profiles.txt", sb_profiles.toString());
		}
		sb_overlap.append("\n").append(sb_overlap_fraction.toString());
		sb_overlap.append(sb_count.toString()).append("\n").append(sb_offset.toString());
		System.out.println(sb_overlap.toString());
		CommonUtils.writeFile(outPrefix+"_spacing_tables.txt", sb_overlap.toString());
		CommonUtils.writeFile(outPrefix+"_spacing_lists.txt", sb_strong_spacings.toString());
	}
	
	private ProfileStats getProfileStats(int[] profile, boolean self){
		int[] tmp=null;
		if (self){
			tmp = profile.clone();
			tmp[profile_range]=0;
		}
		else
			tmp = profile;
		Pair<Integer, TreeSet<Integer>> found = StatUtil.findMax(tmp);
		ProfileStats stats = new ProfileStats();
		stats.max = found.car();
		stats.max_idx = found.cdr();
		return stats;
	}
	
	private class ProfileStats{
		int max;
		TreeSet<Integer> max_idx;
	}

	/**
	 * Plot the aligned binding sites given the anchor, target and sorting sites<br>
	 * output matlab data file, anchor site coordinates and corresponding sequence plot
	 */
	private void plotAlignedSites(){
		String fas[] = anchor_string.split(":");
		String anchorLabel = fas[0];
		int minRange = Integer.parseInt(fas[1]);		// plotting range
		int maxRange = Integer.parseInt(fas[2]);
		if (minRange > maxRange){
			System.err.println("The plotting range is not valid: "+anchor_string);
			return;
		}

		String fts[] = target_string.split(":");
		String targetLabel = fts[0];
		int minOffset = Integer.parseInt(fts[1]);		// target spacing selection range
		int maxOffset = Integer.parseInt(fts[2]);
		if (minOffset > maxOffset){
			System.err.println("The target offset range is not valid: "+target_string);
			return;
		}
	
		String anchorTargetPairName = "a"+anchorLabel +"_t"+ targetLabel;	
		System.out.println("Anchor: "+anchorLabel+ "  Target: "+targetLabel+"  Sortby: "+sort_string);

		ArrayList<String> lines = CommonUtils.readTextFile(cluster_file);
		ArrayList<BMCluster> clusters = new ArrayList<BMCluster>();
		SpacingProfile pf = new SpacingProfile();		// spacing histogram to help find the correct position pair of anchor and target
		for (String l:lines){	// each line is a cluster merged from nearby TFBS
			if (l.startsWith("#"))
				continue;
			
			BMCluster bmc = new BMCluster();
			
			String[] f = l.split("\t");
			Region r = Region.fromString(genome, f[0]);
			bmc.cluster_id = Integer.parseInt(f[2]);
			String[] positions = f[7].split(" ");
			bmc.anchoredSequence = f[8];
			bmc.anchor = new Point(r.getGenome(), r.getChrom(), r.getStart());
			
			// parse all sites
			ArrayList<RSite> sites = new ArrayList<RSite>();
			bmc.sites = sites;
			for (String ps: positions){
				String[] bs = ps.split(":");
				RSite s = new RSite();
				s.id = Integer.parseInt(bs[0].substring(1));
				s.type = bs[0].charAt(0);
				s.label = bs[0];
				switch (s.type){
				case 'B':
					s.pos = Integer.parseInt(bs[1]);
					s.strand = bs[3].charAt(0);
					break;
				case 'M':
				case 'K':
					s.pos = Integer.parseInt(bs[1]);
					if (s.pos<0){				// for PWM and kmer, the strand is encoded as negative pos value
						s.pos =  -s.pos; 
						s.strand = '-';
					}
					else
						s.strand = '+';
					break;
				}
				sites.add(s);
			}
			
			// compute all pairwise spacings between anchor and target across all regions, make spacing histogram
			// anchor is at 0 position, then count the instances of target at each spacing, separate same/opposite strand
			for (RSite anchor: sites){
				String ancStr = anchor.label;
				if (!anchorLabel.equals(ancStr))
					continue;
				for (RSite target: sites){
					String tarStr = target.label;
					if (!targetLabel.equals(tarStr))
						continue;
					Pair<Character, Integer> p = getProfileIndex(anchor, target);
					if (p==null)
						continue;
					switch(p.car()){
					case 's':
						pf.profile_same[p.cdr()]+=1;
						break;
					case 'd':
						pf.profile_diff[p.cdr()]+=1;
						break;
					case 'u':
						pf.profile_unknown[p.cdr()]+=1;
					}
				}	// each site as target
			} // each site as anchor
			clusters.add(bmc);
		}// for each line
			
		// for each cluster, determine the anchor site, and the sorting index (i.e. pos to anchor)
		// (give priority to the sites that shows strong spacing constraint when there are multiple possible combinations)
		// then update strand, sequence and offset based on anchor site
		ArrayList<BMCluster> outClusters = new ArrayList<BMCluster>();
		for (BMCluster c:clusters){
			ArrayList<RSite> anchorSites = new ArrayList<RSite>();
			ArrayList<RSite> targetSites = new ArrayList<RSite>();
			ArrayList<RSite> sortbySites = new ArrayList<RSite>();
			for (RSite s: c.sites){
				if (anchorLabel.equals(s.label))
					anchorSites.add(s);
				if (targetLabel.equals(s.label))
					targetSites.add(s);
				if (sort_string.equals(s.label))
					sortbySites.add(s);
			}
			
			// test all anchor-target paring, give priority to the sites that shows strong spacing constraint when there are multiple possible combinations)
			int bestCount=-1;
			RSite bestAnchor = null;
			int bestOffset = 999;
			int bestSort =999;
			for (RSite as: anchorSites){
				for(RSite ts: targetSites){
					Pair<Character, Integer> p = getProfileIndex(as, ts);
					if (p==null)
						continue;
					int count = -1;
					switch(p.car()){
					case 's':
						count = pf.profile_same[p.cdr()];
						break;
					case 'd':
						count = pf.profile_diff[p.cdr()];
						break;
					case 'u':
						count = pf.profile_unknown[p.cdr()];
					}
					if (count>bestCount){
						bestCount = count;
						bestAnchor = as;
						bestOffset = getProfileIndex(as, ts).cdr()-profile_range;
						if (bestOffset>=minOffset && bestOffset<=maxOffset){
							if (sortbySites.isEmpty())
								bestSort = 999;			// if the sorting TF is not found in this region, move the region to the end
							else{
								Pair<Character, Integer> sort_idx = getProfileIndex(as, sortbySites.get(0)); // just pick one if having 1+ sorting site
								if (sort_idx!=null)
									bestSort = sort_idx.cdr()-profile_range;
							}
						}
					}					
				}
			}
			if (bestAnchor==null)// Skip if no anchor binding site is found 
				continue;
			else{
				c.updateAnchorSite(bestAnchor);
				c.sort_idx = bestSort;
				if (bestOffset>=minOffset && bestOffset<=maxOffset)
					outClusters.add(c);
			}
//			System.out.println(String.format("Count:\t%d\tOffset: %d\tPos:%d", bestCount, c.sort_idx, bestAnchor.pos));
		}// each cluster
		
		if (outClusters.isEmpty()){
			System.err.println("Anchor "+anchor_string+" matched no sites!");
			System.exit(-1);
		}
		else{
			System.out.println(anchorTargetPairName+" found "+outClusters.size()+" regions.");
		}
		Collections.sort(outClusters);
				
		StringBuilder sb = new StringBuilder();
		StringBuilder sb_coords = new StringBuilder();
		sb.append(outClusters.get(0).sites.get(0).getMatlabHeader());
		for (int i=0;i<outClusters.size();i++){
			BMCluster bmc = outClusters.get(i);
			for (RSite s: bmc.sites){
				sb.append(s.toMatlabData(i));				
			}
			sb_coords.append(bmc.getAnchorStrandedPoint()).append("\t").append(bmc.cluster_id).append("\n");
		}
		CommonUtils.writeFile(cluster_file.replace("txt", anchorTargetPairName+"_matlab.txt"), sb.toString());
		CommonUtils.writeFile(cluster_file.replace("txt", anchorTargetPairName+"_coords.txt"), sb_coords.toString());
		
		String[] ss = new String[outClusters.size()];
		for (int i=0;i<outClusters.size();i++){
			BMCluster c = outClusters.get(i);
			String anchoredSequence = c.anchoredSequence;
			int len = anchoredSequence.length();
			int anchor = c.anchorSite.pos;
			if (c.anchorSite.strand=='-'){
				anchoredSequence = SequenceUtils.reverseComplement(anchoredSequence);
				anchor = len - anchor - 1;
			}
			int left = anchor + minRange;
			int leftPadding = 0;
			if (left<0){
				leftPadding = -left;
				left=0;
			}
			int right = anchor + maxRange+1;  // add 1 because substring() is end-exclusive
			int rightPadding = 0;
			if (right>len){
				rightPadding = right-len;
				right=len;
			}
			ss[i] = CommonUtils.padding(leftPadding, 'N')+anchoredSequence.substring(left, right)+CommonUtils.padding(rightPadding, 'N');
		}
		CommonUtils.visualizeSequences(ss, width, height, new File(cluster_file.replace("txt", anchorTargetPairName+"_seqPlot.png")));
	}


	/** compute the relative index of target site given the anchor site as 0 position
	 * 
	 */
	private Pair<Character, Integer> getProfileIndex(RSite anchor, RSite target){
		int offset = target.pos - anchor.pos;
		char type = 'u';
		if (offset>profile_range || offset <-profile_range)			// skip if out of range
			return null;
		if(anchor.strand=='-'){
			offset = -offset;
			offset += profile_range;		// shift to get positive array idx
			if (target.strand=='+')
				type = 'd';
			if (target.strand=='-')
				type = 's';
		}
		if(anchor.strand=='+'){
			offset += profile_range;		// shift to get positive array idx
			if (target.strand=='+')
				type = 's';
			if (target.strand=='-')
				type = 'd';
		}
		if(anchor.strand=='*')
			offset += profile_range;
		return new Pair<Character, Integer>(type, offset);
	}
	
	/** 
	 * RSite: data structure of a binding/motif site relative to an anchor point<br>
	 * The position of the site is relative to an anchor point (stranded)
	 */
	private class RSite implements Cloneable{
		char type = ' ';		// type=0:binding, 1:pwm, 2:kmer
		int id=-1;			// id in that type
		String label="";
		int pos=999;		// position relative to the anchor site
		char strand;		// binding site or motif match strand on the original DNA sequence
		public String toString(){
			return String.format("%d:%s:%s%d ", pos, strand, type, id);
		}
		/** print out this site with a yValue (e.g. the sorted cluster id) specified by user
		 * 
		 */
		public String toMatlabData(int yValue){
			return String.format("%d\t%d\t%d\t%d\t%d\n", type=='B'?0:type=='M'?1:2, id, pos, yValue, (int)strand);
		}
		public String getMatlabHeader(){
			return "B-M-K\tSiteID\tPos\tSort\tStrand\n";
		}
		protected Object clone(){
			RSite clone=null;
			try{
				clone =(RSite)super.clone();
			}
			catch (CloneNotSupportedException e){
				System.out.println(e.toString());
			}
		    return clone;
		}
	}
	/** 
	 * BMCluster: Binding and Motif Cluster<br>
	 * Data structure to store a cluster of binding/motif sites
	 */
	private class BMCluster implements Comparable<BMCluster> {
		private RSite anchorSite;
		int cluster_id = -1;
		int sort_idx = 999;
		Point anchor = null;
		String anchoredSequence = null;
		ArrayList<RSite> sites = null;
		
		/** Comparable default method, sort by increasing offset */
		public int compareTo(BMCluster f) {
			int diff = sort_idx-f.sort_idx;
			return diff>0?1:diff==0?0:-1;
		}
		/** 
		 * Update the anchor site for this cluster of sites<br>
		 * Also update the position of sites as relative positions to this anchor site
		 * @param as
		 */
		void updateAnchorSite(RSite as){
			anchorSite = (RSite)as.clone();
			int pos = as.pos;
			anchor = new Point(anchor.getGenome(), anchor.getChrom(), anchor.getLocation()+pos);
			for (RSite s: sites){
				s.pos -= pos;
				if (as.strand=='-')
					s.pos = -s.pos;
			}
		}
		public String getAnchorStrandedPoint(){
			if (anchorSite.strand=='*')
				return anchor.toString();
			else
				return new StrandedPoint(anchor, anchorSite.strand).toString();
		}
	}
	
	private class SpacingProfile{
		int[] profile_same = new int[profile_range*2+1];		// the anchor site and the subject site are on the same strand (or same orientation of motifs)
		int[] profile_diff = new int[profile_range*2+1];		// the anchor site and the subject site are on the opposite strand (or opposite orientation of motifs)
		int[] profile_unknown = new int[profile_range*2+1];	// unknown, because some GEM binding site does not have a k-mer to assgin strand
	}
	
	private void computeTfbsSpacingDistribution(){
		// classify sites by chrom
		TreeMap<String, ArrayList<Site>> chrom2sites = new TreeMap<String, ArrayList<Site>>();
		for (ArrayList<Site> sites:all_sites){
			for (Site s:sites){
				String chr = s.bs.getChrom();
				if (!chrom2sites.containsKey(chr))
					chrom2sites.put(chr, new ArrayList<Site>());
				chrom2sites.get(chr).add(s);
			}
		}
		
		// sort sites and compute TFBS spacings
		int[] counts = new int[2001];
		for (String chr: chrom2sites.keySet()){
			ArrayList<Site> sites = chrom2sites.get(chr);
			Collections.sort(sites);
//			int previousSpacing = Integer.MAX_VALUE;
			for (int i=1;i<sites.size();i++){	
				int spacing = sites.get(i).bs.getLocation()-sites.get(i-1).bs.getLocation();
//				int minSpacing = Math.min(spacing, previousSpacing);
//				previousSpacing = spacing;
				if (spacing<=2000)
					counts[spacing]++;
			}
		}
		
		// output spacing distributions
		StringBuilder sb = new StringBuilder();
		for (int i=0;i<counts.length;i++){
			sb.append(i+"\t"+counts[i]).append("\n");
		}
		CommonUtils.writeFile("0_BS_spacing_histrogram."+outPrefix+".txt", sb.toString());
		
		sb = new StringBuilder();
		for (int mLen=0;mLen<2500;mLen+=5){
			int sum=0;
			int sum_length=0;
			for (String chr: chrom2sites.keySet()){
				ArrayList<Site> sites = chrom2sites.get(chr);
				ArrayList<Region> rs = new ArrayList<Region>();
				for (Site s:sites){
					rs.add(s.bs.expand(mLen));
				}
				rs = Region.mergeRegions(rs);;
				sum += rs.size();
				for (Region r:rs)
					sum_length+=r.getWidth()-2*mLen;		// subtract the expanded length, border with TF sites
			}
			sb.append(2*mLen+"\t"+sum+"\t"+sum_length).append("\n");
			System.err.print(2*mLen+" ");
		}
		CommonUtils.writeFile("0_BS_mergeLength_stats."+outPrefix+".txt", sb.toString());		
		
		// each site - nearest sites from all other TF data, distance distribution (Yan ... Taipale, 2013, Cell, Cohesin Memory)
		ArrayList<TreeMap<String, int[]>> TF_chrom_coord = new ArrayList<TreeMap<String, int[]>>();
		for (ArrayList<Site> sites:all_sites){
			TreeMap<String, ArrayList<Integer>> chrom_coord = new TreeMap<String, ArrayList<Integer>>();
			for (Site s:sites){
				String chr = s.bs.getChrom();
				if (!chrom_coord.containsKey(chr))
					chrom_coord.put(chr, new ArrayList<Integer>());
				chrom_coord.get(chr).add(s.bs.getLocation());
			}
			TreeMap<String, int[]> chrom_coord2 = new TreeMap<String, int[]>();
			for (String key: chrom_coord.keySet()){
				ArrayList<Integer> coords = chrom_coord.get(key);
				int[] coords2 = new int[coords.size()];
				for (int i=0;i<coords.size();i++)
					coords2[i]=coords.get(i);
				chrom_coord2.put(key, coords2);
			}
			TF_chrom_coord.add(chrom_coord2);
		}
		TreeMap<Integer, Integer> distanceHistogram = new TreeMap<Integer, Integer>();
		for (int i=0;i<all_sites.size();i++){
			ArrayList<Site> sites = all_sites.get(i);
			for (Site s: sites){
				for (int j=0;j<TF_chrom_coord.size();j++){		// each TF, not including TF itself
					if (i==j)
						continue;
					int distance = 0;
					String chr = s.bs.getChrom();
					int coord = s.bs.getLocation();
					TreeMap<String, int[]> chrom_coord = TF_chrom_coord.get(j);
					if (chrom_coord.containsKey(chr)){
						int[] coords = chrom_coord.get(chr);
						int idx = Arrays.binarySearch(coords, coord);
						if (idx>=0){	// found the same coord
							distance = 0;
						}
						else{
							idx = -idx-1;
							if (idx==coords.length)	// all less
								distance = Math.abs(coord-coords[idx-1]);
							else if (idx==0)
								distance = Math.abs(coord-coords[idx]);
							else
								distance = Math.min(Math.abs(coord-coords[idx-1]),Math.abs(coord-coords[idx]));
						}
						if (!distanceHistogram.containsKey(distance)){
							distanceHistogram.put(distance, 0);
						}
						distanceHistogram.put(distance, distanceHistogram.get(distance)+1);
					}
				}
			}
		}
		sb = new StringBuilder();
		for (int d: distanceHistogram.keySet())
			sb.append(d).append("\t").append(distanceHistogram.get(d)).append("\n");
		CommonUtils.writeFile("0_BS_pairwiseTF_distanceHistogram."+outPrefix+".txt", sb.toString());
	}
	
	private void mergedTSS(){
		
		if (tss_file != null){
			ArrayList<Point> tsss = new ArrayList<Point>();
			ArrayList<String> text = CommonUtils.readTextFile(tss_file);
			for (String t: text){
				String[] f = t.split("\t");
				tsss.add(Point.fromString(genome, f[0]));
			}
			
			// classify tss by chrom
			TreeMap<String, ArrayList<Point>> chrom2sites = new TreeMap<String, ArrayList<Point>>();
			for (Point s:tsss){
				String chr = s.getChrom();
				if (!chrom2sites.containsKey(chr))
					chrom2sites.put(chr, new ArrayList<Point>());
				chrom2sites.get(chr).add(s);
			}
			
			// sort tss and form clusters
			ArrayList<ArrayList<Point>> clusters = new ArrayList<ArrayList<Point>>();		
			ArrayList<Point> cluster = new ArrayList<Point>();
			for (String chr: chrom2sites.keySet()){
				ArrayList<Point> sites = chrom2sites.get(chr);
				Collections.sort(sites);

				cluster.add(sites.get(0));
				for (int i=1;i<sites.size();i++){
					Point s = sites.get(i);
					Point p = cluster.get(cluster.size()-1);		//previous
					if (s.getLocation()-p.getLocation()<distance)
						cluster.add(s);
					else{
						cluster.trimToSize();
						clusters.add(cluster);
						cluster = new ArrayList<Point>();
						cluster.add(s);
					}					
				}
				// finish the chromosome
				cluster.trimToSize();
				clusters.add(cluster);
				cluster = new ArrayList<Point>();
			}
			// finish all the sites
			cluster.trimToSize();
			clusters.add(cluster);

			// merge multi-tss in a cluster into one point
			ArrayList<Region> merged = new ArrayList<Region>();			
			for (ArrayList<Point> c:clusters){
				if (c.isEmpty())
					continue;
				Region r = null;
				if (c.size()==1)
					r = c.get(0).expand(0);
				else
					r = new Region(genome, c.get(0).getChrom(), c.get(0).getLocation(), c.get(c.size()-1).getLocation());
				merged.add(r);
			}
			// output
			StringBuilder sb = new StringBuilder();
			for (Region r:merged)
				sb.append(r.getMidpoint().toString()+"\t"+r.getWidth()+"\n");
			CommonUtils.writeFile("0_mergedTSS.txt", sb.toString());
		}
	}
	private void printTssSignal(){

		final int MAXREAD = 1000000;
		
		if (tss_file != null){
			all_TSS = new ArrayList<Point>();
			ArrayList<String> text = CommonUtils.readTextFile(tss_file);
			for (String t: text){
				String[] f = t.split("\t");
				all_TSS.add(Point.fromString(genome, f[0]));
			}
			all_TSS.trimToSize();
		}
		int[][]signals = new int[all_TSS.size()][expts.size()];
		for (int i=0;i<expts.size();i++){
			String readdb_name = readdb_names.get(i);
			List<ChipSeqLocator> rdbexpts = new ArrayList<ChipSeqLocator>();
			String[] pieces = readdb_name.trim().split(";");
            if (pieces.length == 2) {
            	rdbexpts.add(new ChipSeqLocator(pieces[0], pieces[1]));
            } else if (pieces.length == 3) {
            	rdbexpts.add(new ChipSeqLocator(pieces[0], pieces[1], pieces[2]));
            } else {
                throw new RuntimeException("Couldn't parse a ChipSeqLocator from " + readdb_name);
            }
            DeepSeqExpt ip = new DeepSeqExpt(genome, rdbexpts, "readdb", -1);
            ReadCache ipCache = new ReadCache(genome, expts.get(i), null, null);
            
			// cache sorted start positions and counts of all positions
			long tic = System.currentTimeMillis();
			System.err.print("Loading "+ipCache.getName()+" data from ReadDB ... \t");
			List<String> chroms = genome.getChromList();
			if (dev){
				chroms = new ArrayList<String>();
				chroms.add("19");
			}
			// load  data into cache.
			for (String chrom: chroms ){
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
			ipCache.displayStats();
			System.out.println(CommonUtils.timeElapsed(tic));
            
			// now get the data from the cache
            for (int j=0;j<all_TSS.size();j++){
            	Region region = all_TSS.get(j).expand(signal_radius);
            	List<StrandedBase> bases = ipCache.getStrandedBases(region, '+');
            	bases.addAll(ipCache.getStrandedBases(region, '-'));
            	signals[j][i] = (int)StrandedBase.countBaseHits(bases);
            }
		}
		StringBuilder sb = new StringBuilder("#TSS\t");
		for (int i=0;i<expts.size();i++){
			sb.append(expts.get(i)).append("\t");
		}
		CommonUtils.replaceEnd(sb, '\n');
		for (int j=0;j<all_TSS.size();j++){
			sb.append(all_TSS.get(j).toString()).append("\t");
			for (int i=0;i<expts.size();i++){
				sb.append(signals[j][i]).append("\t");
			}
			CommonUtils.replaceEnd(sb, '\n');
		}
		CommonUtils.writeFile("0_TSS_signals."+outPrefix+".txt", sb.toString());
	}
	
	
	private void loadClusterAndTSS(){
		if (tss_file != null){
			all_TSS = new ArrayList<Point>();
			ArrayList<String> text = CommonUtils.readTextFile(tss_file);
			for (String t: text){
				String[] f = t.split("\t");
				all_TSS.add(Point.fromString(genome, f[0]));
			}
		}
		
		if (cluster_file != null){
			all_clusters = new ArrayList<Cluster>();
			ArrayList<String> text = CommonUtils.readTextFile(cluster_file);
			for (String t: text){
				if (t.startsWith("#"))
					continue;
				String[] f = t.split("\t");
				Cluster c = new Cluster();
				all_clusters.add(c);
				c.region = Region.fromString(genome, f[0]);
				String[] f_id = f[4].split(",");
				for (String s: f_id)
					c.TFIDs.add(Integer.parseInt(s));
				String[] f_signal = f[5].split(",");
				for (String s: f_signal)
					c.TF_Signals.add(Double.parseDouble(s));
				String[] f_motif = f[6].split(",");
				for (String s: f_motif)
					c.TF_hasMotifs.add(Integer.parseInt(s)==1);
			}
		}
		
	}
	private void loadClusterAndTssSignals(){
		if (tss_signal_file != null){
			all_TSS = new ArrayList<Point>();
			ArrayList<String> text = CommonUtils.readTextFile(tss_signal_file);
			for (String t: text){
				if (t.charAt(0)=='#')
					continue;
				String[] f = t.split("\t");
				all_TSS.add(Point.fromString(genome, f[0]));
			}
			all_TSS.trimToSize();
			
			String[] fs = text.get(0).split("\t");
			tss_signals = new int[all_TSS.size()][fs.length-1];
			int idx = 0;
			for (String t: text){
				if (t.charAt(0)=='#')
					continue;
				String[] f = t.split("\t");
				for (int i=0;i<f.length-1;i++)
					tss_signals[idx][i]=Integer.parseInt(f[i+1]);
				idx++;
			}
		}
		
		if (cluster_file != null){
			all_clusters = new ArrayList<Cluster>();
			ArrayList<String> text = CommonUtils.readTextFile(cluster_file);
			for (String t: text){
				if (t.startsWith("#"))
					continue;
				String[] f = t.split("\t");
				Cluster c = new Cluster();
				all_clusters.add(c);
				c.region = Region.fromString(genome, f[0]);
				String[] f_id = f[4].split(",");
				for (String s: f_id)
					c.TFIDs.add(Integer.parseInt(s));
				String[] f_signal = f[5].split(",");
				for (String s: f_signal)
					c.TF_Signals.add(Double.parseDouble(s));
				String[] f_motif = f[7].split(",");
				for (String s: f_motif)
					c.TF_hasMotifs.add(s.charAt(0)!='*');
			}
		}
	}

	private void computeCorrelations(){

		StringBuilder sb = new StringBuilder();
		for (Cluster c: all_clusters){
			
			Point anchor = c.region.getMidpoint();
			ArrayList<Point> targets = CommonUtils.getPointsWithinWindow(all_TSS, anchor, range);
			if (exclude_range!=0){
				ArrayList<Point> nearNeighbors = CommonUtils.getPointsWithinWindow(all_TSS, anchor, exclude_range);
				targets.removeAll(nearNeighbors);
			}
			ArrayList<Site_target_corr> list = new ArrayList<Site_target_corr>();
			
			// get unique TF IDs for consideration
			HashSet<Integer> TF_IDs = new HashSet<Integer>();
			for (int i=0;i<c.TF_hasMotifs.size();i++){
				if ( !useDirectBindingOnly || c.TF_hasMotifs.get(i) ){
					TF_IDs.add(c.TFIDs.get(i));	
				}			
			}
			Integer[] TFIDs = new Integer[TF_IDs.size()];
			TF_IDs.toArray(TFIDs);
			if (TFIDs.length<=2)
				continue;
			
			// get signals at anchor site
			List<Double> signals = new ArrayList<Double>();
			for (int i=0;i<TFIDs.length;i++){
				int total = 0;
				for (int j=0;j<c.TFIDs.size();j++)
					if (c.TFIDs.get(j)==TFIDs[i])
						total += c.TF_Signals.get(j);
				signals.add(i, (double)total);
			}
			
			StringBuilder sb1 = new StringBuilder();
			for (int id : TFIDs)
				sb1.append(tf_names.get(id)).append("\t");
			CommonUtils.replaceEnd(sb1, '\t');
			sb.append("#=================================\n#chr"+c.region.toString()+"\t|\t");
			sb.append(sb1.toString()+"|\t");

			for (double s: signals){
				sb.append(s).append("\t");
			}
			CommonUtils.replaceEnd(sb, '\n');
			
			for (Point p:targets){				
				// get corresponding signals at the target sites
				List<Double> target_signals = new ArrayList<Double>();
				int index = Collections.binarySearch(all_TSS, p);
				for (int i=0;i<TFIDs.length;i++){
					target_signals.add(i, (double)tss_signals[index][TFIDs[i]]);
				}
				
				// compute correlation
				double corr = CorrelationSimilarity.computeSimilarity2(signals, target_signals);
				Site_target_corr t = new Site_target_corr();
				t.target = p;
				t.signals = signals;
				t.target_signals = target_signals;
				t.corr = corr;
				if (t.corr>=cutoff)
					list.add(t);
			}
			
			Collections.sort(list);
			for (Site_target_corr t:list){
				sb.append("chr"+anchor.toString()+"-"+t.target.getLocation()+"\t"+t.target.offset(anchor)+"\t"
						+t.signals.size()+"\t"+String.format("%.2f", t.corr)).append("\t|\t");
				for (double s: t.target_signals){
					sb.append(s).append("\t");
				}
				CommonUtils.replaceEnd(sb, '\n');
			}
		}
		System.out.println(sb.toString());
		CommonUtils.writeFile("0_tfbs2target."+outPrefix+".txt", sb.toString());
	}

	private void computeCorrelations_db(){
		// prepare connections to readdb
		ArrayList<DeepSeqExpt> chipseqs = new ArrayList<DeepSeqExpt>();
		for (String readdb_str : readdb_names){
			List<ChipSeqLocator> rdbexpts = new ArrayList<ChipSeqLocator>();
			String[] pieces = readdb_str.trim().split(";");
            if (pieces.length == 2) {
            	rdbexpts.add(new ChipSeqLocator(pieces[0], pieces[1]));
            } else if (pieces.length == 3) {
            	rdbexpts.add(new ChipSeqLocator(pieces[0], pieces[1], pieces[2]));
            } else {
                throw new RuntimeException("Couldn't parse a ChipSeqLocator from " + readdb_str);
            }
            chipseqs.add( new DeepSeqExpt(genome, rdbexpts, "readdb", -1) );
		}
		
		for (Cluster c: all_clusters){
			
			Point anchor = c.region.getMidpoint();
			ArrayList<Point> targets = CommonUtils.getPointsWithinWindow(all_TSS, anchor, range);
			ArrayList<Site_target_corr> list = new ArrayList<Site_target_corr>();
			
			// get unique TF IDs for consideration
			HashSet<Integer> TF_IDs = new HashSet<Integer>();
			for (int i=0;i<c.TF_hasMotifs.size();i++){
				if ( !useDirectBindingOnly || c.TF_hasMotifs.get(i) ){
					TF_IDs.add(c.TFIDs.get(i));	
				}			
			}
			Integer[] TFIDs = new Integer[TF_IDs.size()];
			TF_IDs.toArray(TFIDs);
			if (TFIDs.length<=2)
				continue;
			
			// get signals at anchor site
			List<Double> signals = new ArrayList<Double>();
			for (int i=0;i<TFIDs.length;i++){
				int total = 0;
				for (int j=0;j<c.TFIDs.size();j++)
					if (c.TFIDs.get(j)==TFIDs[i])
						total += c.TF_Signals.get(j);
				signals.add(i, (double)total);
			}
			
			StringBuilder sb1 = new StringBuilder();
			for (int id : TFIDs)
				sb1.append(tf_names.get(id)).append("\t");
			System.out.println("=================================\nchr"+c.region.toString());
			System.out.println(sb1.toString());
			
			for (Point p:targets){				
				// get corresponding signals at the target sites
				List<Double> target_signals = new ArrayList<Double>();
				for (int i=0;i<TFIDs.length;i++){
					target_signals.add(i, (double)chipseqs.get(TFIDs[i]).countHits(p.expand(signal_radius)));
				}
				
				// compute correlation
				double corr = CorrelationSimilarity.computeSimilarity2(signals, target_signals);
				Site_target_corr t = new Site_target_corr();
				t.target = p;
				t.signals = signals;
				t.target_signals = target_signals;
				t.corr = corr;
				if (t.corr>=cutoff)
					list.add(t);
			}
			
			Collections.sort(list);
			for (Site_target_corr t:list){

				System.out.println("chr"+anchor.toString()+"-"+t.target.getLocation()+"\t"+t.target.offset(anchor)+"\t"
						+t.signals.size()+"\t"+String.format("%.2f", t.corr));
				StringBuilder sb = new StringBuilder();
				for (double s: t.signals){
					sb.append(s).append("\t");
				}
				sb.append("\n");
				for (double s: t.target_signals){
					sb.append(s).append("\t");
				}
				System.out.println(sb.toString());
			}
		}
		
		// clean up
		for (DeepSeqExpt e: chipseqs){
			e.closeLoaders();
		}
	}
	
	class Site_target_corr implements Comparable<Site_target_corr>{
		Point target;
		List<Double> signals;
		List<Double> target_signals;
		double corr;
		public int compareTo(Site_target_corr t) {					// descending corr
			return corr>t.corr?-1:corr==t.corr?0:1;
		}
	}
}

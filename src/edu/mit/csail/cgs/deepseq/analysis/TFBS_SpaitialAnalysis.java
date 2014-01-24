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
	private final int TARGET_WIDTH = 250;
	private final int PROFILE_RANGE = 100;
	Genome genome=null;
	ArrayList<String> expts = new ArrayList<String>();
	ArrayList<String> names = new ArrayList<String>();
	ArrayList<String> readdb_names = new ArrayList<String>();
	ArrayList<String> indirect_tf_expts = new ArrayList<String>();
	/** Mapping of TFSS(sequence specific) GEM expt ids (from direct binding id to indirect binding id */
	HashMap<Integer, Integer> directid2indirectid = null;
	
	ArrayList<WeightMatrix> pwms = new ArrayList<WeightMatrix>();
	ArrayList<String> kmers = new ArrayList<String>();
	ArrayList<ArrayList<Site>> all_sites = new ArrayList<ArrayList<Site>>();
	ArrayList<ArrayList<Site>> all_indirect_sites = new ArrayList<ArrayList<Site>>();
	ArrayList<Point> all_TSS;
	ArrayList<Cluster> all_clusters;
	int[][]tss_signals = null;
	
	double gc = 0.42;//mouse		gc=0.41 for human
	int distance = 50;		// distance between TFBS within a cluster
	int range = 1000;		// the range around anchor site to search for targets
	int cluster_motif_padding = 100;  // the padding distance added to the cluster range for motif searching
	int exclude_range = 500;		// the range around anchor site to exclude same site
	int min_site = 1;				// the minimum number of sites for the cluster to be printed out
	int width = 3;
	int height = 3;
	double wm_factor = 0.6;	// PWM threshold, as fraction of max score
	double cutoff = 0.3;	// corr score cutoff
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
	String outPrefix = "out";
	String tfss_file;				// Sequence-specific TFs
	String tss_file;
	String pwm_file;
	String kmer_file;
	String tss_signal_file;
	String cluster_file;
	String exclude_sites_file;
	String anchor_string;	// the id of TF to anchor the sites/regions/sequences
	String sort_string;		// id:left:right, the id of TF to sort the sites/regions/sequences according to its offset from anchor
	
	// command line option:  (the folder contains GEM result folders) 
	// --dir C:\Data\workspace\gse\TFBS_clusters --species "Mus musculus;mm9" --r 2 --pwm_factor 0.6 --expts expt_list.txt [--no_cache --old_format] 
	public static void main(String[] args) {
		TFBS_SpaitialAnalysis mtb = new TFBS_SpaitialAnalysis(args);
		int round = Args.parseInteger(args, "r", 2);
		int type = Args.parseInteger(args, "type", 0);
		ArrayList<ArrayList<Site>> clusters=null;
		switch(type){
		case 0:
			mtb.loadEventsAndMotifs();
			clusters = mtb.mergeTfbsClusters();
			mtb.outputTFBSclusters(clusters);
			break;
		case 1:
			mtb.loadClusterAndTssSignals();
			mtb.computeCorrelations();
			break;
		case 2:
			mtb.loadEventsAndMotifs();
			mtb.computeTfbsSpacingDistribution();
			break;
		case 3:		// to print all the binding sites and motif positions in the clusters for downstream spacing/grammar analysis
			mtb.loadEventsAndMotifs2(round);
			clusters = mtb.mergeTfbsClusters();
			mtb.outputBindingAndMotifSites(clusters);
			break;
		case 4:		// spacing histogram
			mtb.printSpacingHistrograms();
			break;
		case 5:		// anchor/sorted sites for matlab plot and sequences
			mtb.plotAlignedSites();
			break;
		case 11:	// old code
			mtb.loadClusterAndTSS();
			mtb.computeCorrelations_db();
			break;		
		case -1:
			mtb.mergedTSS();
			break;
		case -2:
			mtb.printTssSignal();
			break;
		}
		
	}
	
	public TFBS_SpaitialAnalysis(String[] args){
				
	    try {
	    	Pair<Organism, Genome> pair = Args.parseGenome(args);
	        if(pair != null) {
	            genome = pair.cdr();
	        } else {
	            String genomeString = Args.parseString(args,"g",null);		// text file with chrom lengths
	            if(genomeString != null){
	                genome = new Genome("Genome", new File(genomeString), true);
	            } else{
	                genome=null;
	            }
	        }
	    } catch (NotFoundException e) {
	      e.printStackTrace();
	    }

		Set<String> flags = Args.parseFlags(args);
		outPrefix = Args.parseString(args, "out", outPrefix);
		zero_or_one = flags.contains("zoo");
		dev = flags.contains("dev");
		oldFormat = flags.contains("old_format");
		no_gem_pwm = flags.contains("no_gem_pwm");
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
		names = new ArrayList<String>();
		String info_file = Args.parseString(args, "info", null);
		if (info_file!=null){
			ArrayList<String> info = CommonUtils.readTextFile(info_file);
			for (String txt: info){
				if (!txt.equals("")){
					String[] f = txt.split("\t");
					expts.add(f[0]);
					names.add(f[3]);
					readdb_names.add(f[5]);
				}
			}
		}
		
		anchor_string = Args.parseString(args, "anchor", anchor_string);	// the id of TF/PWM/Kmer to anchor the sites/regions/sequences
		sort_string = Args.parseString(args, "sort", sort_string);			// the id of TF/PWM/Kmer to sort the sites/regions/sequences

		width = Args.parseInteger(args, "width", width);
		height = Args.parseInteger(args, "height", height);
		
		tfss_file = Args.parseString(args, "tfss_file", null);
		if (tfss_file!=null){
			directid2indirectid = new HashMap<Integer, Integer>();
			indirect_tf_expts = CommonUtils.readTextFile(tfss_file);
			for (int i=0;i<indirect_tf_expts.size();i++){
				String e = indirect_tf_expts.get(i);
				int idx = expts.indexOf(e);
				if (idx>=0){		// add fake expt and iNames for the indirect TFSS sites
					expts.add(e);
					names.add("i_"+names.get(idx));
					readdb_names.add(readdb_names.get(idx));
					directid2indirectid.put(idx, expts.size()-1);
				}
			}
		}
		
		tss_file = Args.parseString(args, "tss", null);
		tss_signal_file = Args.parseString(args, "tss_signal", null);
		cluster_file = Args.parseString(args, "cluster", null);
		pwm_file = Args.parseString(args, "pwms", null);				// additional pwms
		kmer_file = Args.parseString(args, "kmers", null);
		exclude_sites_file = Args.parseString(args, "ex", null);
		distance = Args.parseInteger(args, "distance", distance);
		range = Args.parseInteger(args, "range", range);
		cluster_motif_padding = Args.parseInteger(args, "cluster_motif_padding", cluster_motif_padding);
		exclude_range = Args.parseInteger(args, "exclude", exclude_range);
		min_site = Args.parseInteger(args, "min_site", min_site);
		wm_factor = Args.parseDouble(args, "pwm_factor", wm_factor);
		cutoff = Args.parseDouble(args, "cutoff", cutoff);
		gc = Args.parseDouble(args, "gc", gc);
		String genome_dir = Args.parseString(args, "genome", null);
		seqgen = new SequenceGenerator<Region>();
		if (genome_dir!=null)
			seqgen.setGenomePath(genome_dir);
		seqgen.useCache(!flags.contains("no_cache"));
	}
	
	private void loadEventsAndMotifs(){
		ArrayList<Region> ex_regions = new ArrayList<Region>();
		if(exclude_sites_file!=null){
			ex_regions = CommonUtils.loadRegionFile(exclude_sites_file, genome);
		}
		for (int tf=0;tf<names.size();tf++){
			if (names.get(tf).startsWith("i_"))		// names start with i_ are artificially created id for TFSS indirect binding
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
				System.out.println(expt+" does not have valid GPS/GEM event call file.");
				System.exit(1);
			}
		}
		all_sites.addAll(all_indirect_sites);
	}
	/**
	 * This is different from the loadEventsAndMotifs method, in that 
	 * it only loads the event positions, but not nearby motif match positions,
	 * it loads GEM motif associated with the events by default, but
	 * it also loads additional motifs specified by --pwms, or --ksms
	 * @param round
	 */
	private void loadEventsAndMotifs2(int round){
		ArrayList<Region> ex_regions = new ArrayList<Region>();
		if(exclude_sites_file!=null){
			ex_regions = CommonUtils.loadRegionFile(exclude_sites_file, genome);
		}
		for (int tf=0;tf<names.size();tf++){
			String expt = expts.get(tf);
			File dir2= new File(dir, expt);
			dir2= new File(dir2, expt+"_outputs");
			System.err.print(String.format("TF#%d: loading %s", tf, expt));
			
			// load binding event files 
			File gpsFile = new File(dir2, expt+"_"+ (round>=2?round:1) + "_GEM_events.txt");
			String filePath = gpsFile.getAbsolutePath();
			try{
				List<GPSPeak> gpsPeaks = GPSParser.parseGPSOutput(filePath, genome);
				ArrayList<Site> sites = new ArrayList<Site>();
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
		
		// load additional pwms
		if (pwm_file!=null){
			pwms.addAll( CommonUtils.loadPWMs_PFM_file(pwm_file, gc) );
		}
		
		// load kmers
		if (kmer_file!=null){
			kmers.addAll( CommonUtils.readTextFile(kmer_file));
		}
		
		System.out.println(String.format("In total, loaded %d GEM files, %d kmers and %d PWMs.", all_sites.size(), kmers.size(), pwms.size()));
		
		StringBuilder sb = new StringBuilder();
		for (int i=0;i<expts.size();i++){
			sb.append("B"+i+"\t"+expts.get(i)+"\n");
		}
		for (int i=0;i<pwms.size();i++){
			sb.append("M"+i+"\t"+pwms.get(i).getName()+"\n");
		}
		for (int i=0;i<kmers.size();i++){
			sb.append("K"+i+"\t"+kmers.get(i)+"\n");
		}
		CommonUtils.writeFile("0_BS_Motif_clusters."+outPrefix+"_keys.txt", sb.toString());
	}
	
	class Site implements Comparable<Site>{
		int tf_id;
		int event_id;		// original event id in the list of this TF binding data
		Point bs;
		int id;				// global id of all BS
		double signal;
		char motifStrand;								// motif match strand
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
	 * Output all the binding sites in the clusters, for topic modeling analysis or clustering analysis
	 * @param clusters
	 */
	private void outputTFBSclusters(ArrayList<ArrayList<Site>> clusters){
		if (print_full_format){
			StringBuilder sb = new StringBuilder();
			sb.append("#Region\tLength\t#Sites\tTFs\tTFIDs\tSignals\tPos\tMotifs\t#Motif\n");
			for (ArrayList<Site> c:clusters){
				int numSite = c.size();
				Region r = new Region(genome, c.get(0).bs.getChrom(), c.get(0).bs.getLocation(), c.get(numSite-1).bs.getLocation());
				StringBuilder sb_tfs = new StringBuilder();
				StringBuilder sb_tfids = new StringBuilder();
				StringBuilder sb_tf_signals = new StringBuilder();
				StringBuilder sb_tf_positions = new StringBuilder();
				StringBuilder sb_tf_motifs = new StringBuilder();
				int totalMotifs = 0;
				for (Site s:c){
					sb_tfs.append(names.get(s.tf_id)).append(",");
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
			for (ArrayList<Site> c:clusters){
				int numSite = c.size();
				Region r = new Region(genome, c.get(0).bs.getChrom(), c.get(0).bs.getLocation(), c.get(numSite-1).bs.getLocation());
				StringBuilder sb_tfs = new StringBuilder().append(r.toString()).append("\t").append(numSite).append("\t");
				for (Site s:c){
					sb_tfs.append(names.get(s.tf_id)).append(" ");
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
			int[] factorSiteCount = new int[expts.size()];

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
			int[] factorSiteCount = new int[expts.size()];

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
			double[] factorReadCount = new double[expts.size()];

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
			for (int i=0;i<names.size();i++)
				sb.append(names.get(i)).append("\n");
			CommonUtils.writeFile("0_BS_clusters."+outPrefix+".Dictioinary.txt", sb.toString());
		}
		
		if (print_matrix){	// Print region-tfSiteCount matrix, can be used for clustering analysis
			StringBuilder sb = new StringBuilder();
			int[] factorSiteCount = new int[expts.size()];
			sb.append("#Region    ").append("\t");
			for (int i=0;i<names.size();i++)
				sb.append(names.get(i)).append("\t");
			CommonUtils.replaceEnd(sb, '\n');
			
			for (ArrayList<Site> c :clusters){
				int numSite = c.size();
				Region r = new Region(genome, c.get(0).bs.getChrom(), c.get(0).bs.getLocation(), c.get(numSite-1).bs.getLocation());
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
			double[] factorReadCount = new double[expts.size()];
			sb.append("#Region    ").append("\t");
			for (int i=0;i<names.size();i++)
				sb.append(names.get(i)).append("\t");
			CommonUtils.replaceEnd(sb, '\n');
			
			for (ArrayList<Site> c :clusters){
				int numSite = c.size();
				Region r = new Region(genome, c.get(0).bs.getChrom(), c.get(0).bs.getLocation(), c.get(numSite-1).bs.getLocation());
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
		
		for (int id=0;id<clusters.size();id++){
			ArrayList<Site> bindingSites = clusters.get(id);
			int numSite = bindingSites.size();
			Region r = new Region(genome, bindingSites.get(0).bs.getChrom(), bindingSites.get(0).bs.getLocation(), bindingSites.get(numSite-1).bs.getLocation());	
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
			// scan motif matches in the cluster region sequence
			String seq = seqgen.execute(region).toUpperCase();	
			for (int i=0;i<pwms.size();i++){
				ArrayList<Integer> matchPos = CommonUtils.getAllPWMHit(seq, wmLens.get(i), scorers.get(i), wmThresholds.get(i));
				for (int p:matchPos){
					pwmMatchIds.add(i);
					pwmMatchPos.add(p);
				}
			}
			
			// K-mer matches
			ArrayList<Integer> kmerMatchIds = new ArrayList<Integer>();
			ArrayList<Integer> kmerMatchPos = new ArrayList<Integer>();
			for (int i=0;i<kmers.size();i++){
				String kmer = kmers.get(i);
				ArrayList<Integer> matchPos = CommonUtils.getAllKmerHit(seq, kmer);
				for (int p:matchPos){
					kmerMatchIds.add(i);
					kmerMatchPos.add(p);
				}
			}
			sb.append(String.format("%s\t%s\t%d\t%d\t%d\t%d\t%d\t", region.toString(), r.toString(),
					id, r.getWidth(), numSite, pwmMatchIds.size(), kmerMatchIds.size()));
			for (int i=0;i<bindingIds.size();i++)
				sb.append(String.format("B%d:%d:%d:%s:%.1f ", bindingIds.get(i), bindingPos.get(i), eventIds.get(i), strands.get(i), bindingStrength.get(i)));
			for (int i=0;i<pwmMatchIds.size();i++)
				sb.append(String.format("M%d:%d ", pwmMatchIds.get(i), pwmMatchPos.get(i)));
			for (int i=0;i<kmerMatchIds.size();i++)
				sb.append(String.format("K%d:%d ", kmerMatchIds.get(i), kmerMatchPos.get(i)));
			sb.append("\t").append(seq);
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
		ArrayList<String> lines = CommonUtils.readTextFile(cluster_file);
		ArrayList<BMCluster> clusters = new ArrayList<BMCluster>();
		TreeMap<String, TreeMap<String, SpacingProfile>> profiles = new TreeMap<String, TreeMap<String, SpacingProfile>>();
		for (String l:lines){	// each line is a cluster merged from nearby TFBS
			if (l.startsWith("#"))
				continue;
			
			BMCluster bmc = new BMCluster();
			
			String[] f = l.split("\t");
			Region r = Region.fromString(genome, f[0]);
			bmc.cluster_id = Integer.parseInt(f[2]);
			String[] positions = f[7].split(" ");
			bmc.anchoredSequence = f[8];		//TODO: use strand and shift
			bmc.anchor = new Point(r.getGenome(), r.getChrom(), r.getStart());
			
			// parse all sites
			ArrayList<RSite> sites = new ArrayList<RSite>();
			bmc.sites = sites;
			for (String ps: positions){
				String[] bs = ps.split(":");
				RSite s = new RSite();
				s.id = Integer.parseInt(bs[0].substring(1));
				s.type = bs[0].charAt(0);
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
			for (RSite anchor: sites){
				String ancStr = anchor.type+""+anchor.id;
				if (!profiles.containsKey(ancStr))
					profiles.put(ancStr, new TreeMap<String,SpacingProfile>());
				TreeMap<String,SpacingProfile> profiles_anchor = profiles.get(ancStr);
				for (RSite target: sites){
					String tarStr = target.type+""+target.id;
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
				}	// each site as target
			} // each site as anchor
			clusters.add(bmc);
		}// for each line
					
		StringBuilder sb_count = new StringBuilder("Anchor"+"\t"); 		// the count of tallest bar of spacings
		StringBuilder sb_offset = new StringBuilder("Anchor"+"\t");		// the offset of tallest bar of spacings
		
		for (String ancStr: profiles.keySet()){
			sb_count.append(ancStr+"\t");
			sb_offset.append(ancStr+"\t");
		}
		CommonUtils.replaceEnd(sb_count, '\n');
		CommonUtils.replaceEnd(sb_offset, '\n');
		
		for (String ancStr: profiles.keySet()){
			
			StringBuilder sb_profiles = new StringBuilder(ancStr+"\t");			
			for (int i=-PROFILE_RANGE;i<=PROFILE_RANGE;i++)
				sb_profiles.append(i+"\t");
			CommonUtils.replaceEnd(sb_profiles, '\n');
			
			TreeMap<String,SpacingProfile> profiles_anchor = profiles.get(ancStr);
			sb_count.append(ancStr+"\t");
			sb_offset.append(ancStr+"\t");
			for (String tarStr: profiles.keySet()){
				SpacingProfile pf = profiles_anchor.get(tarStr);
				int[] tmp=null;
				if (ancStr.equals(tarStr)){
					tmp = pf.profile_same.clone();
					tmp[100]=0;
				}
				else
					tmp = pf.profile_same;
				Pair<Integer, TreeSet<Integer>> max_same = StatUtil.findMax(tmp);
				if (ancStr.equals(tarStr)){
					tmp = pf.profile_diff.clone();
					tmp[100]=0;
				}
				else
					tmp = pf.profile_diff;
				Pair<Integer, TreeSet<Integer>> max_diff = StatUtil.findMax(tmp);
				if (ancStr.equals(tarStr)){
					tmp = pf.profile_unknown.clone();
					tmp[100]=0;
				}
				else
					tmp = pf.profile_unknown;
				Pair<Integer, TreeSet<Integer>> max_unknown = StatUtil.findMax(tmp);
				
				Pair<Integer, TreeSet<Integer>> max_all=null;
				if (max_same.car()>=max_diff.car())
					max_all = max_same;
				else 
					max_all = max_diff;
				if (max_all.car()<max_diff.car())
					max_all = max_unknown;
				sb_count.append(max_all.car()+"\t");
				sb_offset.append((max_all.cdr().first()-PROFILE_RANGE)+"\t");
				
				sb_profiles.append(tarStr+"_s\t").append(CommonUtils.arrayToString(pf.profile_same)).append("\n");
				sb_profiles.append(tarStr+"_d\t").append(CommonUtils.arrayToString(pf.profile_diff)).append("\n");
				sb_profiles.append(tarStr+"_u\t").append(CommonUtils.arrayToString(pf.profile_unknown)).append("\n");
			}
			CommonUtils.replaceEnd(sb_count, '\n');	
			CommonUtils.replaceEnd(sb_offset, '\n');	
			
			CommonUtils.writeFile(ancStr+"_profiles.txt", sb_profiles.toString());
		}
		System.out.println("The following is similar to the pairwise spacing matrix in GEM paper.");
		System.out.println(sb_count.toString());
		System.out.println(sb_offset.toString());
	}

	/**
	 * Plot the aligned binding sites given the anchor and sorting site label<br>
	 * output matlab data file, anchor site coordinates and corresponding sequence plot
	 */
	private void plotAlignedSites(){
		String fss[] = sort_string.split(":");
		String sortby = fss[0];
		int minOffset = Integer.parseInt(fss[1]);
		int maxOffset = Integer.parseInt(fss[2]);
		String analysisName = "_a"+anchor_string +"_s"+ sortby +"_";	
		System.out.println("Anchor: "+anchor_string+"  Sortby: "+sortby);

		ArrayList<String> lines = CommonUtils.readTextFile(cluster_file);
		ArrayList<BMCluster> clusters = new ArrayList<BMCluster>();
		SpacingProfile pf = new SpacingProfile();
		for (String l:lines){	// each line is a cluster merged from nearby TFBS
			if (l.startsWith("#"))
				continue;
			
			BMCluster bmc = new BMCluster();
			
			String[] f = l.split("\t");
			Region r = Region.fromString(genome, f[0]);
			bmc.cluster_id = Integer.parseInt(f[2]);
			String[] positions = f[7].split(" ");
			bmc.anchoredSequence = f[8];		//TODO: use strand and shift
			bmc.anchor = new Point(r.getGenome(), r.getChrom(), r.getStart());
			
			// parse all sites
			ArrayList<RSite> sites = new ArrayList<RSite>();
			bmc.sites = sites;
			for (String ps: positions){
				String[] bs = ps.split(":");
				RSite s = new RSite();
				s.id = Integer.parseInt(bs[0].substring(1));
				s.type = bs[0].charAt(0);
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
			
			for (RSite anchor: sites){
				String ancStr = anchor.type+""+anchor.id;
				if (!anchor_string.equals(ancStr))
					continue;
				for (RSite target: sites){
					String tarStr = target.type+""+target.id;
					if (!sort_string.equals(tarStr))
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
			ArrayList<RSite> sortbySites = new ArrayList<RSite>();
			for (RSite s: c.sites){
				if (anchor_string.equals(s.type+""+s.id))
					anchorSites.add(s);
				if (sortby.equals(s.type+""+s.id))
					sortbySites.add(s);
			}
			int bestCount=-1;
			RSite bestAnchor = null;
			int bestOffset = 0;
			for (RSite as: anchorSites){
				for(RSite ss: sortbySites){
					Pair<Character, Integer> p = getProfileIndex(as, ss);
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
						bestOffset = p.cdr()-PROFILE_RANGE;
					}					
				}
			}
			if (bestAnchor==null)// Skip if no anchor binding site is found 
				continue;
			else{
				c.updateAnchorSite(bestAnchor);
				c.sort_idx = bestOffset;

				if (bestOffset>=minOffset && bestOffset<=maxOffset)
					outClusters.add(c);
			}
//			System.out.println(String.format("Count:\t%d\tOffset: %d\tPos:%d", bestCount, c.sort_idx, bestAnchor.pos));

		}// each cluster
		
		if (outClusters.isEmpty()){
			System.err.println("Anchor "+anchor_string+" matched no sites!");
			System.exit(-1);
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
		CommonUtils.writeFile(cluster_file.replace(".txt", analysisName+"matlab.txt"), sb.toString());
		CommonUtils.writeFile(cluster_file.replace(".txt", analysisName+"coords.txt"), sb_coords.toString());
		
		String[] ss = new String[outClusters.size()];
		for (int i=0;i<outClusters.size();i++){
			BMCluster c = outClusters.get(i);
			int len = c.anchoredSequence.length();
			int left = c.anchorSite.pos + minOffset;
			int leftPadding = 0;
			if (left<0){
				leftPadding = -left;
				left=0;
			}
			int right = c.anchorSite.pos + maxOffset+1;  // add 1 because substring() is end-exclusive
			int rightPadding = 0;
			if (right>len){
				rightPadding = right-len;
				right=len;
			}
			String seq = CommonUtils.padding(leftPadding, 'N')+c.anchoredSequence.substring(left, right)+CommonUtils.padding(rightPadding, 'N');
			if (c.anchorSite.strand=='-')
				seq = SequenceUtils.reverseComplement(seq);
			ss[i] = seq;
		}
		CommonUtils.visualizeSequences(ss, width, height, new File(cluster_file.replace(".txt", analysisName+"sequences.png")));
	}


	/** compute the relative index of target site given the anchor site as 0 position
	 * 
	 */
	private Pair<Character, Integer> getProfileIndex(RSite anchor, RSite target){
		int offset = target.pos - anchor.pos;
		char type = 'u';
		if (offset>PROFILE_RANGE || offset <-PROFILE_RANGE)			// skip if out of range
			return null;
		if(anchor.strand=='-'){
			offset = -offset;
			offset += PROFILE_RANGE;		// shift to get positive array idx
			if (target.strand=='+')
				type = 'd';
			if (target.strand=='-')
				type = 's';
		}
		if(anchor.strand=='+'){
			offset += PROFILE_RANGE;		// shift to get positive array idx
			if (target.strand=='+')
				type = 's';
			if (target.strand=='-')
				type = 'd';
		}
		if(anchor.strand=='*')
			offset += PROFILE_RANGE;
		return new Pair<Character, Integer>(type, offset);
	}
	
	/** 
	 * RSite: relative site<br>
	 * The position of the site is relative to an anchor point (stranded)
	 */
	private class RSite implements Cloneable{
		char type = ' ';		// type=0:binding, 1:pwm, 2:kmer
		int id=-1;			// id in that type
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
		int[] profile_same = new int[PROFILE_RANGE*2+1];		// the anchor site and the subject site are on the same strand (or same orientation of motifs)
		int[] profile_diff = new int[PROFILE_RANGE*2+1];		// the anchor site and the subject site are on the opposite strand (or opposite orientation of motifs)
		int[] profile_unknown = new int[PROFILE_RANGE*2+1];	// unknown, because some GEM binding site does not have a k-mer to assgin strand
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
            ReadCache ipCache = new ReadCache(genome, expts.get(i));
            
			// cache sorted start positions and counts of all positions
			long tic = System.currentTimeMillis();
			System.err.print("Loading "+ipCache.getName()+" data from ReadDB ... \t");
			List<String> chroms = genome.getChromList();
			if (dev){
				chroms = new ArrayList<String>();
				chroms.add("19");
			}
			for (String chrom: chroms ){
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
			ipCache.displayStats();
			System.out.println(CommonUtils.timeElapsed(tic));
            
			// now get the data from the cache
            for (int j=0;j<all_TSS.size();j++){
            	Region region = all_TSS.get(j).expand(TARGET_WIDTH);
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
				sb1.append(names.get(id)).append("\t");
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
				sb1.append(names.get(id)).append("\t");
			System.out.println("=================================\nchr"+c.region.toString());
			System.out.println(sb1.toString());
			
			for (Point p:targets){				
				// get corresponding signals at the target sites
				List<Double> target_signals = new ArrayList<Double>();
				for (int i=0;i<TFIDs.length;i++){
					target_signals.add(i, (double)chipseqs.get(TFIDs[i]).countHits(p.expand(TARGET_WIDTH)));
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

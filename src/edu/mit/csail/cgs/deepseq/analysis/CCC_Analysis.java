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
import java.util.regex.Matcher;
import java.util.regex.Pattern;

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
import edu.mit.csail.cgs.deepseq.analysis.BindingSpacing_GeneStructure.Site;
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
import edu.mit.csail.cgs.utils.stats.StatUtil;
/**
 * Compute the spatial relationship among multiple TFs' binding sites<br>
 * Find the clusters of multiTF binding and the relative positions of each TF sites
 * @author Yuchun
 *
 */
public class CCC_Analysis {
	private final int TARGET_WIDTH = 250;
	Genome genome=null;
	String[] args;
	ArrayList<String> expts = new ArrayList<String>();
	ArrayList<String> names = new ArrayList<String>();
	ArrayList<String> readdb_names = new ArrayList<String>();
	
	ArrayList<WeightMatrix> pwms = new ArrayList<WeightMatrix>();
	ArrayList<ArrayList<Tss>> all_sites = new ArrayList<ArrayList<Tss>>();
	ArrayList<Point> all_TSS;
	ArrayList<Cluster> all_clusters;
	int[][]tss_signals = null;
	
	int tss_radias = 2000;		// the range around TSS site to search for targets
	int region_padding = 100;	// padding to add to the region, set as ChIP-Seq event call spatial resolution (100bp)
	int exclude_range = 500;		// the range around anchor site to exclude same site
	int min_site = 1;				// the minimum number of sites for the cluster to be printed out
	double wm_factor = 0.6;	// PWM threshold, as fraction of max score
	double cutoff = 0.3;	// corr score cutoff
	File dir;
	private SequenceGenerator<Region> seqgen;
	boolean dev = false;
	boolean zero_or_one = false;	// for each TF, zero or one site per cluster, no multiple sites
	String outPrefix = "out";
	String tss_file;
	String tss_signal_file;
	String cluster_file;
	String exclude_sites_file;

	boolean useDirectBindingOnly = false;
	String gisFile;
	String umassFile;
	String germFile;
	String gencodeFile;
	
	private ArrayList<Interaction> interactions;
	private ArrayList<Tss> tss;
	// command line option:  (the folder contains GEM result folders) 
	// --dir C:\Data\workspace\gse\TFBS_clusters --species "Mus musculus;mm9" --r 2 --pwm_factor 0.6 --expts expt_list.txt [--no_cache --old_format] 
	public static void main(String[] args) {
////		String GIS_PAIR_RE = "chr\\w:\\d\\.\\.\\d\\-chr\\w:\\d\\.\\.\\d,\\d";
//		String GIS_PAIR_RE = "chr\\w:\\d..\\d";
//		Pattern GIS_PATTERN = Pattern.compile(GIS_PAIR_RE);
////		Matcher regionmatcher = GIS_PATTERN.matcher("chr1:3267037..3270061-chr17:41463804..41467791,12");
//		Matcher regionmatcher = GIS_PATTERN.matcher("chr1:3267037..3270061");
//		if (regionmatcher.find()) {
//			System.out.println("Found");
//			System.out.println(regionmatcher.groupCount());
//			if (regionmatcher.groupCount() != 7) {
//				return;
//			}
//		}
		
		CCC_Analysis mtb = new CCC_Analysis(args);
		int type = Args.parseInteger(args, "type", 0);
		switch(type){
		case 0:		// map regulatory regions enriched in topics to target genes (listed by each topic)
			mtb.printTopic2genes();
			break;
		case 1:		// map regulatory regions to target genes (listed by each region)
			mtb.printGenesByRegions();
			break;
		}		
	}
	
	public CCC_Analysis(String[] args){
				
	    try {
	    	Pair<Organism, Genome> pair = Args.parseGenome(args);
	    	if(pair==null){
	    	  System.err.println("No genome provided; provide a Gifford lab DB genome name");
	    	  System.exit(1);
	    	}else{
	    		genome = pair.cdr();
	    	}
	    } catch (NotFoundException e) {
	      e.printStackTrace();
	    }
	    this.args = args;
	    
		Set<String> flags = Args.parseFlags(args);
		outPrefix = Args.parseString(args, "out", outPrefix);
		
		gisFile = Args.parseString(args, "gis", gisFile);
		if (gisFile!=null)
			interactions = loadInteractions(gisFile, "gis");
		else{
			umassFile = Args.parseString(args, "umass", umassFile);
			germFile = Args.parseString(args, "germ", germFile);
		}
		
		gencodeFile = Args.parseString(args, "gencode", gencodeFile);
		
		tss_signal_file = Args.parseString(args, "tss_signal", null);
		cluster_file = Args.parseString(args, "cluster", null);		
		
		
		zero_or_one = flags.contains("zoo");
		dev = flags.contains("dev");
		dir = new File(Args.parseString(args, "dir", "."));
		expts = new ArrayList<String>();
		names = new ArrayList<String>();

		tss_file = Args.parseString(args, "tss", null);

		exclude_sites_file = Args.parseString(args, "ex", null);
		tss_radias = Args.parseInteger(args, "tss_radias", tss_radias);
		region_padding = Args.parseInteger(args, "region_padding", region_padding);
		exclude_range = Args.parseInteger(args, "exclude", exclude_range);
		min_site = Args.parseInteger(args, "min_site", min_site);
		wm_factor = Args.parseDouble(args, "pwm_factor", wm_factor);
		cutoff = Args.parseDouble(args, "cutoff", cutoff);
		seqgen = new SequenceGenerator<Region>();
		seqgen.useCache(!flags.contains("no_cache"));
	}
	
	private ArrayList<Interaction> loadInteractions(String file, String format){
		ArrayList<Interaction> data = new ArrayList<Interaction>();
		ArrayList<String> txt = CommonUtils.readTextFile(file);
		String s = null;
		for (int i=0;i<txt.size();i++){
			if (txt.get(i).startsWith("#"))
				continue;
			if (format.equals("gis")){
				String[] f = txt.get(i).split("\t");
				if (f[3].equals(s))		// skip if it is the second of a pair of lines
					continue;
				s = f[3];		// chr1:3267037..3270061-chr17:41463804..41467791,12				
				Interaction it = new Interaction();
				String[] f2 = s.split(",");
				it.strength = Integer.parseInt(f2[1]);
				String[] f3 = f2[0].split("-");	
				String[] f4 = f3[0].split(":");
				String[] f5 = f4[1].split("\\.\\.");	// . need to be escaped in regExpr
				it.anchor1 = new Region(genome, f4[0].replace("chr", ""), Integer.parseInt(f5[0]), Integer.parseInt(f5[1]));
				String[] f6 = f3[1].split(":");
				String[] f7 = f6[1].split("\\.\\.");
				it.anchor2 = new Region(genome, f6[0].replace("chr", ""), Integer.parseInt(f7[0]), Integer.parseInt(f7[1]));
				data.add(it);
			}			
		}
		data.trimToSize();
		return data;
	}
	
	private void loadTSS(){
		System.out.println("Loading GENCODE annotations ... ");
		ArrayList<String> texts = CommonUtils.readTextFile(gencodeFile);
		tss = new ArrayList<Tss>();					// -1
		for (int i=0;i<texts.size();i++){
			String t = texts.get(i);
			if (t.startsWith("#"))
				continue;
			String f[] = t.split("\t");
			
			// UCSC ucsc_hgTables.txt
			String chr = f[1].replace("chr", "");
			char strand = f[2].charAt(0);
			tss.add(new Tss(i, new Point(genome, chr, Integer.parseInt(f[strand=='+'?3:4])), strand, f[10])); // Gene Symbol
			
			// GENCODE v17
//			String chr = f[2].replace("chr", "");
//			char strand = f[3].charAt(0);
//			tss.add(new Tss(i, new Point(genome, chr, Integer.parseInt(f[strand=='+'?4:5])), strand, f[12])); // Gene Symbol
		}
		tss.trimToSize();
		Collections.sort(tss);
	}
	
	private void printTopic2genes(){
		// load all regoins
		ArrayList<Region> rs = CommonUtils.loadCgsRegionFile(Args.parseString(args, "regions", null), genome);
		
		// loadTSS
		loadTSS();
		ArrayList<Region> tssRegions = new ArrayList<Region>();
		for (Tss t: tss){	// tss has been sorted when it is loaded
			tssRegions.add(t.point.expand(tss_radias));
		}
		tssRegions.trimToSize(); 
		Collections.sort(tssRegions, new Comparator<Region>() {
            public int compare(Region o1, Region o2) {
                return o1.compareToStrict(o2);
            }
        });
		
		// Load HDP assignment file to match topic to regions
//		d w z t
//		131779 30 11 1
//		131779 104 11 0
		ArrayList<String> texts = CommonUtils.readTextFile(Args.parseString(args, "hdp_assignments", null));
		TreeMap<Integer, HashSet<Integer>> topic2region = new TreeMap<Integer, HashSet<Integer>>();
		
		// associate topic to region
		TreeMap<Integer, ArrayList<Integer>> region2topics = new TreeMap<Integer, ArrayList<Integer>>();
		for (String t: texts){
			if (t.startsWith("d"))
				continue;
			String[] f = t.split(" ");
			int topicId = Integer.parseInt(f[2]);
			int regionId = Integer.parseInt(f[0]);
			if (region2topics.containsKey(regionId))
				region2topics.get(regionId).add(topicId);
			else{
				ArrayList<Integer> tids = new ArrayList<Integer>();
				tids.add(topicId);
				region2topics.put(regionId, tids);
			}
		}
		
		for (int regionId:region2topics.keySet()){
			ArrayList<Integer> tids = region2topics.get(regionId);
			Pair<int[], int[]> sorted = StatUtil.sortByOccurences(tids); 
			int[] elements = sorted.car();
			int[] counts = sorted.cdr();
			double meanPlusStd = StatUtil.mean(counts)+StatUtil.std(counts);
//			int maxCount = counts[counts.length-1];
			for (int i=counts.length-1;i>=0;i--){
				if (counts[i]<=meanPlusStd)			// only count highly-used topics
					break;
				int topicId = elements[i];
				if (topic2region.containsKey(topicId))
					topic2region.get(topicId).add(regionId);
				else{
					HashSet<Integer> rids = new HashSet<Integer>();
					rids.add(regionId);
					topic2region.put(topicId, rids);
				}
			}
		}
		
		StringBuilder sb = new StringBuilder();
		for (int topicId:topic2region.keySet()){
			HashSet<Integer> rids = topic2region.get(topicId);
			ArrayList<Region> topicRegions = new ArrayList<Region>();
			for (int rid:rids)
				topicRegions.add(rs.get(rid));
			sb.append(mapRegions2genes(topicRegions, tssRegions, topicId));
		}
		CommonUtils.writeFile(outPrefix+"_tfbs2genes_byTopic.txt", sb.toString());
	}
	
	private void printGenesByRegions(){
		// load all regoins
		ArrayList<Region> rs = CommonUtils.loadCgsRegionFile(Args.parseString(args, "regions", null), genome);
		// loadTSS
		loadTSS();
		ArrayList<Region> tssRegions = new ArrayList<Region>();
		for (Tss t: tss){	// tss has been sorted when it is loaded
			tssRegions.add(t.point.expand(tss_radias));
		}
		tssRegions.trimToSize(); 
		Collections.sort(tssRegions, new Comparator<Region>() {
            public int compare(Region o1, Region o2) {
                return o1.compareToStrict(o2);
            }
        });		
		String s = mapRegions2genes(rs, tssRegions, -1);
		CommonUtils.writeFile(outPrefix+"_tfbs2genes_byRegion.txt", s);
	}
	
	
	private String mapRegions2genes(ArrayList<Region> unexpandedRegions, ArrayList<Region> tssRegions, int topicId){
		ArrayList<Region> regions = new ArrayList<Region>(unexpandedRegions.size());
		for (Region r:unexpandedRegions)
			regions.add(r.expand(region_padding, region_padding));		
		
		// Map regions to overlapping TSS regions
		ArrayList<TreeSet<String>> tssGenes = new ArrayList<TreeSet<String>>();
		for (Region r:regions)
			tssGenes.add(null);
		tssGenes.trimToSize();
		
		for (int j=0;j<regions.size();j++){
			Region r = regions.get(j);
			int idx = Collections.binarySearch(tssRegions, r, new Comparator<Region>() {
	            public int compare(Region o1, Region o2) {
	                return o1.compareToStrict(o2);
	            }
	        });
			if (idx<0){	// not match, point idx to the next index
				idx = Math.min(-idx - 1, tssRegions.size()-1); 
			}
			TreeSet<String> set = new TreeSet<String>();
			for (int i=idx;i<tssRegions.size();i++){
				if (tssRegions.get(i).overlaps(r)){
					set.add(tss.get(i).geneSymbol);					
				}
				else
					break;
			}
			if (!set.isEmpty())
				tssGenes.add(j, set);
		}
		
		// Map regions to ChIA-PET linked (two-way) TSS?
		ArrayList<TreeSet<String>> chiapetGenes = new ArrayList<TreeSet<String>>();
		for (Region r:regions)
			chiapetGenes.add(null);
		chiapetGenes.trimToSize();
		
		for (Interaction ir: interactions){
			Region anchor1 = ir.anchor1;
			Region anchor2 = ir.anchor2;
			// one anchor -- TSS of gene
			ArrayList<Integer> idx1 = getOverlapRegionIndices(anchor1, tssRegions);
			TreeSet<String> geneSymbols = new TreeSet<String>();
			for (int t: idx1){
				geneSymbols.add(tss.get(t).geneSymbol);		
			}
			// the other anchor -- binding region
			ArrayList<Integer> idx2 = getOverlapRegionIndices(anchor2, regions);
			for (int r: idx2){
				TreeSet<String> storedGeneSymbols = chiapetGenes.get(r);
				if (storedGeneSymbols==null)
					chiapetGenes.set(r, geneSymbols);
				else
					storedGeneSymbols.addAll(geneSymbols);
			}
			
			idx1 = getOverlapRegionIndices(anchor2, tssRegions);
			geneSymbols.clear();
			for (int t: idx1){
				geneSymbols.add(tss.get(t).geneSymbol);		
			}
			idx2 = getOverlapRegionIndices(anchor2, regions);
			for (int r: idx2){
				TreeSet<String> storedGeneSymbols = chiapetGenes.get(r);
				if (storedGeneSymbols==null)
					chiapetGenes.set(r, geneSymbols);
				else
					storedGeneSymbols.addAll(geneSymbols);
			}
		}
		
		// output results
		if (topicId==-1){	// return the genes <-- by regions
			StringBuilder sb = new StringBuilder();
			sb.append("#Region\tTSS_Genes\tChIA-PET_Genes\tAll_Genes\n");
			for (int i=0;i<unexpandedRegions.size();i++){
				Region r = unexpandedRegions.get(i);
				TreeSet<String> tssGeneSymbols = tssGenes.get(i);
				TreeSet<String> cpGeneSymbols = chiapetGenes.get(i);
				TreeSet<String> allGeneSymbols = new TreeSet<String>();
				StringBuilder gsb = new StringBuilder();
				if (tssGeneSymbols!=null){	
					allGeneSymbols.addAll(tssGeneSymbols);
					for (String g:tssGeneSymbols)
						gsb.append(g).append(" ");
				}
				StringBuilder cpsb = new StringBuilder();
				if (cpGeneSymbols!=null){
					allGeneSymbols.addAll(cpGeneSymbols);
					for (String g:cpGeneSymbols)
						cpsb.append(g).append(" ");
				}
				StringBuilder allsb = new StringBuilder();
				for (String g:allGeneSymbols) 
					allsb.append(g).append(" ");
				sb.append(String.format("%s\t%s\t%s\t%s\n", r.toString(), gsb.toString(), cpsb.toString(), allsb.toString()));
			}
			return sb.toString();
		}
		else{		// return the genes <-- regions <-- by topic
			StringBuilder sb = new StringBuilder();
			TreeSet<String> tssGeneSymbols = new TreeSet<String>();
			TreeSet<String> cpGeneSymbols = new TreeSet<String>();
			TreeSet<String> allGeneSymbols = new TreeSet<String>();
			for (int i=0;i<unexpandedRegions.size();i++){
				TreeSet<String> tsg = tssGenes.get(i);
				if (tsg!=null)
					tssGeneSymbols.addAll(tsg);
				TreeSet<String> cpg = chiapetGenes.get(i);
				if (cpg!=null)
					cpGeneSymbols.addAll(cpg);			
			}
			allGeneSymbols.addAll(tssGeneSymbols);
			allGeneSymbols.addAll(cpGeneSymbols);
			
			StringBuilder gsb = new StringBuilder();
			if (tssGeneSymbols!=null){	
				for (String g:tssGeneSymbols)
					gsb.append(g).append(" ");
			}
			StringBuilder cpsb = new StringBuilder();
			if (cpGeneSymbols!=null){
				for (String g:cpGeneSymbols)
					cpsb.append(g).append(" ");
			}
			StringBuilder allsb = new StringBuilder();
			if (allGeneSymbols!=null){
				for (String g:allGeneSymbols) 
					allsb.append(g).append(" ");
			}
			sb.append(String.format("%s\t%s\n", (topicId+1)+"_TSS_genes", gsb.toString()));
			sb.append(String.format("%s\t%s\n", (topicId+1)+"_ChIA-PET_genes", cpsb.toString()));
			sb.append(String.format("%s\t%s\n", (topicId+1)+"_ALL_genes", allsb.toString()));
			
			return sb.toString();
		}
	}
	
	private ArrayList<Integer> getOverlapRegionIndices(Region source, ArrayList<Region> targets){
		int idx = Collections.binarySearch(targets, source, new Comparator<Region>() {
            public int compare(Region o1, Region o2) {
                return o1.compareToStrict(o2);
            }
        });
		if (idx<0){	// not match, point idx to the next index
			idx = Math.min(-idx - 1, targets.size()-1); 
		}
		ArrayList<Integer> results = new ArrayList<Integer>();
		for (int i=idx;i<targets.size();i++){
			if (targets.get(i).overlaps(source)){
				results.add(i);					
			}
			else
				break;
		}
		return results;
	}
	
	class Interaction{
		Region anchor1;
		Region anchor2;
		double strength;
		double pvalue;
		
		public int compareByanchor1(Interaction i) {
			return anchor1.compareTo(i.anchor1);
		}
		public int compareByanchor2(Interaction i) {
			return anchor2.compareTo(i.anchor2);
		}
		
	}
	
	class Tss implements Comparable<Tss>{
		Point point;
		int id;
		char strand;
		String geneSymbol;
		public Tss(int id, Point point, char strand, String geneSymbol) {
			this.id = id;
			this.point = point;
			this.strand = strand;
			this.geneSymbol = geneSymbol;
			// TODO Auto-generated constructor stub
		}
		public int compareTo(Tss tss) {					// ascending coordinate
			return(point.compareTo(tss.point));
		}
	}

	class Cluster{
		Region region;
		ArrayList<Integer> TFIDs = new ArrayList<Integer>();
		ArrayList<Double> TF_Signals = new ArrayList<Double>();
		ArrayList<Boolean> TF_hasMotifs = new ArrayList<Boolean>();
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
			ArrayList<Point> targets = CommonUtils.getPointsWithinWindow(all_TSS, anchor, tss_radias);
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
			ArrayList<Integer> targets = CommonUtils.getPointsIdxWithinWindow(all_TSS, anchor.expand(tss_radias));
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
			
			for (int idx:targets){				
				// get corresponding signals at the target sites
				Point p = all_TSS.get(idx);
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

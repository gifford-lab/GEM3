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
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.utils.stats.StatUtil;
/**
 * Compute the spatial relationship among multiple TFs' binding sites<br>
 * Find the clusters of multiTF binding and the relative positions of each TF sites
 * @author Yuchun
 *
 */
public class RegionAnnotator {
	String[] args = null;
	Genome genome=null;
	int padding = 0;  // the padding distance added to the region range
	String outPrefix = "out";
	String region_file;			
	String annotation_file;
	// --species "Homo sapiens;hg19"  --region_file region_test.txt --anno_file anno_encode_regions.txt --out k562_hdp
	public static void main(String[] args) {
		RegionAnnotator ra = new RegionAnnotator(args);
		int type = Args.parseInteger(args, "type", 1);
		switch(type){
		case 0:  	// print out the overlap percentage of each query region for specified annotated regions
			ra.overlap_regions_with_annotations();
			break;
		case 1:		// print out the sorted and aligned anchor coordinates for making line plots
			ra.sort_and_anchor_regions();
			break;
		case 2:
			ra.assign_gene_by_proximity();
			break;
		case 3: 	// print all pairwise distances between distal enhancers of each gene
			ra.compute_enhancer_distances();
			break;
		}
		
	}
	
	public RegionAnnotator(String[] args){
		this.args = args;		
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
		outPrefix = Args.parseString(args, "out", outPrefix);
	}
	
	/**
	 * Overlap regions with a list of annotations, report the length of overlaps
	 */
	private void overlap_regions_with_annotations(){
		ArrayList<Region> queryRegions = CommonUtils.loadCgsRegionFile(Args.parseString(args, "region_file", region_file), genome);
		String anno_file = Args.parseString(args, "anno_file", null);
		if (anno_file!=null){
			ArrayList<String> ids = new ArrayList<String>();
			ArrayList<ArrayList<Integer>> stats = new ArrayList<ArrayList<Integer>>();
			ArrayList<String> lines = CommonUtils.readTextFile(anno_file);
			for (String line: lines){
				String[] f = line.split("\t");
				String id=f[0];
				String name=f[1];
				String type=f[2];
				if (type.equalsIgnoreCase("SS")){	// single state
					ArrayList<Region> annotatedRegions = CommonUtils.load_BED_regions(genome, name).car();
					System.out.println(id);
					ids.add(id);
					ArrayList<Integer> stat = Region.computeBatchOverlapLength(queryRegions, annotatedRegions);
					stats.add(stat);
				}
				if (type.equalsIgnoreCase("MS")){	// multiple states, e.g. HMM segmentation
					Pair<ArrayList<Region>, ArrayList<String>> beds = CommonUtils.load_BED_regions(genome, name);
					ArrayList<Region> regions = beds.car();
					ArrayList<String> annos = beds.cdr();
					TreeSet<String> labels = new TreeSet<String> ();
					for (String l:annos)
						labels.add(l);
					for (String l:labels){
						System.out.println(id+"-"+l);
						ids.add(id+"-"+l);
						ArrayList<Region> annotatedRegions = new ArrayList<Region>();
						for (int i=0;i<annos.size();i++){
							if (annos.get(i).equals(l))
								annotatedRegions.add(regions.get(i));
						}
						Collections.sort(annotatedRegions);
						ArrayList<Integer> stat = Region.computeBatchOverlapLength(queryRegions, annotatedRegions);						
						stats.add(stat);
					}
				}
			}
			
			//output
			StringBuilder sb = new StringBuilder();
			sb.append("#Region\tLength\t");
			for (String id:ids){
				sb.append(id+"\t");
			}
			CommonUtils.replaceEnd(sb, '\n');
			for (int i=0;i<queryRegions.size();i++){
				sb.append(queryRegions.get(i).toString()).append("\t");
				sb.append(queryRegions.get(i).getWidth()).append("\t");
				for (ArrayList<Integer> stat: stats){
					sb.append(String.format("%d\t",stat.get(i)));
				}
				CommonUtils.replaceEnd(sb, '\n');
			}
			CommonUtils.writeFile(outPrefix+"_anno_stats.txt", sb.toString());
		}
	}
	
	/** sort the regions based on the specified data signal,  optionally, use GEM file of a TF to anchor the regions. 
	 * Useful for sorting regions when making line plot of reads
	 * If a region does not contained TF binding site, 
	* either skip that region, or just keep the same region middle point
	 * 
	 */
	private void sort_and_anchor_regions(){
		String regionFileFormat = Args.parseString(args, "rf", "CGS");  
		ArrayList<Region> queryRegions = null;
		if (regionFileFormat.equalsIgnoreCase("BED")){
			Pair<ArrayList<Region>, ArrayList<String>> pair = CommonUtils.load_BED_regions(genome, Args.parseString(args, "region_file", region_file));
			queryRegions=pair.car();
		}
		else
			queryRegions=CommonUtils.loadCgsRegionFile(Args.parseString(args, "region_file", region_file), genome);

		String anchor_GEM_file = Args.parseString(args, "anchor_GEM_file", null); // TF GEM file to anchor the region
		if (anchor_GEM_file !=null){
			List<GPSPeak> gpsPeaks = null;
			ArrayList<Point> sites = new ArrayList<Point>();
			try{
				gpsPeaks = GPSParser.parseGPSOutput(anchor_GEM_file, genome);
			}
			catch (IOException e){
				System.out.println(anchor_GEM_file+" is not a valid GPS/GEM event file.");
				System.exit(1);
			}
			Collections.sort(gpsPeaks);
			for (GPSPeak p: gpsPeaks){
				sites.add((Point)p);
			}
			Set<String> flags = Args.parseFlags(args);
			boolean discardRegion = flags.contains("discardRegion");
			ArrayList<Region> toRemove = new ArrayList<Region>();
			for (int i=0;i<queryRegions.size();i++){
				Region r = queryRegions.get(i);
				ArrayList<Integer> ids = CommonUtils.getPointsWithinWindow(sites,r);
				if (ids.isEmpty()){
					if(discardRegion){
						toRemove.add(r);
						continue;
					}
				}
				else{
					int selected_id = ids.get(0);
					if (ids.size()>1){
						double[] signals = new double[ids.size()];
						for (int j=0;j<ids.size();j++){
							int id =ids.get(j);
							signals[j]=gpsPeaks.get(id).getStrength();
						}
						Pair<Double, TreeSet<Integer>> max = StatUtil.findMax(signals);
						selected_id = ids.get(max.cdr().first());		// if tie, just use first site
					}
					queryRegions.set(i, gpsPeaks.get(selected_id).expand(r.getWidth()/2));  // set the midPoint to the TF location
				}
				
			}
			queryRegions.removeAll(toRemove);
		}

		int window = Args.parseInteger(args, "win", -1);
		
		DeepSeqExpt chipSeqExpt = null;
		String chipSeqFile = Args.parseString(args, "expt", null);
		if (chipSeqFile!=null){
			String fileFormat = Args.parseString(args, "f", "BED");  
			List<File> expts = new ArrayList<File>();
			expts.add(new File(chipSeqFile));
			chipSeqExpt = new DeepSeqExpt(genome, expts, false, fileFormat, -1);
		}
		else{
			List<ChipSeqLocator> rdbexpts = Args.parseChipSeq(args,"rdb");
			chipSeqExpt = new DeepSeqExpt(genome, rdbexpts, "readdb", -1);
		}
		int[] signals = new int[queryRegions.size()];
		for (int i=0;i<queryRegions.size();i++){
			signals[i] = chipSeqExpt.countHits(queryRegions.get(i).getMidpoint().expand(window/2));
		}
		int[] sortedIdx = StatUtil.findSort(signals);
		StringBuilder sb = new StringBuilder();
		for (int i=queryRegions.size()-1;i>=0;i--){			// descending order
			Region r = queryRegions.get(sortedIdx[i]);
			sb.append(r.getMidpoint().toString()).append("\t")
			.append(r.toString()).append("\t").append(r.getWidth()).append("\n");
		}
		CommonUtils.writeFile(outPrefix+".coords_regions.txt", sb.toString());
	}

	/**
	 * Given a list of coordinates, assign them to the nearest genes, optionally constrained with TAD annotations
	 */
	private void assign_gene_by_proximity(){
		if(Args.parseFlags(args).contains("help")){
			String msg = "assign_gene_by_proximity()\nUsage example:\n--species \"Mus musculus;mm10\"  --type 2 --genes mm10.refseq.ucsc_hgTables.txt --coords IL_k27ac.p300_sorted.coords.txt --tad TAD.mm10.bed\n"+
					"code for isInTAD, 0: nearest TSS, but no TSS found inside the TAD; 1: absolute nearest TSS, in the same TAD; 2: not absolute nearest, but nearest in TAD;  -9: coord is not in any TAD.";
			System.err.println(msg);
			return;
		}
		String coords_file = Args.parseString(args, "coords", null);
		ArrayList<Point> coords = CommonUtils.loadCgsPointFile(coords_file, genome);
		String tad_file = Args.parseString(args, "tad", null);
		String tad_name = Args.parseString(args, "tad_name", "no_TAD");
		ArrayList<Region> tads = new ArrayList<Region>();
		if (tad_file!=null){
			tads = CommonUtils.load_BED_regions(genome, tad_file).car();
			Collections.sort(tads);
		}
		HashMap<String, ArrayList<Region>> chr2tads = new HashMap<String, ArrayList<Region>>();
		for (Region r:tads){
			String chr = r.getChrom();
			if (!chr2tads.containsKey(chr))
				chr2tads.put(chr, new ArrayList<Region>());
			chr2tads.get(chr).add(r);
		}
		
		// load refSeq gene annotation
		int tssRange = Args.parseInteger(args, "tss_range", -1);
		ArrayList<String> texts = CommonUtils.readTextFile(Args.parseString(args, "genes", null));
		TreeMap<StrandedPoint, TreeSet<String>> tss2genes = new TreeMap<StrandedPoint, TreeSet<String>>();
		for (int i=0;i<texts.size();i++){
			String t = texts.get(i);
			if (t.startsWith("#"))
				continue;
			String f[] = t.split("\t");
			String chr = f[2].replace("chr", "");
			char strand = f[3].charAt(0);
			StrandedPoint tss = new StrandedPoint(genome, chr, Integer.parseInt(f[strand=='+'?4:5]), strand);
			String symbol = f[12];
			if (!tss2genes.containsKey(tss))
				tss2genes.put(tss, new TreeSet<String>());
			tss2genes.get(tss).add(symbol);
		}
		HashMap<String, ArrayList<StrandedPoint>> chr2tsss = new HashMap<String, ArrayList<StrandedPoint>>();
		for (StrandedPoint tss: tss2genes.keySet()){
			String chr = tss.getChrom();
			if (!chr2tsss.containsKey(chr))
				chr2tsss.put(chr, new ArrayList<StrandedPoint>());
			chr2tsss.get(chr).add(tss);			
		}
		// for each TSS, add nearby genes (within tssRange) into the TSS
		if (tssRange!=-1){
			for (String chr: chr2tsss.keySet()){
				ArrayList<StrandedPoint> tsss = chr2tsss.get(chr);
				for (int idx=0;idx<tsss.size();idx++){
					StrandedPoint tss = tsss.get(idx);
					for (int j=idx+1; j<tsss.size();j++){
						if (tss.distance(tsss.get(j))<=tssRange)
							tss2genes.get(tss).addAll(tss2genes.get(tsss.get(j)));
					}
					for (int j=idx-1; j>=0;j--){
						if (tss.distance(tsss.get(j))<=tssRange)
							tss2genes.get(tss).addAll(tss2genes.get(tsss.get(j)));
					}
				}
			}
		}
		
		StringBuilder sb = new StringBuilder();
		sb.append("#Region\tGene\tMidReg\tTSS\tdistance\tisInTAD\n");
		for (Point p:coords){			
			String chr = p.getChrom();
			ArrayList<StrandedPoint> tsss_in_chr = chr2tsss.get(chr);
			int idx_tss = Collections.binarySearch(tsss_in_chr, p);
			Point nearestTSS = null;
			if (idx_tss<0){
				idx_tss = -(idx_tss+1);
				// compare distance with points before and after the query point
				if (idx_tss==0)
					nearestTSS = tsss_in_chr.get(idx_tss);
				else{
					if (idx_tss<tsss_in_chr.size() &&
							p.distance(tsss_in_chr.get(idx_tss-1)) > p.distance(tsss_in_chr.get(idx_tss))){
						nearestTSS = tsss_in_chr.get(idx_tss);
					}
					else
						nearestTSS = tsss_in_chr.get(idx_tss-1);
				}
			}
			else // exactly on TSS
				nearestTSS = tsss_in_chr.get(idx_tss);
			
			// if with TAD constraint
			if (tad_file!=null){
				int idx = Collections.binarySearch(chr2tads.get(chr), p.expand(0));
				Region tad = null;
				if (idx<0){
					idx = -(idx+1) -1;  // insert point - 1 ==> Previous object
					tad = chr2tads.get(chr).get(idx);
					if (!tad.contains(p)){
//						System.err.println(String.format("Point %s is not within any TAD!", p.toString()));
						for (String g:tss2genes.get(nearestTSS))
							sb.append(p.expand(500).toString()+"\t"+g+"\t"+p.toString()+"\t"+nearestTSS.toString()+"\t"+p.distance(nearestTSS)+"\t"+-9+"\n");
					}
					else{
						// now tad contains the enhancer coord
						if (tad.contains(nearestTSS)){
							// the nearest TSS is within the TAD
							for (String g:tss2genes.get(nearestTSS))
								sb.append(p.expand(500).toString()+"\t"+g+"\t"+p.toString()+"\t"+nearestTSS.toString()+"\t"+p.distance(nearestTSS)+"\t"+1+"\n");
						}
						else{ // the nearest TSS is not within the TAD
							// do we have other TSS in TAD?
							int idx_tad_start = Collections.binarySearch(tsss_in_chr, new Point(genome, chr, tad.getStart()));
							if (idx_tad_start<0)
								idx_tad_start = -(idx_tad_start+1);
							int idx_tad_end = Collections.binarySearch(tsss_in_chr, new Point(genome, chr, tad.getEnd()));
							if (idx_tad_end<0)
								idx_tad_end = -(idx_tad_end+1)-1;
							if (idx_tad_start<=idx_tad_end){	// found TSSs within TAD
								// either the start or the end should be the nearest in TAD
								if (p.distance(tsss_in_chr.get(idx_tad_start)) > p.distance(tsss_in_chr.get(idx_tad_end))){
									nearestTSS = tsss_in_chr.get(idx_tad_end);
								}
								else
									nearestTSS = tsss_in_chr.get(idx_tad_start);
							}
							if (tad.contains(nearestTSS)){
								// the new nearest TSS is within the TAD
								for (String g:tss2genes.get(nearestTSS))
									sb.append(p.expand(500).toString()+"\t"+g+"\t"+p.toString()+"\t"+nearestTSS.toString()+"\t"+p.distance(nearestTSS)+"\t"+2+"\n");
							}
							else{ 							
								for (String g:tss2genes.get(nearestTSS))
	//								sb.append(String.format("%s\t%s\t%s\t%s\t%d\t%d\n", p.expand(500).toString(), g, p.toString(), nearestTSS.toString(), p.distance(nearestTSS), 0));
									sb.append(p.expand(500).toString()+"\t"+g+"\t"+p.toString()+"\t"+nearestTSS.toString()+"\t"+p.distance(nearestTSS)+"\t"+0+"\n");
							}
						}
					}
				}
				else{
					System.err.println(String.format("Point %s matches TAD!", chr2tads.get(chr).get(idx).toString()));
				}
			}
			else{ // if without TAD constraint
				for (String g:tss2genes.get(nearestTSS))
					sb.append(p.expand(500).toString()+"\t"+g+"\t"+p.toString()+"\t"+nearestTSS.toString()+"\t"+p.distance(nearestTSS)+"\t"+0+"\n");
			}
		} // for each coord
		
//		System.out.print(sb.toString());
		CommonUtils.writeFile(coords_file.replace(".txt", (tssRange!=-1?(".tss"+tssRange):"") 
				+ "."+tad_name + ".geneAssignments.txt"), sb.toString());

	}
	
	private void compute_enhancer_distances(){
		ArrayList<String> lines = CommonUtils.readTextFile(Args.parseString(args, "g2e", null));
		String prev = "";
		ArrayList<Point> ps = new ArrayList<Point>();
		for (String l:lines){
			String f[] = l.split("\t");
			if (prev.equals(f[0])){		// same gene
				ps.add(Point.fromString(genome, f[1]));
			}
			else{	// different gene, if prev only 1 region, ignore, if more regions, compute distance
				if (ps.size()>1){		// more than 1 enhancer
					ArrayList<Integer> distances = new ArrayList<Integer>();
					for (int i=0;i<ps.size();i++)
						for (int j=i+1;j<ps.size();j++){
							int d = ps.get(i).distance(ps.get(j));
							distances.add(d);
							System.out.println(d);
						}
				}
				prev = f[0];
				ps.clear();
				ps.add(Point.fromString(genome, f[1]));
			}
		}
		if (ps.size()>1){		// more than 1 enhancer
			ArrayList<Integer> distances = new ArrayList<Integer>();
			for (int i=0;i<ps.size();i++)
				for (int j=i+1;j<ps.size();j++){
					int d = ps.get(i).distance(ps.get(j));
					distances.add(d);
					System.out.println(d);
				}
		}
	}
	
}

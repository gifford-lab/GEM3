package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
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

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedPoint;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.DeepSeqExpt;
import edu.mit.csail.cgs.deepseq.discovery.KPPMixture;
import edu.mit.csail.cgs.deepseq.features.ComponentFeature;
import edu.mit.csail.cgs.deepseq.features.Feature;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.deepseq.utilities.ReadCache;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSPeak;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public class CID {
	Genome genome;
	Set<String> flags;
	String[] args;
	
	boolean isDev=false;
	boolean use_1_end_reads=false;	// for estimating anchor point position, not for MICC quantification
	
	int min_span = 4000;		// min PET span to exclude self-ligation reads
	int dc = 300;				// distance_cutoff for density clustering
	int read_1d_merge_dist = 2000;
	int max_cluster_merge_dist = 5000;
	int distance_factor = 40;
	int distance_base = 200;
	int tss_merge_dist = 500;
	int tss_radius = 2000;
	int chiapet_radius = 2000;
	int span_anchor_ratio = 15;		// span / anchor_width 
	int min_pet_count = 2; 		// minimum number of PET count to be called as an interaction
	int max_pet_count = 10000;	// max number of PET count for a region pair to run density clustering, split if larger
	int numQuantile = 100;
	int micc_min_pet = 1;

	TreeMap<Region, InteractionCall> r2it = new TreeMap<Region, InteractionCall>();
	String fileName = null;
	
	public CID(String[] args) {
		genome = CommonUtils.parseGenome(args);
		flags = Args.parseFlags(args);
		this.args = args;
		isDev = flags.contains("dev");
		use_1_end_reads = flags.contains("1end");
		
		dc = Args.parseInteger(args, "dc", dc);
		read_1d_merge_dist = Args.parseInteger(args, "read_merge_dist", read_1d_merge_dist);
		tss_merge_dist = Args.parseInteger(args, "tss_merge_dist", tss_merge_dist);
		max_cluster_merge_dist = Args.parseInteger(args, "max_cluster_merge_dist", max_cluster_merge_dist);
		distance_factor = Args.parseInteger(args, "distance_factor", distance_factor);
		tss_radius = Args.parseInteger(args, "tss_radius", tss_radius);
		chiapet_radius = Args.parseInteger(args, "chiapet_radius", chiapet_radius);
		numQuantile = Args.parseInteger(args, "num_span_quantile", 100);
		min_span = Args.parseInteger(args, "min_span", 4000);
		micc_min_pet = Args.parseInteger(args, "micc", micc_min_pet);
	}

	public static void main(String args[]) {
		CID analysis = new CID(args);
		int type = Args.parseInteger(args, "type", 3);

		switch (type) {
		case 0:
			analysis.cleanUpOverlaps();
			analysis.countGenesPerRegion();
			analysis.StatsTAD();
			break;
		case 1: // count distal read pairs per gene (old: step1)
			analysis.countReadPairs();
			break;
		case 2: // count distal read pairs per gene (old: step2)
			analysis.clusterDistalReads();
			break;
		case 3: // region (1D merged-read) based clustering
			analysis.findAllInteractions();
			break;
		case 30: // region (1D merged-read) based clustering
//			analysis.findAllInteractionsMerging();
			break;
		case 4: // find gene-based dense cluster of read pairs
			postProcessing(args);
			break;
		case 5:
			annotateRegions(args);
			break;
		case 7:
			annotateTADs(args);
			break;
		case 6: // merged-TSS based clustering
			getPetLength(args);
			break;
		}
	}

	/**
	 * Clean up data. Merge overlap distal regions if they are connected to the
	 * same TSS. Optionally remove distal regions that overlap with the
	 * connected TSS
	 */
	void cleanUpOverlaps() {
		ArrayList<String> texts = CommonUtils.readTextFile(fileName);
		// use tssString as the key to ensure uniqueness, tss Point objects have
		// different references even if from same tss
		TreeMap<String, ArrayList<Region>> tss2distalRegions = new TreeMap<String, ArrayList<Region>>();
		for (String line : texts) {
			String f[] = line.trim().split("\t");
			InteractionCall it = new InteractionCall();
			it.tssString = f[6];
			it.tss = Point.fromString(genome, f[6]);
			it.geneID = f[7];
			it.geneSymbol = f[8];
			if (f.length <= 13) {
				it.distal = Region.fromString(genome, f[9]);
				it.pvalue = Double.parseDouble(f[11]);
			} else {
				it.distal = Region.fromString(genome, f[13]);
				it.pvalue = Double.parseDouble(f[15]);
			}

			// skip interactions that have distal regions containing the TSS?
			if (flags.contains("rm_self")) {
				if (it.distal.contains(it.tss))
					continue;
			}
			// if duplicate, shrink 1bp to make it uniquie,
			// if connect to same TSS, it will be merged later
			while (r2it.containsKey(it.distal))
				it.distal = it.distal.expand(-1, -1);
			r2it.put(it.distal, it);
			if (!tss2distalRegions.containsKey(it.tssString))
				tss2distalRegions.put(it.tssString, new ArrayList<Region>());
			tss2distalRegions.get(it.tssString).add(it.distal);
		}

		// for each tss, merge overlapping distal regions (<1kb)
		for (String tss : tss2distalRegions.keySet()) {
			ArrayList<Region> regions = tss2distalRegions.get(tss);
			ArrayList<Region> mergedRegions = new ArrayList<Region>();
			Collections.sort(regions);
			Region previous = regions.get(0);
			ArrayList<Region> previousRegions = new ArrayList<Region>();
			previousRegions.add(previous);

			for (int i = 1; i < regions.size(); i++) {
				Region region = regions.get(i);
				// if overlaps with previous region, combine the regions, take
				// the best p-values
				if (previous.overlaps(region)) {
					previous = previous.combine(region);
					previousRegions.add(region);
				} else { // not overlap any more, update, then move to next one
					mergedRegions.add(previous);
					// merge overlapping regions, update interactions
					if (previousRegions.size() > 1) { // merged
						InteractionCall it = null;
						double bestPvalue = 1;
						for (Region r : previousRegions) {
							it = r2it.get(r);
							r2it.remove(r); // remove old one
							bestPvalue = Math.min(bestPvalue, it.pvalue);
						}
						it.distal = previous; // previous has been merged
						it.pvalue = bestPvalue;
						r2it.put(previous, it); // add merged region
					}
					previousRegions.clear();
					previous = region;
					previousRegions.add(previous);
				}
			}
			mergedRegions.add(previous);
			if (previousRegions.size() > 1) { // merged
				InteractionCall it = null;
				// merge overlapping regions, update interactions
				double bestPvalue = 1;
				for (Region r : previousRegions) {
					it = r2it.get(r);
					r2it.remove(r);
					bestPvalue = Math.min(bestPvalue, it.pvalue);
				}
				it.distal = previous;
				it.pvalue = bestPvalue;
				r2it.put(previous, it);
			}
			mergedRegions.trimToSize();
		}

		// print out cleaned up data
		if (flags.contains("print_merged")) {
			StringBuilder sb = new StringBuilder();
			for (Region r : r2it.keySet()) {
				InteractionCall it = r2it.get(r);
				sb.append(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%.2e", it.distal.toBED(),
						it.tss.expand(2000).toBED(), it.tss.toString(), it.geneID, it.geneSymbol, it.distal.toString(),
						it.distal.getWidth(), it.pvalue)).append("\n");
			}
			CommonUtils.writeFile(
					fileName.replace(".bedpe", (flags.contains("rm_self") ? ".rmSelf" : "") + ".mergeDistal.bedpe"),
					sb.toString());
		}
	}

	/**
	 * Count number of genes connected by each region (overlapping regions are
	 * merged).
	 * 
	 * @param fileName
	 */
	void countGenesPerRegion() {
		// load RefSeq gene annotation
		int tssRange = Args.parseInteger(args, "tss_range", 100);
		ArrayList<String> texts = CommonUtils.readTextFile(Args.parseString(args, "genes", null));
		TreeMap<StrandedPoint, TreeSet<String>> tss2genes = new TreeMap<StrandedPoint, TreeSet<String>>();
		for (int i = 0; i < texts.size(); i++) {
			String t = texts.get(i);
			if (t.startsWith("#"))
				continue;
			String f[] = t.split("\t");
			String chr = f[2].replace("chr", "");
			char strand = f[3].charAt(0);
			StrandedPoint tss = new StrandedPoint(genome, chr, Integer.parseInt(f[strand == '+' ? 4 : 5]), strand);
			String symbol = f[12];
			if (!tss2genes.containsKey(tss))
				tss2genes.put(tss, new TreeSet<String>());
			tss2genes.get(tss).add(symbol);
		}
		/** all TSS annotatons organized by chroms */
		HashMap<String, ArrayList<StrandedPoint>> chr2tsss = new HashMap<String, ArrayList<StrandedPoint>>();
		for (StrandedPoint tss : tss2genes.keySet()) {
			String chr = tss.getChrom();
			if (!chr2tsss.containsKey(chr))
				chr2tsss.put(chr, new ArrayList<StrandedPoint>());
			chr2tsss.get(chr).add(tss);
		}
		for (String chr : chr2tsss.keySet()) {
			ArrayList<StrandedPoint> tsss = chr2tsss.get(chr);
			Collections.sort(tsss);
			// print inter-TSS distances
			// for (int i=0;i<tsss.size()-1;i++){
			// System.out.println(tsss.get(i+1).distance(tsss.get(i)));
			// }
		}

		// merge nearby regions, collect genes within tssRange of the TSSs.
		TreeSet<String> allGenes = new TreeSet<String>();
		ArrayList<Region> rs = new ArrayList<Region>();
		// if coords are provided, only report the subset of distal regions that
		// overlaps
		String coords_file = Args.parseString(args, "coords", null);
		String coords_name = Args.parseString(args, "coords_name", null);
		if (coords_file != null) {
			ArrayList<Point> coords = CommonUtils.loadCgsPointFile(coords_file, genome);
			for (Region r : r2it.keySet()) {
				for (Point p : coords) {
					if (!r.getChrom().equalsIgnoreCase(p.getChrom()))
						continue;
					if (r.distance(p) <= 500) { // add distal anchor if overlap
												// with enhancer coords
						rs.add(r);
						r2it.get(r).overlapCoords.add(p);
					}
				}
			}

			// print the enhancer_coords and the tss_geneSymbol pairs
			StringBuilder sb1 = new StringBuilder();
			sb1.append("#Coord\tTSS\tGene\n");
			for (Region r : rs) {
				InteractionCall it = r2it.get(r);
				Point tss = it.tss;
				ArrayList<StrandedPoint> tsss = chr2tsss.get(tss.getChrom());
				ArrayList<StrandedPoint> tss_linked = new ArrayList<StrandedPoint>();
				if (tsss == null) {
					continue;
				}
				int idx = Collections.binarySearch(tsss, tss);
				if (idx < 0)
					idx = -(idx + 1); // insert point
				for (int j = idx; j < tsss.size(); j++) {
					if (tss.distance(tsss.get(j)) > tssRange)
						break;
					else
						tss_linked.add(tsss.get(j));
				}
				for (int j = idx - 1; j >= 0; j--) {
					if (tss.distance(tsss.get(j)) > tssRange)
						break;
					else
						tss_linked.add(tsss.get(j));
				}
				for (StrandedPoint t : tss_linked) {
					for (Point p : it.overlapCoords) {
						TreeSet<String> genes = tss2genes.get(t);
						for (String g : genes)
							sb1.append(p.toString()).append("\t").append(t.toString()).append("\t").append(g)
									.append("\n");
					}
				}
			}
			CommonUtils.writeFile(fileName.replace(".bedpe", (flags.contains("rm_self") ? ".rmSelf" : "")
					+ (coords_file != null ? ("." + coords_name) : "") + ".tss" + tssRange + ".coord2genes.txt"),
					sb1.toString());

		} else
			rs.addAll(r2it.keySet()); // add all distal regions, rs is sorted,
										// r2it is a TreeMap

		Region merged = rs.get(0);
		ArrayList<Integer> ids = new ArrayList<Integer>();
		ids.add(0);
		StringBuilder sb = new StringBuilder("#Merged_region\twidth\tr_id\tgenes\tcount\n");

		for (int i = 1; i < rs.size(); i++) {
			Region r = rs.get(i);
			if (merged.getChrom().equals(r.getChrom())
					&& merged.distance(r) < Args.parseInteger(args, "distance", 1000)) {
				merged = merged.combine(r);
				ids.add(i);
			} else {
				sb.append(merged.toString()).append("\t").append(merged.getWidth()).append("\t");

				// for each merged region, get all TSS, find refSeq genes within
				// 100bp window
				TreeSet<String> genes = new TreeSet<String>();
				for (int id : ids) {
					InteractionCall it = r2it.get(rs.get(id));
					sb.append(id).append(",");
					Point tss = it.tss;
					ArrayList<StrandedPoint> tsss = chr2tsss.get(tss.getChrom());
					if (tsss == null) {
						genes.add(it.geneSymbol);
						continue;
					}
					int idx = Collections.binarySearch(tsss, tss);
					if (idx < 0)
						idx = -(idx + 1); // insert point
					for (int j = idx; j < tsss.size(); j++) {
						if (tss.distance(tsss.get(j)) > tssRange)
							break;
						else
							genes.addAll(tss2genes.get(tsss.get(j)));
					}
					for (int j = idx - 1; j >= 0; j--) {
						if (tss.distance(tsss.get(j)) > tssRange)
							break;
						else
							genes.addAll(tss2genes.get(tsss.get(j)));
					}
				}
				CommonUtils.replaceEnd(sb, '\t');
				for (String s : genes)
					sb.append(s).append(",");
				CommonUtils.replaceEnd(sb, '\t');
				sb.append(genes.size()).append("\n");
				ids.clear();
				ids.add(i);
				allGenes.addAll(genes);

				merged = r; // setup for next merge
			}
		}
		// finish the last merge
		sb.append(merged.toString()).append("\t").append(merged.getWidth()).append("\t");
		TreeSet<String> genes = new TreeSet<String>();
		for (int id : ids) {
			sb.append(id).append(",");
			Point tss = r2it.get(rs.get(id)).tss;
			ArrayList<StrandedPoint> tsss = chr2tsss.get(tss.getChrom());
			if (tsss == null) {
				genes.add(r2it.get(rs.get(id)).geneSymbol);
				continue;
			}
			int idx = Collections.binarySearch(tsss, tss);
			if (idx < 0)
				idx = -(idx + 1); // insert point
			for (int j = idx; j < tsss.size(); j++) {
				if (tss.distance(tsss.get(j)) > tssRange)
					break;
				else
					genes.addAll(tss2genes.get(tsss.get(j)));
			}
			for (int j = idx - 1; j >= 0; j--) {
				if (tss.distance(tsss.get(j)) > tssRange)
					break;
				else
					genes.addAll(tss2genes.get(tsss.get(j)));
			}
		}
		CommonUtils.replaceEnd(sb, '\t');
		if (genes.isEmpty())
			sb.append("None").append(",");
		for (String s : genes)
			sb.append(s).append(",");
		CommonUtils.replaceEnd(sb, '\t');
		sb.append(genes.size()).append("\n");
		allGenes.addAll(genes);

		// System.out.println(sb.toString());
		CommonUtils.writeFile(fileName.replace(".bedpe", (flags.contains("rm_self") ? ".rmSelf" : "")
				+ (coords_file != null ? ("." + coords_name) : "") + ".tss" + tssRange + ".per_region_count.txt"),
				sb.toString());

		// print out all linked genes
		sb = new StringBuilder();
		for (String g : allGenes)
			sb.append(g).append("\n");
		CommonUtils.writeFile(
				fileName.replace(".bedpe", (flags.contains("rm_self") ? ".rmSelf" : "")
						+ (coords_file != null ? ("." + coords_name) : "") + ".tss" + tssRange + ".geneSymbols.txt"),
				sb.toString());
	}

	void StatsTAD() {
		String tad_file = Args.parseString(args, "tad", null);
		ArrayList<Region> tads = CommonUtils.load_BED_regions(genome, tad_file).car();
		Collections.sort(tads);
		ArrayList<InteractionCall> itNonTAD = new ArrayList<InteractionCall>();
		ArrayList<InteractionCall> itSameTAD = new ArrayList<InteractionCall>();
		ArrayList<InteractionCall> itCrossTAD = new ArrayList<InteractionCall>();

		for (Region r : r2it.keySet()) {
			Region mid = r.getMidpoint().expand(0);
			int idx = Collections.binarySearch(tads, mid);
			Region tad = null;
			if (idx < 0) {
				idx = -(idx + 1) - 1; // insert point - 1 ==> Previous object
				tad = tads.get(idx);
				if (!tad.contains(mid)) {
					// System.err.println(String.format("Point %s is not within
					// any TAD!", p.toString()));
					itNonTAD.add(r2it.get(r));
				} else {// now tad contains the distal coord
					if (tad.contains(r2it.get(r).tss))
						itSameTAD.add(r2it.get(r));
					else { // find the tad that TSS is in
						Region tss = r2it.get(r).tss.expand(0);
						idx = Collections.binarySearch(tads, tss);
						if (idx < 0) {
							idx = -(idx + 1) - 1; // insert point - 1 ==>
													// Previous object
							tad = tads.get(idx);
							if (!tad.contains(tss)) // TSS is not in a TAD
								itNonTAD.add(r2it.get(r));
							else // in TAD, must be another TAD
								itCrossTAD.add(r2it.get(r));
						}
					}
				}
			}
		}

		System.out.println(String.format("In same TAD:\t %d\nCross TAD:\t %d\nNot in TAD:\t %d\n", itSameTAD.size(),
				itCrossTAD.size(), itNonTAD.size()));

		StringBuilder sb = new StringBuilder(
				"#Interaction\tdistance\tdistal\ttss\tSymbol\tgeneID\tp_-lg10\tTAD_status\n");
		for (InteractionCall it : itSameTAD)
			sb.append(it.toString()).append("\t1\n");
		for (InteractionCall it : itCrossTAD)
			sb.append(it.toString()).append("\t2\n");
		for (InteractionCall it : itNonTAD)
			sb.append(it.toString()).append("\t0\n");
		CommonUtils.writeFile(
				fileName.replace(".bedpe", (flags.contains("rm_self") ? ".rmSelf" : "") + ".StatsTAD.txt"),
				sb.toString());
	}

	private void countReadPairs() {
		long tic = System.currentTimeMillis();
		HashSet<String> geneSet = new HashSet<String>();
		String gString = Args.parseString(args, "genes", null);
		if (gString == null) {
			String gFile = Args.parseString(args, "gene_file", null);
			ArrayList<String> lines = CommonUtils.readTextFile(gFile);
			for (String g : lines)
				geneSet.add(g.trim());
		} else {
			String genes[] = Args.parseString(args, "genes", null).split(",");
			for (String g : genes)
				geneSet.add(g.trim());
		}

		// load refSeq gene annotation
		int tssRadius = Args.parseInteger(args, "tss_range", 10001) / 2;
		int chiapetRadius = Args.parseInteger(args, "chiapet_radius", 2000);
		ArrayList<String> gene_annots = CommonUtils.readTextFile(Args.parseString(args, "gene_anno", null));
		TreeMap<String, TreeSet<StrandedPoint>> gene2tss = new TreeMap<String, TreeSet<StrandedPoint>>();
		for (int i = 0; i < gene_annots.size(); i++) {
			String t = gene_annots.get(i);
			if (t.startsWith("#"))
				continue;
			String f[] = t.split("\t");
			String symbol = f[12];
			if (!geneSet.contains(symbol))
				continue;
			String chr = f[2].replace("chr", "");
			char strand = f[3].charAt(0);
			StrandedPoint tss = new StrandedPoint(genome, chr, Integer.parseInt(f[strand == '+' ? 4 : 5]), strand);
			if (!gene2tss.containsKey(symbol))
				gene2tss.put(symbol, new TreeSet<StrandedPoint>());
			gene2tss.get(symbol).add(tss);
		}

		// load read pairs
		// for now, just loop through every read pair. To run faster, should
		// sort by each end and use binary search to find tss end overlaps
		ArrayList<String> read_pairs = CommonUtils.readTextFile(Args.parseString(args, "read_pair", null));
		HashMap<String, Pair<ArrayList<Point>, ArrayList<Point>>> chr2reads = new HashMap<String, Pair<ArrayList<Point>, ArrayList<Point>>>();
		for (String s : read_pairs) {
			String[] f = s.split("\t");
			Point r1 = Point.fromString(genome, f[0]);
			String r1Chrom = r1.getChrom();
			Point r2 = Point.fromString(genome, f[1]);
			if (!r1Chrom.equals(r2.getChrom())) // skip if not from the same
												// chromosome
				continue;
			if (!chr2reads.containsKey(r1Chrom)) {
				ArrayList<Point> r1s = new ArrayList<Point>();
				ArrayList<Point> r2s = new ArrayList<Point>();
				chr2reads.put(r1Chrom, new Pair<ArrayList<Point>, ArrayList<Point>>(r1s, r2s));
			}
			Pair<ArrayList<Point>, ArrayList<Point>> reads = chr2reads.get(r1Chrom);
			reads.car().add(r1);
			reads.cdr().add(r2);
		}

		System.out.println("Loaded ChIA-PET read pairs: " + CommonUtils.timeElapsed(tic));
		System.out.println();

		// load TF sites
		ArrayList<String> tfs = CommonUtils.readTextFile(Args.parseString(args, "tf_sites", null));
		ArrayList<List<GPSPeak>> allPeaks = new ArrayList<List<GPSPeak>>();
		for (int i = 0; i < tfs.size(); i++) {
			try {
				allPeaks.add(GPSParser.parseGPSOutput(tfs.get(i), genome));
				System.out.println("Loaded " + tfs.get(i));
			} catch (IOException e) {
				System.out.println(tfs.get(i) + " does not have a valid GPS/GEM event call file.");
				e.printStackTrace(System.err);
				System.exit(1);
			}
		}

		// load histone mark regions
		ArrayList<String> hms = CommonUtils.readTextFile(Args.parseString(args, "regions", null));
		ArrayList<List<Region>> allRegions = new ArrayList<List<Region>>();
		for (int i = 0; i < hms.size(); i++) {
			allRegions.add(CommonUtils.load_BED_regions(genome, hms.get(i)).car());
			System.out.println("Loaded " + hms.get(i));
		}
		System.out.println();
		// TreeMap<String, ArrayList<Integer>> gene2distances = new
		// TreeMap<String, ArrayList<Integer>>();
		// compute distance for each gene
		ArrayList<String> geneList = new ArrayList<String>();
		geneList.addAll(gene2tss.keySet());
		StringBuilder sb = new StringBuilder();

		for (int id = 0; id < geneList.size(); id++) {
			String g = geneList.get(id);
			System.out.print(g + " ");
			ArrayList<Integer> distances = new ArrayList<Integer>();
			ArrayList<ArrayList<Integer>> isTfBounds = new ArrayList<ArrayList<Integer>>();
			TreeSet<StrandedPoint> coords = gene2tss.get(g);
			// if gene has multiple TSSs, use the center position
			int count = coords.size();
			StrandedPoint centerPoint = null;
			for (StrandedPoint p : coords) {
				if (count < coords.size() / 2)
					break;
				else {
					centerPoint = p;
					count--;
				}
			}

			// if one end of the read pair is near TSS, compute the offset of
			// the other end
			boolean isMinus = centerPoint.getStrand() == '-';
			Pair<ArrayList<Point>, ArrayList<Point>> reads = chr2reads.get(centerPoint.getChrom());
			if (reads == null)
				continue;
			ArrayList<Point> read1s = reads.car();
			ArrayList<Point> read2s = reads.cdr();
			for (int i = 0; i < read1s.size(); i++) {
				int offset_p1 = read1s.get(i).offset(centerPoint);
				int offset_p2 = read2s.get(i).offset(centerPoint);
				int dist_p1 = Math.abs(offset_p1);
				int dist_p2 = Math.abs(offset_p2);
				// only add distance to the list if one read is within
				// TSS_Radius, the other read is outside of TSS_Radius
				if (dist_p1 < tssRadius) {
					if (dist_p2 > tssRadius) {
						distances.add(isMinus ? -offset_p2 : offset_p2);
						ArrayList<Integer> isBound = new ArrayList<Integer>();
						Point p = read2s.get(i);
						for (int j = 0; j < allPeaks.size(); j++) {
							List<GPSPeak> peaks = allPeaks.get(j);
							int bound = 0;
							for (GPSPeak gps : peaks) {
								if (gps.getChrom().equals(p.getChrom()) && gps.distance(p) <= chiapetRadius) {
									bound = 1;
									break;
								}
							}
							isBound.add(bound);
						}
						for (int j = 0; j < allRegions.size(); j++) {
							List<Region> rs = allRegions.get(j);
							int bound = 0;
							for (Region r : rs) {
								// if the region r contains point p, or the
								// distance between midPoint of r and p is less
								// than ChIAPET_radias
								if (r.getChrom().equals(p.getChrom())
										&& (r.contains(p) || r.getMidpoint().distance(p) <= chiapetRadius)) {
									bound = 1;
									break;
								}
							}
							isBound.add(bound);
						}
						isTfBounds.add(isBound);
					}
				} else {
					if (dist_p2 < tssRadius) {
						distances.add(isMinus ? -offset_p1 : offset_p1);
						ArrayList<Integer> isBound = new ArrayList<Integer>();
						Point p = read1s.get(i);
						for (int j = 0; j < allPeaks.size(); j++) {
							List<GPSPeak> peaks = allPeaks.get(j);
							int bound = 0;
							for (GPSPeak gps : peaks) {
								if (gps.getChrom().equals(p.getChrom()) && gps.distance(p) <= chiapetRadius) {
									bound = 1;
									break;
								}
							}
							isBound.add(bound);
						}
						for (int j = 0; j < allRegions.size(); j++) {
							List<Region> rs = allRegions.get(j);
							int bound = 0;
							for (Region r : rs) {
								// if the region r contains point p, or the
								// distance between midPoint of r and p is less
								// than ChIAPET_radias
								if (r.getChrom().equals(p.getChrom())
										&& (r.contains(p) || r.getMidpoint().distance(p) <= chiapetRadius)) {
									bound = 1;
									break;
								}
							}
							isBound.add(bound);
						}
						isTfBounds.add(isBound);
					}
				}
			}

			if (!distances.isEmpty()) {
				// gene2distances.put(g, distances);
				for (int i = 0; i < distances.size(); i++) {
					sb.append(g).append("\t").append(centerPoint.toString()).append("\t").append(id);
					sb.append("\t").append(distances.get(i));
					for (int b : isTfBounds.get(i))
						sb.append("\t").append(b);
					sb.append("\n");
				}
			}
		} // for each gene
		CommonUtils.writeFile("all_genes.distal_offsets.txt", sb.toString());

		System.out.println("\n\n" + CommonUtils.timeElapsed(tic));
	}

	private void clusterDistalReads() {
		int tss_exclude = Args.parseInteger(args, "tss_exclude", 8000);
		int step = Args.parseInteger(args, "merge_dist", 1500);
		int minRead = Args.parseInteger(args, "min_count", 2);

		// load data
		ArrayList<String> lines = CommonUtils.readTextFile(Args.parseString(args, "tss_reads", null));
		TSSwithReads tss = new TSSwithReads();
		tss.symbol = "---";
		ArrayList<TSSwithReads> allTss = new ArrayList<TSSwithReads>();
		for (String l : lines) { // each line is a distal read
			String f[] = l.split("\t");
			int offset = Integer.parseInt(f[3]);
			if (Math.abs(offset) < tss_exclude) // skip the read if it is within
												// TSS exclude region
				continue;
			if (!f[0].equals(tss.symbol)) { // a new gene
				tss = new TSSwithReads();
				tss.symbol = f[0];
				tss.coord = StrandedPoint.fromString(genome, f[1]);
				tss.id = Integer.parseInt(f[2]);
				tss.reads = new TreeMap<Integer, ArrayList<Boolean>>();
				allTss.add(tss);
			}
			ArrayList<Boolean> isBound = new ArrayList<Boolean>();
			for (int i = 4; i < f.length; i++) {
				isBound.add(f[i].equals("1"));
			}
			tss.reads.put(offset, isBound);
		}

		lines = CommonUtils.readTextFile(Args.parseString(args, "germ", null));
		HashMap<Point, ArrayList<Point>> germTss2distals = new HashMap<Point, ArrayList<Point>>();
		ArrayList<Point> germTss = new ArrayList<Point>();
		for (String l : lines) { // each line is a call
			String f[] = l.split("\t");
			Point t = new Region(genome, f[3].replace("chr", ""), Integer.parseInt(f[4]), Integer.parseInt(f[5]))
					.getMidpoint();
			if (!germTss2distals.containsKey(t))
				germTss2distals.put(t, new ArrayList<Point>());
			germTss2distals.get(t)
					.add(new Region(genome, f[0].replace("chr", ""), Integer.parseInt(f[1]), Integer.parseInt(f[2]))
							.getMidpoint());
		}
		germTss.addAll(germTss2distals.keySet());
		Collections.sort(germTss);

		lines = CommonUtils.readTextFile(Args.parseString(args, "mango", null));
		HashMap<Point, ArrayList<Point>> a2bs = new HashMap<Point, ArrayList<Point>>();
		HashMap<Point, ArrayList<Point>> b2as = new HashMap<Point, ArrayList<Point>>();
		for (String l : lines) { // each line is a call
			String f[] = l.split("\t");
			Point a = new Region(genome, f[0].replace("chr", ""), Integer.parseInt(f[1]), Integer.parseInt(f[2]))
					.getMidpoint();
			Point b = new Region(genome, f[3].replace("chr", ""), Integer.parseInt(f[4]), Integer.parseInt(f[5]))
					.getMidpoint();
			if (!a2bs.containsKey(a))
				a2bs.put(a, new ArrayList<Point>());
			a2bs.get(a).add(b);
			if (!b2as.containsKey(b))
				b2as.put(b, new ArrayList<Point>());
			b2as.get(b).add(a);
		}
		ArrayList<Point> aPoints = new ArrayList<Point>();
		aPoints.addAll(a2bs.keySet());
		Collections.sort(aPoints);
		ArrayList<Point> bPoints = new ArrayList<Point>();
		bPoints.addAll(b2as.keySet());
		Collections.sort(bPoints);

		// cluster the reads
		for (TSSwithReads t : allTss) {
			ArrayList<Integer> cluster = new ArrayList<Integer>();
			for (int offset : t.reads.keySet()) {
				if (cluster.isEmpty() || offset - cluster.get(cluster.size() - 1) < step) {
					cluster.add(offset);
				} else { // have a large distance, finish old cluster, create
							// new cluster
					if (cluster.size() >= minRead) { // at least 2 reads
						int median = cluster.get(cluster.size() / 2);
						Point tssPoint = t.coord;
						Point distalPoint = new Point(genome, t.coord.getChrom(),
								t.coord.getLocation() + (t.coord.getStrand() == '+' ? median : -median));
						Region distalRegion = null;
						if (t.coord.getStrand() == '+') {
							distalRegion = new Region(genome, t.coord.getChrom(),
									t.coord.getLocation() + cluster.get(0),
									t.coord.getLocation() + cluster.get(cluster.size() - 1));
						} else {
							distalRegion = new Region(genome, t.coord.getChrom(),
									t.coord.getLocation() - cluster.get(cluster.size() - 1),
									t.coord.getLocation() - cluster.get(0));
						}
						// print result if the read cluster is not in the tss
						// exclusion range
						System.out.print(String.format("%s\t%s\t%s\t%s\t%d\t%d\t%d\t", t.symbol,
								t.coord.getLocationString(), distalRegion.getLocationString(),
								distalPoint.getLocationString(), median, cluster.size(), distalRegion.getWidth()));

						// print binding overlap information
						int count = t.reads.get(cluster.get(0)).size();
						for (int c = 0; c < count; c++) {
							boolean isBound = false;
							for (int clusterOffset : cluster) {
								isBound = isBound || t.reads.get(clusterOffset).get(c);
							}
							System.out.print(isBound ? "1\t" : "0\t");
						}

						// print ChIA-PET call overlap info
						Point tssLeft = new Point(genome, t.coord.getChrom(), t.coord.getLocation() - 2000);
						Point tssRight = new Point(genome, t.coord.getChrom(), t.coord.getLocation() + 2000);

						// GERM
						int index = Collections.binarySearch(germTss, tssLeft);
						if (index < 0) // if key not found
							index = -(index + 1);
						int indexRight = Collections.binarySearch(germTss, tssRight);
						if (indexRight < 0) // if key not found
							indexRight = -(indexRight + 1);
						// if key match found, continue to search (
						// binarySearch() give undefined index with multiple
						// matches)
						boolean isOverlapped = false;
						indexRange: for (int i = index - 1; i <= indexRight + 2; i++) {
							if (i < 0 || i >= germTss.size())
								continue;
							try {
								Point tt = germTss.get(i);
								if (tt.distance(tssPoint) <= 2000) {
									if (!germTss2distals.containsKey(tt))
										continue;
									for (Point d : germTss2distals.get(tt)) {
										// System.out.print(tt.getLocationString()+"\t"+d.getLocationString());
										if (d.distance(distalPoint) <= 2000) {
											isOverlapped = true;
											// System.out.println("\tHIT");
											break indexRange;
										}
										// else
										// System.out.println();
									}
								}
							} catch (IllegalArgumentException e) { // ignore
							}
						}
						System.out.print(isOverlapped ? "1\t" : "0\t");

						// Mango
						index = Collections.binarySearch(aPoints, tssLeft);
						if (index < 0) // if key not found
							index = -(index + 1);
						indexRight = Collections.binarySearch(aPoints, tssRight);
						if (indexRight < 0) // if key not found
							indexRight = -(indexRight + 1);
						// if key match found, continue to search (
						// binarySearch() give undefined index with multiple
						// matches)
						isOverlapped = false;
						indexA: for (int i = index - 1; i <= indexRight + 2; i++) {
							if (i < 0 || i >= aPoints.size())
								continue;
							try {
								Point a = aPoints.get(i);
								if (a.distance(tssPoint) <= 2000) {
									if (!a2bs.containsKey(a))
										continue;
									for (Point b : a2bs.get(a)) {
										if (b.distance(distalPoint) <= 2000) {
											isOverlapped = true;
											break indexA;
										}
									}
								}
							} catch (IllegalArgumentException e) { // ignore
							}
						}
						if (isOverlapped)
							System.out.print("1\t");
						else {
							index = Collections.binarySearch(bPoints, tssLeft);
							if (index < 0) // if key not found
								index = -(index + 1);
							indexRight = Collections.binarySearch(bPoints, tssRight);
							if (indexRight < 0) // if key not found
								indexRight = -(indexRight + 1);
							// if key match found, continue to search (
							// binarySearch() give undefined index with multiple
							// matches)
							isOverlapped = false;
							indexB: for (int i = index - 1; i <= indexRight + 2; i++) {
								if (i < 0 || i >= bPoints.size())
									continue;
								try {
									Point b = bPoints.get(i);
									if (b.distance(tssPoint) <= 2000) {
										if (!b2as.containsKey(b))
											continue;
										for (Point a : b2as.get(b)) {
											if (a.distance(distalPoint) <= 2000) {
												isOverlapped = true;
												break indexB;
											}
										}
									}
								} catch (IllegalArgumentException e) { // ignore
								}
							}
							System.out.print(isOverlapped ? "1\t" : "0\t");
						}

						System.out.println();
					}
					cluster.clear();
					cluster.add(offset);
					continue;
				}
			}
		}
	}

	private int span2mergingDist(int span){
		return Math.min(max_cluster_merge_dist,
				Math.max(dc, distance_base+span/distance_factor));
	}
	
	// find interactions by 2D density clustering
	private void findAllInteractions() {
		long tic0 = System.currentTimeMillis();
		String outName = Args.parseString(args, "out", "CID");

		System.out.println(String.format("Options: --g \"%s\" --data \"%s\" --out \"%s\" --dc %d --read_merge_dist %d --distance_factor %d --max_cluster_merge_dist %d --min_span %d\n", 
				Args.parseString(args, "g", null), Args.parseString(args, "data", null), Args.parseString(args, "out", "Result"),
				dc, read_1d_merge_dist, distance_factor, max_cluster_merge_dist, min_span));

		// load read pairs
		// PETS: same chromosome, and longer than min_distance distance; 
		// the left read is required to be lower than the right read; if not, flip 

		// sort by each end so that we can search to find matches or overlaps
		System.out.println("Running CID on "+outName);
		System.out.println("\nLoading ChIA-PET read pairs ... ");
		
		String exChroms = Args.parseString(args, "ex", "M");
		HashSet<String> excludedChroms = new HashSet<String>();
		String exf[] = exChroms.split(",");
		for (String s: exf) {
			excludedChroms.add(s.trim().replace("chr", ""));
		}
		
//		ArrayList<Integer> dist_minus_plus = new ArrayList<Integer>();
//		ArrayList<Integer> dist_plus_minus = new ArrayList<Integer>();
//		ArrayList<Integer> dist_minus_minus = new ArrayList<Integer>();
//		ArrayList<Integer> dist_plus_plus = new ArrayList<Integer>();

		ArrayList<String> read_pair_lines = CommonUtils.readTextFile(Args.parseString(args, "data", null));
		String[] f1 = read_pair_lines.get(0).split("\t");
		boolean isBEDPE = f1.length >= 6;
		if (isBEDPE){
			System.out.println("\nDetected input data to be BEDPE format!");
			if (f1.length<10){
				System.err.println("Wrong BEDPE format. The columns 9 and 10 should be the strand information of the two read ends.");
				System.exit(-1);
			}
		}
		// store single ends
		// intra-chrom long PETs: always keep
		// inter-chrom PETs: keep if considering inter-chrom interactions
		// intra-chrom short PETs: i.e. self-ligation, keep when using 1d reads to estimate anchor points
		// single-end-reads: i.e. unpaired, keep only if use_1_end_reads flag is true, but may be useful when using 1d reads to estimate anchor points
		ArrayList<StrandedPoint> reads = new ArrayList<StrandedPoint>(); 
		
		// all PETs sorted by the low end
		ArrayList<ReadPair> low = new ArrayList<ReadPair>(); 
		ArrayList<ReadPair> high = new ArrayList<ReadPair>(); // sort high end
		StrandedPoint tmp1 = null;
		for (String s : read_pair_lines) {
			String[] f = s.split("\t");
			StrandedPoint r1;
			StrandedPoint r2;
			if (!isBEDPE){		// cgsPoints 
				r1 = StrandedPoint.fromString(genome, f[0]);
				r2 = StrandedPoint.fromString(genome, f[1]);
				if (excludedChroms.contains(r1.getChrom()) || excludedChroms.contains(r2.getChrom()))
					continue;
			}
			else{	// BEDPE format
				if (!use_1_end_reads && (f[0].charAt(0)=='*' || f[3].charAt(0)=='*'))
					continue;
				
				char strand1 = f[8].charAt(0);
//				r1 = new StrandedPoint(genome, f[0].replace("chr", ""), (Integer.parseInt(f[1])+Integer.parseInt(f[2]))/2, strand1);
				r1 = new StrandedPoint(genome, f[0].replace("chr", ""), strand1=='+'?Integer.parseInt(f[1]):Integer.parseInt(f[2]), strand1);
				char strand2 = f[9].charAt(0);
//				r2 = new StrandedPoint(genome, f[3].replace("chr", ""), (Integer.parseInt(f[4])+Integer.parseInt(f[5]))/2, strand2);
				r2 = new StrandedPoint(genome, f[3].replace("chr", ""), strand1=='+'?Integer.parseInt(f[4]):Integer.parseInt(f[5]), strand2);
				if (excludedChroms.contains(r1.getChrom()) || excludedChroms.contains(r2.getChrom()))
					continue;
				// if not both ends are aligned properly, skip, 
				// but if one of the read is mapped, add the read to the single-end reads object 
				if (r1.getChrom().equals("*")){
					// add read2 as single-end if mapped
					if (!r2.getChrom().equals("*"))	
						reads.add(r2);
					continue;
				}
				if (r2.getChrom().equals("*")){
					// add read1 as single end if mapped
					if (!r1.getChrom().equals("*")) 
						reads.add(r1);
					continue;
				}
			}
			
			// TODO: change next line if predicting inter-chrom interactions
			// r1 and r2 should be on the same chromosome for PETs
			if (!r1.getChrom().equals(r2.getChrom())) 
				continue;
			
			// TODO: should we treat inter-chrom reads as single-end reads???
			reads.add(r1);
			reads.add(r2);
			
			int dist = r1.distance(r2);
			if (dist < min_span)
				continue;
			if (r1.getLocation() > r2.getLocation()){	// r1 should be lower than r2
				tmp1 = r1;
				r2 = r1;
				r1 = tmp1;
			}
//			// count PETs by strand-orientation
//			if (r1.getStrand() == '-') {
//				if (r2.getStrand() == '+')
//					dist_minus_plus.add(dist);
//				else if (r2.getStrand() == '-')
//					dist_minus_minus.add(dist);
//			} else if (r1.getStrand() == '+') {
//				if (r2.getStrand() == '+')
//					dist_plus_plus.add(dist);
//				else if (r2.getStrand() == '-')
//					dist_plus_minus.add(dist);
//			}

			ReadPair rp = new ReadPair();
			rp.r1 = r1;
			rp.r2 = r2;
			low.add(rp);
			ReadPair rp2 = new ReadPair();
			rp2.r1 = r1;
			rp2.r2 = r2;
			high.add(rp2);
		}

		low.trimToSize();
		high.trimToSize();
		// sort by low end read1: default
		Collections.sort(low);
		Collections.sort(high, new Comparator<ReadPair>() {		// sort by high end read2
			public int compare(ReadPair o1, ReadPair o2) {
				return o1.compareRead2(o2);
			}
		});

		ArrayList<Point> lowEnds = new ArrayList<Point>();
		for (ReadPair r : low)
			lowEnds.add(r.r1);
		lowEnds.trimToSize();
		ArrayList<Point> highEnds = new ArrayList<Point>();
		for (ReadPair r : high)
			highEnds.add(r.r2);
		highEnds.trimToSize();

		reads.trimToSize();
		Collections.sort(reads);

		System.out.println("\nLoaded total single reads = " + reads.size() + ", filtered PETs =" + highEnds.size()
		+ " : " + CommonUtils.timeElapsed(tic0));

		// a hack to print out PETs that support a list of BEDPE-format loops
		String loop_file = Args.parseString(args, "loops", null);
		if (loop_file != null) {
			String outprefix = loop_file.replace(".bedpe", "");
			int binSize = Args.parseInteger(args, "bin", 100);
//			System.out.println(loop_file);
			ArrayList<String> lines = CommonUtils.readTextFile(loop_file);
			int count = 0;
			for (String anchorString : lines) {
				// System.out.println(anchorString);
				String[] f = anchorString.split("\t");
				Region region1 = new Region(genome, f[0].replace("chr", ""), Integer.parseInt(f[1]),
						Integer.parseInt(f[2]));
				Region region2 = new Region(genome, f[3].replace("chr", ""), Integer.parseInt(f[4]),
						Integer.parseInt(f[5]));
//				if (region1.overlaps(region2)){		// if overlap, merge
					region1 = new Region(genome, region1.getChrom(), region1.getStart(), region2.getEnd());
					region2 = region1;
//				}
				int r1start = region1.getStart();
				int r2start = region2.getStart();
				int bin1 = region1.getWidth()/binSize;		// j:read1:x
				int bin2 = region2.getWidth()/binSize;		// i:read2:y
				int[][] map = new int[bin2][bin1];
				for (int i=0;i<bin2;i++)			// init to 0
					for (int j=0;j<bin1;j++)	
						map[i][j] = 0;
				ArrayList<Integer> idx = CommonUtils.getPointsIdxWithinWindow(lowEnds, region1);
				for (int id : idx) {
					ReadPair rp = low.get(id);
					if (region2.contains(rp.r2)){
//						System.out.println(rp + "\t" + (rp.r1.getLocation() - r1start) + "\t"
//								+ (rp.r2.getLocation() - r2start));
						int j = (rp.r1.getLocation() - r1start) / binSize;	// j:read1:x
						if (j>=bin1)
							j=bin1-1;
						int i = (rp.r2.getLocation() - r2start) / binSize;	// i:read2:y
						if (i>=bin2)
							i=bin2-1;
						map[i][j] = map[i][j] +1;
					}
				}
				StringBuilder sb = new StringBuilder();
				sb.append("> ").append(region1.toString()).append(" ").append(region2.toString())
				.append(" ").append(f[6]).append(" ").append(region2.getEnd()-r1start).append("\n");
				sb.append(r2start-r1start).append("\t");
				for (int j=0;j<bin1;j++)
					sb.append((j+1)*binSize).append("\t");
				CommonUtils.replaceEnd(sb, '\n');
				for (int i=0;i<bin2;i++)	
					sb.append((i+1)*binSize).append("\t").append(CommonUtils.arrayToString(map[i])).append("\n");
				CommonUtils.writeFile(outprefix+"."+count+".map.txt", sb.toString());
				count++;
			}
			System.exit(0);
		}

//		// 	to get only a subset of data for testing
//		String test_loops = Args.parseString(args, "subset_loops", null);
//		if (test_loops != null) {
//			ArrayList<StrandedPoint> reads2 = new ArrayList<StrandedPoint>(); 
//			// all PET sorted by the low end
//			ArrayList<ReadPair> low2 = new ArrayList<ReadPair>(); 
//
//			ArrayList<String> lines = CommonUtils.readTextFile(test_loops);
//			for (String anchorString : lines) {
//				// System.out.println(anchorString);
//				String[] f = anchorString.split("\t");
//				Region region1 = new Region(genome, f[0].replace("chr", ""), Integer.parseInt(f[1]),
//						Integer.parseInt(f[2]));
//				Region region2 = new Region(genome, f[3].replace("chr", ""), Integer.parseInt(f[4]),
//						Integer.parseInt(f[5]));
//				if (region1.overlaps(region2)){		// if overlap, merge
//					region1 = new Region(genome, region1.getChrom(), region1.getStart(), region2.getEnd());
//					region2 = region1;
//				}
//				ArrayList<Integer> idx = CommonUtils.getPointsIdxWithinWindow(lowEnds, region1);
//				for (int id : idx) {
//					ReadPair rp = low.get(id);
//					if (region2.contains(rp.r2)){
//						low2.add(rp);
//					}
//				}
//				ArrayList<Integer> idx2 = CommonUtils.getPointsIdxWithinWindow(reads, region1);
//				for (int id : idx2) 
//					reads2.add(reads.get(id));
//				idx2 = CommonUtils.getPointsIdxWithinWindow(reads, region2);
//				for (int id : idx2) 
//					reads2.add(reads.get(id));
//			}
//			
//			// re-set the data structures and sort
//			low.clear();
//			low.addAll(low2);low2 = null;
//			low.trimToSize();
//			high.clear();
//			high.addAll(low);
//			high.trimToSize();
//			Collections.sort(low);
//			Collections.sort(high, new Comparator<ReadPair>() {		// sort by high end read2
//				public int compare(ReadPair o1, ReadPair o2) {
//					return o1.compareRead2(o2);
//				}
//			});
//
//			lowEnds.clear();
//			for (ReadPair r : low)
//				lowEnds.add(r.r1);
//			lowEnds.trimToSize();
//			highEnds.clear();
//			for (ReadPair r : high)
//				highEnds.add(r.r2);
//			highEnds.trimToSize();
//
//			reads.clear();
//			reads.addAll(reads2);
//			reads.trimToSize();
//			Collections.sort(reads);	
//		}	// get subset of data
		
		if (low.isEmpty()){
			System.err.println("\nNo read pair was loaded. Exit here!\n");
			System.exit(-1);
		}
		
//		/** 
//		 * compute self-ligation fraction for minus-plus PETs
//		 */
//		ArrayList<Integer> dist_other = new ArrayList<Integer>();
//		dist_other.addAll(dist_plus_plus);
//		dist_plus_plus = null;
//		dist_other.addAll(dist_plus_minus);
//		dist_plus_minus = null;
//		dist_other.addAll(dist_minus_minus);
//		dist_minus_minus = null;
//		dist_other.trimToSize();
//		Collections.sort(dist_other);
//		dist_minus_plus.trimToSize();
//		Collections.sort(dist_minus_plus);
//
//		int step = dist_other.size() / numQuantile;
//		// [e0 e1), end exclusive
//		ArrayList<Integer> binEdges = new ArrayList<Integer>(); 
//		// first idx for the number that is equal or larger than edge
//		ArrayList<Integer> indexes_other = new ArrayList<Integer>();
//		ArrayList<Integer> indexes_minus_plus = new ArrayList<Integer>(); 
//		for (int i = 0; i <= numQuantile; i++) {
//			int edge = dist_other.get(i * step);
//			binEdges.add(edge);
//			indexes_other.add(CommonUtils.findKey(dist_other, edge));
//			indexes_minus_plus.add(CommonUtils.findKey(dist_minus_plus, edge));
//		}
//		ArrayList<Double> mpNonSelfFraction = new ArrayList<Double>();
//		for (int i = 0; i < binEdges.size() - 1; i++) {
//			double mpNonSelfCount = (indexes_other.get(i + 1) - indexes_other.get(i)) / 3.0;
//			int mpCount = indexes_minus_plus.get(i + 1) - indexes_minus_plus.get(i);
//			double frac = mpNonSelfCount / mpCount;
//			if (frac > 1) {
//				mpNonSelfFraction.add(1.0);
//				break;
//			} else
//				mpNonSelfFraction.add(frac);
//		}
//		indexes_other = null;
//		indexes_minus_plus = null;
//		for (int i = binEdges.size() - 1; i >= mpNonSelfFraction.size(); i--)
//			binEdges.remove(i);
//		int maxEdge = binEdges.get(binEdges.size() - 1);
//
//		// output the distance-fraction table
//		StringBuilder dfsb = new StringBuilder();
//		for (int i = 0; i < binEdges.size(); i++)
//			dfsb.append(binEdges.get(i)).append("\t").append(mpNonSelfFraction.get(i)).append("\n");
//		CommonUtils.writeFile(outName + ".minusPlusFraction.txt", dfsb.toString());

//		System.out.println("\nAnalyzed strand-orientation of PETs: " + CommonUtils.timeElapsed(tic0));

		/*************************************************************************
		 *  run GPS to calling peaks (for improving spatial accuracy)
		 *************************************************************************/
		List<Feature> peaks = null;
		if (flags.contains("gps")) {
			boolean run_gem = false;
		    	if (Args.parseInteger(args,"k", -1)!=-1 || Args.parseInteger(args,"k_min", -1)!=-1 || Args.parseInteger(args,"kmin", -1)!=-1
		    			|| Args.parseString(args, "seed", null)!=null)
		    		run_gem = true;
		    	else
		    		System.err.println("Warning: no options (--k, --k_min & --k_max, or --seed) to run motif discovery. GPS will be run.");

			System.out.println("\nRunning "+(run_gem?"GEM":"GPS")+" on single-end reads: " + CommonUtils.timeElapsed(tic0));

			// prepare single end read data and run GEM			
			KPPMixture mixture = new KPPMixture(genome, prepareGEMData(reads), args);
	        int round = 0;
	        outName = mixture.getOutName();		// get the new path with GEM_output folder
			mixture.setOutName(outName+"_"+round);
//			mixture.plotAllReadDistributions(mixture.getAllModels(), outName+"_"+round);  // for testing
	        peaks = mixture.execute();
	        mixture.printFeatures(round);
	        mixture.printFilteredFeatures(round);
            mixture.printInsignificantFeatures(round);
	        mixture.releaseMemory();
	        round++;
	        boolean not_update_model = false;
	        int minLeft = Args.parseInteger(args,"d_l", 600);
	        int minRight = Args.parseInteger(args,"d_r", 500);
	        boolean constant_model_range = Args.parseFlags(args).contains("constant_model_range");
	        if (!not_update_model){
		        if (!constant_model_range){
		            Pair<Integer, Integer> newEnds = mixture.getModel().getNewEnds(minLeft, minRight);
		            mixture.updateBindingModel(newEnds.car(), newEnds.cdr(), outName+"_"+round);
		        }
		        else
		            mixture.updateBindingModel(-mixture.getModel().getMin(), mixture.getModel().getMax(), outName+"_"+round);
	        }
	        
	        mixture.runKMAC();	
			mixture.setOutName(outName+"_"+round);
	        peaks = mixture.execute();
	        peaks.addAll(mixture.getInsignificantFeatures());
	        mixture.printFeatures(round);
	        mixture.printFilteredFeatures(round);
            mixture.printInsignificantFeatures(round);
	        mixture.releaseMemory();
	        round++;
	        if (!not_update_model){
		        if (!constant_model_range){
		            Pair<Integer, Integer> newEnds = mixture.getModel().getNewEnds(minLeft, minRight);
		            mixture.updateBindingModel(newEnds.car(), newEnds.cdr(), outName+"_"+round);
		        }
		        else
		            mixture.updateBindingModel(-mixture.getModel().getMin(), mixture.getModel().getMax(), outName+"_"+round);
	        }
	        mixture.plotAllReadDistributions(mixture.getAllModels(), outName+"_"+round);
	        System.out.println("\nDone running "+(run_gem?"GEM":"GPS")+": " + CommonUtils.timeElapsed(tic0));
	        outName = Args.parseString(args, "out", "CID");		// set outName back to the original
		}	// if running GPS
		
		
		/***********************************************************************
		 * One dimension read segmentation (similar to GEM code)
		 ***********************************************************************/
		// only consider PETs here, but not single-end-mapped reads, because the goal is to partition PETs
		// compared with reads, pets1d do NOT include cross-chrom reads, single-ended reads, and self-ligation reads.
		// it is used for 1d segmentation and MICC anchor region quantification.
//		ArrayList<Point> pets1d = new ArrayList<Point>();
//		for (ReadPair r : low){			// low and high holds the same data, just sorted differently 
//			pets1d.add(r.r1);
//			pets1d.add(r.r2);
//		}
//		pets1d.trimToSize();
//		Collections.sort(pets1d);

		// TODO: use 1D cross correlation to determine the distance to shift
		ArrayList<Region> rs0 = new ArrayList<Region>();
//		ArrayList<Point> summits = new ArrayList<Point>();
		// cut the pooled reads into independent regions
		int start0 = 0;
		int minCount = 3;
		for (int i = 1; i < lowEnds.size(); i++) {
			Point p0 = lowEnds.get(i-1);
			Point p1 = lowEnds.get(i);
			// not same chorm, or a large enough gap to cut
			if ((!p0.getChrom().equals(p1.getChrom())) || p1.getLocation() - p0.getLocation() > read_1d_merge_dist) { 
				// only select region with read count larger than minimum count
				int count = i - start0;
				if (count >= minCount) {
					Region r = new Region(genome, p0.getChrom(), lowEnds.get(start0).getLocation(),
							lowEnds.get(i - 1).getLocation());
					rs0.add(r);
					ArrayList<Point> ps = new ArrayList<Point>();
					for (int j = start0; j < i; j++)
						ps.add(lowEnds.get(j));
//					int maxCount = 0;
//					int maxIdx = -1;
//					for (int j = 0; j < ps.size(); j++) {
//						Point mid = ps.get(j);
//						int c = CommonUtils.getPointsWithinWindow(ps, mid, read_1d_merge_dist).size();
//						if (c > maxCount) {
//							maxCount = c;
//							maxIdx = start0 + j;
//						}
//					}
//					summits.add(reads.get(maxIdx));
				}
				start0 = i;
			}
		}
		// the last region
		int count = lowEnds.size() - start0;
		if (count >= minCount) {
			Region r = new Region(genome, lowEnds.get(start0).getChrom(), lowEnds.get(start0).getLocation(),
					lowEnds.get(lowEnds.size() - 1).getLocation());
			rs0.add(r);
			ArrayList<Point> ps = new ArrayList<Point>();
			for (int j = start0; j < lowEnds.size(); j++)
				ps.add(lowEnds.get(j));
//			int maxCount = 0;
//			int maxIdx = -1;
//			for (int j = 0; j < ps.size(); j++) {
//				Point mid = ps.get(j);
//				int c = CommonUtils.getPointsWithinWindow(ps, mid, read_1d_merge_dist).size();
//				if (c > maxCount) {
//					maxCount = c;
//					maxIdx = start0 + j;
//				}
//			}
//			summits.add(reads.get(maxIdx));
		}
		
		System.out.println("\nSegment left ends of PETs into " + rs0.size() + " regions, " + CommonUtils.timeElapsed(tic0));

		
		/**************************
		 * Load other data for annotations
		 **************************/
		// load gene annotation
		ArrayList<String> lines = CommonUtils.readTextFile(Args.parseString(args, "gene_anno", null));
		ArrayList<Point> allTSS = new ArrayList<Point>();
		HashMap<StrandedPoint, ArrayList<String>> tss2geneSymbols = new HashMap<StrandedPoint, ArrayList<String>>();
		for (int i = 0; i < lines.size(); i++) {
			String t = lines.get(i);
			if (t.startsWith("#"))
				continue;
			String f[] = t.split("\t");
			String chr = f[2].replace("chr", "");
			char strand = f[3].charAt(0);
			StrandedPoint tss = new StrandedPoint(genome, chr, Integer.parseInt(f[strand == '+' ? 4 : 5]), strand);
			allTSS.add(tss);
			if (!tss2geneSymbols.containsKey(tss))
				tss2geneSymbols.put(tss, new ArrayList<String>());
			tss2geneSymbols.get(tss).add(f[12]);
		}
		allTSS.trimToSize();
		Collections.sort(allTSS);

		
		if (flags.contains("print_cluster")){
			StringBuilder sbDensityDetails = new StringBuilder();
			ReadPairCluster rpc = new ReadPairCluster();
			for (ReadPair rp: high)
				rpc.addReadPair(rp);
			rpc.update(false);
			sbDensityDetails.append(String.format("# %s:%d-%d\n", rpc.leftRegion.getChrom(), rpc.r1min, rpc.r2max));
			densityClustering(rpc, span2mergingDist(rpc.getLoopRegionWidth()), sbDensityDetails);
			CommonUtils.writeFile(String.format("%s.cluster.density.txt", outName), sbDensityDetails.toString());
			System.out.println(String.format("\nClustering results have been written to %s.cluster.density.txt", outName));
			System.exit(0);
		}		

		
		/***********************************************************
		 * find dense PET cluster for each 1D clustered region
		 ***********************************************************/
		ArrayList<Interaction> interactions = new ArrayList<Interaction>();
		HashSet<ReadPair> usedPETs = new HashSet<ReadPair>();
		ArrayList<ReadPairCluster> clustersCalled = new ArrayList<ReadPairCluster>();
		
		for (int j = 0; j < rs0.size(); j++) { // for all regions
			Region region = rs0.get(j);
//			if (region.contains(new Point(region.getGenome(), "13", 23562420)))	//13:23562420-23564665
//				j +=0;
//			
			// get the PETs with read1 in the region, sort and merge by read2
			ArrayList<Integer> idx = CommonUtils.getPointsIdxWithinWindow(lowEnds, region);
			if (idx.size() > 1) {
				ArrayList<ReadPair> rps = new ArrayList<ReadPair>();
				for (int i : idx) {
					rps.add(low.get(i));
				}
				Collections.sort(rps, new Comparator<ReadPair>() {
					public int compare(ReadPair o1, ReadPair o2) {
						return o1.compareRead2(o2);
					}
				});
				ArrayList<ReadPairCluster> rpcs = new ArrayList<ReadPairCluster>();
				int current = -100000;
				ReadPairCluster c = new ReadPairCluster();
				for (ReadPair rp : rps) {
					// a big gap
					if (rp.r2.getLocation() - current > max_cluster_merge_dist) {  // merge all possible PETs
						if (c.pets.size() >= min_pet_count) {
							c.update(false);
							rpcs.add(c);
						}
						c = new ReadPairCluster();
					}
					c.addReadPair(rp);
					current = rp.r2.getLocation();
				}
				if (c.pets.size() >= min_pet_count) { // finish up the last cluster
					c.update(false);
					rpcs.add(c);
				}
				if (rpcs.isEmpty())
					continue;
				
//				ArrayList<ReadPairCluster> rpcs2 = splitRecursively(rpcs, true, true, true);
//				if (rpcs2 != null && !rpcs2.isEmpty()) {
//					rpcs = rpcs2;
//					rpcs2 = null;
//				}
				ArrayList<ReadPairCluster> rpcs2 = new ArrayList<ReadPairCluster>();
				for (ReadPairCluster rpc: rpcs){
					ArrayList<ReadPairCluster> splited = splitRecursively(rpc, true, true);
					if (!splited.isEmpty())
						rpcs2.addAll(splited);
				}
				if (rpcs2 != null && !rpcs2.isEmpty()) {
					rpcs = rpcs2;
					rpcs2 = null;
				}
				else
					continue;
				// after the recursive splitting, all the rpcs should be 3kb away from each other, 
				// with each PET having at least 1 PET within less than 3kb
				// therefore, the density clustering and merging operate on independent rpcs
				
				// density clustering of PETs 
				int r1min=Integer.MAX_VALUE;
				int r2max=0;
				for (ReadPairCluster rpc: rpcs){
					if (r1min>rpc.r1min)
						r1min = rpc.r1min;
					if (r2max<rpc.r2max)
						r2max = rpc.r2max;
				}
				StringBuilder sbDensityDetails = null;
				if (flags.contains("print_cluster")) {
					sbDensityDetails = new StringBuilder();
					sbDensityDetails.append(String.format("# %s:%d-%d\n", region.getChrom(), r1min, r2max));
				}
				// Density clustering, start with largest possible span, i.e. the whole width
				ArrayList<ReadPairCluster> tmp = new ArrayList<ReadPairCluster>();
				for (ReadPairCluster rpc: rpcs){
					int span = rpc.getLoopRegionWidth();
					int dc = this.span2mergingDist(span);
					// TODO: remove singleton using dc (if it is not the max)
					if (dc<this.max_cluster_merge_dist){
						
					}
					tmp.addAll(densityClustering(rpc, span2mergingDist(span), sbDensityDetails));
				}
				rpcs = tmp;
				if (rpcs.isEmpty())
					continue;
				if (flags.contains("print_cluster")) 
					CommonUtils.writeFile(String.format("%s.cluster.%d.txt", outName, j), sbDensityDetails.toString());
				
				// merge nearby clusters
				Collections.sort(rpcs);
				for (int i = 0; i < rpcs.size(); i++) {
					ReadPairCluster c1 = rpcs.get(i);
//					if (c1.leftRegion.overlaps(new Region(genome, "13", 23562420, 23564665)))	//13:23562420-23564665
//						i+=0;
					ArrayList<ReadPairCluster> toRemoveClusters = new ArrayList<ReadPairCluster>();
					int c1Span = c1.span;
					if (c1.r1width*span_anchor_ratio>c1Span || c1.r2width*span_anchor_ratio>c1Span)
						continue;	// if c1 anchors are too wide, skip merging c1
					for (int jj = i+1; jj < rpcs.size(); jj++) {
						ReadPairCluster c2 = rpcs.get(jj);
						int c2Span = c2.span;
						if (c2.r1width*span_anchor_ratio>c2Span || c2.r2width*span_anchor_ratio>c2Span)
							continue;	// if c2 anchors are too wide, skip merging c2
						// cluster_merge_dist is dependent on the PET span distance
						int dist = Math.min(c1Span, c2Span);
						if ((Math.max(c1.r1max, c2.r1max)-Math.min(c1.r1min, c2.r1min))*span_anchor_ratio>dist || 
								(Math.max(c1.r2max, c2.r2max)-Math.min(c1.r2min, c2.r2min))*span_anchor_ratio>dist )
							continue;	// if merged anchors are too wide, skip merging
						int cluster_merge_dist = span2mergingDist(dist);
						if (c1.leftPoint.distance(c2.leftPoint) > cluster_merge_dist*2 || 
								c1.rightPoint.distance(c2.rightPoint) > cluster_merge_dist*2 ||
								c1.leftRegion.distance(c2.leftRegion) > cluster_merge_dist ||
								c1.rightRegion.distance(c2.rightRegion) > cluster_merge_dist) 
							continue;
						// if close enough, simply merge c2 to c1
						toRemoveClusters.add(c2);
						boolean toSetAnchorPoints = false;
						Point left = null; 
						Point right  = null;
						if (c1.pets.size()==c2.pets.size()){		// case 1
							toSetAnchorPoints = true;
						}
						else if (c1.pets.size()<c2.pets.size()){	// case 2
							left = c2.leftPoint;
							right = c2.rightPoint;
						} 
//						else if (c1.pets.size()>c2.pets.size()) : 	// case 3: Do nothing
							
						for (ReadPair rp2 : c2.pets)
							c1.addReadPair(rp2);
						c1.update(toSetAnchorPoints);
						
						if (!toSetAnchorPoints && left!=null){		// only with case 2, set anchors to c2
							c1.leftPoint = left;
							c1.rightPoint = right;
							c1.span = left.distance(right);
						}
							
						c1Span = c1.span;
//						if (isDev)
//							System.err.println(String.format("Merged %s with %s to %s - %s.", c1_old, c2.toString(), c1.toString(), 
//									c1.getLoopRegionString(2000)));
						if (c1.r1width*span_anchor_ratio>c1Span || c1.r2width*span_anchor_ratio>c1Span)
							break;		// if c1 anchors are too wide, stop merging c1
					} // for each pair of nearby clusters
					if (!toRemoveClusters.isEmpty()){
						rpcs.removeAll(toRemoveClusters);
						toRemoveClusters.clear();
					}
					toRemoveClusters = null;
				}
				clustersCalled.addAll(rpcs);
			} // all regions with long PETs
		} // loop over all regions


		// refresh the PETs again because some PET1 might not be
		// included but are within the cluster_merge_dist range
		HashSet<ReadPair> unusedPET2 = new HashSet<ReadPair>();
		unusedPET2.addAll(low);
		for (ReadPairCluster cc : clustersCalled) {
			unusedPET2.removeAll(cc.pets);
		}
		
		for (ReadPairCluster cc : clustersCalled) {
//			if (cc.leftRegion.contains(new Point(genome, "13", 23562486)) && cc.rightRegion.contains(new Point(genome, "13", 23746753))) {	//13:23562486 13:23746753
//				System.out.println(cc.pets.size());
//				System.out.println(cc.leftRegion.toString()+cc.rightRegion.toString());
//				for (ReadPair rp: cc.pets) 
//					System.out.println(rp.toString());
//				dc+=0;
//			}
//			int c = cc.pets.size();
			ArrayList<Integer> idx2 = CommonUtils.getPointsIdxWithinWindow(lowEnds, cc.leftRegion);
//			cc.pets.clear();
			for (int i : idx2) {
				ReadPair rp = low.get(i);
				if (cc.rightRegion.contains(rp.r2) && unusedPET2.contains(rp))
					cc.addReadPair(rp);
			}
			cc.update(false);
//			if (cc.pets.size()!=c)
//				System.out.println(String.format("%d --> %d: %s -- %s", c, cc.pets.size(), cc.leftRegion.toString(), cc.rightRegion.toString()));
//			if (cc.leftRegion.contains(new Point(genome, "13", 23562486)) && cc.rightRegion.contains(new Point(genome, "13", 23746753))) {	//13:23562486 13:23746753
//				System.out.println(cc.pets.size());
//				System.out.println(cc.leftRegion.toString()+cc.rightRegion.toString());
//				for (ReadPair rp: cc.pets) 
//					System.out.println(rp.toString());
//				dc+=0;
//			}
		}
		
//		ArrayList<Region> rs = new ArrayList<Region>();
//		for (ReadPairCluster cc : clustersCalled) {
//			rs.add(cc.leftRegion);
//			rs.add(cc.rightRegion);
//		}
//		System.out.println(rs.size());
//		rs = Region.mergeRegions(rs);
//		System.out.println(rs.size());
//		int ll=0;
//		for (Region r:rs)
//			ll+=r.getWidth();
//		System.out.println(genome.getChromLength(rs.get(0).getChrom()));
//		System.out.println(ll);
		
		// Get GPS peak calls
		int interPeakDistance = Args.parseInteger(args, "ipd", 1000);
		ArrayList<Point> hits = new ArrayList<Point>();
		if (peaks!=null) {		// have run GPS
			List<GPSPeak> ps = new ArrayList<GPSPeak>();
			for (Feature f: peaks) {
				ComponentFeature cf = (ComponentFeature)f;
				Point p = cf.getPeak();
				ps.add(new GPSPeak(p.getGenome(), p.getChrom(), p.getLocation(), cf.getTotalEventStrength()));
			}
			Collections.sort(ps);
			if (ps.size()>1) {
				ArrayList<GPSPeak> toRemove = new ArrayList<GPSPeak>();
				for (int i=0;i<ps.size();i++) {
					GPSPeak p1 = ps.get(i);
					for (int j=i+1;j<ps.size();j++) {
						GPSPeak p2 = ps.get(j);
						if (p1.getChrom().equals(p2.getChrom()) && p1.distance(p2)<interPeakDistance) {	
							// too close (arbitrary interPeakDistance bp), merge to the stronger peak
							if (p1.getStrength()<p2.getStrength()) {
								toRemove.add(p1);
								break;		// move to next p1
							}
							else {
								toRemove.add(p2);
							}
						}
						else {
							i=j-1;		// j will be the next p1
							break;
						}
					}
				}
				ps.removeAll(toRemove);
			}
			hits.addAll(ps);
		}
		else {		// NOT run GPS, load GPS calls, for easy testing
			String fileString = Args.parseString(args, "GPS", null);
			if (fileString!=null) {
				File gpsFile = new File(fileString);
				try {
					List<GPSPeak> ps = GPSParser.parseGPSOutput(gpsFile.getAbsolutePath(), genome);
					Collections.sort(ps);
					if (ps.size()>1) {
						ArrayList<GPSPeak> toRemove = new ArrayList<GPSPeak>();
						for (int i=0;i<ps.size();i++) {
							GPSPeak p1 = ps.get(i);
							for (int j=i+1;j<ps.size();j++) {
								GPSPeak p2 = ps.get(j);
								if (p1.getChrom().equals(p2.getChrom()) && p1.distance(p2)<interPeakDistance) {	// too close (arbitrary interPeakDistance bp), merge to the stronger peak
									if (p1.getStrength()<p2.getStrength()) {
										toRemove.add(p1);
										break;		// move to next p1
									}
									else {
										toRemove.add(p2);
									}
								}
								else {
									i=j-1;		// j will be the next p1
									break;
								}
							}
						}
						ps.removeAll(toRemove);
					}
					hits.addAll(ps);
				}
				catch(Exception e) {}
			}
		}
		
		// re-assign PETs 
		ArrayList<ReadPairCluster> clustersAssigned = new ArrayList<ReadPairCluster>();
		for (ReadPairCluster cc : clustersCalled) {
//			if (cc.leftRegion.overlaps(new Region(genome, "13", 23562420, 23564665)))	//13:23562420-23564665
//				dc+=0;
			int mergeDist = span2mergingDist(cc.span);
			ArrayList<Integer> id1 = CommonUtils.getPointsIdxWithinWindow(hits, cc.leftRegion);
			ArrayList<Integer> id2 = CommonUtils.getPointsIdxWithinWindow(hits, cc.rightRegion);
			ArrayList<Pair<Point,Point>> peakPairs = new ArrayList<Pair<Point,Point>>();
			for (int i1: id1) {
				for (int i2: id2) {
					Point p1 = hits.get(i1);
					Point p2 = hits.get(i2);
					int offset = p2.getLocation()-p1.getLocation();
					if (Math.abs(offset)<this.min_span)
						continue;
					if (offset<0) {		// make sure p1 < p2
						Point tmp = p1;
						p1 = p2;
						p2 = tmp;
					}
					peakPairs.add(new Pair<Point,Point>(p1,p2));
				}
			}
			if (peakPairs.isEmpty()){			// No peak found, just keep the old one
				clustersAssigned.add(cc);
				continue;
			}
				
//			System.out.println("\n"+cc.toString()+"\t"+cc.getLoopRegionString(cc.span/20)+"\t"+cc.leftRegion.toString()+"\t"+cc.rightRegion.toString());
			ReadPairCluster[] rpcs = new ReadPairCluster[peakPairs.size()];
//			System.out.print("    Peaks #="+peakPairs.size()+": ");
			for (int i=0;i<peakPairs.size();i++) {
//				System.out.print(String.format("[%d] %s; ", i, peakPairs.get(i).toString()));
				Pair<Point, Point> p = peakPairs.get(i);
				if (p.car().getLocation()==23562486 && p.cdr().getLocation()==23746753)
					dc+=0;
				rpcs[i] = new ReadPairCluster();
			}
//			System.out.println();
			// initial assignment
			for (ReadPair rp: cc.pets) {
//				System.out.print("\t"+ rp.toString());
				int min = Integer.MAX_VALUE;
				int minId = -1;
				for (int i=0;i<peakPairs.size();i++) {
					Pair<Point,Point> pp = peakPairs.get(i);
//					System.out.print(pp.toString());
					int dist  = Math.max(pp.car().distance(rp.r1), pp.cdr().distance(rp.r2));		// Chebyshev distance
////					System.out.print(" "+dist+"| ");
					if (dist<min) {		//  && dist<mergeDist
						minId = i;
						min = dist;
					}
				}
//				System.out.println(" --> "+min+"["+minId+"]");
				if (minId>-1)
					rpcs[minId].addReadPair(rp);
			}
			// secondary assignment
			for (int i=0;i<rpcs.length;i++) {
				int petCount = rpcs[i].pets.size();
				if (petCount==0) 
					continue;
				if (petCount==1) {
					ReadPair rp = rpcs[i].pets.get(0);
					int min = Integer.MAX_VALUE;
					int minId = -1;
					for (int j=0;j<peakPairs.size();j++) {
						if (i==j || rpcs[j].pets.isEmpty())	// exclude itself and empty cluster centers
							continue;
						Pair<Point,Point> pp = peakPairs.get(j);
						int dist  = Math.max(pp.car().distance(rp.r1), pp.cdr().distance(rp.r2));		// Chebyshev distance
						if (dist<min) {		//  && dist<mergeDist
							minId = j;
							min = dist;
						}
					}
					if (minId>-1)
						rpcs[minId].addReadPair(rp);
					rpcs[i].pets.clear();
				}
			}
			// update clusters
			int rpcCount=0;
			for (int i=0;i<rpcs.length;i++) {
				int petCount = rpcs[i].pets.size();
				if (petCount<2)
					continue;
				rpcCount++;
				rpcs[i].update(false);
				// set the peak positions as the anchor points
				rpcs[i].leftPoint = peakPairs.get(i).car();
				rpcs[i].rightPoint = peakPairs.get(i).cdr();
				rpcs[i].d_c = cc.d_c;
				rpcs[i].span = rpcs[i].leftPoint.distance(rpcs[i].rightPoint);
			}
			// merge RP clusters
			if (rpcCount>1) {
//				dc += 0;
				for (int i=0;i<rpcs.length;i++) {
					ReadPairCluster rpci = rpcs[i];
					for (int j=i+1;j<rpcs.length;j++) {
						ReadPairCluster rpcj = rpcs[j];
						if (rpci.pets.isEmpty() || rpcj.pets.isEmpty())
							continue;
						// if nearby, or PET2 (if not merged, PET2 will be deleted), merge
						if (rpci.distance(rpcj)<mergeDist || 
								(rpci.pets.size()==2 && rpci.span>dc) || 
								(rpcj.pets.size()==2 && rpcj.span>dc)) {	
							ReadPairCluster donor=null, receiver=null;
							if (rpci.pets.size()>rpcj.pets.size()) {
								donor = rpcj;
								receiver = rpci;
							}
							else if (rpci.pets.size()<rpcj.pets.size()){
								donor = rpci;
								receiver = rpcj;
							} else {		// equal size, take the stronger peak-pair
								double iStrength = ((GPSPeak)peakPairs.get(i).car()).getStrength()+((GPSPeak)peakPairs.get(i).cdr()).getStrength();
								double jStrength = ((GPSPeak)peakPairs.get(j).car()).getStrength()+((GPSPeak)peakPairs.get(j).cdr()).getStrength();
								if (iStrength>=jStrength){
									donor = rpcj;
									receiver = rpci;
								}
								else{
									donor = rpci;
									receiver = rpcj;
								}
							}
							for(ReadPair rp: donor.pets)
								receiver.addReadPair(rp);
							receiver.update(false);
							donor.pets.clear();

						}
					}
				}
//				dc += 0;
			}
			for (int i=0;i<rpcs.length;i++) {
				if (rpcs[i].pets.size()>=2) {
					clustersAssigned.add(rpcs[i]);
//					System.out.println("    "+ rpcs[i].toString());
				}
			}
//			dc += 0;
		}
		clustersCalled.clear();
		clustersCalled.addAll(clustersAssigned);
		
		// finalize PET clusters
		usedPETs.clear();
		for (ReadPairCluster cc : clustersCalled) {
			ArrayList<ReadPair> pets = cc.pets;
			if (pets.size() < min_pet_count)
				continue;
			// skip PET2 that have anchor wider than dc, or 1D read count is also 2 (not significant by MICC)
			if (pets.size()==2){
				if (cc.r1width>dc && cc.r2width>dc)	
					continue;
				if (flags.contains("strict")){		// more stringent calls
					int l = CommonUtils.getPointsIdxWithinWindow(lowEnds, cc.leftRegion.expand(dc, dc)).size();
					int r = CommonUtils.getPointsIdxWithinWindow(highEnds, cc.rightRegion.expand(dc, dc)).size();
					if (l==2 || r==2 || (l<=4 && r<=4))
						continue;
				}
			}

			// mark all PETs in the cluster as used (PET2+)
			// to get real PET1 (no PET1 from the m-p adjustment)
			usedPETs.addAll(pets);	
			
			int totalCount = pets.size();
			int minusPlusCount = 0;
			int adjustedCount = totalCount;
			ReadPairCluster rpc = cc;
			
//			// new PET cluster with adjustment
//			ReadPairCluster rpc = new ReadPairCluster(); 
//			if (flags.contains("mp_adjust")){
//				// count minus-plus PETs to adjust the PET counts
//				ArrayList<ReadPair> mpRPs = new ArrayList<ReadPair>();
//				for (ReadPair rp : pets)
//					if (rp.r1.getStrand() == '-' && rp.r2.getStrand() == '+')
//						mpRPs.add(rp);
//				minusPlusCount = mpRPs.size();
//				pets.removeAll(mpRPs);
//				adjustedCount = -1;
//				Collections.sort(mpRPs);		// sort by low end read1: default
//				if (pets.isEmpty()) { //  with only minus-plus PETs
//					int dist = cc.span;
//					if (dist >= maxEdge)
//						adjustedCount = minusPlusCount;
//					else {
//						int index = Collections.binarySearch(binEdges, dist);
//						if (index < 0) // if key not found
//							index = -(index + 1);
//						adjustedCount = (int) (minusPlusCount * mpNonSelfFraction.get(index));
//					}
//					// add the adjusted m-p PETs in the middle
//					int midIndexMP = mpRPs.size() / 2 - adjustedCount /2;
//					int endIndex = midIndexMP + adjustedCount;
//					for (int k = midIndexMP; k < endIndex; k++)
//						rpc.addReadPair(mpRPs.get(k));
//					rpc.update();
//				} else {
//					for (ReadPair rp : pets)
//						rpc.addReadPair(rp);
//					rpc.update();
//					int dist = rpc.span;
//					if (dist >= maxEdge)
//						adjustedCount = totalCount;
//					else {
//						int index = Collections.binarySearch(binEdges, dist);
//						if (index < 0) // if key not found
//							index = -(index + 1);
//						adjustedCount = totalCount - minusPlusCount
//								+ (int) (minusPlusCount * mpNonSelfFraction.get(index));
//					}
//					// add the adjusted m-p PETs in the middle
//					int extra = adjustedCount - (totalCount - minusPlusCount);
//					int midIndexMP = mpRPs.size() / 2 - extra /2;
//					int endIndex = midIndexMP + extra;
//					for (int k = midIndexMP; k < endIndex; k++)
//						rpc.addReadPair(mpRPs.get(k));
//					rpc.update();
//				}
//				if (adjustedCount < min)
//					continue;
//			}
//			else	// not adjustment, rpc is the cluster cc
//				rpc = cc;
								
			Interaction it = new Interaction();
			interactions.add(it);
			it.count = totalCount;
			it.count2 = totalCount - minusPlusCount;
			it.adjustedCount = adjustedCount;
			it.d_c = rpc.d_c;

			// add gene annotations
			pets = rpc.pets;
			it.leftRegion = rpc.leftRegion;
			it.leftPoint = rpc.leftPoint;
			ArrayList<Integer> ts = CommonUtils.getPointsIdxWithinWindow(allTSS,
					it.leftRegion.expand(tss_radius, tss_radius));
			StringBuilder tsb = new StringBuilder();
			TreeSet<String> gSymbols = new TreeSet<String>();
			for (int i : ts) {
				gSymbols.addAll(tss2geneSymbols.get(allTSS.get(i)));
			}
			for (String s : gSymbols)
				tsb.append(s).append(",");
			if (tsb.length() == 0)
				it.leftLabel = "nonTSS";
			else
				it.leftLabel = tsb.toString();

			it.rightRegion = rpc.rightRegion;
			it.rightPoint = rpc.rightPoint;
			ts = CommonUtils.getPointsIdxWithinWindow(allTSS, it.rightRegion.expand(tss_radius, tss_radius));
			tsb = new StringBuilder();
			gSymbols.clear();
			for (int i : ts) {
				gSymbols.addAll(tss2geneSymbols.get(allTSS.get(i)));
			}
			for (String s : gSymbols)
				tsb.append(s).append(",");
			if (tsb.length() == 0)
				it.rightLabel = "nonTSS";
			else
				it.rightLabel = tsb.toString();
		}  // all clusters
		
		interactions.trimToSize();
		Collections.sort(interactions);		
		System.out.println("\nCalled " + interactions.size() + " PET clusters, " + CommonUtils.timeElapsed(tic0));

		// mark PET1 (after removing the used PET2+)
		low.removeAll(usedPETs);
		low.trimToSize();
		high.clear();
		high = null;
		System.out.println("\nClustered PETs n=" + usedPETs.size() + "\nSingle PETs n=" + low.size());

		/******************************
		 * Annotate and report
		 *******************************/
		annotateInteractionCallsAndOutput(interactions, lowEnds, highEnds, low, tic0);
		
		System.out.println("\nDone: " + CommonUtils.timeElapsed(tic0));
	}
	
	/** Density Clustering, start with a rough clustering to reduce to smaller clusters, then refine with more accurate span distances */	
	private ArrayList<ReadPairCluster> densityClustering(ReadPairCluster cc, int d_c, StringBuilder sb){
		ArrayList<ReadPairCluster> results = new ArrayList<ReadPairCluster>();
		ArrayList<ReadPair> pets = cc.pets;
		if (d_c==-1){
			d_c = span2mergingDist(cc.span);
			if (cc.r1width<d_c && cc.r2width<d_c){
				cc.d_c = d_c;
				results.add(cc);
				return results;
			}
		}
		
//		Region left = Region.fromString(genome, "17:41380000-41385000");
//		Region right = Region.fromString(genome, "17:41443000-41447000");
////		if (cc.leftRegion.overlaps(left) && cc.rightRegion.overlaps(right))
////			dc += 0;
		int count = pets.size();
		long tic=-1;
		if (count>5000 && isDev){
			System.err.print(String.format("Warning: #PETs=%s; %s; used ", cc.toString(), 
					cc.getLoopRegionString(cc.getLoopRegionWidth()/20)));
			tic = System.currentTimeMillis();
		}
		int s = cc.r1min;
		// Binning the data if there are too many data points. Now it works fine for each point in its own bin.
		PetBin bins[] = new PetBin[count];
		for (int i=0;i<count;i++){
			PetBin bin = new PetBin(pets.get(i).r1.getLocation()-s, pets.get(i).r2.getLocation()-s, i);
//			if (bin.r2==41444438-s)
//				dc += 0;
			bin.addId(i);
			bins[i] = bin;
		}
		
//			System.out.println(String.format("Region pair %d-%d and %d-%d, span=%d, dc=%d", cc.r1min, cc.r1max, cc.r2min, cc.r2max, span, distance_cutoff));
		// distance matrix
		int[][] dist = new int[count][count];
		int maxDist = 0;
		for (int i=0;i<count;i++){
			dist[i][i] = 0;
			for (int j=i+1;j<count;j++){
				dist[i][j] = Math.max(Math.abs(bins[i].r1-bins[j].r1), Math.abs(bins[i].r2-bins[j].r2)); // Chebyshev distance
				dist[j][i] = dist[i][j];
				if (maxDist<dist[i][j])
					maxDist = dist[i][j];
			}
		}
		// density
		for (int i=0;i<count;i++){
			PetBin b = bins[i];
			for (int j=0;j<count;j++)
				if (dist[i][j] <= d_c)
					b.addNbBin(bins[j]);
			b.density = b.nbBins.size();
		}
		Arrays.sort(bins);
		
		//if tie for highest density, find the median position as center
		if (bins[0].density==bins[1].density){
			ArrayList<PetBin> topBins = new ArrayList<PetBin>();
			ArrayList<Integer> r1s = new ArrayList<Integer>();
			ArrayList<Integer> r2s = new ArrayList<Integer>();
			int densityTop = bins[0].density;
			for (int i=0;i<bins.length;i++){
				PetBin pb = bins[i];
				if (pb.density == densityTop){
					topBins.add(pb);
					r1s.add(pb.r1);
					r2s.add(pb.r2);
				}
				else
					break;
			}
			Collections.sort(r1s);
			Collections.sort(r2s);
			int median1 = r1s.get(r1s.size()/2);
			int median2 = r2s.get(r2s.size()/2);
			int minDist = Integer.MAX_VALUE;
			int minIdx = -1;
			for (int i=0;i<topBins.size();i++){
				PetBin pb = topBins.get(i);
				int distance = Math.abs(pb.r1-median1)+Math.abs(pb.r2-median2);
				if (distance<minDist){
					minDist = distance;
					minIdx = i;
				}
			}
			PetBin tmp = topBins.get(0);
			bins[0] = bins[minIdx];		// move the most median point to the top as cluster center
			bins[minIdx] = tmp;
		}
		
		// delta and gamma
		bins[0].delta = maxDist;
		bins[0].gamma = bins[0].delta * bins[0].density;
		for (int i=1;i<count;i++){		// idx of sorted bin
			int minDist = maxDist;
			int ii = bins[i].binId;		// original idx of this bin
			for (int j=0;j<i;j++){		// for points with higher density
				int jj = bins[j].binId;
				if (minDist>dist[ii][jj])
					minDist = dist[ii][jj];
			}
			bins[i].delta = minDist;
			bins[i].gamma = bins[i].delta * bins[i].density;
		}
		
		// sort by gamma
		Arrays.sort(bins, new Comparator<PetBin>() {		// sort by descending gamma
			public int compare(PetBin o1, PetBin o2) {
				return o1.compareByGamma(o2);
			}
		});
		
		// assign clusters
		ArrayList<PetBin> centers = new ArrayList<PetBin>();
		ArrayList<PetBin> singletons = new ArrayList<PetBin>();
		for (int i=0;i<count;i++){
			PetBin b = bins[i];
			if (b.clusterBinId == -1) {		// unassigned
				if (b.density<=1){
					singletons.add(b);
//					System.err.println("Density<=1");
					continue;
				}
				// find all higher-density neighbor points (within d_c or other distance)
				ArrayList<PetBin> higherDensityNeighbors = new ArrayList<PetBin>();
				for (int j=0; j<count;j++){		
					PetBin hb = bins[j];
					// the neighbor test below can be replaced with other distance criteria
					// it is a neighbor, higher density, and assigned
					if (b.nbBins.contains(hb) && hb.clusterBinId!=-1 && 
							(hb.density>b.density || (hb.density==b.density && hb.delta>=b.delta)) )	
						higherDensityNeighbors.add(hb);
				}
				// find the nearest higher-density point
				int clusterId = -1;
				int minDist = Integer.MAX_VALUE;
				int bid = b.binId;
				for (PetBin hnb:higherDensityNeighbors){
					if (dist[bid][hnb.binId] < minDist){
						minDist = dist[bid][hnb.binId];
						clusterId = hnb.clusterBinId;
					}
				}
				// if the nearest higher-density point (within d_c or other distance) has been assigned,
		        // follow the same assignment, otherwise, start a new cluster
				if (clusterId != -1)	// nearest higher point has been assigned
					b.clusterBinId = clusterId;
				else{
					// this is a new cluster center
					b.clusterBinId = b.binId;
					centers.add(b);
//						for (PetBin nb:b.nbBins)
//							if (nb.clusterBinId == -1)		// if the cluster neighbors has not been assigned
//								nb.clusterBinId = b.binId;
				}
			}
		} // assign clusters

		// re-assign singletons (cluster with only 1 point)
		eachCenter: for (int j=0;j<centers.size();j++){
			PetBin b = centers.get(j);
			int bid = b.binId;
			if (bid != b.clusterBinId)
				continue;
			int clusterSize = 0;
			for (int i=0;i<count;i++){
				if (bid == bins[i].clusterBinId){
					clusterSize++;
					if (clusterSize>=2)
						continue eachCenter;	// skip if more than 2 PETs
				}
			}
			
			// b is a singleton
			// assign it to the nearest higher-density point within dc
			int clusterId = -1;
			int minDist = d_c;
//			for (int i=0;i<count;i++){
//				if (bid!=bins[i].binId){
//					int dist_ji = dist[bid][bins[i].binId];
//					if ( dist_ji < minDist){
//						minDist = dist_ji;
//						clusterId = bins[i].clusterBinId;
//					}
//				}
//			}
			for (PetBin nb: b.nbBins){	// check neighbors of b, assign b to same cluster as its NN
				if (bid!=nb.binId){
					int dist_ji = dist[bid][nb.binId];
					if ( dist_ji < minDist){
						minDist = dist_ji;
						clusterId = nb.clusterBinId;
					}
				}
			}
			if (clusterId != -1)
				b.clusterBinId = clusterId;
			else
				singletons.add(b);
		}
		centers.removeAll(singletons);
		
		if (sb!=null)
			sb.append("Read1\tRead2\tClusterId\tDensity\tDelta\tGamma\tPetId\tClusterCenterId\tDist\n");
		for (int j=0;j<centers.size();j++){
			PetBin b = centers.get(j);
			int bid = b.binId;
			ReadPairCluster rpc = new ReadPairCluster();
			for (int i=0;i<count;i++){
				PetBin m = bins[i];
				if (m.clusterBinId == bid){		// cluster member
					ReadPair rp = pets.get(m.binId);
					rpc.addReadPair(rp);
					if (bid==m.binId){		// cluster center
						rpc.leftPoint = rp.r1;
						rpc.rightPoint = rp.r2;
					}
					if (sb!=null)
						sb.append(String.format("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", rp.r1.getLocation(), rp.r2.getLocation(), j+1,
							m.density,m.delta, m.gamma, m.binId,m.clusterBinId, d_c));
				}
			}
			if (rpc.pets.size()>=2){
				rpc.update(true);
				rpc.d_c = d_c;
				results.add(rpc);
			}
		}
		if (sb!=null){		// print out singletons
			for (PetBin m : singletons){
				ReadPair rp = pets.get(m.binId);
				sb.append(String.format("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", rp.r1.getLocation(), rp.r2.getLocation(), -1,
						m.density,m.delta, m.gamma, m.binId,-1, d_c));
			}
		}
		if (tic!=-1)
			System.err.println(CommonUtils.timeElapsed(tic));
		
		// fine-tune clusters based on the updated span and anchor width
		ArrayList<ReadPairCluster> newResults = new ArrayList<ReadPairCluster>();
		for (ReadPairCluster rpc: results){
			int dcNew = span2mergingDist(rpc.span);
			if (d_c > dcNew)		// fine tune RPCs with their own d_c
				newResults.addAll(densityClustering(rpc, dcNew, sb));
			// if anchor width is too large and span is short, reduce d_c
			else if (rpc.span < span_anchor_ratio * Math.max(rpc.r1width, rpc.r2width) && rpc.span<100000)
				newResults.addAll(densityClustering(rpc, (int)(d_c*0.75), sb));
			else
				newResults.add(rpc);
		}
		return newResults;
	}


	/**
	 * Annotate after CID interaction calling and then print output files
	 */
	private void annotateInteractionCallsAndOutput(ArrayList<Interaction> interactions, ArrayList<Point> lowEnds, ArrayList<Point> highEnds, 
			ArrayList<ReadPair> singleton_pets, long tic0){

		// load TF sites
		String tfs_file = Args.parseString(args, "tf_sites", null);
		ArrayList<ArrayList<Point>> allPeaks = new ArrayList<ArrayList<Point>>();
		if (tfs_file != null) {
			ArrayList<String> tfs = CommonUtils.readTextFile(tfs_file);
			for (int i = 0; i < tfs.size(); i++) {
				try {
					ArrayList<Point> ps = new ArrayList<Point>();
					ps.addAll(GPSParser.parseGPSOutput(tfs.get(i), genome));
					ps.trimToSize();
					Collections.sort(ps);
					allPeaks.add(ps);
					System.out.println("Loaded " + tfs.get(i));
				} catch (IOException e) {
					System.out.println(tfs.get(i) + " does not have a valid GPS/GEM event call file.");
					e.printStackTrace(System.err);
					System.exit(1);
				}
			}
			allPeaks.trimToSize();
		}

		// load histone mark or DHS, SE regions
		System.out.println();
		String hms_file = Args.parseString(args, "regions", null);
		ArrayList<List<Region>> allRegions = new ArrayList<List<Region>>();
		if (hms_file != null) {
			ArrayList<String> hms = CommonUtils.readTextFile(hms_file);
			for (int i = 0; i < hms.size(); i++) {
				allRegions.add(CommonUtils.load_BED_regions(genome, hms.get(i)).car());
				System.out.println("Loaded " + hms.get(i));
			}
			allRegions.trimToSize();
		}

		// load other Interaction calls
		String germ_file = Args.parseString(args, "germ", null);
		ArrayList<Point> tPoints = new ArrayList<Point>();
		HashMap<Point, ArrayList<Point>> t2ds = new HashMap<Point, ArrayList<Point>>();
		if (germ_file != null) {
			ArrayList<String> lines = CommonUtils.readTextFile(germ_file);
			for (String l : lines) { // each line is a call
				String f[] = l.split("\t");
				Point t = new Region(genome, f[3].replace("chr", ""), Integer.parseInt(f[4]), Integer.parseInt(f[5]))
						.getMidpoint();
				Point d = new Region(genome, f[0].replace("chr", ""), Integer.parseInt(f[1]), Integer.parseInt(f[2]))
						.getMidpoint();
				if (t.getLocation() > d.getLocation()) { // make sure t < d
					Point tmp = t;
					t = d;
					d = tmp;
				}
				if (!t2ds.containsKey(t))
					t2ds.put(t, new ArrayList<Point>());
				t2ds.get(t).add(d);
			}
			tPoints.addAll(t2ds.keySet());
			tPoints.trimToSize();
			Collections.sort(tPoints);
		}

		String mango_file = Args.parseString(args, "mango", null);
		HashMap<Point, ArrayList<Point>> a2bs = new HashMap<Point, ArrayList<Point>>();
		ArrayList<Point> aPoints = new ArrayList<Point>();
		if (mango_file != null) {
			ArrayList<String> lines = CommonUtils.readTextFile(mango_file);
			for (String l : lines) { // each line is a call
				String f[] = l.split("\t");
				Point a = new Region(genome, f[0].replace("chr", ""), Integer.parseInt(f[1]), Integer.parseInt(f[2]))
						.getMidpoint();
				Point b = new Region(genome, f[3].replace("chr", ""), Integer.parseInt(f[4]), Integer.parseInt(f[5]))
						.getMidpoint();
				if (a.getLocation() > b.getLocation()) { // make sure a < b
					Point tmp = a;
					a = b;
					b = tmp;
				}

				if (!a2bs.containsKey(a))
					a2bs.put(a, new ArrayList<Point>());
				a2bs.get(a).add(b);
			}
			aPoints.addAll(a2bs.keySet());
			aPoints.trimToSize();
			Collections.sort(aPoints);
		}

		System.out.println("\nAnnotate and report, " + CommonUtils.timeElapsed(tic0));
		// report the interactions and annotations
		StringBuilder sb = new StringBuilder();
		CommonUtils.writeFile(Args.parseString(args, "out", "Result") + ".readClusters.txt", "");
		for (Interaction it : interactions) {
			sb.append(it.toString()).append("\t");

			// annotate the proximal and distal anchors with TF and HM and
			// regions
			ArrayList<Integer> isOverlapped = new ArrayList<Integer>();
			// left ancor
			int radius = it.leftRegion.getWidth() / 2 + chiapet_radius;
			for (ArrayList<Point> ps : allPeaks) {
				ArrayList<Integer> p = CommonUtils.getPointsIdxWithinWindow(ps, it.leftPoint.expand(radius));
				isOverlapped.add(p.size());
			}
			for (List<Region> rs : allRegions) {
				isOverlapped.add(CommonUtils.getRegionsOverlapsWindow(rs, it.leftRegion, chiapet_radius).size());
			}
			// right anchor
			radius = it.rightRegion.getWidth() / 2 + chiapet_radius;
			for (ArrayList<Point> ps : allPeaks) {
				ArrayList<Integer> p = CommonUtils.getPointsIdxWithinWindow(ps, it.rightPoint.expand(radius));
				isOverlapped.add(p.size());
			}
			for (List<Region> rs : allRegions) {
				isOverlapped.add(CommonUtils.getRegionsOverlapsWindow(rs, it.rightRegion, chiapet_radius).size());
			}
			// print out TF and region overlaps
			if (hms_file != null || tfs_file != null) {
				for (int b : isOverlapped)
					sb.append(b).append("\t");
			}
			// print ChIA-PET call overlap info
			if (germ_file != null || mango_file != null) {
				int leftHalfWidth = it.leftRegion.getWidth() / 2 + chiapet_radius;
				int rightHalfWidth = it.rightRegion.getWidth() / 2 + chiapet_radius;
				Point leftAnchorLeft = new Point(genome, it.leftPoint.getChrom(),
						it.leftPoint.getLocation() - chiapet_radius);
				Point leftAnchorRight = new Point(genome, it.leftPoint.getChrom(),
						it.leftPoint.getLocation() + chiapet_radius);

				// GERM
				int index = Collections.binarySearch(tPoints, leftAnchorLeft);
				if (index < 0) // if key not found
					index = -(index + 1);
				int indexRight = Collections.binarySearch(tPoints, leftAnchorRight);
				if (indexRight < 0) // if key not found
					indexRight = -(indexRight + 1);
				// if key match found, continue to search ( binarySearch() give
				// undefined index with multiple matches)
				boolean isGermOverlapped = false;
				indexRange: for (int i = index - 1; i <= indexRight + 2; i++) {
					if (i < 0 || i >= tPoints.size())
						continue;
					try {
						Point tt = tPoints.get(i);
						if (tt.distance(it.leftPoint) <= leftHalfWidth) {
							if (!t2ds.containsKey(tt))
								continue;
							for (Point d : t2ds.get(tt)) {
								if (d.distance(it.rightPoint) <= rightHalfWidth) {
									isGermOverlapped = true;
									break indexRange;
								}
							}
						}
					} catch (IllegalArgumentException e) { // ignore
					}
				}
				if (isGermOverlapped)
					sb.append("1\t");
				else
					sb.append("0\t");

				// Mango
				// aPoints and bPoints are the midPoint of the two anchorRegions
				index = Collections.binarySearch(aPoints, leftAnchorLeft);
				if (index < 0) // if key not found
					index = -(index + 1);
				indexRight = Collections.binarySearch(aPoints, leftAnchorRight);
				if (indexRight < 0) // if key not found
					indexRight = -(indexRight + 1);
				// if key match found, continue to search ( binarySearch() give
				// undefined index with multiple matches)
				boolean isMangoOverlapped = false;
				indexA: for (int i = index - 1; i <= indexRight + 2; i++) {
					if (i < 0 || i >= aPoints.size())
						continue;
					try {
						Point a = aPoints.get(i);
						if (a.distance(it.leftPoint) <= leftHalfWidth) {
							if (!a2bs.containsKey(a))
								continue;
							for (Point b : a2bs.get(a)) {
								if (b.distance(it.rightPoint) <= rightHalfWidth) {
									isMangoOverlapped = true;
									break indexA;
								}
							}
						}
					} catch (IllegalArgumentException e) { // ignore
					}
				}
				if (isMangoOverlapped)
					sb.append("1\t");
				else
					sb.append("0\t");
			}
			CommonUtils.replaceEnd(sb, '\n');
			
			if (sb.length()>10000000){
				CommonUtils.appendFile(Args.parseString(args, "out", "Result") + ".readClusters.txt", sb.toString());
				sb=new StringBuilder();
			}
		} // each interaction
		CommonUtils.appendFile(Args.parseString(args, "out", "Result") + ".readClusters.txt", sb.toString());

		
		if (isDev){
			StringBuilder sbSprout = new StringBuilder(); //.append("coordA\tcoordB\tcount\tid\n");
			StringBuilder sbLoop = new StringBuilder();
			int id=0;
			for (Interaction it: interactions){
				sbSprout.append(it.toSproutString()).append("\t").append(it.leftRegion.toString()).append("\t")
				.append(it.rightRegion.toString()).append("\t").append(it.toLoopString()).append("\t").append(it.getSpan()).append("\t").append(id).append("\n");
				sbLoop.append(it.toLoopString()).append("\t").append(id).append("\t").append(it.adjustedCount).append("\n");
				id++;
			}
			CommonUtils.writeFile(Args.parseString(args, "out", "Result") + ".unfiltered.sprout_anchors.txt", sbSprout.toString());
			CommonUtils.writeFile(Args.parseString(args, "out", "Result") + ".unfiltered.loopRegions.txt", sbLoop.toString());
		}
		
		/** 
		 * output BEDPE format
		 */
		// HERE we need to also include PET1 for MICC and ChiaSig analysis
		if (micc_min_pet==1){
			for (ReadPair rp : singleton_pets) {
				Interaction it = new Interaction();
				interactions.add(it);
				it.leftPoint = rp.r1;
				it.rightPoint = rp.r2;
				it.leftRegion = new Region(rp.r1.getGenome(), rp.r1.getChrom(), rp.r1.getLocation(), rp.r1.getLocation());
				it.rightRegion = new Region(rp.r2.getGenome(), rp.r2.getChrom(), rp.r2.getLocation(), rp.r2.getLocation());
				it.count = 1;
				it.count2 = 1;
				it.adjustedCount = 1;
			}
			singleton_pets.clear();
			singleton_pets = null;
		}

		sb = new StringBuilder();
		String filename = Args.parseString(args, "out", "Result") + ".bedpe";
		CommonUtils.writeFile(filename, sb.toString());
		for (Interaction it : interactions) {
			// Note: the BED coordinate is centered on the anchor Point, 
			// but the read counts are from the anchor region + padding
			Region leftLocal = it.leftRegion.expand(dc, dc);	//TODO: what is the best bp to expand
			Region rightLocal = it.rightRegion.expand(dc, dc);
			if (isDev)
				sb.append(String.format("%s\t%s\t%d\t%d\t%d\t%d\n", 
						it.leftPoint.expand(leftLocal.getWidth()/2).toBED(), 
						it.rightPoint.expand(rightLocal.getWidth()/2).toBED(), it.adjustedCount,
						CommonUtils.getPointsIdxWithinWindow(lowEnds, leftLocal).size(), 
						CommonUtils.getPointsIdxWithinWindow(highEnds, rightLocal).size(), it.d_c));
			else
				sb.append(String.format("%s\t%s\t%d\t%d\t%d\n", 
						it.leftPoint.expand(leftLocal.getWidth()/2).toBED(), 
						it.rightPoint.expand(rightLocal.getWidth()/2).toBED(), it.adjustedCount,
						CommonUtils.getPointsIdxWithinWindow(lowEnds, leftLocal).size(), 
						CommonUtils.getPointsIdxWithinWindow(highEnds, rightLocal).size()));
			
			if (sb.length()>10000000){
				CommonUtils.appendFile(filename, sb.toString());
				sb=new StringBuilder();
			}
		}
		CommonUtils.appendFile(filename, sb.toString());
	}

	/**
	 * split read pair cluster recursively <br>
	 * at gaps larger than cluster_merge_dist, on both ends alternately
	 * because splitting at one end may remove some PETs and introduce gaps at
	 * the other end
	 */
	ArrayList<ReadPairCluster> splitRecursively(ReadPairCluster cc, boolean toSplitLeftAnchor, boolean isFirstSplit) {
		ArrayList<ReadPairCluster> rpcs = new ArrayList<ReadPairCluster>();
		
		// a point may be used by multiple PETs
		HashMap<Point, ArrayList<ReadPair>> map = new HashMap<Point, ArrayList<ReadPair>>();
		ArrayList<Point> points = new ArrayList<Point>();
		for (ReadPair rp : cc.pets) {
			Point t = toSplitLeftAnchor ? rp.r1 : rp.r2;
			if (!map.containsKey(t))
				map.put(t, new ArrayList<ReadPair>());
			map.get(t).add(rp);
		}
		points.addAll(map.keySet());
		points.trimToSize();
		Collections.sort(points);
		
		ArrayList<Integer> splitIndices = new ArrayList<Integer>();
		int skip=1;		// how many points to skip
		while(true){
			for (int i=skip;i<points.size();i++){
				if (points.get(i).distance(points.get(i-skip)) > max_cluster_merge_dist){		//TODO: change if inter-chrom
					splitIndices.add(i);
				}
			}
			if (!splitIndices.isEmpty() || cc.pets.size() < max_pet_count)
				break;
			else
				skip++;
		}
		if(splitIndices.isEmpty()){
			if (isFirstSplit)
				return splitRecursively(cc, !toSplitLeftAnchor, false);
			else{	// return the original rpc
				rpcs.add(cc);
				return rpcs;
			}
		}
		else{	// splited, split at the other end
			ArrayList<ReadPairCluster> rpcs2 = new ArrayList<ReadPairCluster>();
			int start = 0;
			splitIndices.add(points.size());
			for (int idx: splitIndices){
				ReadPairCluster c = new ReadPairCluster();
				for (int i=start; i<idx; i++){
					for (ReadPair rp : map.get(points.get(i)))
						c.addReadPair(rp);
				}
				start = idx;
				if (c.pets.size() < min_pet_count)
					continue;
				c.update(false);
				ArrayList<ReadPairCluster> splited = splitRecursively(c, !toSplitLeftAnchor, false);
				if (!splited.isEmpty())
					rpcs2.addAll(splited);
			}
			return rpcs2;
		}
		
//		int countSplit = 0;
//		int curr = -100000;
//		ReadPairCluster c = new ReadPairCluster();
//		countSplit--; // first split is not real, subtract count here
//		for (Point p : points) {
//			if (p.getLocation() - curr > max_cluster_merge_dist) { // a big gap
//				countSplit++;
//				if (c.pets.size() >= min)
//					rpcs.add(c);
//				c = new ReadPairCluster();
//			}
//			for (ReadPair rp : map.get(p))
//				c.addReadPair(rp);
//			c.update(false);
//			curr = p.getLocation();
//		}
//		if (c.pets.size() >= min){ // finish up the last cluster
//			rpcs.add(c);
//		}
//
//		// decide to split at the other end, or return
//		if (countSplit > 0 || isFirstSplit) {
//			ArrayList<ReadPairCluster> rpcs2 = new ArrayList<ReadPairCluster>();
//			for (ReadPairCluster rpc: rpcs){
//				ArrayList<ReadPairCluster> splited = splitRecursively(rpc, !toSplitLeftAnchor, false);
//				if (!splited.isEmpty())
//					rpcs2.addAll(splited);
//			}
//			return rpcs2;
//		}
//		else 
//			return rpcs;
	}
	
//	/**
//	 * split read pair cluster recursively <br>
//	 * at gaps larger than cluster_merge_dist, on both ends alternately
//	 * because splitting at one end may remove some PETs and introduce gaps at
//	 * the other end
//	 */
//	ArrayList<ReadPairCluster> splitRecursively(ArrayList<ReadPairCluster> rpcs, boolean toSplitLeftAnchor, 
//			boolean toUseMaxDistance, boolean isFirstSplit) {
//		if (rpcs.isEmpty())
//			return null;
//	
//		int countSplit = 0;
//		ArrayList<ReadPairCluster> rpcs2 = new ArrayList<ReadPairCluster>();
//		for (ReadPairCluster cc : rpcs) {
//			HashMap<Point, ArrayList<ReadPair>> map = new HashMap<Point, ArrayList<ReadPair>>();
//			ArrayList<Point> splitPoints = new ArrayList<Point>();
//			HashSet<Point> tmp = new HashSet<Point>();
//			for (ReadPair rp : cc.pets) {
//				Point t = toSplitLeftAnchor ? rp.r1 : rp.r2;
//				tmp.add(t);
//				if (!map.containsKey(t))
//					map.put(t, new ArrayList<ReadPair>());
//				map.get(t).add(rp);
//			}
//			splitPoints.addAll(tmp);
//			Collections.sort(splitPoints);
//			int curr = -100000;
//			ReadPairCluster c = new ReadPairCluster();
//			countSplit--; // first split is not real, subtract count here
//			for (Point p : splitPoints) {
//				if (p.getLocation() - curr > max_cluster_merge_dist) { // a big gap
//					countSplit++;
//					if (c.pets.size() >= min)
//						rpcs2.add(c);
//					c = new ReadPairCluster();
//				}
//				for (ReadPair rp : map.get(p))
//					c.addReadPair(rp);
//				c.update(false);
//				curr = p.getLocation();
//			}
//			if (c.pets.size() >= min){ // finish up the last cluster
//				rpcs2.add(c);
//			}
//		}
//		if (countSplit > 0 || isFirstSplit) {
//			// split at the other end
//			ArrayList<ReadPairCluster> rpcs3 = splitRecursively(rpcs2, !toSplitLeftAnchor, toUseMaxDistance, false); 
//			return (rpcs3==null || rpcs3.isEmpty()) ? null : rpcs3;
//		} else if (rpcs2.isEmpty())
//			return null;
//		else 
//			return rpcs2;
//	}

	/**
	 * Overlap distal anchors of CID interaction calls with some annotation as regions.
	 * @param args
	 */
	private static void annotateRegions(String[] args) {
		Genome genome = CommonUtils.parseGenome(args);
		String cpcFile = Args.parseString(args, "cpc", null);
		String bed1File = Args.parseString(args, "bed1", null);
		String bed2File = Args.parseString(args, "bed2", null);
		int win = Args.parseInteger(args, "win", 1500);

		Pair<ArrayList<Region>, ArrayList<String>> tmp = CommonUtils.load_BED_regions(genome, bed1File);
		ArrayList<Region> r1s = tmp.car();
		ArrayList<String> s1s = tmp.cdr();

		tmp = CommonUtils.load_BED_regions(genome, bed2File);
		ArrayList<Region> r2s = tmp.car();
		ArrayList<String> s2s = tmp.cdr();

		ArrayList<String> lines = CommonUtils.readTextFile(cpcFile);
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < lines.size(); i++) {
			String t = lines.get(i);
			String f[] = t.split("\t");
			// field 5 (i.e. 6th): distal anchor region. the input cpc file has been swapped anchors if 6th column is the TSS
			Region distal = Region.fromString(genome, f[5]);	
			ArrayList<Integer> idx1 = CommonUtils.getRegionIdxOverlapsWindow(r1s, distal, win);
			ArrayList<Integer> idx2 = CommonUtils.getRegionIdxOverlapsWindow(r2s, distal, win);
			for (int j = 0; j < f.length; j++)
				sb.append(f[j]).append("\t");
			
			for (int id : idx1)
				sb.append(s1s.get(id)).append(",");
			if (idx1.isEmpty())
				sb.append("NULL");
			sb.append("\t");
			for (int id : idx2)
				sb.append(s2s.get(id)).append(",");
			if (idx2.isEmpty())
				sb.append("NULL");
			sb.append("\t").append(idx1.size()).append("\t").append(idx2.size()).append("\t");
			// field 10 (i.e. 11th): ChIA-PET PET count; the last field: ChIA-PET q value
			sb.append(f[10]).append("\t").append(f[f.length - 1]);
			sb.append("\n");
		}
		CommonUtils.writeFile(cpcFile.replace("txt", "") + "annotated.txt", sb.toString());
	}

	private ArrayList<Pair<ReadCache,ReadCache>> prepareGEMData(ArrayList<StrandedPoint> reads) {
		ArrayList<Pair<ReadCache,ReadCache>> expts = new ArrayList<Pair<ReadCache,ReadCache>>();
		ReadCache ipCache = new ReadCache(genome, "IP", null, null);
		int coord = reads.get(0).getLocation();
		int plus =0;
		int minus =0;
		String chrom = reads.get(0).getChrom();
		ArrayList<Integer> coordPlus = new ArrayList<Integer>();
		ArrayList<Integer> coordMinus = new ArrayList<Integer>();
		ArrayList<Float> countPlus = new ArrayList<Float>();
		ArrayList<Float> countMinus = new ArrayList<Float>();
//		System.out.println("Total # reads = " + reads.size());
		for (StrandedPoint p: reads) {
			String chr = p.getChrom();
//			if (!chr.equals("13"))
//				continue;
			if (chrom.equals(chr)) {
				if (p.getLocation()==coord) {
					if(p.getStrand()=='+')
						plus++;
					else
						minus++;
				}
				else {
					if (plus>0) {
						coordPlus.add(coord);
						countPlus.add((float)plus);
					}
					if (minus>0) {
						coordMinus.add(coord);
						countMinus.add((float)minus);
					}	
					coord = p.getLocation();
					plus = 0;
					minus = 0;
					if(p.getStrand()=='+')
						plus++;
					else
						minus++;
				}
			}
			else {	// a new chrom
				// finished previous chrom and init for the new chrom
				if (plus>0) {
					coordPlus.add(coord);
					countPlus.add((float)plus);
				}
				if (minus>0) {
					coordMinus.add(coord);
					countMinus.add((float)minus);
				}
				coord = p.getLocation();
				plus = 0;
				minus = 0;
				if(p.getStrand()=='+')
					plus++;
				else
					minus++;
				ipCache.addHits(chrom, '+', coordPlus, countPlus);
				ipCache.addHits(chrom, '-', coordMinus, countMinus);
				coordPlus = new ArrayList<Integer>();
				coordMinus = new ArrayList<Integer>();
				countPlus = new ArrayList<Float>();
				countMinus = new ArrayList<Float>();
				chrom = chr;
			}
		}
		// finish every thing
		if (plus>0) {
			coordPlus.add(coord);
			countPlus.add((float)plus);
		}
		if (minus>0) {
			coordMinus.add(coord);
			countMinus.add((float)minus);
		}
		ipCache.addHits(chrom, '+', coordPlus, countPlus);
		ipCache.addHits(chrom, '-', coordMinus, countMinus);

		ipCache.populateArrays(true);
//		ipCache.displayStats();
		
		expts.add(new Pair<ReadCache,ReadCache>(ipCache, null));
		return expts;
	}
	/**
	 * Overlap CPC interaction calls with some annotation as regions.
	 * @param args
	 */
	private static void annotateTADs(String[] args) {
		Genome genome = CommonUtils.parseGenome(args);
		String cpcFile = Args.parseString(args, "cpc", null);
		String tadFile = Args.parseString(args, "tad", null);
		int win = Args.parseInteger(args, "win", 1500);

		Pair<ArrayList<Region>, ArrayList<String>> tmp = CommonUtils.load_BED_regions(genome, tadFile);
		ArrayList<Region> r1s = tmp.car();
		ArrayList<String> s1s = tmp.cdr();

		ArrayList<String> lines = CommonUtils.readTextFile(cpcFile);
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < lines.size(); i++) {
			String t = lines.get(i);
			String f[] = t.split("\t");
			// distal anchor region, the input cpc file has swapped anchors if 6th column is the TSS
			Region distal = Region.fromString(genome, f[5]);	
			ArrayList<Integer> idx1 = CommonUtils.getRegionIdxOverlapsWindow(r1s, distal, win);
			for (int j = 0; j < f.length; j++)
				sb.append(f[j]).append("\t");
			for (int id : idx1)
				sb.append(s1s.get(id)).append(",");
			if (idx1.isEmpty())
				sb.append("NULL");
			sb.append("\t").append(idx1.size()).append("\t");
			sb.append(f[10]).append("\t").append(f[f.length - 1]);
			sb.append("\n");
		}
		CommonUtils.writeFile(cpcFile.replace("txt", "") + "annotated.txt", sb.toString());
	}
	
	private static void postProcessing(String[] args) {
		String cpcFile = Args.parseString(args, "cpc", null);
		ArrayList<String> lines = CommonUtils.readTextFile(cpcFile);
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < lines.size(); i++) {
			String t = lines.get(i);
			String f[] = t.split("\t");
			String g1s[] = f[0].trim().split(",");
			String g2s[] = f[3].trim().split(",");
			for (String g1 : g1s) {
				for (String g2 : g2s) {
					sb.append(g1).append("\t").append(f[1]).append("\t").append(f[2]).append("\t").append(g2);
					for (int j = 4; j < f.length; j++)
						sb.append("\t").append(f[j]);
					sb.append("\n");
				}
			}
		}
		CommonUtils.writeFile(cpcFile.replace("txt", "") + "per_gene.txt", sb.toString());
	}

	private static void getPetLength(String[] args) {
		Genome genome = CommonUtils.parseGenome(args);
		ArrayList<String> read_pairs = CommonUtils.readTextFile(Args.parseString(args, "read_pair", null));
		boolean isBEDPE = Args.parseString(args, "format", "rp").equalsIgnoreCase("bedpe");
		StringBuilder sb_minus_plus = new StringBuilder();
		StringBuilder sb_plus_minus = new StringBuilder();
		StringBuilder sb_minus_minus = new StringBuilder();
		StringBuilder sb_plus_plus = new StringBuilder();
		StrandedPoint tmp = null;
//		int count=0;
		for (String s : read_pairs) {
//			System.out.print((count++)+" ");
			String[] f = s.split("\t");
			StrandedPoint r1;
			StrandedPoint r2;
			if (!isBEDPE){
				r1 = StrandedPoint.fromString(genome, f[0]);
				r2 = StrandedPoint.fromString(genome, f[1]);
			}
			else{
				char strand1 = f[8].charAt(0);
				r1 = new StrandedPoint(genome, f[0].replace("chr", ""), (Integer.parseInt(f[1])+Integer.parseInt(f[2]))/2, strand1);
				char strand2 = f[9].charAt(0);
				r2 = new StrandedPoint(genome, f[3].replace("chr", ""), (Integer.parseInt(f[4])+Integer.parseInt(f[5]))/2, strand2);
				// if not both ends are aligned properly, skip
				if (r1.getChrom().equals("*") || r2.getChrom().equals("*"))	
					continue;
			}
			// r1 and r2 should be on the same chromosome
			if (!r1.getChrom().equals(r2.getChrom())) 
				continue;
			if (r1.getLocation() > r2.getLocation()){
				tmp = r1;
				r2 = r1;
				r1 = tmp;
			}
			int dist = r1.distance(r2);
			if (r1.getStrand() == '-') {
				if (r2.getStrand() == '+')
					sb_minus_plus.append(dist).append("\n");
				else if (r2.getStrand() == '-')
					sb_minus_minus.append(dist).append("\n");
			} else if (r1.getStrand() == '+') {
				if (r2.getStrand() == '+')
					sb_plus_plus.append(dist).append("\n");
				else if (r2.getStrand() == '-')
					sb_plus_minus.append(dist).append("\n");
			}
		}
		CommonUtils.writeFile(Args.parseString(args, "out", "Result") + ".minusPlus.length.txt",
				sb_minus_plus.toString());
		CommonUtils.writeFile(Args.parseString(args, "out", "Result") + ".plusMinus.length.txt",
				sb_plus_minus.toString());
		CommonUtils.writeFile(Args.parseString(args, "out", "Result") + ".minusMinus.length.txt",
				sb_minus_minus.toString());
		CommonUtils.writeFile(Args.parseString(args, "out", "Result") + ".plusPlus.length.txt",
				sb_plus_plus.toString());
		System.exit(0);
	}
	
	private class PetBin implements Comparable<PetBin> {
		int binId;
		int r1;
		int r2;
		ArrayList<Integer> pids = new ArrayList<Integer>();		// the original PET ids
		int density = 0;
		int delta;
		int gamma;
		int clusterBinId;
		ArrayList<PetBin> nbBins = new ArrayList<PetBin>();		// BinId of neighbors
		
		PetBin(int r1, int r2, int binId){
			this.r1 = r1;
			this.r2 = r2;
			this.binId = binId;
			clusterBinId = -1;		// init
		}
		void addId(int pid){
			pids.add(pid);
		}
		void addNbBin(PetBin nb){
			nbBins.add(nb);
		}
		public int compareTo(PetBin b) {	// default, descending density
			return this.density==b.density ? 0 : (this.density<b.density?1:-1);
		}
		public int compareByGamma(PetBin b) {	// default, descending density
			return this.gamma==b.gamma ? 0 : (this.gamma<b.gamma?1:-1);
		}
		public String toString(){
			return String.format("id=%d d=%d delta=%d gamma=%d cid=%d", binId, density, delta, gamma, clusterBinId);
		}
	}
	
	class TSSwithReads {
		String symbol;
		int id;
		StrandedPoint coord;
		TreeMap<Integer, ArrayList<Boolean>> reads; // distal read offset -->
													// binary binding indicator
	}

	/** r1 should have lower coordinate than r2 */
	class ReadPair implements Comparable<ReadPair> {
		StrandedPoint r1;
		StrandedPoint r2;

		int compareRead1(ReadPair i) {
			return r1.compareTo(i.r1);
		}

		int compareRead2(ReadPair i) {
			return r2.compareTo(i.r2);
		}

		@Override
		public int compareTo(ReadPair rp) {
			return compareRead1(rp);
		}

		public String toString() {
			return r1.toString() + "--" + r2.toString();
		}
		
		int distance(ReadPair rp) {
			return Math.max(r1.distance(rp.r1), r2.distance(rp.r2));
		}
	}

	class ReadPairCluster implements Comparable<ReadPairCluster> {
		int r1min = Integer.MAX_VALUE;
		int r1max = -1;
		int r2min = Integer.MAX_VALUE;
		int r2max = -1;
//		ReadPair centerPET;
		int r1width = -1;
		int r2width = -1;
		Point leftPoint = null;
		Point rightPoint = null;
		Region leftRegion = null;
		Region rightRegion = null;		
		/** The span is defined as the distance between the median positions of the two anchors */
		int span = -1;
		private ArrayList<ReadPair> pets = new ArrayList<ReadPair>();
		int d_c = -1;

		/**
		 * Update the stats of the RPC. <br>
		 * This method is usually called after a set of RPs are added.<br>
		 * Anchor regions and their widths will be set regardless the value of toSetAnchorPoints
		 */
		void update(boolean toSetAnchorPoints){
			r1width = r1max - r1min;
			r2width = r2max - r2min;
			Point r1 = pets.get(0).r1;
			Point r2 = pets.get(0).r2;
			if (toSetAnchorPoints){
				sortByRead1();
				int size = pets.size();
				if (size % 2 == 1)
					leftPoint = (Point) pets.get(size/2).r1;
				else
					leftPoint = new Point(r1.getGenome(), r1.getChrom(), 
							(pets.get(size/2).r1.getLocation()+pets.get(size/2-1).r1.getLocation())/2);
				sortByRead2();
				if (size % 2 == 1)
					rightPoint = (Point) pets.get(size/2).r2;
				else
					rightPoint = new Point(r2.getGenome(), r2.getChrom(), 
							(pets.get(size/2).r2.getLocation()+pets.get(size/2-1).r2.getLocation())/2);
				span = leftPoint.distance(rightPoint);
			}
			leftRegion = new Region(r1.getGenome(), r1.getChrom(), r1min, r1max);
			rightRegion = new Region(r2.getGenome(), r2.getChrom(), r2min, r2max);
		}
		
		void addReadPair(ReadPair rp) {
			if (r1min > rp.r1.getLocation())
				r1min = rp.r1.getLocation();
			if (r2min > rp.r2.getLocation())
				r2min = rp.r2.getLocation();
			if (r1max < rp.r1.getLocation())
				r1max = rp.r1.getLocation();
			if (r2max < rp.r2.getLocation())
				r2max = rp.r2.getLocation();
			pets.add(rp);
		}

		public int compareToByMidPoints(ReadPairCluster rpc) {
			int offset  = this.leftPoint.offset(rpc.leftPoint);
			if (offset>0)
				return 1;
			else if (offset<0)
				return -1;
			else
				return 0;
		}
		
		@Override
		public int compareTo(ReadPairCluster rpc) {
			if (r1min<rpc.r1min)
				return -1;
			else if (r1min>rpc.r1min)
				return 1;
			else if (r2max<rpc.r2max)
				return -1;
			else if (r2max>rpc.r2max)
				return 1;
			else
				return 0;
		}
		void sortByRead1() {
			Collections.sort(pets, new Comparator<ReadPair>() {
				public int compare(ReadPair o1, ReadPair o2) {
					return o1.compareRead1(o2);
				}
			});
		}

		void sortByRead2() {
			Collections.sort(pets, new Comparator<ReadPair>() {
				public int compare(ReadPair o1, ReadPair o2) {
					return o1.compareRead2(o2);
				}
			});
		}
		
		String getLoopRegionString (int padding) {
			return String.format("%s:%d-%d", pets.get(0).r1.getChrom(), r1min-padding, r2max+padding);
		}
		
		int getLoopRegionWidth() {
			return r2max - r1min;
		}
		
		public int distance(ReadPairCluster rpc) {
			return Math.max(leftPoint.distance(rpc.leftPoint), rightPoint.distance(rpc.rightPoint));
		}
		
		public String toString2() {
			StringBuilder sb = new StringBuilder();
			sb.append(pets.size()).append("=<");
			for (ReadPair rp : pets)
				sb.append(rp.toString()).append(",");
			CommonUtils.replaceEnd(sb, '>');
			return sb.toString();
		}

		public String toString0() {
			StringBuilder sb = new StringBuilder();
			sb.append(pets.size()).append("=<");
			sb.append(pets.get(0).r1.getChrom()).append(":").append(r1min).append("-").append(r1max);
			sb.append("==").append(r2min).append("-").append(r2max).append(">");
			return sb.toString();
		}
		
		public String toString() {
			StringBuilder sb = new StringBuilder();
			sb.append(pets.size()).append(": <");
			sb.append(r1max - r1min).append("--").append((r2max+r2min-r1max-r1min)/2);
			sb.append("--").append(r2max - r2min).append(">");
			return sb.toString();
		}
	}

	class Interaction0 {
		Point tss;
		String geneSymbol;
		Region tssRegion;
		int geneID;
		Point distalPoint;
		boolean isTfAnchord;
		Region distalRegion;
		int count;
		int indirectCount;
		double density;

		// double pvalue;
		public String toString() {
			// return String.format("%d %.1f\t< %s %s -- %s >", count, density,
			// geneSymbol, tssRegion, distalRegion);
			int dist = distalPoint.offset(tss);
			int padding = Math.abs(dist / 20);
			return String.format("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%.1f", geneSymbol,
					(tss instanceof StrandedPoint) ? (StrandedPoint) tss : tss, tssRegion, distalPoint, distalRegion,
					tss.getChrom() + ":"
							+ (Math.min(Math.min(tssRegion.getStart(), distalRegion.getStart()), tss.getLocation())
									- padding)
							+ "-"
							+ (Math.max(Math.max(tssRegion.getEnd(), distalRegion.getEnd()), tss.getLocation())
									+ padding),
					tssRegion.getWidth(), distalRegion.getWidth(), dist, count, indirectCount, density);
		}
	}

	/**
	 * Interaction object <br>
	 * 2016/09/01, cluster PETs without gene annotation, therefore, no more tss
	 * or distal distinction<br>
	 * The anchor regions will overlap gene annotation to label gene
	 * information, but only after the clustering step
	 * 
	 * @author yguo
	 *
	 */
	class Interaction implements Comparable<Interaction> {
		Point leftPoint;
		Region leftRegion;
		String leftLabel;
		Point rightPoint;
		Region rightRegion;
		String rightLabel;
		/** All read counts regardless PET orientation */
		int count;
		/** NON-minus-plus PET counts */
		int count2;
		/** PET counts after adjusting for the minus-plus fraction */
		int adjustedCount;
		int d_c = -1;
		double density;

		public int getSpan(){
			return leftPoint.distance(rightPoint);
		}
		public String toString() {
			// return String.format("%d %.1f\t< %s %s -- %s >", count, density,
			// geneSymbol, tssRegion, distalRegion);
			int dist = rightPoint.offset(leftPoint);
			return String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d", 
					leftLabel, new Point(leftPoint).toString(), leftRegion, rightLabel, new Point(rightPoint).toString(), rightRegion,
					toLoopString(), leftRegion.getWidth(), rightRegion.getWidth(), dist, adjustedCount, count, count2);
		}
		public String toSproutString() {
			return String.format("%s\t%s\t%d", new Point(leftPoint).toString(), new Point(rightPoint).toString(),  adjustedCount);
		}
		public String toLoopString(){
			int dist = rightPoint.offset(leftPoint);
			int padding = Math.abs(dist / 20);
			return leftPoint.getChrom() + ":"
					+ Math.max(1, (Math.min(Math.min(leftRegion.getStart(), rightRegion.getStart()),
							leftPoint.getLocation()) - padding))
					+ "-"
					+ (Math.max(Math.max(leftRegion.getEnd(), rightRegion.getEnd()), leftPoint.getLocation())
							+ padding); // whole it loop region
		}
		@Override
		public int compareTo(Interaction it) {		// descending PET count
			if (adjustedCount > it.adjustedCount)
				return -1;
			else if (adjustedCount < it.adjustedCount)
				return 1;
			else
				return 0;
		}
	}

	class InteractionCall {
		Point tss;
		String tssString;
		String geneSymbol;
		String geneID;
		Region distal;
		double pvalue;
		TreeSet<Point> overlapCoords = new TreeSet<Point>();

		public String toString() {
			int start, end;
			if (tss.getLocation() < distal.getMidpoint().getLocation()) { // if
																			// TSS
																			// is
																			// upstream
				start = tss.getLocation();
				end = distal.getEnd();
			} else {
				start = distal.getStart();
				end = tss.getLocation();
			}
			Region it = new Region(genome, tss.getChrom(), start, end).expand(2000, 2000);
			return String.format("%s\t%d\t%s\t%s\t%s\t%s\t%.2f", it, distal.getMidpoint().distance(tss),
					distal.toString(), tssString, geneSymbol, geneID, -Math.log10(pvalue));
		}
	}
}

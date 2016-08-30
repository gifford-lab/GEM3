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
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.StrandedBase;
import edu.mit.csail.cgs.deepseq.analysis.CCC_Analysis.Interaction;
import edu.mit.csail.cgs.deepseq.analysis.CCC_Analysis.Tss;
import edu.mit.csail.cgs.deepseq.analysis.TFBS_SpaitialAnalysis.Site;
import edu.mit.csail.cgs.deepseq.discovery.kmer.Kmer;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSPeak;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public class ChIAPET_analysis {
	Genome genome;
	Set<String> flags;
	String[] args;
	int read_merge_dist = 500;
	int tss_merge_dist = 500;
	int max_cluster_merge_dist = 3000;
	int distance_factor = 3;
	int self_exclude = 8000;
	int tss_radius = 2000;
	int chiapet_radius = 2000;
	double overlap_ratio = 0.8;
	
	TreeMap<Region, InteractionCall> r2it = new TreeMap<Region, InteractionCall>();
	String fileName = null;
	
	public ChIAPET_analysis(String[] args){
	    genome = CommonUtils.parseGenome(args);
	    
		flags = Args.parseFlags(args);
		this.args = args;
		
		fileName = Args.parseString(args, "bedpe", null);	
		
		read_merge_dist = Args.parseInteger(args, "read_merge_dist", 500);
		tss_merge_dist = Args.parseInteger(args, "tss_merge_dist", 500);
		max_cluster_merge_dist = Args.parseInteger(args, "max_cluster_merge_dist", 3000);
		distance_factor = Args.parseInteger(args, "distance_factor", 3);
		self_exclude = Args.parseInteger(args, "self_exclude", 8000);
		tss_radius = Args.parseInteger(args, "tss_radius", 2000);
		chiapet_radius = Args.parseInteger(args, "chiapet_radius", 2000);
		overlap_ratio = Args.parseDouble(args, "overlap_ratio", 0.8);
		
	}
	
	public static void main(String args[]){
		ChIAPET_analysis analysis = new ChIAPET_analysis(args);
		int type = Args.parseInteger(args, "type", 0);

		switch(type){
		case 0:
			analysis.cleanUpOverlaps();
			analysis.countGenesPerRegion();
			analysis.StatsTAD();
			break;
		case 1:		// count distal read pairs per gene (old: step1)
			analysis.countReadPairs();
			break;
		case 2:		// count distal read pairs per gene (old: step2)
			analysis.clusterDistalReads();
			break;
		case 4:		// find gene-based dense cluster of read pairs
			analysis.findDenseClusters();
			break;
		case 3:		// merged-TSS based clustering
			analysis.findTssInteractions();
			break;
		}
		
	}
	/**
	 * Clean up data. Merge overlap distal regions if they are connected to the same TSS. Optionally remove distal regions that overlap with the connected TSS
	 */
	void cleanUpOverlaps(){
		ArrayList<String> texts = CommonUtils.readTextFile(fileName);
		// use tssString as the key to ensure uniqueness, tss Point objects have different references even if from same tss
		TreeMap<String, ArrayList<Region>> tss2distalRegions = new TreeMap<String, ArrayList<Region>>();
		for (String line:texts){
			String f[] = line.trim().split("\t");
			InteractionCall it = new InteractionCall();
			it.tssString = f[6];
			it.tss = Point.fromString(genome, f[6]);
			it.geneID = f[7];
			it.geneSymbol = f[8];
			if (f.length<=13){
				it.distal = Region.fromString(genome, f[9]);
				it.pvalue = Double.parseDouble(f[11]);
			}
			else{
				it.distal = Region.fromString(genome, f[13]);
				it.pvalue = Double.parseDouble(f[15]);
			}
			
			// skip interactions that have distal regions containing the TSS?
			if (flags.contains("rm_self")){
				if (it.distal.contains(it.tss))
					continue;
			}

			while (r2it.containsKey(it.distal))
				it.distal = it.distal.expand(-1, -1);	// if duplicate, shrink 1bp to make it uniquie, if connect to same TSS, it will be merged later
			r2it.put(it.distal, it);
			if (!tss2distalRegions.containsKey(it.tssString))
				tss2distalRegions.put(it.tssString, new ArrayList<Region>() );
			tss2distalRegions.get(it.tssString).add(it.distal);
		}
		
		// for each tss, merge overlapping distal regions (<1kb)
		for (String tss: tss2distalRegions.keySet()){
			ArrayList<Region> regions = tss2distalRegions.get(tss);
			ArrayList<Region> mergedRegions = new ArrayList<Region>();
			Collections.sort(regions);
			Region previous = regions.get(0);
			ArrayList<Region> previousRegions = new ArrayList<Region>();
			previousRegions.add(previous);
			
			for (int i = 1; i < regions.size(); i++) {
	          Region region = regions.get(i);
				// if overlaps with previous region, combine the regions, take the best p-values
				if (previous.overlaps(region)){
					previous = previous.combine(region);
					previousRegions.add(region);
				} 
				else{	// not overlap any more, update, then move to next one
					mergedRegions.add(previous);
					// merge overlapping regions, update interactions
					if (previousRegions.size()>1){	// merged
						InteractionCall it=null;
						double bestPvalue = 1;
						for (Region r:previousRegions){
							it = r2it.get(r);
							r2it.remove(r);			// remove old one
							bestPvalue = Math.min(bestPvalue, it.pvalue);
						}
						it.distal = previous;		// previous has been merged
						it.pvalue = bestPvalue;
						r2it.put(previous, it);		// add merged region
					}
					previousRegions.clear();
					previous = region;
					previousRegions.add(previous);
				}
			}
			mergedRegions.add(previous);
			if (previousRegions.size()>1){	// merged
				InteractionCall it=null;
				// merge overlapping regions, update interactions
				double bestPvalue = 1;
				for (Region r:previousRegions){
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
		if (flags.contains("print_merged")){
			StringBuilder sb = new StringBuilder();
			for (Region r:r2it.keySet()){
				InteractionCall it = r2it.get(r);
				sb.append(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%.2e", it.distal.toBED(), it.tss.expand(2000).toBED(), it.tss.toString(), it.geneID, it.geneSymbol, it.distal.toString(), it.distal.getWidth(), it.pvalue)).append("\n");
			}
			CommonUtils.writeFile(fileName.replace(".bedpe", (flags.contains("rm_self")?".rmSelf":"") + ".mergeDistal.bedpe"), sb.toString());
		}
	}
	
	/**
	 * Count number of genes connected by each region (overlapping regions are merged).
	 * @param fileName
	 */
	void countGenesPerRegion(){
		// load RefSeq gene annotation
		int tssRange = Args.parseInteger(args, "tss_range", 100);
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
		/** all TSS annotatons organized by chroms */
		HashMap<String, ArrayList<StrandedPoint>> chr2tsss = new HashMap<String, ArrayList<StrandedPoint>>();
		for (StrandedPoint tss: tss2genes.keySet()){
			String chr = tss.getChrom();
			if (!chr2tsss.containsKey(chr))
				chr2tsss.put(chr, new ArrayList<StrandedPoint>());
			chr2tsss.get(chr).add(tss);			
		}
		for (String chr:chr2tsss.keySet()){
			ArrayList<StrandedPoint> tsss = chr2tsss.get(chr);
			Collections.sort(tsss);
			// print inter-TSS distances
//			for (int i=0;i<tsss.size()-1;i++){
//				System.out.println(tsss.get(i+1).distance(tsss.get(i)));
//			}
		}
		
		
		// merge nearby regions, collect genes within tssRange of the TSSs.
		TreeSet<String> allGenes = new TreeSet<String>();
		ArrayList<Region> rs = new ArrayList<Region>();
		// if coords are provided, only report the subset of distal regions that overlaps
		String coords_file = Args.parseString(args, "coords", null);
		String coords_name = Args.parseString(args, "coords_name", null);
		if (coords_file!=null){
			ArrayList<Point> coords = CommonUtils.loadCgsPointFile(coords_file, genome);
			for (Region r:r2it.keySet()){
				for (Point p: coords){
					if (!r.getChrom().equalsIgnoreCase(p.getChrom()))
						continue;
					if (r.distance(p)<=500)	{	// add distal anchor if overlap with enhancer coords
						rs.add(r);
						r2it.get(r).overlapCoords.add(p);
					}
				}
			}
			
			// print the enhancer_coords and the tss_geneSymbol pairs
			StringBuilder sb1 = new StringBuilder();
			sb1.append("#Coord\tTSS\tGene\n");
			for (Region r: rs){
				InteractionCall it = r2it.get(r);
				Point tss = it.tss;
				ArrayList<StrandedPoint> tsss = chr2tsss.get(tss.getChrom());
				ArrayList<StrandedPoint> tss_linked = new ArrayList<StrandedPoint> ();
				if (tsss==null){
					continue;
				}
				int idx = Collections.binarySearch(tsss, tss);
				if (idx<0)
					idx = -(idx+1);  // insert point 
				for (int j=idx; j<tsss.size();j++){
					if (tss.distance(tsss.get(j))>tssRange)
						break;
					else
						tss_linked.add(tsss.get(j));
				}
				for (int j=idx-1; j>=0;j--){
					if (tss.distance(tsss.get(j))>tssRange)
						break;
					else
						tss_linked.add(tsss.get(j));
				}
				for (StrandedPoint t: tss_linked){
					for (Point p: it.overlapCoords){
						TreeSet<String> genes = tss2genes.get(t);
						for (String g: genes)
							sb1.append(p.toString()).append("\t").append(t.toString()).append("\t").append(g).append("\n");
					}
				}
			}
			CommonUtils.writeFile(fileName.replace(".bedpe", (flags.contains("rm_self")?".rmSelf":"")
					 + (coords_file!=null?("."+coords_name):"") + ".tss" + tssRange + ".coord2genes.txt"), sb1.toString());	

		}
		else		
			rs.addAll(r2it.keySet());		// add all distal regions, rs is sorted, r2it is a TreeMap
		
		Region merged = rs.get(0);
		ArrayList<Integer> ids = new ArrayList<Integer>();
		ids.add(0);
		StringBuilder sb = new StringBuilder("#Merged_region\twidth\tr_id\tgenes\tcount\n");

		for (int i=1;i<rs.size();i++){
			Region r = rs.get(i);
			if (merged.getChrom().equals(r.getChrom()) && merged.distance(r)<Args.parseInteger(args, "distance", 1000)){
				merged = merged.combine(r);
				ids.add(i);
			}
			else{
				sb.append(merged.toString()).append("\t").append(merged.getWidth()).append("\t");
				
				// for each merged region, get all TSS, find refSeq genes within 100bp window
				TreeSet<String> genes = new TreeSet<String>();
				for (int id:ids){
					InteractionCall it = r2it.get(rs.get(id));
					sb.append(id).append(",");
					Point tss = it.tss;
					ArrayList<StrandedPoint> tsss = chr2tsss.get(tss.getChrom());
					if (tsss==null){
						genes.add(it.geneSymbol);
						continue;
					}
					int idx = Collections.binarySearch(tsss, tss);
					if (idx<0)
						idx = -(idx+1);  // insert point 
					for (int j=idx; j<tsss.size();j++){
						if (tss.distance(tsss.get(j))>tssRange)
							break;
						else
							genes.addAll(tss2genes.get(tsss.get(j)));
					}
					for (int j=idx-1; j>=0;j--){
						if (tss.distance(tsss.get(j))>tssRange)
							break;
						else
							genes.addAll(tss2genes.get(tsss.get(j)));
					}	
				}
				CommonUtils.replaceEnd(sb, '\t');
				for (String s:genes)
					sb.append(s).append(",");
				CommonUtils.replaceEnd(sb, '\t');
				sb.append(genes.size()).append("\n");
				ids.clear();
				ids.add(i);
				allGenes.addAll(genes);
				
				merged = r;		// setup for next merge
			}
		}
		// finish the last merge
		sb.append(merged.toString()).append("\t").append(merged.getWidth()).append("\t");
		TreeSet<String> genes = new TreeSet<String>();
		for (int id:ids){
			sb.append(id).append(",");
			Point tss = r2it.get(rs.get(id)).tss;
			ArrayList<StrandedPoint> tsss = chr2tsss.get(tss.getChrom());
			if (tsss==null){
				genes.add(r2it.get(rs.get(id)).geneSymbol);
				continue;
			}
			int idx = Collections.binarySearch(tsss, tss);
			if (idx<0)
				idx = -(idx+1);  // insert point 
			for (int j=idx; j<tsss.size();j++){
				if (tss.distance(tsss.get(j))>tssRange)
					break;
				else
					genes.addAll(tss2genes.get(tsss.get(j)));
			}
			for (int j=idx-1; j>=0;j--){
				if (tss.distance(tsss.get(j))>tssRange)
					break;
				else
					genes.addAll(tss2genes.get(tsss.get(j)));
			}	
		}
		CommonUtils.replaceEnd(sb, '\t');
		if (genes.isEmpty())
			sb.append("None").append(",");
		for (String s:genes)
			sb.append(s).append(",");
		CommonUtils.replaceEnd(sb, '\t');
		sb.append(genes.size()).append("\n");
		allGenes.addAll(genes);
		
//		System.out.println(sb.toString());
		CommonUtils.writeFile(fileName.replace(".bedpe", (flags.contains("rm_self")?".rmSelf":"")
				 + (coords_file!=null?("."+coords_name):"") + ".tss" + tssRange + ".per_region_count.txt"), sb.toString());
		
		// print out all linked genes
		sb = new StringBuilder();
		for (String g: allGenes)
			sb.append(g).append("\n");
		CommonUtils.writeFile(fileName.replace(".bedpe", (flags.contains("rm_self")?".rmSelf":"")
				 + (coords_file!=null?("."+coords_name):"") + ".tss" + tssRange + ".geneSymbols.txt"), sb.toString());	
	}
	
	void StatsTAD(){
		String tad_file = Args.parseString(args, "tad", null);
		ArrayList<Region> tads = CommonUtils.load_BED_regions(genome, tad_file).car();
		Collections.sort(tads);
		ArrayList<InteractionCall> itNonTAD = new ArrayList<InteractionCall>();
		ArrayList<InteractionCall> itSameTAD = new ArrayList<InteractionCall>();
		ArrayList<InteractionCall> itCrossTAD = new ArrayList<InteractionCall>();
		
		for (Region r: r2it.keySet()){
			Region mid = r.getMidpoint().expand(0);
			int idx = Collections.binarySearch(tads, mid);
			Region tad = null;
			if (idx<0){
				idx = -(idx+1) -1;  // insert point - 1 ==> Previous object
				tad = tads.get(idx);
				if (!tad.contains(mid)){
//					System.err.println(String.format("Point %s is not within any TAD!", p.toString()));
					itNonTAD.add(r2it.get(r));
				}
				else{// now tad contains the distal coord
					if (tad.contains(r2it.get(r).tss))
						itSameTAD.add(r2it.get(r));
					else{	// find the tad that TSS is in
						Region tss = r2it.get(r).tss.expand(0);
						idx = Collections.binarySearch(tads, tss);
						if (idx<0){
							idx = -(idx+1) -1;  // insert point - 1 ==> Previous object
							tad = tads.get(idx);
							if (!tad.contains(tss))	// TSS is not in a TAD
								itNonTAD.add(r2it.get(r));
							else	// in TAD, must be another TAD
								itCrossTAD.add(r2it.get(r));
						}
					}
				}
			}
		}

		System.out.println(String.format("In same TAD:\t %d\nCross TAD:\t %d\nNot in TAD:\t %d\n", itSameTAD.size(), itCrossTAD.size(), itNonTAD.size()));

		StringBuilder sb = new StringBuilder("#Interaction\tdistance\tdistal\ttss\tSymbol\tgeneID\tp_-lg10\tTAD_status\n");
		for (InteractionCall it: itSameTAD)
			sb.append(it.toString()).append("\t1\n");
		for (InteractionCall it: itCrossTAD)
			sb.append(it.toString()).append("\t2\n");
		for (InteractionCall it: itNonTAD)
			sb.append(it.toString()).append("\t0\n");
		CommonUtils.writeFile(fileName.replace(".bedpe", (flags.contains("rm_self")?".rmSelf":"") + ".StatsTAD.txt"), sb.toString());
		
	}
	
	
	private void countReadPairs(){
		long tic = System.currentTimeMillis();
		HashSet<String> geneSet = new HashSet<String>();
		String gString = Args.parseString(args, "genes", null);
		if (gString==null){
			String gFile = Args.parseString(args, "gene_file", null);
			ArrayList<String> lines = CommonUtils.readTextFile(gFile);
			for (String g:lines)
				geneSet.add(g.trim());
		}
		else{
			String genes[] = Args.parseString(args, "genes", null).split(",");
			for (String g:genes)
				geneSet.add(g.trim());
		}
		
		// load refSeq gene annotation
		int tssRadius = Args.parseInteger(args, "tss_range", 10001)/2;
		int chiapetRadius = Args.parseInteger(args, "chiapet_radius", 2000);
		ArrayList<String> gene_annots = CommonUtils.readTextFile(Args.parseString(args, "gene_anno", null));
		TreeMap<String, TreeSet<StrandedPoint>> gene2tss = new TreeMap<String, TreeSet<StrandedPoint>>();
		for (int i=0;i<gene_annots.size();i++){
			String t = gene_annots.get(i);
			if (t.startsWith("#"))
				continue;
			String f[] = t.split("\t");
			String symbol = f[12];
			if (!geneSet.contains(symbol))
				continue;
			String chr = f[2].replace("chr", "");
			char strand = f[3].charAt(0);
			StrandedPoint tss = new StrandedPoint(genome, chr, Integer.parseInt(f[strand=='+'?4:5]), strand);
			if (!gene2tss.containsKey(symbol))
				gene2tss.put(symbol, new TreeSet<StrandedPoint>());
			gene2tss.get(symbol).add(tss);
		}
		
		// load read pairs
		// for now, just loop through every read pair. To run faster, should sort by each end and use binary search to find tss end overlaps
		ArrayList<String> read_pairs = CommonUtils.readTextFile(Args.parseString(args, "read_pair", null));
		HashMap<String, Pair<ArrayList<Point>, ArrayList<Point>>> chr2reads = new HashMap<String, Pair<ArrayList<Point>, ArrayList<Point>>> ();
		for (String s: read_pairs){
			String[] f = s.split("\t");
			Point r1 = Point.fromString(genome, f[0]);
			String r1Chrom=r1.getChrom();
			Point r2 = Point.fromString(genome, f[1]);
			if (!r1Chrom.equals(r2.getChrom()))		// skip if not from the same chromosome
				continue;
			if (!chr2reads.containsKey(r1Chrom)){
				ArrayList<Point> r1s=new ArrayList<Point>();
				ArrayList<Point> r2s=new ArrayList<Point>();
				chr2reads.put(r1Chrom, new Pair<ArrayList<Point>, ArrayList<Point>>(r1s, r2s));
			}
			Pair<ArrayList<Point>, ArrayList<Point>> reads = chr2reads.get(r1Chrom);
			reads.car().add(r1);
			reads.cdr().add(r2);
		}

		System.out.println("Loaded ChIA-PET read pairs: "+CommonUtils.timeElapsed(tic));
		System.out.println();
		
		// load TF sites
		ArrayList<String> tfs = CommonUtils.readTextFile(Args.parseString(args, "tf_sites", null));
		ArrayList<List<GPSPeak>> allPeaks = new ArrayList<List<GPSPeak>>();
		for (int i=0;i<tfs.size();i++){
			try{
				allPeaks.add(GPSParser.parseGPSOutput(tfs.get(i), genome));
				System.out.println("Loaded "+tfs.get(i));
			}
			catch (IOException e){
				System.out.println(tfs.get(i)+" does not have a valid GPS/GEM event call file.");
				e.printStackTrace(System.err);
				System.exit(1);
			}
		}

		// load histone mark regions
		ArrayList<String> hms = CommonUtils.readTextFile(Args.parseString(args, "regions", null));
		ArrayList<List<Region>> allRegions = new ArrayList<List<Region>>();
		for (int i=0;i<hms.size();i++){
			allRegions.add(CommonUtils.load_BED_regions(genome, hms.get(i)).car());
			System.out.println("Loaded "+hms.get(i));
		}
		System.out.println();
//		TreeMap<String, ArrayList<Integer>> gene2distances = new TreeMap<String, ArrayList<Integer>>();		
		// compute distance for each gene
		ArrayList<String> geneList = new ArrayList<String>();
		geneList.addAll(gene2tss.keySet());
		StringBuilder sb = new StringBuilder();
		
		for (int id=0;id<geneList.size();id++){
			String g = geneList.get(id);
			System.out.print(g+" ");
			ArrayList<Integer> distances = new ArrayList<Integer> (); 
			ArrayList<ArrayList<Integer>> isTfBounds = new ArrayList<ArrayList<Integer>>();
			TreeSet<StrandedPoint> coords = gene2tss.get(g);
			// if gene has multiple TSSs, use the center position
			int count = coords.size();
			StrandedPoint centerPoint = null;
			for (StrandedPoint p:coords){
				if (count<coords.size()/2)
					break;
				else{
					centerPoint = p;
					count--;
				}
			}
			
			// if one end of the read pair is near TSS, compute the offset of the other end
			boolean isMinus = centerPoint.getStrand()=='-';
			Pair<ArrayList<Point>, ArrayList<Point>> reads = chr2reads.get(centerPoint.getChrom());
			if (reads==null)
				continue;
			ArrayList<Point> read1s = reads.car();
			ArrayList<Point> read2s = reads.cdr();
			for (int i=0;i<read1s.size();i++){
				int offset_p1 = read1s.get(i).offset(centerPoint);
				int offset_p2 = read2s.get(i).offset(centerPoint);
				int dist_p1 = Math.abs(offset_p1);
				int dist_p2 = Math.abs(offset_p2);
				// only add distance to the list if one read is within TSS_Radius, the other read is outside of TSS_Radius
				if (dist_p1<tssRadius){
					if (dist_p2>tssRadius){
						distances.add(isMinus?-offset_p2:offset_p2);
						ArrayList<Integer> isBound = new ArrayList<Integer>();
						Point p = read2s.get(i);
						for (int j=0;j<allPeaks.size();j++){
							List<GPSPeak> peaks = allPeaks.get(j);
							int bound = 0;
							for (GPSPeak gps: peaks){
								if (gps.getChrom().equals(p.getChrom()) && gps.distance(p)<=chiapetRadius){
									bound = 1;
									break;
								}
							}
							isBound.add(bound);
						}
						for (int j=0;j<allRegions.size();j++){
							List<Region> rs = allRegions.get(j);
							int bound = 0;
							for (Region r: rs){
								// if the region r contains point p, or the distance between midPoint of r and p is less than ChIAPET_radias
								if (r.getChrom().equals(p.getChrom()) && (r.contains(p) || r.getMidpoint().distance(p)<=chiapetRadius)){
									bound = 1;
									break;
								}
							}
							isBound.add(bound);
						}
						isTfBounds.add(isBound);
					}
				}
				else{
					if (dist_p2<tssRadius){
						distances.add(isMinus?-offset_p1:offset_p1);
						ArrayList<Integer> isBound = new ArrayList<Integer>();
						Point p = read1s.get(i);
						for (int j=0;j<allPeaks.size();j++){
							List<GPSPeak> peaks = allPeaks.get(j);
							int bound = 0;
							for (GPSPeak gps: peaks){
								if (gps.getChrom().equals(p.getChrom()) && gps.distance(p)<=chiapetRadius){
									bound = 1;
									break;
								}
							}
							isBound.add(bound);
						}
						for (int j=0;j<allRegions.size();j++){
							List<Region> rs = allRegions.get(j);
							int bound = 0;
							for (Region r: rs){
								// if the region r contains point p, or the distance between midPoint of r and p is less than ChIAPET_radias
								if (r.getChrom().equals(p.getChrom()) && (r.contains(p) || r.getMidpoint().distance(p)<=chiapetRadius)){
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
			
			if (!distances.isEmpty()){
//				gene2distances.put(g, distances);
				for (int i=0;i<distances.size();i++){
					sb.append(g).append("\t").append(centerPoint.toString()).append("\t").append(id);
					sb.append("\t").append(distances.get(i));
					for (int b: isTfBounds.get(i))
						sb.append("\t").append(b);
					sb.append("\n");
				}
			}
		}  // for each gene
		CommonUtils.writeFile("all_genes.distal_offsets.txt", sb.toString());
		
		System.out.println("\n\n"+CommonUtils.timeElapsed(tic));
	}
	
	private void clusterDistalReads(){
			int tss_exclude = Args.parseInteger(args, "tss_exclude", 8000);
			int step = Args.parseInteger(args, "merge_dist", 1500);
			int minRead = Args.parseInteger(args, "min_count", 2);
			
			// load data
			ArrayList<String> lines = CommonUtils.readTextFile(Args.parseString(args, "tss_reads", null));
			TSSwithReads tss = new TSSwithReads();
			tss.symbol = "---";
			ArrayList<TSSwithReads> allTss = new ArrayList<TSSwithReads>();
			for (String l: lines){		// each line is a distal read
				String f[] = l.split("\t");
				int offset = Integer.parseInt(f[3]);
				if (Math.abs(offset)<tss_exclude)		// skip the read if it is within TSS exclude region						
					continue;
				if (!f[0].equals(tss.symbol)){	// a new gene
					tss = new TSSwithReads();
					tss.symbol = f[0];
					tss.coord = StrandedPoint.fromString(genome, f[1]);
					tss.id = Integer.parseInt(f[2]);
					tss.reads = new TreeMap<Integer, ArrayList<Boolean>>();
					allTss.add(tss);				
				}
				ArrayList<Boolean> isBound = new ArrayList<Boolean>();
				for (int i=4;i<f.length;i++){
					isBound.add(f[i].equals("1"));
				}
				tss.reads.put(offset, isBound);
			}
	
			lines = CommonUtils.readTextFile(Args.parseString(args, "germ", null));
			HashMap<Point,ArrayList<Point>> germTss2distals = new HashMap<Point,ArrayList<Point>>();
			ArrayList<Point> germTss = new ArrayList<Point>();
			for (String l: lines){		// each line is a call
				String f[] = l.split("\t");
				Point t = new Region(genome, f[3].replace("chr", ""), Integer.parseInt(f[4]), Integer.parseInt(f[5])).getMidpoint();
				if (!germTss2distals.containsKey(t))
					germTss2distals.put(t, new ArrayList<Point>());
				germTss2distals.get(t).add(new Region(genome, f[0].replace("chr", ""), Integer.parseInt(f[1]), Integer.parseInt(f[2])).getMidpoint());
			}
			germTss.addAll(germTss2distals.keySet());
			Collections.sort(germTss);
			
			lines = CommonUtils.readTextFile(Args.parseString(args, "mango", null));
			HashMap<Point, ArrayList<Point>> a2bs = new HashMap<Point, ArrayList<Point>>();
			HashMap<Point, ArrayList<Point>> b2as = new HashMap<Point, ArrayList<Point>>();
			for (String l: lines){		// each line is a call
				String f[] = l.split("\t");
				Point a = new Region(genome, f[0].replace("chr", ""), Integer.parseInt(f[1]), Integer.parseInt(f[2])).getMidpoint();
				Point b = new Region(genome, f[3].replace("chr", ""), Integer.parseInt(f[4]), Integer.parseInt(f[5])).getMidpoint();
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
			for (TSSwithReads t:allTss){
				ArrayList<Integer> cluster = new ArrayList<Integer>();
				for (int offset:t.reads.keySet()){
					if (cluster.isEmpty() || offset-cluster.get(cluster.size()-1)<step){
						cluster.add(offset);
					}
					else{		// have a large distance, finish old cluster, create new cluster
						if (cluster.size()>=minRead){	// at least 2 reads
							int median = cluster.get(cluster.size()/2);
							Point tssPoint = t.coord;
							Point distalPoint = new Point(genome, t.coord.getChrom(),t.coord.getLocation()+(t.coord.getStrand()=='+'?median:-median));
							Region distalRegion = null;
							if (t.coord.getStrand()=='+'){
								distalRegion = new Region(genome, t.coord.getChrom(),t.coord.getLocation()+cluster.get(0), 
										t.coord.getLocation()+cluster.get(cluster.size()-1));
							}
							else{
								distalRegion = new Region(genome, t.coord.getChrom(),t.coord.getLocation()-cluster.get(cluster.size()-1), 
										t.coord.getLocation()-cluster.get(0));
							}
							// print result if the read cluster is not in the tss exclusion range
							System.out.print(String.format("%s\t%s\t%s\t%s\t%d\t%d\t%d\t", 
									t.symbol, t.coord.getLocationString(), distalRegion.getLocationString(), distalPoint.getLocationString(), median, cluster.size(), 
									distalRegion.getWidth()));
							
							// print binding overlap information
							int count=t.reads.get(cluster.get(0)).size();
							for (int c=0;c<count;c++){
								boolean isBound = false;
								for (int clusterOffset: cluster){
									isBound = isBound || t.reads.get(clusterOffset).get(c);
								}
								System.out.print(isBound?"1\t":"0\t");
							}
							
							// print ChIA-PET call overlap info
							Point tssLeft = new Point(genome, t.coord.getChrom(),t.coord.getLocation()-2000);
							Point tssRight = new Point(genome, t.coord.getChrom(),t.coord.getLocation()+2000);
							
							// GERM
							int index = Collections.binarySearch(germTss,  tssLeft);
							if( index < 0 )  							// if key not found
								index = -(index+1); 
							int indexRight = Collections.binarySearch(germTss,  tssRight);
							if( indexRight < 0 )  							// if key not found
								indexRight = -(indexRight+1); 
							// if key match found, continue to search ( binarySearch() give undefined index with multiple matches)
							boolean isOverlapped = false;
							indexRange: for (int i=index-1;i<=indexRight+2;i++){
								if (i<0 || i>=germTss.size())
									continue;
								try{
									Point tt = germTss.get(i);
									if (tt.distance(tssPoint)<=2000){ 
										if (!germTss2distals.containsKey(tt))
											continue;
										for (Point d:germTss2distals.get(tt)){
	//										System.out.print(tt.getLocationString()+"\t"+d.getLocationString());
											if (d.distance(distalPoint)<=2000){
												isOverlapped = true;
	//											System.out.println("\tHIT");
												break indexRange;
											}
	//										else
	//											System.out.println();
										}
									}
								}
								catch (IllegalArgumentException e){	// ignore								
								}								
							}
							System.out.print(isOverlapped?"1\t":"0\t");
							
							// Mango
							index = Collections.binarySearch(aPoints,  tssLeft);
							if( index < 0 )  							// if key not found
								index = -(index+1); 
							indexRight = Collections.binarySearch(aPoints,  tssRight);
							if( indexRight < 0 )  							// if key not found
								indexRight = -(indexRight+1); 
							// if key match found, continue to search ( binarySearch() give undefined index with multiple matches)
							isOverlapped = false;
							indexA: for (int i=index-1;i<=indexRight+2;i++){
								if (i<0 || i>=aPoints.size())
									continue;
								try{
									Point a = aPoints.get(i);
									if (a.distance(tssPoint)<=2000){ 
										if (!a2bs.containsKey(a))
											continue;
										for (Point b:a2bs.get(a)){
											if (b.distance(distalPoint)<=2000){
												isOverlapped = true;
												break indexA;
											}
										}
									}
								}
								catch (IllegalArgumentException e){	// ignore								
								}
							}
							if (isOverlapped)
								System.out.print("1\t");
							else{
								index = Collections.binarySearch(bPoints,  tssLeft);
								if( index < 0 )  							// if key not found
									index = -(index+1); 
								indexRight = Collections.binarySearch(bPoints,  tssRight);
								if( indexRight < 0 )  							// if key not found
									indexRight = -(indexRight+1); 
								// if key match found, continue to search ( binarySearch() give undefined index with multiple matches)
								isOverlapped = false;
								indexB: for (int i=index-1;i<=indexRight+2;i++){
									if (i<0 || i>=bPoints.size())
										continue;
									try{
										Point b = bPoints.get(i);
										if (b.distance(tssPoint)<=2000){ 
											if (!b2as.containsKey(b))
												continue;
											for (Point a:b2as.get(b)){
												if (a.distance(distalPoint)<=2000){
													isOverlapped = true;
													break indexB;
												}
											}
										}
									}
									catch (IllegalArgumentException e){	// ignore								
									}								
								}
								System.out.print(isOverlapped?"1\t":"0\t");
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

	private void findDenseClusters(){
		long tic0 = System.currentTimeMillis();
		long tic = System.currentTimeMillis();
		int read_merge_dist = Args.parseInteger(args, "read_merge_dist", 500);
		int max_cluster_merge_dist = Args.parseInteger(args, "max_cluster_merge_dist", 3000);
		int distance_factor = Args.parseInteger(args, "distance_factor", 3);
		int tss_exclude = Args.parseInteger(args, "tss_exclude", 8000);
		int tss_radius = Args.parseInteger(args, "tss_radius", 2000);
		int chiapet_radius = Args.parseInteger(args, "chiapet_radius", 2000);
		double overlap_ratio = Args.parseDouble(args, "overlap_ratio", 0.8);

		// load the genes to find interactions
		// geneSet can be supplied as a string on command line, or a file.
		// if neither is supplied, all the genes in the gene_annotation is used.
		HashSet<String> geneSet = new HashSet<String>();
		String gString = Args.parseString(args, "genes", null);
		if (gString==null){
			String gFile = Args.parseString(args, "gene_file", null);
			if (gFile!=null){
				ArrayList<String> lines = CommonUtils.readTextFile(gFile);
				for (String g:lines)
					geneSet.add(g.trim());
			}
		}
		else{
			String genes[] = Args.parseString(args, "genes", "").split(",");
			if (!genes[0].equals("")){
				for (String g:genes)
					geneSet.add(g.trim());
			}
		}
		
		// load refSeq gene annotation
		ArrayList<String> gene_annots = CommonUtils.readTextFile(Args.parseString(args, "gene_anno", null));
		TreeMap<String, TreeSet<StrandedPoint>> gene2tss = new TreeMap<String, TreeSet<StrandedPoint>>();
		for (int i=0;i<gene_annots.size();i++){
			String t = gene_annots.get(i);
			if (t.startsWith("#"))
				continue;
			String f[] = t.split("\t");
			String symbol = f[12];
			if (!geneSet.isEmpty()){
				if (!geneSet.contains(symbol))
					continue;
			}
			String chr = f[2].replace("chr", "");
			char strand = f[3].charAt(0);
			StrandedPoint tss = new StrandedPoint(genome, chr, Integer.parseInt(f[strand=='+'?4:5]), strand);
			if (!gene2tss.containsKey(symbol))
				gene2tss.put(symbol, new TreeSet<StrandedPoint>());
			gene2tss.get(symbol).add(tss);
		}
		
		
		// load read pairs
		// only use read pairs on the same chromosome
		// the left read is required to be lower than the right read, if not, flip the two
		
		// sort by each end so that we can use binary search to find matches or overlaps
		System.out.println("Loading ChIA-PET read pairs: "+CommonUtils.timeElapsed(tic));

		ArrayList<String> read_pairs = CommonUtils.readTextFile(Args.parseString(args, "read_pair", null));
		ArrayList<ReadPair> low = new ArrayList<ReadPair> ();	// all read pairs sorted by the low end
		ArrayList<ReadPair> high = new ArrayList<ReadPair> ();	// sorted by high end
		for (String s: read_pairs){
			String[] f = s.split("\t");
			Point r1 = Point.fromString(genome, f[0]);
			String r1Chrom=r1.getChrom();
			Point r2 = Point.fromString(genome, f[1]);
			// TODO: change next line if prediction cross-chrom interactions
			if (!r1Chrom.equals(r2.getChrom()))		// r1 and r2 should be on the same chromosome
				continue;
			ReadPair rp = new ReadPair();
			if (r1.compareTo(r2)<0){		// r1 should be lower than r2
				rp.r1 = r1;
				rp.r2 = r2;
			}
			else{
				rp.r1 = r2;
				rp.r2 = r1;
			}
			low.add(rp);
			ReadPair rp2 = new ReadPair();
			rp2.r1 = rp.r1;
			rp2.r2 = rp.r2;
			high.add(rp2);
		}
		low.trimToSize();
		high.trimToSize();
		Collections.sort(low, new Comparator<ReadPair>(){
            public int compare(ReadPair o1, ReadPair o2) {
                return o1.compareRead1(o2);
            }
        });
		Collections.sort(high, new Comparator<ReadPair>(){
            public int compare(ReadPair o1, ReadPair o2) {
                return o1.compareRead2(o2);
            }
        });
		ArrayList<Point> lowEnds = new ArrayList<Point>();
		for (ReadPair r:low)
			lowEnds.add(r.r1);
		lowEnds.trimToSize();
		ArrayList<Point> highEnds = new ArrayList<Point>();
		for (ReadPair r:high)
			highEnds.add(r.r2);
		highEnds.trimToSize();
		
		System.out.println("Loaded "+ highEnds.size() +" ChIA-PET read pairs: "+CommonUtils.timeElapsed(tic));
		System.out.println();
		
		
		// load TF sites
		ArrayList<String> tfs = CommonUtils.readTextFile(Args.parseString(args, "tf_sites", null));
		ArrayList<ArrayList<Point>> allPeaks = new ArrayList<ArrayList<Point>>();
		for (int i=0;i<tfs.size();i++){
			try{
				ArrayList<Point> ps = new ArrayList<Point>();
				ps.addAll(GPSParser.parseGPSOutput(tfs.get(i), genome));
				ps.trimToSize();
				Collections.sort(ps);
				allPeaks.add(ps);
				System.out.println("Loaded "+tfs.get(i));
			}
			catch (IOException e){
				System.out.println(tfs.get(i)+" does not have a valid GPS/GEM event call file.");
				e.printStackTrace(System.err);
				System.exit(1);
			}
		}
		allPeaks.trimToSize();

		// load histone mark or DHS, SE regions
		ArrayList<String> hms = CommonUtils.readTextFile(Args.parseString(args, "regions", null));
		ArrayList<List<Region>> allRegions = new ArrayList<List<Region>>();
		for (int i=0;i<hms.size();i++){
			allRegions.add(CommonUtils.load_BED_regions(genome, hms.get(i)).car());
			System.out.println("Loaded "+hms.get(i));
		}
		allRegions.trimToSize();
		System.out.println();
	
		// load other Interaction calls
		ArrayList<String> lines = CommonUtils.readTextFile(Args.parseString(args, "germ", null));
		HashMap<Point,ArrayList<Point>> germTss2distals = new HashMap<Point,ArrayList<Point>>();
		ArrayList<Point> germTss = new ArrayList<Point>();
		for (String l: lines){		// each line is a call
			String f[] = l.split("\t");
			Point t = new Region(genome, f[3].replace("chr", ""), Integer.parseInt(f[4]), Integer.parseInt(f[5])).getMidpoint();
			if (!germTss2distals.containsKey(t))
				germTss2distals.put(t, new ArrayList<Point>());
			germTss2distals.get(t).add(new Region(genome, f[0].replace("chr", ""), Integer.parseInt(f[1]), Integer.parseInt(f[2])).getMidpoint());
		}
		germTss.addAll(germTss2distals.keySet());
		germTss.trimToSize();
		Collections.sort(germTss);
		
		lines = CommonUtils.readTextFile(Args.parseString(args, "mango", null));
		HashMap<Point, ArrayList<Point>> a2bs = new HashMap<Point, ArrayList<Point>>();
		HashMap<Point, ArrayList<Point>> b2as = new HashMap<Point, ArrayList<Point>>();
		for (String l: lines){		// each line is a call
			String f[] = l.split("\t");
			Point a = new Region(genome, f[0].replace("chr", ""), Integer.parseInt(f[1]), Integer.parseInt(f[2])).getMidpoint();
			Point b = new Region(genome, f[3].replace("chr", ""), Integer.parseInt(f[4]), Integer.parseInt(f[5])).getMidpoint();
			if (!a2bs.containsKey(a))
				a2bs.put(a, new ArrayList<Point>());
			a2bs.get(a).add(b);
			if (!b2as.containsKey(b))
				b2as.put(b, new ArrayList<Point>());
			b2as.get(b).add(a);
		}
		ArrayList<Point> aPoints = new ArrayList<Point>();
		aPoints.addAll(a2bs.keySet());
		aPoints.trimToSize();
		Collections.sort(aPoints);
		ArrayList<Point> bPoints = new ArrayList<Point>();
		bPoints.addAll(b2as.keySet());
		bPoints.trimToSize();
		Collections.sort(bPoints);
		System.out.println("Loaded all the data, "+CommonUtils.timeElapsed(tic0));
		
		// find dense read cluster for each gene
		// Two big chunks of codes for 1) Distal--TSS or 2) TSS--Distal, according to their coordinates
		ArrayList<String> geneList = new ArrayList<String>();
		geneList.addAll(gene2tss.keySet());
		geneList.trimToSize();
		ArrayList<Interaction> interactions = new ArrayList<Interaction>();
		tic = System.currentTimeMillis();
		
		for (int id=0;id<geneList.size();id++){
			String g = geneList.get(id);
			System.out.print(g+" ");
			TreeSet<StrandedPoint> coords = gene2tss.get(g);
			// if a gene has multiple TSSs, use the center position
			int count = coords.size();
			Point centerPoint = null;
			for (StrandedPoint p:coords){
				if (count<coords.size()/2)
					break;
				else{
					centerPoint = p;
					count--;
				}
			}
			
		/** 1) For TSS at higher coordinates, distal anchors are at lower coordinates */
			
			// get the distal ends, merge nearby read pairs
			Region tssRegion = centerPoint.expand(tss_radius);
			Region excludeRegion = centerPoint.expand(tss_exclude);
			int exStart = excludeRegion.getStart();
			int exEnd = excludeRegion.getEnd();
			ArrayList<Integer> idx = CommonUtils.getPointsWithinWindow(highEnds, tssRegion);
			if (!idx.isEmpty()){
//				System.out.println("\n"+g+"\tCluster reads, "+CommonUtils.timeElapsed(tic));
//				tic = System.currentTimeMillis();
				
				ArrayList<ReadPair> rps = new ArrayList<ReadPair> ();
				for (int i: idx){
					ReadPair rp = high.get(i);
					int pos = rp.r1.getLocation();
					if(pos<exStart || pos>exEnd)
						rps.add(rp);
				}
				Collections.sort(rps, new Comparator<ReadPair>(){
		            public int compare(ReadPair o1, ReadPair o2) {
		                return o1.compareRead1(o2);
		            }
		        });
				ArrayList<ReadPairCluster> rpcs = new ArrayList<ReadPairCluster>();
				int current = -100000;
				ReadPairCluster c = new ReadPairCluster();
				// merge reads
				for (ReadPair rp: rps){
					if (rp.r1.getLocation()-current>read_merge_dist){	// a big gap
						if (c.reads.size()>=2){
							rpcs.add(c);
						}
						c = new ReadPairCluster();
					}
					c.addReadPair(rp);
					current = rp.r1.getLocation();
				}
				if (c.reads.size()>=2){		// finish up the last cluster
					rpcs.add(c);
				}
				
				// test whether to merge clusters
				// if two nearby clusters are within the cluster_merge_dist, make a new cluster that covers both distal regions
				// also try to expand the TSS regions a little to see if we can merge more nearby reads
	
//				System.out.println("Merge clusters, "+CommonUtils.timeElapsed(tic));
//				tic = System.currentTimeMillis();
				
	//			System.out.println("\nDistal region at lower coord than TSS\nBefore merging, number of clusters = "+rpcs.size());
				ArrayList<ReadPairCluster> toRemoveClusters = new ArrayList<ReadPairCluster>();
				for (int i=1; i<rpcs.size();i++){
					ReadPairCluster c1=rpcs.get(i-1);
					ReadPairCluster c2=rpcs.get(i);
					int dist = Math.min(c1.r2min+c1.r2max-c1.r1min-c1.r1max, c2.r2min+c2.r2max-c2.r1min-c2.r1max)/2;
					int cluster_merge_dist = Math.min(max_cluster_merge_dist, read_merge_dist + (int)Math.sqrt(dist) * distance_factor);
					double d = Math.min(c1.getDensity(cluster_merge_dist), c2.getDensity(cluster_merge_dist));
					if (c2.r1min - c1.r1max < cluster_merge_dist){
						Region distalRegionMerged = new Region(centerPoint.getGenome(), centerPoint.getChrom(), c1.r1min, c2.r1max);
						idx = CommonUtils.getPointsWithinWindow(lowEnds, distalRegionMerged);
						// TSS can expand at the high end, but the low end is dependent on the distal read positions
						int tssStart = tssRegion.getMidpoint().getLocation()-tss_exclude - Math.max(c1.r1max, c2.r1max);
						tssStart = Math.min(Math.max(tssStart,0), cluster_merge_dist);
						Region tssExpanded = tssRegion.expand(tssStart, cluster_merge_dist);
						int start = tssExpanded.getStart();
						int end = tssExpanded.getEnd();
						rps = new ArrayList<ReadPair> ();
						for (int ii: idx){
							ReadPair rp2 = low.get(ii);
							int pos = rp2.r2.getLocation();
							if(pos>=start && pos<=end)
								rps.add(rp2);
						}
						Collections.sort(rps, new Comparator<ReadPair>(){
				            public int compare(ReadPair o1, ReadPair o2) {
				                return o1.compareRead2(o2);
				            }
				        });
						// coord2 are the read2 positions of those pairs that read1 is in the merged region
						ArrayList<Integer> coord2 = new ArrayList<Integer>();	
						for (ReadPair rp: rps)
							coord2.add(rp.r2.getLocation());
						int idxMin = Collections.binarySearch(coord2, Math.min(c1.r2min, c2.r2min));
						if (idxMin<0){
							System.out.println("c1.r2min, c2.r2min: " + c1.r2min + "," + c2.r2min);
							idxMin = Math.max(0, -idxMin-1);
						}
						int idxMax = Collections.binarySearch(coord2, Math.max(c1.r2max, c2.r2max));
						if (idxMax<0){
							System.out.println("c1.r2max, c2.r2max: " + c1.r2max + "," + c2.r2max);
							idxMax = Math.min(-idxMax-1, rps.size()-1);
						}
						ReadPairCluster cNew = new ReadPairCluster();
						for (int ii=idxMin; ii<=idxMax;ii++){
							cNew.addReadPair(rps.get(ii));
						}
						double best_density = cNew.getDensity(cluster_merge_dist);
						int best_idxMin = idxMin;
						int best_idxMax = idxMax;
						// expand to the lower end first
						for (int ii=idxMin-1; ii>=0; ii--){
							cNew.addReadPair(rps.get(ii));
							double new_density = cNew.getDensity(cluster_merge_dist);
	//						System.out.println("New density "+new_density);
							if (new_density>best_density){
	//							System.out.println("Better density "+best_density+"-->"+new_density);
								best_idxMin = ii;
								best_density = new_density;
							}
						}
						// update to the best set of readpairs after expanding the lower end
						cNew = new ReadPairCluster();		
						for (int ii=best_idxMin; ii<=best_idxMax;ii++){
							cNew.addReadPair(rps.get(ii));
						}
						
						// expand to the higher end
	//					System.out.println("Best density "+best_density);
						for (int ii=idxMax+1; ii<rps.size(); ii++){
							cNew.addReadPair(rps.get(ii));
							double new_density = cNew.getDensity(cluster_merge_dist);
	//						System.out.println("New density "+new_density);
							if (new_density>best_density){
	//							System.out.println("Better density "+best_density+"-->"+new_density);
								best_idxMax = ii;
								best_density = new_density;
							}
						}
						if (best_density>d){	// better density, merge the regions
	//						System.out.println("Merged "+c1.toString()+" and \n\t"+c2.toString()+" to "+cNew.toString());
							// update to the best set of readpairs
							cNew = new ReadPairCluster();
							for (int ii=best_idxMin; ii<=best_idxMax;ii++){
								cNew.addReadPair(rps.get(ii));
							}
							rpcs.set(i, cNew);
							toRemoveClusters.add(c1);
						}
					}
				}	// for each pair of nearby clusters
				rpcs.removeAll(toRemoveClusters);
	//			System.out.println("After merging,  number of clusters = "+rpcs.size());
	
				// report
				for (ReadPairCluster cc: rpcs){
					Interaction it = new Interaction();
					interactions.add(it);
					it.geneSymbol = g;
					it.geneID = id;
					it.tss = centerPoint;
					it.tssRegion = new Region(centerPoint.getGenome(), centerPoint.getChrom(), cc.r2min, cc.r2max);
					it.distalRegion = new Region(centerPoint.getGenome(), centerPoint.getChrom(), cc.r1min, cc.r1max);
					it.distalPoint = it.distalRegion.getMidpoint();
					it.count = cc.reads.size();
					int cluster_merge_dist = Math.min(max_cluster_merge_dist, 
							read_merge_dist + (int)Math.sqrt(Math.abs(it.tss.offset(it.distalPoint))) * distance_factor);
					it.density = cc.getDensity(cluster_merge_dist);
				}
				rpcs = null;
//				System.gc();

//				System.out.println("Done merging clusters, "+CommonUtils.timeElapsed(tic));
//				tic = System.currentTimeMillis();
			}
			
		/** 2) For TSS at lower coordinates */
			
			// get the distal ends, merge nearby read pairs
			idx = CommonUtils.getPointsWithinWindow(lowEnds, tssRegion);
			if (!idx.isEmpty()){
				ArrayList<ReadPair> rps = new ArrayList<ReadPair> ();
				for (int i: idx){
					ReadPair rp = low.get(i);
					int pos = rp.r2.getLocation();
					if(pos<exStart || pos>exEnd)
						rps.add(rp);
				}
				Collections.sort(rps, new Comparator<ReadPair>(){
		            public int compare(ReadPair o1, ReadPair o2) {
		                return o1.compareRead2(o2);
		            }
		        });
				ArrayList<ReadPairCluster> rpcs = new ArrayList<ReadPairCluster>();
				int current = -100000;
				ReadPairCluster c = new ReadPairCluster();
				for (ReadPair rp: rps){
					if (rp.r2.getLocation()-current>read_merge_dist){	// a big gap
						if (c.reads.size()>=2){
							rpcs.add(c);
						}
						c = new ReadPairCluster();
					}
					c.addReadPair(rp);
					current = rp.r2.getLocation();
				}
				if (c.reads.size()>=2){		// finish up the last cluster
					rpcs.add(c);
				}
				
				// test whether to merge clusters
//				System.out.println("Start merge clusters, "+CommonUtils.timeElapsed(tic));
//				tic = System.currentTimeMillis();
				
				ArrayList<ReadPairCluster> toRemoveClusters = new ArrayList<ReadPairCluster>();
				for (int i=1; i<rpcs.size();i++){
					ReadPairCluster c1=rpcs.get(i-1);
					ReadPairCluster c2=rpcs.get(i);
					// cluster_merge_dist is dependent on the distance between two anchor regions
					int dist = Math.min(c1.r2min+c1.r2max-c1.r1min-c1.r1max, c2.r2min+c2.r2max-c2.r1min-c2.r1max)/2;
					int cluster_merge_dist = Math.min(max_cluster_merge_dist, read_merge_dist + (int)Math.sqrt(dist) * distance_factor);
					double d = Math.min(c1.getDensity(cluster_merge_dist), c2.getDensity(cluster_merge_dist));
					if (c2.r2min - c1.r2max < cluster_merge_dist){
//						System.out.println("Close enough, "+CommonUtils.timeElapsed(tic));
//						tic = System.currentTimeMillis();
						Region rMerged = new Region(centerPoint.getGenome(), centerPoint.getChrom(), c1.r2min, c2.r2max);
						idx = CommonUtils.getPointsWithinWindow(highEnds, rMerged);
//						System.out.println("Got points, "+CommonUtils.timeElapsed(tic));
//						tic = System.currentTimeMillis();
						// TSS can expand at the lower end, but the high end is dependent on the distal read positions
						int tssEnd = Math.min(c1.r2min, c2.r2min) - (tssRegion.getMidpoint().getLocation()+tss_exclude);
						tssEnd = Math.min(Math.max(tssEnd,0), cluster_merge_dist);
						Region tssExpanded = tssRegion.expand(cluster_merge_dist, tssEnd);	
						int start = tssExpanded.getStart();
						int end = tssExpanded.getEnd();
						rps = new ArrayList<ReadPair> ();
						for (int ii: idx){
							ReadPair rp2 = high.get(ii);
							int pos = rp2.r1.getLocation();
							if(pos>=start && pos<=end)
								rps.add(rp2);
						}
						Collections.sort(rps, new Comparator<ReadPair>(){
				            public int compare(ReadPair o1, ReadPair o2) {
				                return o1.compareRead1(o2);
				            }
				        });
//						System.out.println("Sorted read1, "+CommonUtils.timeElapsed(tic));
//						tic = System.currentTimeMillis();
						// coord1 are the read1 positions of those pairs that read2 is in the merged region
						ArrayList<Integer> coord1 = new ArrayList<Integer>();	
						for (ReadPair rp: rps)
							coord1.add(rp.r1.getLocation());
						int idxMin = Collections.binarySearch(coord1, Math.min(c1.r1min, c2.r1min));
						if (idxMin<0){
							System.out.println("c1.r1min, c2.r1min: " + c1.r1min + "," + c2.r1min);
							idxMin = Math.max(0, -idxMin-1);
						}
						int idxMax = Collections.binarySearch(coord1, Math.max(c1.r1max, c2.r1max));
						if (idxMax<0){
							System.out.println("c1.r1max, c2.r1max: " + c1.r1max + "," + c2.r1max);
							idxMax = Math.min(-idxMax-1, rps.size()-1);
						}
						ReadPairCluster cNew = new ReadPairCluster();
						for (int ii=idxMin; ii<=idxMax;ii++){
							cNew.addReadPair(rps.get(ii));
						}
						double best_density = cNew.getDensity(cluster_merge_dist);
						int best_idxMin = idxMin;
						int best_idxMax = idxMax;

//						System.out.println("To expand lower TSS, "+CommonUtils.timeElapsed(tic));
//						tic = System.currentTimeMillis();
						
						// expand to the lower end first
						for (int ii=idxMin-1; ii>=0; ii--){
							cNew.addReadPair(rps.get(ii));
							double new_density = cNew.getDensity(cluster_merge_dist);
	//						System.out.println("New density "+new_density);
							if (new_density>best_density){
	//							System.out.println("Better density "+best_density+"-->"+new_density);
								best_idxMin = ii;
								best_density = new_density;
							}
						}
						// update to the best set of readpairs after expanding the lower end
						cNew = new ReadPairCluster();		
						for (int ii=best_idxMin; ii<=best_idxMax;ii++){
							cNew.addReadPair(rps.get(ii));
						}
						
						// expand to the higher end
//						System.out.println("To expand higher TSS, "+CommonUtils.timeElapsed(tic));
//						tic = System.currentTimeMillis();
	//					System.out.println("Best density "+best_density);
						for (int ii=idxMax+1; ii<rps.size(); ii++){
							cNew.addReadPair(rps.get(ii));
							double new_density = cNew.getDensity(cluster_merge_dist);
	//						System.out.println("New density "+new_density);
							if (new_density>best_density){
	//							System.out.println("Better density "+best_density+"-->"+new_density);
								best_idxMax = ii;
								best_density = new_density;
							}
						}
						if (best_density>d){	// better density, merge the regions
	//						System.out.println("Merged "+c1.toString()+" and \n\t"+c2.toString()+" to "+cNew.toString());
							// update to the best set of readpairs
							cNew = new ReadPairCluster();
							for (int ii=best_idxMin; ii<=best_idxMax;ii++){
								cNew.addReadPair(rps.get(ii));
							}
							rpcs.set(i, cNew);
							toRemoveClusters.add(c1);
						}
					}
				}	// for each pair of nearby clusters
				rpcs.removeAll(toRemoveClusters);
	//			System.out.println("After merging,  number of clusters = "+rpcs.size());

//				System.out.println("Got clusters, "+CommonUtils.timeElapsed(tic));
//				tic = System.currentTimeMillis();
				for (ReadPairCluster cc: rpcs){
					Interaction it = new Interaction();
					interactions.add(it);
					it.geneSymbol = g;
					it.geneID = id;
					it.tss = centerPoint;
					it.tssRegion = new Region(centerPoint.getGenome(), centerPoint.getChrom(), cc.r1min, cc.r1max);
					it.distalRegion = new Region(centerPoint.getGenome(), centerPoint.getChrom(), cc.r2min, cc.r2max);
					it.distalPoint = it.distalRegion.getMidpoint();
					it.count = cc.reads.size();
					int cluster_merge_dist = Math.min(max_cluster_merge_dist, 
							read_merge_dist + (int)Math.sqrt(Math.abs(it.tss.offset(it.distalPoint))) * distance_factor);
					it.density = cc.getDensity(cluster_merge_dist);
				}
				interactions.trimToSize();
				rpcs = null;
//				System.gc();
//				System.out.println(String.format("Done merging clusters, n=%d, %s, total %s", 
//						interactions.size(), CommonUtils.timeElapsed(tic), CommonUtils.timeElapsed(tic0)));
//				tic = System.currentTimeMillis();
			}  
		}// for each gene

		// consolidate interaction anchors (for nearby TSSs)
		System.out.println("\n\nConsolidate nearby Tss interactions, "+CommonUtils.timeElapsed(tic0));
		TreeMap<Point, ArrayList<Interaction>> tss2it = new TreeMap<Point, ArrayList<Interaction>>();
		for (Interaction it: interactions){
			if (!tss2it.containsKey(it.tss))
				tss2it.put(it.tss, new ArrayList<Interaction>());
			tss2it.get(it.tss).add(it);
		}
		ArrayList<Point> tssList = new ArrayList<Point>();
		tssList.addAll(tss2it.keySet());
		for (int i=1;i<tssList.size();i++){
			Point p1 = tssList.get(i-1);
			Point p2 = tssList.get(i);
			if (p2.getChrom().equals(p1.getChrom()) && p2.offset(p1)<max_cluster_merge_dist){
				ArrayList<Interaction> its1 = tss2it.get(p1);
				ArrayList<Interaction> its2 = tss2it.get(p2);
				for (Interaction it1:its1){
					int start1 = it1.distalRegion.getStart();
					int end1 = it1.distalRegion.getEnd();
					for (Interaction it2:its2){
						if (it1.distalRegion.equals(it2.distalRegion))
							continue;
						int start2 = it2.distalRegion.getStart();
						int end2 = it2.distalRegion.getEnd();
//						int overlapWidth = Math.min(end1, end2)-Math.max(start2, start1);
						int overlapWidth = it2.distalRegion.getOverlapSize(it1.distalRegion);
						// if the two distal regions are highly overlapped
						if (overlapWidth>it1.distalRegion.getWidth()*overlap_ratio && overlapWidth>it2.distalRegion.getWidth()*overlap_ratio){
							//TODO: it would be more accurate to update the read pair count with the new distal region
							Region r = new Region(it1.distalRegion.getGenome(), it1.distalRegion.getChrom(), 
									Math.min(start2, start1), Math.max(end1, end2));
							it1.distalRegion  = r;
							it2.distalRegion  = r;
							Point mid = r.getMidpoint();
							if (mid.distance(it1.distalPoint)<mid.distance(it2.distalPoint))
								it2.distalPoint = it1.distalPoint;
							else
								it1.distalPoint = it2.distalPoint;
						}
					}
				}
			}		
		}

		
		// for each TSS and its distal anchor, count how many read pairs lead to other distal anchors
		// it is like finding the other two edges of the triangle. The count is the sum (over all distal anchors) 
		// of min connecting read count (btw tss-distal2 and distal-distal2).
		System.out.println("\nCount indirect / triangle read pairs, "+CommonUtils.timeElapsed(tic0));
		HashMap<String, ArrayList<Interaction>> gene2it = new HashMap<String, ArrayList<Interaction>>();
		for (Interaction it: interactions){
			if (!gene2it.containsKey(it.geneSymbol))
				gene2it.put(it.geneSymbol, new ArrayList<Interaction>());
			gene2it.get(it.geneSymbol).add(it);
		}
		for (String gene: gene2it.keySet()){
			ArrayList<Interaction> its = gene2it.get(gene);
			for (Interaction it: its){
//				System.out.println(it);
				ArrayList<Integer> indirectCounts = new ArrayList<Integer>();
				Region distal = it.distalRegion.expand(chiapet_radius, chiapet_radius);
				for (Interaction it2: its){
					if (it==it2)
						continue;
					int start = it2.distalRegion.getStart()-chiapet_radius;
					int end = it2.distalRegion.getEnd()+chiapet_radius;
					if (distal.getStart() > start){		// if it1 has higher coord, select by high and then low
						ArrayList<Integer> idx = CommonUtils.getPointsWithinWindow(highEnds, distal);
						ArrayList<ReadPair> rps = new ArrayList<ReadPair> ();
						for (int ii: idx){
							ReadPair rp2 = high.get(ii);
							int pos = rp2.r1.getLocation();
							if(pos>=start && pos<=end)
								rps.add(rp2);
						}
						indirectCounts.add(Math.min(rps.size(), it2.count));		// distal2--distal , it2_count
					}
					else{		// select by low and then high
						ArrayList<Integer> idx = CommonUtils.getPointsWithinWindow(lowEnds, distal);
						ArrayList<ReadPair> rps = new ArrayList<ReadPair> ();
						for (int ii: idx){
							ReadPair rp2 = low.get(ii);
							int pos = rp2.r2.getLocation();
							if(pos>=start && pos<=end)
								rps.add(rp2);
						}
						indirectCounts.add(Math.min(rps.size(), it2.count));		// distal--distal2 , it2_count
					}
				}
				int sum = 0;
				for (int c: indirectCounts)
					sum += c;
				it.indirectCount = sum;
				
//				System.out.println();
//				for (int c: indirectCounts)
//					System.out.print(c);
//				System.out.println();
			}
		}
	
		
		System.out.println("\n\nAnnotate and report, "+CommonUtils.timeElapsed(tic0));
		// report the interactions and annotations
		// annotate the proximal and distal anchors with TF and HM and regions
		StringBuilder sb = new StringBuilder();
		for (Interaction it: interactions){
			ArrayList<Integer> isOverlapped = new ArrayList<Integer>();
			Point mid = it.distalRegion.getMidpoint();
			int radius = it.distalRegion.getWidth()/2+chiapet_radius;
			for (ArrayList<Point> ps : allPeaks){
				ArrayList<Point> p = CommonUtils.getPointsWithinWindow(ps, mid, radius);
				isOverlapped.add(p.size());
				if (!it.isTfAnchord && !p.isEmpty()){	// if the distal point has not been anchored by a TF site
					it.distalPoint = p.get(p.size()/2);
					it.isTfAnchord = true;
				}
			}
			for (List<Region> rs: allRegions){
				isOverlapped.add(CommonUtils.getRegionsOverlapsWindow(rs, it.distalRegion, chiapet_radius).size());
			}
			// proximal
			for (ArrayList<Point> ps : allPeaks){
				ArrayList<Point> p = CommonUtils.getPointsWithinWindow(ps, it.tss, chiapet_radius);
				isOverlapped.add(p.size());
			}
			Region tssRegion = it.tss.expand(chiapet_radius);
			for (List<Region> rs: allRegions){
				isOverlapped.add(CommonUtils.getRegionsOverlapsWindow(rs, tssRegion, chiapet_radius).size());
			}
			// print out TF and region overlaps
			sb.append(it.toString()).append("\t");
			for (int b: isOverlapped)
				sb.append(b).append("\t");
			
			// print ChIA-PET call overlap info
			Point tssPoint = it.tss;
			int tssHalfWidth = it.tssRegion.getWidth()/2+chiapet_radius;
			int distalHalfWidth = it.distalRegion.getWidth()/2+chiapet_radius;
			Point distalPoint = it.distalPoint;
			Point tssLeft = new Point(genome, tssPoint.getChrom(), tssPoint.getLocation()-2000);
			Point tssRight = new Point(genome, tssPoint.getChrom(), tssPoint.getLocation()+2000);
			
			// GERM
			int index = Collections.binarySearch(germTss,  tssLeft);
			if( index < 0 )  							// if key not found
				index = -(index+1); 
			int indexRight = Collections.binarySearch(germTss,  tssRight);
			if( indexRight < 0 )  							// if key not found
				indexRight = -(indexRight+1); 
			// if key match found, continue to search ( binarySearch() give undefined index with multiple matches)
			boolean isGermOverlapped = false;
			indexRange: for (int i=index-1;i<=indexRight+2;i++){
				if (i<0 || i>=germTss.size())
					continue;
				try{
					Point tt = germTss.get(i);
					if (tt.distance(tssPoint)<=tssHalfWidth){ 
						if (!germTss2distals.containsKey(tt))
							continue;
						for (Point d:germTss2distals.get(tt)){
							if (d.distance(distalPoint)<=distalHalfWidth){
								isGermOverlapped = true;
								break indexRange;
							}
						}
					}
				}
				catch (IllegalArgumentException e){	// ignore								
				}								
			}
			sb.append(isGermOverlapped?"1\t":"0\t");
			
			// Mango
			// aPoints and bPoints are the midPoint of the two anchorRegions
			index = Collections.binarySearch(aPoints,  tssLeft);
			if( index < 0 )  							// if key not found
				index = -(index+1); 
			indexRight = Collections.binarySearch(aPoints,  tssRight);
			if( indexRight < 0 )  							// if key not found
				indexRight = -(indexRight+1); 
			// if key match found, continue to search ( binarySearch() give undefined index with multiple matches)
			boolean isMangoOverlapped = false;
			indexA: for (int i=index-1;i<=indexRight+2;i++){
				if (i<0 || i>=aPoints.size())
					continue;
				try{
					Point a = aPoints.get(i);
					if (a.distance(tssPoint)<=tssHalfWidth){ 
						if (!a2bs.containsKey(a))
							continue;
						for (Point b:a2bs.get(a)){
							if (b.distance(distalPoint)<=distalHalfWidth){
								isMangoOverlapped = true;
								break indexA;
							}
						}
					}
				}
				catch (IllegalArgumentException e){	// ignore								
				}
			}
			if (isMangoOverlapped)
				sb.append("1\t");
			else{
				index = Collections.binarySearch(bPoints,  tssLeft);
				if( index < 0 )  							// if key not found
					index = -(index+1); 
				indexRight = Collections.binarySearch(bPoints,  tssRight);
				if( indexRight < 0 )  							// if key not found
					indexRight = -(indexRight+1); 
				// if key match found, continue to search ( binarySearch() give undefined index with multiple matches)
				isMangoOverlapped = false;
				indexB: for (int i=index-1;i<=indexRight+2;i++){
					if (i<0 || i>=bPoints.size())
						continue;
					try{
						Point b = bPoints.get(i);
						if (b.distance(tssPoint)<=tssHalfWidth){ 
							if (!b2as.containsKey(b))
								continue;
							for (Point a:b2as.get(b)){
								if (a.distance(distalPoint)<=distalHalfWidth){
									isMangoOverlapped = true;
									break indexB;
								}
							}
						}
					}
					catch (IllegalArgumentException e){	// ignore								
					}								
				}
				sb.append(isMangoOverlapped?"1\t":"0\t");
			}
			
			CommonUtils.replaceEnd(sb, '\n');
			
		}
		CommonUtils.writeFile(Args.parseString(args, "out", "Result")+".readClusters.txt", sb.toString());
		
		// output BEDPE format
		sb = new StringBuilder();
		for (Interaction it: interactions){
			Region distalLocal = it.distalRegion.expand(read_merge_dist, read_merge_dist);
			Region tssLocal = it.tssRegion.expand(read_merge_dist, read_merge_dist);
			int distalLocalCount, tssLocalCount;
			if (it.distalPoint.offset(it.tss)<0){	// distal is in lower coord
				distalLocalCount = CommonUtils.getPointsWithinWindow(lowEnds, distalLocal).size();
				tssLocalCount = CommonUtils.getPointsWithinWindow(highEnds, tssLocal).size();
				sb.append(String.format("%s\t%s\t%d\t%d\t%d\n", distalLocal.toBED(), tssLocal.toBED(),
						it.count, distalLocalCount, tssLocalCount));
			}
			else{	// distal is in higher coord
				distalLocalCount = CommonUtils.getPointsWithinWindow(highEnds, distalLocal).size();
				tssLocalCount = CommonUtils.getPointsWithinWindow(lowEnds, tssLocal).size();
				sb.append(String.format("%s\t%s\t%d\t%d\t%d\n", tssLocal.toBED(), distalLocal.toBED(),
						it.count, tssLocalCount, distalLocalCount));
			}
		}
		CommonUtils.writeFile(Args.parseString(args, "out", "Result")+".bedpe", sb.toString());
		
		System.out.println("\n\nDone: "+CommonUtils.timeElapsed(tic0));
	}
	
	private void findTssInteractions(){
		long tic0 = System.currentTimeMillis();
		long tic = System.currentTimeMillis();
		
		boolean simple_merge = !flags.contains("density_merge");

		// load gene annotation
		ArrayList<String> lines = CommonUtils.readTextFile(Args.parseString(args, "gene_anno", null));
		ArrayList<TSS> allTSS = new ArrayList<TSS>();
//		TreeMap<String, TreeSet<StrandedPoint>> gene2tss = new TreeMap<String, TreeSet<StrandedPoint>>();
		for (int i=0;i<lines.size();i++){
			String t = lines.get(i);
			if (t.startsWith("#"))
				continue;
			String f[] = t.split("\t");
			String symbol = f[12];
			String chr = f[2].replace("chr", "");
			char strand = f[3].charAt(0);
			StrandedPoint tss = new StrandedPoint(genome, chr, Integer.parseInt(f[strand=='+'?4:5]), strand);
//			if (!gene2tss.containsKey(symbol))
//				gene2tss.put(symbol, new TreeSet<StrandedPoint>());
//			gene2tss.get(symbol).add(tss);
			allTSS.add(new TSS(symbol, tss,i));
		}
		allTSS.trimToSize();
		Collections.sort(allTSS);
		
		// merge nearby TSSs into TSS_clusters
		ArrayList<ArrayList<TSS>> clusters = new ArrayList<ArrayList<TSS>>();		
		ArrayList<TSS> cluster = new ArrayList<TSS>();
		cluster.add(allTSS.get(0));
		for (int i=1;i<allTSS.size();i++){
			TSS t = allTSS.get(i);
			TSS p = cluster.get(cluster.size()-1);		//previous tss
			if (t.coord.getChrom().equals(p.coord.getChrom()) 
					&& t.coord.distance(p.coord)<tss_merge_dist){
				if (cluster.get(0).coord.distance(t.coord) > tss_radius*2){	// if the cluster may be too big, stop here
					cluster.trimToSize();
					clusters.add(cluster);
					cluster = new ArrayList<TSS>();
					cluster.add(t);
				}
				else
					cluster.add(t);
			}
			else{
				cluster.trimToSize();
				clusters.add(cluster);
				cluster = new ArrayList<TSS>();
				cluster.add(t);
			}
		}
		// finish all the sites
		cluster.trimToSize();
		clusters.add(cluster);
		System.out.println("Loaded "+ allTSS.size() + " TSSs, merged to "+ clusters.size()+", "+CommonUtils.timeElapsed(tic));
		allTSS.clear();
		allTSS = null;
		
		// load read pairs
		// only use read pairs on the same chromosome, and longer than self_exclude distance
		// the left read is required to be lower than the right read; if not, flip the two
		
		// sort by each end so that we can use binary search to find matches or overlaps
		System.out.println("Loading ChIA-PET read pairs: "+CommonUtils.timeElapsed(tic));

		ArrayList<String> read_pairs = CommonUtils.readTextFile(Args.parseString(args, "read_pair", null));
		ArrayList<Point> reads = new ArrayList<Point>();	// store PET as single ends
		ArrayList<ReadPair> low = new ArrayList<ReadPair> ();	// all PET sorted by the low end
		ArrayList<ReadPair> high = new ArrayList<ReadPair> ();	// sorted by high end
		for (String s: read_pairs){
			String[] f = s.split("\t");
			Point r1 = Point.fromString(genome, f[0]);
			String r1Chrom=r1.getChrom();
			Point r2 = Point.fromString(genome, f[1]);
			reads.add(r1);
			reads.add(r2);
			// TODO: change next line if prediction cross-chrom interactions
			if (!r1Chrom.equals(r2.getChrom()))		// r1 and r2 should be on the same chromosome
				continue;
			if (r1.distance(r2)<self_exclude)	// PETs shorter than 8kb are considered self-ligation reads
				continue;
			ReadPair rp = new ReadPair();
			if (r1.compareTo(r2)<0){		// r1 should be lower than r2
				rp.r1 = r1;
				rp.r2 = r2;
			}
			else{
				rp.r1 = r2;
				rp.r2 = r1;
			}
			low.add(rp);
			ReadPair rp2 = new ReadPair();
			rp2.r1 = rp.r1;
			rp2.r2 = rp.r2;
			high.add(rp2);
		}
		low.trimToSize();
		high.trimToSize();
		Collections.sort(low, new Comparator<ReadPair>(){
            public int compare(ReadPair o1, ReadPair o2) {
                return o1.compareRead1(o2);
            }
        });
		Collections.sort(high, new Comparator<ReadPair>(){
            public int compare(ReadPair o1, ReadPair o2) {
                return o1.compareRead2(o2);
            }
        });
				
		ArrayList<Point> lowEnds = new ArrayList<Point>();
		for (ReadPair r:low)
			lowEnds.add(r.r1);
		lowEnds.trimToSize();
		ArrayList<Point> highEnds = new ArrayList<Point>();
		for (ReadPair r:high)
			highEnds.add(r.r2);
		highEnds.trimToSize();

		reads.trimToSize();
		Collections.sort(reads);
		
		System.out.println("Loaded total="+ (reads.size()/2) +", filtered="+ highEnds.size() +" ChIA-PET read pairs: "+CommonUtils.timeElapsed(tic));
		System.out.println();
		
		// one dimension clustering to define anchors (similar to GEM code)
		// TODO: use cross correlation to determine the distance to shift
		ArrayList<Region> rs0 = new ArrayList<Region>();
		ArrayList<Point> summits = new ArrayList<Point>();
		// cut the pooled reads into independent regions
		int start0=0;
		int minCount = 3;
		for (int i=1;i<reads.size();i++){
			Point p0 = reads.get(i-1);
			Point p1 = reads.get(i);
			if ((!p0.getChrom().equals(p1.getChrom())) || p1.getLocation()-p0.getLocation() > read_merge_dist){ // not same chorm, or a large enough gap to cut
				// only select region with read count larger than minimum count
				int count = i-start0;
				if (count >= minCount){
					Region r = new Region(genome, p0.getChrom(), reads.get(start0).getLocation(), reads.get(i-1).getLocation());
					rs0.add(r);
					ArrayList<Point> ps = new ArrayList<Point>();
					for (int j=start0;j<i;j++)
						ps.add(reads.get(j));
					int maxCount = 0;
					int maxIdx = -1;
					for (int j=0;j<ps.size();j++){
						Point mid = ps.get(j);
						int c = CommonUtils.getPointsWithinWindow(ps, mid, read_merge_dist).size();
						if (c>maxCount){
							maxCount = c;
							maxIdx = start0+j;
						}
					}
					summits.add(reads.get(maxIdx));
				}
				start0 = i;
			}
		}
		// the last region
		int count = reads.size()-start0;
		if (count >= minCount){
			Region r = new Region(genome, reads.get(start0).getChrom(), 
					reads.get(start0).getLocation(), reads.get(reads.size()-1).getLocation());
			rs0.add(r);
			ArrayList<Point> ps = new ArrayList<Point>();
			for (int j=start0;j<reads.size();j++)
				ps.add(reads.get(j));
			int maxCount = 0;
			int maxIdx = -1;
			for (int j=0;j<ps.size();j++){
				Point mid = ps.get(j);
				int c = CommonUtils.getPointsWithinWindow(ps, mid, read_merge_dist).size();
				if (c>maxCount){
					maxCount = c;
					maxIdx = start0+j;
				}
			}
			summits.add(reads.get(maxIdx));
		}
		reads.clear();
		reads = null;
		System.out.println("\nMerge all PETs into "+rs0.size()+" regions.");
		
		// load TF sites
		ArrayList<String> tfs = CommonUtils.readTextFile(Args.parseString(args, "tf_sites", null));
		ArrayList<ArrayList<Point>> allPeaks = new ArrayList<ArrayList<Point>>();
		for (int i=0;i<tfs.size();i++){
			try{
				ArrayList<Point> ps = new ArrayList<Point>();
				ps.addAll(GPSParser.parseGPSOutput(tfs.get(i), genome));
				ps.trimToSize();
				Collections.sort(ps);
				allPeaks.add(ps);
				System.out.println("Loaded "+tfs.get(i));
			}
			catch (IOException e){
				System.out.println(tfs.get(i)+" does not have a valid GPS/GEM event call file.");
				e.printStackTrace(System.err);
				System.exit(1);
			}
		}
		allPeaks.trimToSize();

		// load histone mark or DHS, SE regions
		ArrayList<String> hms = CommonUtils.readTextFile(Args.parseString(args, "regions", null));
		ArrayList<List<Region>> allRegions = new ArrayList<List<Region>>();
		for (int i=0;i<hms.size();i++){
			allRegions.add(CommonUtils.load_BED_regions(genome, hms.get(i)).car());
			System.out.println("Loaded "+hms.get(i));
		}
		allRegions.trimToSize();
		System.out.println();
	
		// load other Interaction calls
		lines = CommonUtils.readTextFile(Args.parseString(args, "germ", null));
		HashMap<Point,ArrayList<Point>> germTss2distals = new HashMap<Point,ArrayList<Point>>();
		ArrayList<Point> germTss = new ArrayList<Point>();
		for (String l: lines){		// each line is a call
			String f[] = l.split("\t");
			Point t = new Region(genome, f[3].replace("chr", ""), Integer.parseInt(f[4]), Integer.parseInt(f[5])).getMidpoint();
			if (!germTss2distals.containsKey(t))
				germTss2distals.put(t, new ArrayList<Point>());
			germTss2distals.get(t).add(new Region(genome, f[0].replace("chr", ""), Integer.parseInt(f[1]), Integer.parseInt(f[2])).getMidpoint());
		}
		germTss.addAll(germTss2distals.keySet());
		germTss.trimToSize();
		Collections.sort(germTss);
		
		lines = CommonUtils.readTextFile(Args.parseString(args, "mango", null));
		HashMap<Point, ArrayList<Point>> a2bs = new HashMap<Point, ArrayList<Point>>();
		HashMap<Point, ArrayList<Point>> b2as = new HashMap<Point, ArrayList<Point>>();
		for (String l: lines){		// each line is a call
			String f[] = l.split("\t");
			Point a = new Region(genome, f[0].replace("chr", ""), Integer.parseInt(f[1]), Integer.parseInt(f[2])).getMidpoint();
			Point b = new Region(genome, f[3].replace("chr", ""), Integer.parseInt(f[4]), Integer.parseInt(f[5])).getMidpoint();
			if (!a2bs.containsKey(a))
				a2bs.put(a, new ArrayList<Point>());
			a2bs.get(a).add(b);
			if (!b2as.containsKey(b))
				b2as.put(b, new ArrayList<Point>());
			b2as.get(b).add(a);
		}
		ArrayList<Point> aPoints = new ArrayList<Point>();
		aPoints.addAll(a2bs.keySet());
		aPoints.trimToSize();
		Collections.sort(aPoints);
		ArrayList<Point> bPoints = new ArrayList<Point>();
		bPoints.addAll(b2as.keySet());
		bPoints.trimToSize();
		Collections.sort(bPoints);
		
		System.out.println("Loaded all the data, "+CommonUtils.timeElapsed(tic0));
		
		/********************************************** 
		 * find dense read cluster for each TSS cluster 
		 **********************************************/
		// Two big chunks of codes for 1) Distal--TSS or 2) TSS--Distal, according to the relative coordinates of TSS vs distal
//		ArrayList<String> geneList = new ArrayList<String>();
//		geneList.addAll(gene2tss.keySet());
//		geneList.trimToSize();
		ArrayList<Interaction> interactions = new ArrayList<Interaction>();
		HashSet<ReadPair> usedPETs = new HashSet<ReadPair>();
		
		tic = System.currentTimeMillis();
		
		for (int id=0;id<clusters.size();id++){
			ArrayList<TSS> cTSS = clusters.get(id);
			StringBuilder sb = new StringBuilder();
			for (TSS t: cTSS)
				sb.append(t.symbol).append(",");
			String g = sb.toString();
//			System.out.print(g+" ");
			Point centerPoint = cTSS.get(cTSS.size()/2).coord;	// median position of the TSS cluster
			
			/** 1) For TSS at higher coordinates, distal anchors are at lower coordinates */
			
			// get the distal ends, merge nearby read pairs
			Region tssRegion = centerPoint.expand(tss_radius);
			Region excludeRegion = centerPoint.expand(self_exclude);
			int exStart = excludeRegion.getStart();
			int exEnd = excludeRegion.getEnd();
			ArrayList<Integer> idx = CommonUtils.getPointsWithinWindow(highEnds, tssRegion);
			if (idx.size()>1){
				ArrayList<ReadPair> rps = new ArrayList<ReadPair> ();
				for (int i: idx){
					ReadPair rp = high.get(i);
					int pos = rp.r1.getLocation();
					if(pos<exStart || pos>exEnd)
						rps.add(rp);
				}
				Collections.sort(rps, new Comparator<ReadPair>(){	//sort by read1 (distal end)
		            public int compare(ReadPair o1, ReadPair o2) {
		                return o1.compareRead1(o2);
		            }
		        });
				ArrayList<ReadPairCluster> rpcs = new ArrayList<ReadPairCluster>();
				int current = -100000;
				ReadPairCluster c = new ReadPairCluster();
				// merge reads
				for (ReadPair rp: rps){
					if (rp.r1.getLocation()-current>read_merge_dist){	// a big gap
						if (c.reads.size()>=2){
							rpcs.add(c);
						}
						c = new ReadPairCluster();
					}
					c.addReadPair(rp);
					current = rp.r1.getLocation();
				}
				if (c.reads.size()>=2){		// finish up the last cluster
					rpcs.add(c);
				}
				
				// test whether to merge clusters
				// if two nearby clusters are within the cluster_merge_dist, make a new cluster that covers both distal regions
				// also try to expand the TSS regions a little to see if we can merge more nearby reads
	
//				System.out.println("\nDistal region at lower coord than TSS "+g+"\nBefore merging, number of clusters = "+rpcs.size());
				ArrayList<ReadPairCluster> toRemoveClusters = new ArrayList<ReadPairCluster>();
				for (int i=1; i<rpcs.size();i++){
					ReadPairCluster c1=rpcs.get(i-1);
					ReadPairCluster c2=rpcs.get(i);
					int dist = Math.min(c1.r2min+c1.r2max-c1.r1min-c1.r1max, c2.r2min+c2.r2max-c2.r1min-c2.r1max)/2;
					int cluster_merge_dist = Math.min(max_cluster_merge_dist, read_merge_dist + (int)Math.sqrt(dist) * distance_factor);
					double d = Math.min(c1.getDensity(cluster_merge_dist), c2.getDensity(cluster_merge_dist));
					if (c2.r1min - c1.r1max < cluster_merge_dist){
						if (simple_merge){
							// simply merge c1 to c2
							for (ReadPair rp2: c1.reads)
								c2.addReadPair(rp2);
							toRemoveClusters.add(c1);
							continue;
						}
						
						// more complicated merging based on density
						Region distalRegionMerged = new Region(centerPoint.getGenome(), centerPoint.getChrom(), c1.r1min, c2.r1max);
						idx = CommonUtils.getPointsWithinWindow(lowEnds, distalRegionMerged);
						// TSS can expand at the high end, but the low end is dependent on the distal read positions
						int tssStart = tssRegion.getMidpoint().getLocation()-self_exclude - Math.max(c1.r1max, c2.r1max);
						tssStart = Math.min(Math.max(tssStart,0), cluster_merge_dist);
						Region tssExpanded = tssRegion.expand(tssStart, cluster_merge_dist);
						int start = tssExpanded.getStart();
						int end = tssExpanded.getEnd();
						rps = new ArrayList<ReadPair> ();
						for (int ii: idx){
							ReadPair rp2 = low.get(ii);
							int pos = rp2.r2.getLocation();
							if(pos>=start && pos<=end)
								rps.add(rp2);
						}
						Collections.sort(rps, new Comparator<ReadPair>(){
				            public int compare(ReadPair o1, ReadPair o2) {
				                return o1.compareRead2(o2);
				            }
				        });
						// coord2 are the read2 positions of those pairs that read1 is in the merged region
						ArrayList<Integer> coord2 = new ArrayList<Integer>();	
						for (ReadPair rp: rps)
							coord2.add(rp.r2.getLocation());
						int idxMin = Collections.binarySearch(coord2, Math.min(c1.r2min, c2.r2min));
						if (idxMin<0){
							System.out.println("c1.r2min, c2.r2min: " + c1.r2min + "," + c2.r2min);
							idxMin = Math.max(0, -idxMin-1);
						}
						int idxMax = Collections.binarySearch(coord2, Math.max(c1.r2max, c2.r2max));
						if (idxMax<0){
							System.out.println("c1.r2max, c2.r2max: " + c1.r2max + "," + c2.r2max);
							idxMax = Math.min(-idxMax-1, rps.size()-1);
						}
						ReadPairCluster cNew = new ReadPairCluster();
						for (int ii=idxMin; ii<=idxMax;ii++){
							cNew.addReadPair(rps.get(ii));
						}
						double best_density = cNew.getDensity(cluster_merge_dist);
						int best_idxMin = idxMin;
						int best_idxMax = idxMax;
						// expand to the lower end first
						for (int ii=idxMin-1; ii>=0; ii--){
							cNew.addReadPair(rps.get(ii));
							double new_density = cNew.getDensity(cluster_merge_dist);
	//						System.out.println("New density "+new_density);
							if (new_density>best_density){
	//							System.out.println("Better density "+best_density+"-->"+new_density);
								best_idxMin = ii;
								best_density = new_density;
							}
						}
						// update to the best set of readpairs after expanding the lower end
						cNew = new ReadPairCluster();		
						for (int ii=best_idxMin; ii<=best_idxMax;ii++){
							cNew.addReadPair(rps.get(ii));
						}
						
						// expand to the higher end
	//					System.out.println("Best density "+best_density);
						for (int ii=idxMax+1; ii<rps.size(); ii++){
							cNew.addReadPair(rps.get(ii));
							double new_density = cNew.getDensity(cluster_merge_dist);
	//						System.out.println("New density "+new_density);
							if (new_density>best_density){
	//							System.out.println("Better density "+best_density+"-->"+new_density);
								best_idxMax = ii;
								best_density = new_density;
							}
						}
						if (best_density>d){	// better density, merge the regions
	//						System.out.println("Merged "+c1.toString()+" and \n\t"+c2.toString()+" to "+cNew.toString());
							// update to the best set of readpairs
							cNew = new ReadPairCluster();
							for (int ii=best_idxMin; ii<=best_idxMax;ii++){
								cNew.addReadPair(rps.get(ii));
							}
							rpcs.set(i, cNew);		// replace c2
							toRemoveClusters.add(c1);
						}
						else{
//							System.out.println("Not merged! "+ best_density +"<"+d);
						}
					}
				}	// for each pair of nearby clusters
				rpcs.removeAll(toRemoveClusters);
				toRemoveClusters = null;
//				System.out.println("After merging,  number of clusters = "+rpcs.size());
				
				ArrayList<ReadPairCluster> rpcs2 = splitRecursively(rpcs, false);
				if (rpcs2!=null){
					rpcs = rpcs2;
//					System.out.println(g+"\tAfter merging,  number of clusters = "+rpcs.size());
//					System.out.println("After spliting,  number of clusters = "+rpcs2.size());
					rpcs2 = null;
				}
				
				// report
				for (ReadPairCluster cc: rpcs){
					Interaction it = new Interaction();
					interactions.add(it);
					it.geneSymbol = g;
					it.geneID = id;
					it.tss = centerPoint;
					it.tssRegion = new Region(centerPoint.getGenome(), centerPoint.getChrom(), cc.r2min, cc.r2max);
					it.distalRegion = new Region(centerPoint.getGenome(), centerPoint.getChrom(), cc.r1min, cc.r1max);
					it.distalPoint = it.distalRegion.getMidpoint();
					it.count = cc.reads.size();
					int cluster_merge_dist = Math.min(max_cluster_merge_dist, 
							read_merge_dist + (int)Math.sqrt(Math.abs(it.tss.offset(it.distalPoint))) * distance_factor);
					it.density = cc.getDensity(cluster_merge_dist);
					usedPETs.addAll(cc.reads);
				}
				rpcs = null;
//				System.gc();

//				System.out.println("Done merging clusters, "+CommonUtils.timeElapsed(tic));
//				tic = System.currentTimeMillis();
			} // TSS at higher coordinates
			
			/** 2) For TSS at lower coordinates */
			
			// get the distal ends, merge nearby read pairs
			idx = CommonUtils.getPointsWithinWindow(lowEnds, tssRegion);
			if (idx.size()>1){
				ArrayList<ReadPair> rps = new ArrayList<ReadPair> ();
				for (int i: idx){
					ReadPair rp = low.get(i);
					int pos = rp.r2.getLocation();
					if(pos<exStart || pos>exEnd)
						rps.add(rp);
				}
				Collections.sort(rps, new Comparator<ReadPair>(){
		            public int compare(ReadPair o1, ReadPair o2) {
		                return o1.compareRead2(o2);
		            }
		        });
				ArrayList<ReadPairCluster> rpcs = new ArrayList<ReadPairCluster>();
				int current = -100000;
				ReadPairCluster c = new ReadPairCluster();
				for (ReadPair rp: rps){
					if (rp.r2.getLocation()-current>read_merge_dist){	// a big gap
						if (c.reads.size()>=2){
							rpcs.add(c);
						}
						c = new ReadPairCluster();
					}
					c.addReadPair(rp);
					current = rp.r2.getLocation();
				}
				if (c.reads.size()>=2){		// finish up the last cluster
					rpcs.add(c);
				}
				
				// test whether to merge clusters
//				System.out.println("Start merge clusters, "+CommonUtils.timeElapsed(tic));
//				tic = System.currentTimeMillis();
				
				ArrayList<ReadPairCluster> toRemoveClusters = new ArrayList<ReadPairCluster>();
				for (int i=1; i<rpcs.size();i++){
					ReadPairCluster c1=rpcs.get(i-1);
					ReadPairCluster c2=rpcs.get(i);
					// cluster_merge_dist is dependent on the distance between two anchor regions
					int dist = Math.min(c1.r2min+c1.r2max-c1.r1min-c1.r1max, c2.r2min+c2.r2max-c2.r1min-c2.r1max)/2;
					int cluster_merge_dist = Math.min(max_cluster_merge_dist, read_merge_dist + (int)Math.sqrt(dist) * distance_factor);
					double d = Math.min(c1.getDensity(cluster_merge_dist), c2.getDensity(cluster_merge_dist));
					if (c2.r2min - c1.r2max < cluster_merge_dist){
						if (simple_merge){
							// simply merge c1 to c2
							for (ReadPair rp2: c1.reads)
								c2.addReadPair(rp2);
							toRemoveClusters.add(c1);
							continue;
						}
						
						// complicated merge
						Region rMerged = new Region(centerPoint.getGenome(), centerPoint.getChrom(), c1.r2min, c2.r2max);
						idx = CommonUtils.getPointsWithinWindow(highEnds, rMerged);
						// TSS can expand at the lower end, but the high end is dependent on the distal read positions
						int tssEnd = Math.min(c1.r2min, c2.r2min) - (tssRegion.getMidpoint().getLocation()+self_exclude);
						tssEnd = Math.min(Math.max(tssEnd,0), cluster_merge_dist);
						Region tssExpanded = tssRegion.expand(cluster_merge_dist, tssEnd);	
						int start = tssExpanded.getStart();
						int end = tssExpanded.getEnd();
						rps = new ArrayList<ReadPair> ();
						for (int ii: idx){
							ReadPair rp2 = high.get(ii);
							int pos = rp2.r1.getLocation();
							if(pos>=start && pos<=end)
								rps.add(rp2);
						}
						Collections.sort(rps, new Comparator<ReadPair>(){
				            public int compare(ReadPair o1, ReadPair o2) {
				                return o1.compareRead1(o2);
				            }
				        });
//						System.out.println("Sorted read1, "+CommonUtils.timeElapsed(tic));
//						tic = System.currentTimeMillis();
						// coord1 are the read1 positions of those pairs that read2 is in the merged region
						ArrayList<Integer> coord1 = new ArrayList<Integer>();	
						for (ReadPair rp: rps)
							coord1.add(rp.r1.getLocation());
						int idxMin = Collections.binarySearch(coord1, Math.min(c1.r1min, c2.r1min));
						if (idxMin<0){
							System.out.println("c1.r1min, c2.r1min: " + c1.r1min + "," + c2.r1min);
							idxMin = Math.max(0, -idxMin-1);
						}
						int idxMax = Collections.binarySearch(coord1, Math.max(c1.r1max, c2.r1max));
						if (idxMax<0){
							System.out.println("c1.r1max, c2.r1max: " + c1.r1max + "," + c2.r1max);
							idxMax = Math.min(-idxMax-1, rps.size()-1);
						}
						ReadPairCluster cNew = new ReadPairCluster();
						for (int ii=idxMin; ii<=idxMax;ii++){
							cNew.addReadPair(rps.get(ii));
						}
						double best_density = cNew.getDensity(cluster_merge_dist);
						int best_idxMin = idxMin;
						int best_idxMax = idxMax;

						// expand to the lower end first
						for (int ii=idxMin-1; ii>=0; ii--){
							cNew.addReadPair(rps.get(ii));
							double new_density = cNew.getDensity(cluster_merge_dist);
							if (new_density>best_density){
								best_idxMin = ii;
								best_density = new_density;
							}
						}
						// update to the best set of readpairs after expanding the lower end
						cNew = new ReadPairCluster();		
						for (int ii=best_idxMin; ii<=best_idxMax;ii++){
							cNew.addReadPair(rps.get(ii));
						}
						
						// expand to the higher end
						for (int ii=idxMax+1; ii<rps.size(); ii++){
							cNew.addReadPair(rps.get(ii));
							double new_density = cNew.getDensity(cluster_merge_dist);
							if (new_density>best_density){
								best_idxMax = ii;
								best_density = new_density;
							}
						}
						if (best_density>d){	// better density, merge the regions
							// update to the best set of readpairs
							cNew = new ReadPairCluster();
							for (int ii=best_idxMin; ii<=best_idxMax;ii++){
								cNew.addReadPair(rps.get(ii));
							}
							rpcs.set(i, cNew);
							toRemoveClusters.add(c1);
						}
					}
				}	// for each pair of nearby clusters
				rpcs.removeAll(toRemoveClusters);
	//			System.out.println("After merging,  number of clusters = "+rpcs.size());
				ArrayList<ReadPairCluster> rpcs2 = splitRecursively(rpcs, true);
				if (rpcs2!=null){
					rpcs = rpcs2;
					rpcs2 = null;
				}
				for (ReadPairCluster cc: rpcs){
					Interaction it = new Interaction();
					interactions.add(it);
					it.geneSymbol = g;
					it.geneID = id;
					it.tss = centerPoint;
					it.tssRegion = new Region(centerPoint.getGenome(), centerPoint.getChrom(), cc.r1min, cc.r1max);
					it.distalRegion = new Region(centerPoint.getGenome(), centerPoint.getChrom(), cc.r2min, cc.r2max);
					it.distalPoint = it.distalRegion.getMidpoint();
					it.count = cc.reads.size();
					int cluster_merge_dist = Math.min(max_cluster_merge_dist, 
							read_merge_dist + (int)Math.sqrt(Math.abs(it.tss.offset(it.distalPoint))) * distance_factor);
					it.density = cc.getDensity(cluster_merge_dist);
					usedPETs.addAll(cc.reads);
				}
				rpcs = null;
//				System.out.println(String.format("Done merging clusters, n=%d, %s, total %s", 
//						interactions.size(), CommonUtils.timeElapsed(tic), CommonUtils.timeElapsed(tic0)));
//				tic = System.currentTimeMillis();
			}  // TSS at lower coordinates
		}// for each TSS_cluster
		interactions.trimToSize();
		int itTssCount = interactions.size();
		System.out.println("\n\nTss interactions ="+itTssCount+", "+CommonUtils.timeElapsed(tic0));
		
//		/** consolidate interaction anchors (for nearby TSSs) */
//		System.out.println("\n\nConsolidate nearby Tss interactions, "+CommonUtils.timeElapsed(tic0));
//		TreeMap<Point, ArrayList<Interaction>> tss2it = new TreeMap<Point, ArrayList<Interaction>>();
//		for (Interaction it: interactions){
//			if (!tss2it.containsKey(it.tss))
//				tss2it.put(it.tss, new ArrayList<Interaction>());
//			tss2it.get(it.tss).add(it);
//		}
//		ArrayList<Point> tssList = new ArrayList<Point>();
//		tssList.addAll(tss2it.keySet());
//		for (int i=1;i<tssList.size();i++){
//			Point p1 = tssList.get(i-1);
//			Point p2 = tssList.get(i);
//			if (p2.getChrom().equals(p1.getChrom()) && p2.offset(p1)<max_cluster_merge_dist){
//				ArrayList<Interaction> its1 = tss2it.get(p1);
//				ArrayList<Interaction> its2 = tss2it.get(p2);
//				for (Interaction it1:its1){
//					int start1 = it1.distalRegion.getStart();
//					int end1 = it1.distalRegion.getEnd();
//					for (Interaction it2:its2){
//						if (it1.distalRegion.equals(it2.distalRegion))
//							continue;
//						int start2 = it2.distalRegion.getStart();
//						int end2 = it2.distalRegion.getEnd();
////						int overlapWidth = Math.min(end1, end2)-Math.max(start2, start1);
//						int overlapWidth = it2.distalRegion.getOverlapSize(it1.distalRegion);
//						// if the two distal regions are highly overlapped
//						if (overlapWidth>it1.distalRegion.getWidth()*overlap_ratio && overlapWidth>it2.distalRegion.getWidth()*overlap_ratio){
//							//TODO: it would be more accurate to update the read pair count with the new distal region
//							Region r = new Region(it1.distalRegion.getGenome(), it1.distalRegion.getChrom(), 
//									Math.min(start2, start1), Math.max(end1, end2));
//							it1.distalRegion  = r;
//							it2.distalRegion  = r;
//							Point mid = r.getMidpoint();
//							if (mid.distance(it1.distalPoint)<mid.distance(it2.distalPoint))
//								it2.distalPoint = it1.distalPoint;
//							else
//								it1.distalPoint = it2.distalPoint;
//						}
//					}
//				}
//			}		
//		}
		
		/** call interactions from non-TSS anchors. */
		// We only need to consider nonTSS vs nonTSS, because if any of them involves a TSS, it should have been called already.
		System.out.println("\nCall interactions from non-TSS anchors.");
		// get the non-TSS anchors by removing those that overlap with TSS anchors
		ArrayList<Point> tssPoints = new ArrayList<Point>();
		HashSet<Point> tssPts = new HashSet<Point>();
		for (int id=0;id<clusters.size();id++){
			ArrayList<TSS> cTSS = clusters.get(id);
			Point centerPoint = cTSS.get(cTSS.size()/2).coord;
			tssPts.add(centerPoint);
		}
		tssPoints.addAll(tssPts);
		tssPts.clear();
		tssPts = null;
		tssPoints.trimToSize();
		Collections.sort(tssPoints);
		
		System.out.println("After merging, "+rs0.size()+" regions.");
		ArrayList<Integer> tssIdx = new ArrayList<Integer>();
		for (int i=0;i<summits.size();i++){
			if (!CommonUtils.getPointsWithinWindow(tssPoints, summits.get(i), tss_radius).isEmpty())
				tssIdx.add(i);
		}
		for (int i=tssIdx.size()-1; i>=0; i--){
			summits.remove(tssIdx.get(i).intValue());
			rs0.remove(tssIdx.get(i).intValue());
		}
		System.out.println("Remove TSS regions, "+rs0.size()+" regions.");
		tssIdx.clear();
		tssIdx = null;
		// pre-screen for the regions that overlap with long PETs
		ArrayList<Integer> toRemove = new ArrayList<Integer>();
		for (int i=0;i<rs0.size();i++){
			Region r = rs0.get(i);
			if (CommonUtils.getPointsWithinWindow(highEnds, r).size()<2 
					&& CommonUtils.getPointsWithinWindow(lowEnds, r).size()<2)
				toRemove.add(i);
		}
		for (int i=toRemove.size()-1; i>=0; i--){
			summits.remove(toRemove.get(i).intValue());
			rs0.remove(toRemove.get(i).intValue());
		}
		toRemove = null;
		System.out.println("Remove short PET regions, "+rs0.size()+" regions.");

		for (int j=0;j<rs0.size();j++){	// for all nonTSS regions
			Point centerPoint = summits.get(j);
			Region tssRegion = null;	// pretend to be a TSS
			if (rs0.get(j).getWidth()<tss_radius*2)
				tssRegion = rs0.get(j);
			else
				tssRegion = centerPoint.expand(tss_radius);
			
			// get the distal ends, merge nearby read pairs
			ArrayList<Integer>idx = CommonUtils.getPointsWithinWindow(lowEnds, tssRegion);
			if (idx.size()>1){
				int exStart = centerPoint.getLocation()-self_exclude;
				int exEnd = centerPoint.getLocation()+self_exclude;
				ArrayList<ReadPair> rps = new ArrayList<ReadPair> ();
				for (int i: idx){
					ReadPair rp = low.get(i);
					int pos = rp.r2.getLocation();
					if(pos<exStart || pos>exEnd)
						rps.add(rp);
				}
				Collections.sort(rps, new Comparator<ReadPair>(){
		            public int compare(ReadPair o1, ReadPair o2) {
		                return o1.compareRead2(o2);
		            }
		        });
				ArrayList<ReadPairCluster> rpcs = new ArrayList<ReadPairCluster>();
				int current = -100000;
				ReadPairCluster c = new ReadPairCluster();
				for (ReadPair rp: rps){
					if (rp.r2.getLocation()-current>read_merge_dist){	// a big gap
						if (c.reads.size()>=2){
							rpcs.add(c);
						}
						c = new ReadPairCluster();
					}
					c.addReadPair(rp);
					current = rp.r2.getLocation();
				}
				if (c.reads.size()>=2){		// finish up the last cluster
					rpcs.add(c);
				}
				
				// remove those distal ends that overlap with tssPoints
				ArrayList<ReadPairCluster> toRemoveClusters = new ArrayList<ReadPairCluster>();
				for (ReadPairCluster cc: rpcs){
					Region distal = new Region(centerPoint.getGenome(), centerPoint.getChrom(), cc.r2min, cc.r2max);
					if (!CommonUtils.getPointsWithinWindow(tssPoints, distal).isEmpty())
						toRemoveClusters.add(cc);
				}
				rpcs.removeAll(toRemoveClusters);
				
				// test whether to merge clusters
				toRemoveClusters = new ArrayList<ReadPairCluster>();
				for (int i=1; i<rpcs.size();i++){
					ReadPairCluster c1=rpcs.get(i-1);
					ReadPairCluster c2=rpcs.get(i);
					// cluster_merge_dist is dependent on the distance between two anchor regions
					int dist = Math.min(c1.r2min+c1.r2max-c1.r1min-c1.r1max, c2.r2min+c2.r2max-c2.r1min-c2.r1max)/2;
					int cluster_merge_dist = Math.min(max_cluster_merge_dist, read_merge_dist + (int)Math.sqrt(dist) * distance_factor);
					double d = Math.min(c1.getDensity(cluster_merge_dist), c2.getDensity(cluster_merge_dist));
					if (c2.r2min - c1.r2max < cluster_merge_dist){
						if (simple_merge){
							// simply merge c1 to c2
							for (ReadPair rp2: c1.reads)
								c2.addReadPair(rp2);
							toRemoveClusters.add(c1);
							continue;
						}
						Region rMerged = new Region(centerPoint.getGenome(), centerPoint.getChrom(), c1.r2min, c2.r2max);
						idx = CommonUtils.getPointsWithinWindow(highEnds, rMerged);
						// TSS can expand at the lower end, but the high end is dependent on the distal read positions
						int tssEnd = Math.min(c1.r2min, c2.r2min) - (tssRegion.getMidpoint().getLocation()+self_exclude);
						tssEnd = Math.min(Math.max(tssEnd,0), cluster_merge_dist);
						Region tssExpanded = tssRegion.expand(cluster_merge_dist, tssEnd);	
						int start = tssExpanded.getStart();
						int end = tssExpanded.getEnd();
						rps = new ArrayList<ReadPair> ();
						for (int ii: idx){
							ReadPair rp2 = high.get(ii);
							int pos = rp2.r1.getLocation();
							if(pos>=start && pos<=end)
								rps.add(rp2);
						}
						Collections.sort(rps, new Comparator<ReadPair>(){
				            public int compare(ReadPair o1, ReadPair o2) {
				                return o1.compareRead1(o2);
				            }
				        });
						// coord1 are the read1 positions of those pairs that read2 is in the merged region
						ArrayList<Integer> coord1 = new ArrayList<Integer>();	
						for (ReadPair rp: rps)
							coord1.add(rp.r1.getLocation());
						int idxMin = Collections.binarySearch(coord1, Math.min(c1.r1min, c2.r1min));
						if (idxMin<0){
							System.out.println("c1.r1min, c2.r1min: " + c1.r1min + "," + c2.r1min);
							idxMin = Math.max(0, -idxMin-1);
						}
						int idxMax = Collections.binarySearch(coord1, Math.max(c1.r1max, c2.r1max));
						if (idxMax<0){
							System.out.println("c1.r1max, c2.r1max: " + c1.r1max + "," + c2.r1max);
							idxMax = Math.min(-idxMax-1, rps.size()-1);
						}
						ReadPairCluster cNew = new ReadPairCluster();
						for (int ii=idxMin; ii<=idxMax;ii++){
							cNew.addReadPair(rps.get(ii));
						}
						double best_density = cNew.getDensity(cluster_merge_dist);
						int best_idxMin = idxMin;
						int best_idxMax = idxMax;
						
						// expand to the lower end first
						for (int ii=idxMin-1; ii>=0; ii--){
							cNew.addReadPair(rps.get(ii));
							double new_density = cNew.getDensity(cluster_merge_dist);
							if (new_density>best_density){
								best_idxMin = ii;
								best_density = new_density;
							}
						}
						// update to the best set of readpairs after expanding the lower end
						cNew = new ReadPairCluster();		
						for (int ii=best_idxMin; ii<=best_idxMax;ii++){
							cNew.addReadPair(rps.get(ii));
						}
						
						// expand to the higher end
						for (int ii=idxMax+1; ii<rps.size(); ii++){
							cNew.addReadPair(rps.get(ii));
							double new_density = cNew.getDensity(cluster_merge_dist);
							if (new_density>best_density){
								best_idxMax = ii;
								best_density = new_density;
							}
						}
						if (best_density>d){	// better density, merge the regions
							// update to the best set of readpairs
							cNew = new ReadPairCluster();
							for (int ii=best_idxMin; ii<=best_idxMax;ii++){
								cNew.addReadPair(rps.get(ii));
							}
							rpcs.set(i, cNew);
							toRemoveClusters.add(c1);
						}
					}
				}	// for each pair of nearby clusters
				rpcs.removeAll(toRemoveClusters);
				toRemoveClusters.clear();
				toRemoveClusters=null;
				
				ArrayList<ReadPairCluster> rpcs2 = splitRecursively(rpcs, true);
				if (rpcs2!=null){
					rpcs = rpcs2;
					rpcs2 = null;
				}				
				for (ReadPairCluster cc: rpcs){
					Interaction it = new Interaction();
					interactions.add(it);
					it.geneSymbol = centerPoint.toString();
//					it.geneID = id;
					it.tss = centerPoint;
					it.tssRegion = new Region(centerPoint.getGenome(), centerPoint.getChrom(), cc.r1min, cc.r1max);
					it.distalRegion = new Region(centerPoint.getGenome(), centerPoint.getChrom(), cc.r2min, cc.r2max);
					it.distalPoint = it.distalRegion.getMidpoint();
					it.count = cc.reads.size();
					int cluster_merge_dist = Math.min(max_cluster_merge_dist, 
							read_merge_dist + (int)Math.sqrt(Math.abs(it.tss.offset(it.distalPoint))) * distance_factor);
					it.density = cc.getDensity(cluster_merge_dist);
					usedPETs.addAll(cc.reads);
				}
				rpcs = null;
			}  // nonTSS regions with long PETs
		} // loop over all nonTSS regions
		interactions.trimToSize();
		System.out.println("\n\nNon-Tss interactions ="+(interactions.size()-itTssCount)+
				", total interactions ="+interactions.size()+", "+CommonUtils.timeElapsed(tic0));

		
		/** Indirect counts */
		// for each TSS and its distal anchor, count how many read pairs lead to other distal anchors
		// it is like finding the other two edges of the triangle. The count is the sum (over all distal anchors) 
		// of min connecting read count (btw tss-distal2 and distal-distal2).
		System.out.println("\nCount indirect / triangle read pairs, "+CommonUtils.timeElapsed(tic0));
		HashMap<String, ArrayList<Interaction>> gene2it = new HashMap<String, ArrayList<Interaction>>();
		for (Interaction it: interactions){
			if (!gene2it.containsKey(it.geneSymbol))
				gene2it.put(it.geneSymbol, new ArrayList<Interaction>());
			gene2it.get(it.geneSymbol).add(it);
		}
		for (String gene: gene2it.keySet()){
			ArrayList<Interaction> its = gene2it.get(gene);
			for (Interaction it: its){
//				System.out.println(it);
				ArrayList<Integer> indirectCounts = new ArrayList<Integer>();
				Region distal = it.distalRegion.expand(chiapet_radius, chiapet_radius);
				for (Interaction it2: its){
					if (it==it2)
						continue;
					int start = it2.distalRegion.getStart()-chiapet_radius;
					int end = it2.distalRegion.getEnd()+chiapet_radius;
					if (distal.getStart() > start){		// if it1 has higher coord, select by high and then low
						ArrayList<Integer> idx = CommonUtils.getPointsWithinWindow(highEnds, distal);
						ArrayList<ReadPair> rps = new ArrayList<ReadPair> ();
						for (int ii: idx){
							ReadPair rp2 = high.get(ii);
							int pos = rp2.r1.getLocation();
							if(pos>=start && pos<=end)
								rps.add(rp2);
						}
						indirectCounts.add(Math.min(rps.size(), it2.count));		// distal2--distal , it2_count
					}
					else{		// select by low and then high
						ArrayList<Integer> idx = CommonUtils.getPointsWithinWindow(lowEnds, distal);
						ArrayList<ReadPair> rps = new ArrayList<ReadPair> ();
						for (int ii: idx){
							ReadPair rp2 = low.get(ii);
							int pos = rp2.r2.getLocation();
							if(pos>=start && pos<=end)
								rps.add(rp2);
						}
						indirectCounts.add(Math.min(rps.size(), it2.count));		// distal--distal2 , it2_count
					}
				}
				int sum = 0;
				for (int c: indirectCounts)
					sum += c;
				it.indirectCount = sum;
				
//				System.out.println();
//				for (int c: indirectCounts)
//					System.out.print(c);
//				System.out.println();
			}
		}
		

		// mark PET1 (after removing the PET2+)
		high.removeAll(usedPETs);
		high.trimToSize();
		low.clear();
		low = null;
		
	
		/** Annotate and report */
		System.out.println("\n\nAnnotate and report, "+CommonUtils.timeElapsed(tic0));
		// report the interactions and annotations
		// annotate the proximal and distal anchors with TF and HM and regions
		StringBuilder sb = new StringBuilder();
		for (Interaction it: interactions){
			ArrayList<Integer> isOverlapped = new ArrayList<Integer>();
			Point mid = it.distalRegion.getMidpoint();
			int radius = it.distalRegion.getWidth()/2+chiapet_radius;
			for (ArrayList<Point> ps : allPeaks){
				ArrayList<Point> p = CommonUtils.getPointsWithinWindow(ps, mid, radius);
				isOverlapped.add(p.size());
				if (!it.isTfAnchord && !p.isEmpty()){	// if the distal point has not been anchored by a TF site
					it.distalPoint = p.get(p.size()/2);
					it.isTfAnchord = true;
				}
			}
			for (List<Region> rs: allRegions){
				isOverlapped.add(CommonUtils.getRegionsOverlapsWindow(rs, it.distalRegion, chiapet_radius).size());
			}
			// proximal
			for (ArrayList<Point> ps : allPeaks){
				ArrayList<Point> p = CommonUtils.getPointsWithinWindow(ps, it.tss, chiapet_radius);
				isOverlapped.add(p.size());
			}
			Region tssRegion = it.tss.expand(chiapet_radius);
			for (List<Region> rs: allRegions){
				isOverlapped.add(CommonUtils.getRegionsOverlapsWindow(rs, tssRegion, chiapet_radius).size());
			}
			// print out TF and region overlaps
			sb.append(it.toString()).append("\t");
			for (int b: isOverlapped)
				sb.append(b).append("\t");
			
			// print ChIA-PET call overlap info
			Point tssPoint = it.tss;
			int tssHalfWidth = it.tssRegion.getWidth()/2+chiapet_radius;
			int distalHalfWidth = it.distalRegion.getWidth()/2+chiapet_radius;
			Point distalPoint = it.distalPoint;
			Point tssLeft = new Point(genome, tssPoint.getChrom(), tssPoint.getLocation()-2000);
			Point tssRight = new Point(genome, tssPoint.getChrom(), tssPoint.getLocation()+2000);
			
			// GERM
			int index = Collections.binarySearch(germTss,  tssLeft);
			if( index < 0 )  							// if key not found
				index = -(index+1); 
			int indexRight = Collections.binarySearch(germTss,  tssRight);
			if( indexRight < 0 )  							// if key not found
				indexRight = -(indexRight+1); 
			// if key match found, continue to search ( binarySearch() give undefined index with multiple matches)
			boolean isGermOverlapped = false;
			indexRange: for (int i=index-1;i<=indexRight+2;i++){
				if (i<0 || i>=germTss.size())
					continue;
				try{
					Point tt = germTss.get(i);
					if (tt.distance(tssPoint)<=tssHalfWidth){ 
						if (!germTss2distals.containsKey(tt))
							continue;
						for (Point d:germTss2distals.get(tt)){
							if (d.distance(distalPoint)<=distalHalfWidth){
								isGermOverlapped = true;
								break indexRange;
							}
						}
					}
				}
				catch (IllegalArgumentException e){	// ignore								
				}								
			}
			sb.append(isGermOverlapped?"1\t":"0\t");
			
			// Mango
			// aPoints and bPoints are the midPoint of the two anchorRegions
			index = Collections.binarySearch(aPoints,  tssLeft);
			if( index < 0 )  							// if key not found
				index = -(index+1); 
			indexRight = Collections.binarySearch(aPoints,  tssRight);
			if( indexRight < 0 )  							// if key not found
				indexRight = -(indexRight+1); 
			// if key match found, continue to search ( binarySearch() give undefined index with multiple matches)
			boolean isMangoOverlapped = false;
			indexA: for (int i=index-1;i<=indexRight+2;i++){
				if (i<0 || i>=aPoints.size())
					continue;
				try{
					Point a = aPoints.get(i);
					if (a.distance(tssPoint)<=tssHalfWidth){ 
						if (!a2bs.containsKey(a))
							continue;
						for (Point b:a2bs.get(a)){
							if (b.distance(distalPoint)<=distalHalfWidth){
								isMangoOverlapped = true;
								break indexA;
							}
						}
					}
				}
				catch (IllegalArgumentException e){	// ignore								
				}
			}
			if (isMangoOverlapped)
				sb.append("1\t");
			else{
				index = Collections.binarySearch(bPoints,  tssLeft);
				if( index < 0 )  							// if key not found
					index = -(index+1); 
				indexRight = Collections.binarySearch(bPoints,  tssRight);
				if( indexRight < 0 )  							// if key not found
					indexRight = -(indexRight+1); 
				// if key match found, continue to search ( binarySearch() give undefined index with multiple matches)
				isMangoOverlapped = false;
				indexB: for (int i=index-1;i<=indexRight+2;i++){
					if (i<0 || i>=bPoints.size())
						continue;
					try{
						Point b = bPoints.get(i);
						if (b.distance(tssPoint)<=tssHalfWidth){ 
							if (!b2as.containsKey(b))
								continue;
							for (Point a:b2as.get(b)){
								if (a.distance(distalPoint)<=distalHalfWidth){
									isMangoOverlapped = true;
									break indexB;
								}
							}
						}
					}
					catch (IllegalArgumentException e){	// ignore								
					}								
				}
				sb.append(isMangoOverlapped?"1\t":"0\t");
			}
			
			CommonUtils.replaceEnd(sb, '\n');
			
		}
		CommonUtils.writeFile(Args.parseString(args, "out", "Result")+".readClusters.txt", sb.toString());
		
		// output BEDPE format
		// HERE we need to also include PET1 for MICC and ChiaSig analysis
		System.out.println("Single PETs n=" + high.size());
		for (ReadPair rp: high){
			Interaction it = new Interaction();
			interactions.add(it);
//			it.geneSymbol = centerPoint.toString();
//			it.geneID = id;
			it.tss = rp.r1;
			it.distalPoint = rp.r2;
			it.tssRegion = new Region(rp.r1.getGenome(), rp.r1.getChrom(), rp.r1.getLocation(), rp.r1.getLocation());
			it.distalRegion = new Region(rp.r2.getGenome(), rp.r2.getChrom(), rp.r2.getLocation(), rp.r2.getLocation());
			it.count = 1;
//			int cluster_merge_dist = Math.min(max_cluster_merge_dist, 
//					read_merge_dist + (int)Math.sqrt(Math.abs(it.tss.offset(it.distalPoint))) * distance_factor);
//			it.density = cc.getDensity(cluster_merge_dist);
		}
		high.clear();
		high=null;
		
		sb = new StringBuilder();
		for (Interaction it: interactions){
			Region distalLocal = it.distalRegion.expand(read_merge_dist, read_merge_dist);
			Region tssLocal = it.tssRegion.expand(read_merge_dist, read_merge_dist);
			int distalLocalCount, tssLocalCount;
			if (it.distalPoint.offset(it.tss)<0){	// distal is in lower coord
				distalLocalCount = CommonUtils.getPointsWithinWindow(lowEnds, distalLocal).size();
				tssLocalCount = CommonUtils.getPointsWithinWindow(highEnds, tssLocal).size();
				sb.append(String.format("%s\t%s\t%d\t%d\t%d\n", distalLocal.toBED(), tssLocal.toBED(),
						it.count, distalLocalCount, tssLocalCount));
			}
			else{	// distal is in higher coord
				distalLocalCount = CommonUtils.getPointsWithinWindow(highEnds, distalLocal).size();
				tssLocalCount = CommonUtils.getPointsWithinWindow(lowEnds, tssLocal).size();
				sb.append(String.format("%s\t%s\t%d\t%d\t%d\n", tssLocal.toBED(), distalLocal.toBED(),
						it.count, tssLocalCount, distalLocalCount));
			}
		}
		CommonUtils.writeFile(Args.parseString(args, "out", "Result")+".bedpe", sb.toString());
		
		System.out.println("\n\nDone: "+CommonUtils.timeElapsed(tic0));
	}
	
	/** split read pair cluster recursively <br>
	 *  at gaps larger than cluster_merge_dist, on both ends alternatively 
	 *  because splitting at one end may remove some PETs that introduce gaps at the other end */
	ArrayList<ReadPairCluster> splitRecursively(ArrayList<ReadPairCluster> rpcs, boolean toSplitLeftAnchor){
		if (rpcs.isEmpty())
			return null;
		
		int countSplit=0;
		ArrayList<ReadPairCluster> rpcs2 = new ArrayList<ReadPairCluster>();
		for (ReadPairCluster cc: rpcs){
			HashMap<Point, ArrayList<ReadPair>> map = new HashMap<Point, ArrayList<ReadPair>>();
			ArrayList<Point> splitPoints = new ArrayList<Point>();
			HashSet<Point> tmp = new HashSet<Point>();
			int dist = cc.r2max - cc.r1min;
			int cluster_merge_dist = Math.min(max_cluster_merge_dist, read_merge_dist + (int)Math.sqrt(dist) * distance_factor);
			for (ReadPair rp: cc.reads){
				Point t = toSplitLeftAnchor ? rp.r1 : rp.r2;
				tmp.add(t);
				if (!map.containsKey(t))
					map.put(t, new ArrayList<ReadPair>());
				map.get(t).add(rp);
			}
			splitPoints.addAll(tmp);
			Collections.sort(splitPoints);
			int curr = -100000;
			ReadPairCluster c = new ReadPairCluster();
			countSplit--;		// first split is not real, subtract count here
			for (Point p: splitPoints){
				if (p.getLocation()-curr>cluster_merge_dist){	// a big gap
					countSplit++;
					if (c.reads.size()>=2)
						rpcs2.add(c);
					c = new ReadPairCluster();
				}
				for (ReadPair rp: map.get(p))
					c.addReadPair(rp);
				curr = p.getLocation();
			}
			if (c.reads.size()>=2)		// finish up the last cluster
				rpcs2.add(c);
		}
		if (countSplit>0){
			ArrayList<ReadPairCluster> rpcs3 = splitRecursively(rpcs2, !toSplitLeftAnchor);	// split at the other end
			return rpcs3==null ? rpcs2 : rpcs3;
		}else
			return null;
	}

	private void clusterPETs(){		// NOT USED FOR NOW
		long tic0 = System.currentTimeMillis();
		long tic = System.currentTimeMillis();
		int read_merge_dist = Args.parseInteger(args, "read_merge_dist", 500);
		int max_cluster_merge_dist = Args.parseInteger(args, "max_cluster_merge_dist", 3000);
		int distance_factor = Args.parseInteger(args, "distance_factor", 3);
		int self_exclude = Args.parseInteger(args, "self_exclude", 8000);	// cutoff for self-ligation PET
		int tss_radius = Args.parseInteger(args, "tss_radius", 2000);
		int chiapet_radius = Args.parseInteger(args, "chiapet_radius", 2000);
		double overlap_ratio = Args.parseDouble(args, "overlap_ratio", 0.8);

		
		// load refSeq gene annotation
		ArrayList<String> lines = CommonUtils.readTextFile(Args.parseString(args, "gene_anno", null));
		ArrayList<TSS> allTSS = new ArrayList<TSS>();
		TreeMap<String, TreeSet<StrandedPoint>> gene2tss = new TreeMap<String, TreeSet<StrandedPoint>>();
		for (int i=0;i<lines.size();i++){
			String t = lines.get(i);
			if (t.startsWith("#"))
				continue;
			String f[] = t.split("\t");
			String chr = f[2].replace("chr", "");
			char strand = f[3].charAt(0);
			TSS tss = new TSS(f[12], new StrandedPoint(genome, chr, Integer.parseInt(f[strand=='+'?4:5]), strand),i);
			allTSS.add(tss);
		}
		
		
		// load read pairs
		// only use read pairs on the same chromosome
		// the left read is required to be lower than the right read, if not, flip the two
		
		// sort by each end so that we can use binary search to find matches or overlaps
		System.out.println("Loading ChIA-PET read pairs: "+CommonUtils.timeElapsed(tic));

		ArrayList<String> read_pairs = CommonUtils.readTextFile(Args.parseString(args, "read_pair", null));
		ArrayList<ReadPair> low = new ArrayList<ReadPair> ();	// all read pairs sorted by the low end
		ArrayList<ReadPair> high = new ArrayList<ReadPair> ();	// sorted by high end
		for (String s: read_pairs){
			String[] f = s.split("\t");
			Point r1 = Point.fromString(genome, f[0]);
			String r1Chrom=r1.getChrom();
			Point r2 = Point.fromString(genome, f[1]);
			// TODO: change next line if prediction cross-chrom interactions
			if (!r1Chrom.equals(r2.getChrom()))		// r1 and r2 should be on the same chromosome
				continue;
			if (r1.distance(r2)>self_exclude)	// PET shorter than 8kb are considered self-ligation reads
				continue;
			ReadPair rp = new ReadPair();
			if (r1.compareTo(r2)<0){		// r1 should be lower than r2
				rp.r1 = r1;
				rp.r2 = r2;
			}
			else{
				rp.r1 = r2;
				rp.r2 = r1;
			}
			low.add(rp);
			ReadPair rp2 = new ReadPair();
			rp2.r1 = rp.r1;
			rp2.r2 = rp.r2;
			high.add(rp2);
		}
		low.trimToSize();
		high.trimToSize();
		Collections.sort(low, new Comparator<ReadPair>(){
            public int compare(ReadPair o1, ReadPair o2) {
                return o1.compareRead1(o2);
            }
        });
		Collections.sort(high, new Comparator<ReadPair>(){
            public int compare(ReadPair o1, ReadPair o2) {
                return o1.compareRead2(o2);
            }
        });
		ArrayList<Point> lowEnds = new ArrayList<Point>();
		for (ReadPair r:low)
			lowEnds.add(r.r1);
		lowEnds.trimToSize();
		ArrayList<Point> highEnds = new ArrayList<Point>();
		for (ReadPair r:high)
			highEnds.add(r.r2);
		highEnds.trimToSize();
		
		System.out.println("Loaded "+ highEnds.size() +" ChIA-PET read pairs: "+CommonUtils.timeElapsed(tic));
		System.out.println();
		
		
		// load TF sites
		ArrayList<String> tfs = CommonUtils.readTextFile(Args.parseString(args, "tf_sites", null));
		ArrayList<ArrayList<Point>> allPeaks = new ArrayList<ArrayList<Point>>();
		for (int i=0;i<tfs.size();i++){
			try{
				ArrayList<Point> ps = new ArrayList<Point>();
				ps.addAll(GPSParser.parseGPSOutput(tfs.get(i), genome));
				ps.trimToSize();
				Collections.sort(ps);
				allPeaks.add(ps);
				System.out.println("Loaded "+tfs.get(i));
			}
			catch (IOException e){
				System.out.println(tfs.get(i)+" does not have a valid GPS/GEM event call file.");
				e.printStackTrace(System.err);
				System.exit(1);
			}
		}
		allPeaks.trimToSize();

		// load histone mark or DHS, SE regions
		ArrayList<String> hms = CommonUtils.readTextFile(Args.parseString(args, "regions", null));
		ArrayList<List<Region>> allRegions = new ArrayList<List<Region>>();
		for (int i=0;i<hms.size();i++){
			allRegions.add(CommonUtils.load_BED_regions(genome, hms.get(i)).car());
			System.out.println("Loaded "+hms.get(i));
		}
		allRegions.trimToSize();
		System.out.println();
	
		// load other Interaction calls
		lines = CommonUtils.readTextFile(Args.parseString(args, "germ", null));
		HashMap<Point,ArrayList<Point>> germTss2distals = new HashMap<Point,ArrayList<Point>>();
		ArrayList<Point> germTss = new ArrayList<Point>();
		for (String l: lines){		// each line is a call
			String f[] = l.split("\t");
			Point t = new Region(genome, f[3].replace("chr", ""), Integer.parseInt(f[4]), Integer.parseInt(f[5])).getMidpoint();
			if (!germTss2distals.containsKey(t))
				germTss2distals.put(t, new ArrayList<Point>());
			germTss2distals.get(t).add(new Region(genome, f[0].replace("chr", ""), Integer.parseInt(f[1]), Integer.parseInt(f[2])).getMidpoint());
		}
		germTss.addAll(germTss2distals.keySet());
		germTss.trimToSize();
		Collections.sort(germTss);
		
		lines = CommonUtils.readTextFile(Args.parseString(args, "mango", null));
		HashMap<Point, ArrayList<Point>> a2bs = new HashMap<Point, ArrayList<Point>>();
		HashMap<Point, ArrayList<Point>> b2as = new HashMap<Point, ArrayList<Point>>();
		for (String l: lines){		// each line is a call
			String f[] = l.split("\t");
			Point a = new Region(genome, f[0].replace("chr", ""), Integer.parseInt(f[1]), Integer.parseInt(f[2])).getMidpoint();
			Point b = new Region(genome, f[3].replace("chr", ""), Integer.parseInt(f[4]), Integer.parseInt(f[5])).getMidpoint();
			if (!a2bs.containsKey(a))
				a2bs.put(a, new ArrayList<Point>());
			a2bs.get(a).add(b);
			if (!b2as.containsKey(b))
				b2as.put(b, new ArrayList<Point>());
			b2as.get(b).add(a);
		}
		ArrayList<Point> aPoints = new ArrayList<Point>();
		aPoints.addAll(a2bs.keySet());
		aPoints.trimToSize();
		Collections.sort(aPoints);
		ArrayList<Point> bPoints = new ArrayList<Point>();
		bPoints.addAll(b2as.keySet());
		bPoints.trimToSize();
		Collections.sort(bPoints);
		System.out.println("Loaded all the data, "+CommonUtils.timeElapsed(tic0));
		
		// find dense read cluster for each gene
		// Two big chunks of codes for 1) Distal--TSS or 2) TSS--Distal, according to their coordinates
		ArrayList<String> geneList = new ArrayList<String>();
		geneList.addAll(gene2tss.keySet());
		geneList.trimToSize();
		ArrayList<Interaction> interactions = new ArrayList<Interaction>();
		tic = System.currentTimeMillis();
		
		for (int id=0;id<geneList.size();id++){
			String g = geneList.get(id);
			System.out.print(g+" ");
			TreeSet<StrandedPoint> coords = gene2tss.get(g);
			// if a gene has multiple TSSs, use the center position
			int count = coords.size();
			Point centerPoint = null;
			for (StrandedPoint p:coords){
				if (count<coords.size()/2)
					break;
				else{
					centerPoint = p;
					count--;
				}
			}
			
		/** 1) For TSS at higher coordinates, distal anchors are at lower coordinates */
			
			// get the distal ends, merge nearby read pairs
			Region tssRegion = centerPoint.expand(tss_radius);
			Region excludeRegion = centerPoint.expand(self_exclude);
			int exStart = excludeRegion.getStart();
			int exEnd = excludeRegion.getEnd();
			ArrayList<Integer> idx = CommonUtils.getPointsWithinWindow(highEnds, tssRegion);
			if (!idx.isEmpty()){
//				System.out.println("\n"+g+"\tCluster reads, "+CommonUtils.timeElapsed(tic));
//				tic = System.currentTimeMillis();
				
				ArrayList<ReadPair> rps = new ArrayList<ReadPair> ();
				for (int i: idx){
					ReadPair rp = high.get(i);
					int pos = rp.r1.getLocation();
					if(pos<exStart || pos>exEnd)
						rps.add(rp);
				}
				Collections.sort(rps, new Comparator<ReadPair>(){
		            public int compare(ReadPair o1, ReadPair o2) {
		                return o1.compareRead1(o2);
		            }
		        });
				ArrayList<ReadPairCluster> rpcs = new ArrayList<ReadPairCluster>();
				int current = -100000;
				ReadPairCluster c = new ReadPairCluster();
				// merge reads
				for (ReadPair rp: rps){
					if (rp.r1.getLocation()-current>read_merge_dist){	// a big gap
						if (c.reads.size()>=2){
							rpcs.add(c);
						}
						c = new ReadPairCluster();
					}
					c.addReadPair(rp);
					current = rp.r1.getLocation();
				}
				if (c.reads.size()>=2){		// finish up the last cluster
					rpcs.add(c);
				}
				
				// test whether to merge clusters
				// if two nearby clusters are within the cluster_merge_dist, make a new cluster that covers both distal regions
				// also try to expand the TSS regions a little to see if we can merge more nearby reads
	
//				System.out.println("Merge clusters, "+CommonUtils.timeElapsed(tic));
//				tic = System.currentTimeMillis();
				
	//			System.out.println("\nDistal region at lower coord than TSS\nBefore merging, number of clusters = "+rpcs.size());
				ArrayList<ReadPairCluster> toRemoveClusters = new ArrayList<ReadPairCluster>();
				for (int i=1; i<rpcs.size();i++){
					ReadPairCluster c1=rpcs.get(i-1);
					ReadPairCluster c2=rpcs.get(i);
					int dist = Math.min(c1.r2min+c1.r2max-c1.r1min-c1.r1max, c2.r2min+c2.r2max-c2.r1min-c2.r1max)/2;
					int cluster_merge_dist = Math.min(max_cluster_merge_dist, read_merge_dist + (int)Math.sqrt(dist) * distance_factor);
					double d = Math.min(c1.getDensity(cluster_merge_dist), c2.getDensity(cluster_merge_dist));
					if (c2.r1min - c1.r1max < cluster_merge_dist){
						Region distalRegionMerged = new Region(centerPoint.getGenome(), centerPoint.getChrom(), c1.r1min, c2.r1max);
						idx = CommonUtils.getPointsWithinWindow(lowEnds, distalRegionMerged);
						// TSS can expand at the high end, but the low end is dependent on the distal read positions
						int tssStart = tssRegion.getMidpoint().getLocation()-self_exclude - Math.max(c1.r1max, c2.r1max);
						tssStart = Math.min(Math.max(tssStart,0), cluster_merge_dist);
						Region tssExpanded = tssRegion.expand(tssStart, cluster_merge_dist);
						int start = tssExpanded.getStart();
						int end = tssExpanded.getEnd();
						rps = new ArrayList<ReadPair> ();
						for (int ii: idx){
							ReadPair rp2 = low.get(ii);
							int pos = rp2.r2.getLocation();
							if(pos>=start && pos<=end)
								rps.add(rp2);
						}
						Collections.sort(rps, new Comparator<ReadPair>(){
				            public int compare(ReadPair o1, ReadPair o2) {
				                return o1.compareRead2(o2);
				            }
				        });
						// coord2 are the read2 positions of those pairs that read1 is in the merged region
						ArrayList<Integer> coord2 = new ArrayList<Integer>();	
						for (ReadPair rp: rps)
							coord2.add(rp.r2.getLocation());
						int idxMin = Collections.binarySearch(coord2, Math.min(c1.r2min, c2.r2min));
						if (idxMin<0){
							System.out.println("c1.r2min, c2.r2min: " + c1.r2min + "," + c2.r2min);
							idxMin = Math.max(0, -idxMin-1);
						}
						int idxMax = Collections.binarySearch(coord2, Math.max(c1.r2max, c2.r2max));
						if (idxMax<0){
							System.out.println("c1.r2max, c2.r2max: " + c1.r2max + "," + c2.r2max);
							idxMax = Math.min(-idxMax-1, rps.size()-1);
						}
						ReadPairCluster cNew = new ReadPairCluster();
						for (int ii=idxMin; ii<=idxMax;ii++){
							cNew.addReadPair(rps.get(ii));
						}
						double best_density = cNew.getDensity(cluster_merge_dist);
						int best_idxMin = idxMin;
						int best_idxMax = idxMax;
						// expand to the lower end first
						for (int ii=idxMin-1; ii>=0; ii--){
							cNew.addReadPair(rps.get(ii));
							double new_density = cNew.getDensity(cluster_merge_dist);
	//						System.out.println("New density "+new_density);
							if (new_density>best_density){
	//							System.out.println("Better density "+best_density+"-->"+new_density);
								best_idxMin = ii;
								best_density = new_density;
							}
						}
						// update to the best set of readpairs after expanding the lower end
						cNew = new ReadPairCluster();		
						for (int ii=best_idxMin; ii<=best_idxMax;ii++){
							cNew.addReadPair(rps.get(ii));
						}
						
						// expand to the higher end
	//					System.out.println("Best density "+best_density);
						for (int ii=idxMax+1; ii<rps.size(); ii++){
							cNew.addReadPair(rps.get(ii));
							double new_density = cNew.getDensity(cluster_merge_dist);
	//						System.out.println("New density "+new_density);
							if (new_density>best_density){
	//							System.out.println("Better density "+best_density+"-->"+new_density);
								best_idxMax = ii;
								best_density = new_density;
							}
						}
						if (best_density>d){	// better density, merge the regions
	//						System.out.println("Merged "+c1.toString()+" and \n\t"+c2.toString()+" to "+cNew.toString());
							// update to the best set of readpairs
							cNew = new ReadPairCluster();
							for (int ii=best_idxMin; ii<=best_idxMax;ii++){
								cNew.addReadPair(rps.get(ii));
							}
							rpcs.set(i, cNew);
							toRemoveClusters.add(c1);
						}
					}
				}	// for each pair of nearby clusters
				rpcs.removeAll(toRemoveClusters);
	//			System.out.println("After merging,  number of clusters = "+rpcs.size());
	
				// report
				for (ReadPairCluster cc: rpcs){
					Interaction it = new Interaction();
					interactions.add(it);
					it.geneSymbol = g;
					it.geneID = id;
					it.tss = centerPoint;
					it.tssRegion = new Region(centerPoint.getGenome(), centerPoint.getChrom(), cc.r2min, cc.r2max);
					it.distalRegion = new Region(centerPoint.getGenome(), centerPoint.getChrom(), cc.r1min, cc.r1max);
					it.distalPoint = it.distalRegion.getMidpoint();
					it.count = cc.reads.size();
					int cluster_merge_dist = Math.min(max_cluster_merge_dist, 
							read_merge_dist + (int)Math.sqrt(Math.abs(it.tss.offset(it.distalPoint))) * distance_factor);
					it.density = cc.getDensity(cluster_merge_dist);
				}
				rpcs = null;
//				System.gc();

//				System.out.println("Done merging clusters, "+CommonUtils.timeElapsed(tic));
//				tic = System.currentTimeMillis();
			}
			
		/** 2) For TSS at lower coordinates */
			
			// get the distal ends, merge nearby read pairs
			idx = CommonUtils.getPointsWithinWindow(lowEnds, tssRegion);
			if (!idx.isEmpty()){
				ArrayList<ReadPair> rps = new ArrayList<ReadPair> ();
				for (int i: idx){
					ReadPair rp = low.get(i);
					int pos = rp.r2.getLocation();
					if(pos<exStart || pos>exEnd)
						rps.add(rp);
				}
				Collections.sort(rps, new Comparator<ReadPair>(){
		            public int compare(ReadPair o1, ReadPair o2) {
		                return o1.compareRead2(o2);
		            }
		        });
				ArrayList<ReadPairCluster> rpcs = new ArrayList<ReadPairCluster>();
				int current = -100000;
				ReadPairCluster c = new ReadPairCluster();
				for (ReadPair rp: rps){
					if (rp.r2.getLocation()-current>read_merge_dist){	// a big gap
						if (c.reads.size()>=2){
							rpcs.add(c);
						}
						c = new ReadPairCluster();
					}
					c.addReadPair(rp);
					current = rp.r2.getLocation();
				}
				if (c.reads.size()>=2){		// finish up the last cluster
					rpcs.add(c);
				}
				
				// test whether to merge clusters
//				System.out.println("Start merge clusters, "+CommonUtils.timeElapsed(tic));
//				tic = System.currentTimeMillis();
				
				ArrayList<ReadPairCluster> toRemoveClusters = new ArrayList<ReadPairCluster>();
				for (int i=1; i<rpcs.size();i++){
					ReadPairCluster c1=rpcs.get(i-1);
					ReadPairCluster c2=rpcs.get(i);
					// cluster_merge_dist is dependent on the distance between two anchor regions
					int dist = Math.min(c1.r2min+c1.r2max-c1.r1min-c1.r1max, c2.r2min+c2.r2max-c2.r1min-c2.r1max)/2;
					int cluster_merge_dist = Math.min(max_cluster_merge_dist, read_merge_dist + (int)Math.sqrt(dist) * distance_factor);
					double d = Math.min(c1.getDensity(cluster_merge_dist), c2.getDensity(cluster_merge_dist));
					if (c2.r2min - c1.r2max < cluster_merge_dist){
//						System.out.println("Close enough, "+CommonUtils.timeElapsed(tic));
//						tic = System.currentTimeMillis();
						Region rMerged = new Region(centerPoint.getGenome(), centerPoint.getChrom(), c1.r2min, c2.r2max);
						idx = CommonUtils.getPointsWithinWindow(highEnds, rMerged);
//						System.out.println("Got points, "+CommonUtils.timeElapsed(tic));
//						tic = System.currentTimeMillis();
						// TSS can expand at the lower end, but the high end is dependent on the distal read positions
						int tssEnd = Math.min(c1.r2min, c2.r2min) - (tssRegion.getMidpoint().getLocation()+self_exclude);
						tssEnd = Math.min(Math.max(tssEnd,0), cluster_merge_dist);
						Region tssExpanded = tssRegion.expand(cluster_merge_dist, tssEnd);	
						int start = tssExpanded.getStart();
						int end = tssExpanded.getEnd();
						rps = new ArrayList<ReadPair> ();
						for (int ii: idx){
							ReadPair rp2 = high.get(ii);
							int pos = rp2.r1.getLocation();
							if(pos>=start && pos<=end)
								rps.add(rp2);
						}
						Collections.sort(rps, new Comparator<ReadPair>(){
				            public int compare(ReadPair o1, ReadPair o2) {
				                return o1.compareRead1(o2);
				            }
				        });
//						System.out.println("Sorted read1, "+CommonUtils.timeElapsed(tic));
//						tic = System.currentTimeMillis();
						// coord1 are the read1 positions of those pairs that read2 is in the merged region
						ArrayList<Integer> coord1 = new ArrayList<Integer>();	
						for (ReadPair rp: rps)
							coord1.add(rp.r1.getLocation());
						int idxMin = Collections.binarySearch(coord1, Math.min(c1.r1min, c2.r1min));
						if (idxMin<0){
							System.out.println("c1.r1min, c2.r1min: " + c1.r1min + "," + c2.r1min);
							idxMin = Math.max(0, -idxMin-1);
						}
						int idxMax = Collections.binarySearch(coord1, Math.max(c1.r1max, c2.r1max));
						if (idxMax<0){
							System.out.println("c1.r1max, c2.r1max: " + c1.r1max + "," + c2.r1max);
							idxMax = Math.min(-idxMax-1, rps.size()-1);
						}
						ReadPairCluster cNew = new ReadPairCluster();
						for (int ii=idxMin; ii<=idxMax;ii++){
							cNew.addReadPair(rps.get(ii));
						}
						double best_density = cNew.getDensity(cluster_merge_dist);
						int best_idxMin = idxMin;
						int best_idxMax = idxMax;

//						System.out.println("To expand lower TSS, "+CommonUtils.timeElapsed(tic));
//						tic = System.currentTimeMillis();
						
						// expand to the lower end first
						for (int ii=idxMin-1; ii>=0; ii--){
							cNew.addReadPair(rps.get(ii));
							double new_density = cNew.getDensity(cluster_merge_dist);
	//						System.out.println("New density "+new_density);
							if (new_density>best_density){
	//							System.out.println("Better density "+best_density+"-->"+new_density);
								best_idxMin = ii;
								best_density = new_density;
							}
						}
						// update to the best set of readpairs after expanding the lower end
						cNew = new ReadPairCluster();		
						for (int ii=best_idxMin; ii<=best_idxMax;ii++){
							cNew.addReadPair(rps.get(ii));
						}
						
						// expand to the higher end
//						System.out.println("To expand higher TSS, "+CommonUtils.timeElapsed(tic));
//						tic = System.currentTimeMillis();
	//					System.out.println("Best density "+best_density);
						for (int ii=idxMax+1; ii<rps.size(); ii++){
							cNew.addReadPair(rps.get(ii));
							double new_density = cNew.getDensity(cluster_merge_dist);
	//						System.out.println("New density "+new_density);
							if (new_density>best_density){
	//							System.out.println("Better density "+best_density+"-->"+new_density);
								best_idxMax = ii;
								best_density = new_density;
							}
						}
						if (best_density>d){	// better density, merge the regions
	//						System.out.println("Merged "+c1.toString()+" and \n\t"+c2.toString()+" to "+cNew.toString());
							// update to the best set of readpairs
							cNew = new ReadPairCluster();
							for (int ii=best_idxMin; ii<=best_idxMax;ii++){
								cNew.addReadPair(rps.get(ii));
							}
							rpcs.set(i, cNew);
							toRemoveClusters.add(c1);
						}
					}
				}	// for each pair of nearby clusters
				rpcs.removeAll(toRemoveClusters);
	//			System.out.println("After merging,  number of clusters = "+rpcs.size());

//				System.out.println("Got clusters, "+CommonUtils.timeElapsed(tic));
//				tic = System.currentTimeMillis();
				for (ReadPairCluster cc: rpcs){
					Interaction it = new Interaction();
					interactions.add(it);
					it.geneSymbol = g;
					it.geneID = id;
					it.tss = centerPoint;
					it.tssRegion = new Region(centerPoint.getGenome(), centerPoint.getChrom(), cc.r1min, cc.r1max);
					it.distalRegion = new Region(centerPoint.getGenome(), centerPoint.getChrom(), cc.r2min, cc.r2max);
					it.distalPoint = it.distalRegion.getMidpoint();
					it.count = cc.reads.size();
					int cluster_merge_dist = Math.min(max_cluster_merge_dist, 
							read_merge_dist + (int)Math.sqrt(Math.abs(it.tss.offset(it.distalPoint))) * distance_factor);
					it.density = cc.getDensity(cluster_merge_dist);
				}
				interactions.trimToSize();
				rpcs = null;
//				System.gc();
//				System.out.println(String.format("Done merging clusters, n=%d, %s, total %s", 
//						interactions.size(), CommonUtils.timeElapsed(tic), CommonUtils.timeElapsed(tic0)));
//				tic = System.currentTimeMillis();
			}  
		}// for each gene

		// consolidate interaction anchors (for nearby TSSs)
		System.out.println("\n\nConsolidate nearby Tss interactions, "+CommonUtils.timeElapsed(tic0));
		TreeMap<Point, ArrayList<Interaction>> tss2it = new TreeMap<Point, ArrayList<Interaction>>();
		for (Interaction it: interactions){
			if (!tss2it.containsKey(it.tss))
				tss2it.put(it.tss, new ArrayList<Interaction>());
			tss2it.get(it.tss).add(it);
		}
		ArrayList<Point> tssList = new ArrayList<Point>();
		tssList.addAll(tss2it.keySet());
		for (int i=1;i<tssList.size();i++){
			Point p1 = tssList.get(i-1);
			Point p2 = tssList.get(i);
			if (p2.getChrom().equals(p1.getChrom()) && p2.offset(p1)<max_cluster_merge_dist){
				ArrayList<Interaction> its1 = tss2it.get(p1);
				ArrayList<Interaction> its2 = tss2it.get(p2);
				for (Interaction it1:its1){
					int start1 = it1.distalRegion.getStart();
					int end1 = it1.distalRegion.getEnd();
					for (Interaction it2:its2){
						if (it1.distalRegion.equals(it2.distalRegion))
							continue;
						int start2 = it2.distalRegion.getStart();
						int end2 = it2.distalRegion.getEnd();
//						int overlapWidth = Math.min(end1, end2)-Math.max(start2, start1);
						int overlapWidth = it2.distalRegion.getOverlapSize(it1.distalRegion);
						// if the two distal regions are highly overlapped
						if (overlapWidth>it1.distalRegion.getWidth()*overlap_ratio && overlapWidth>it2.distalRegion.getWidth()*overlap_ratio){
							//TODO: it would be more accurate to update the read pair count with the new distal region
							Region r = new Region(it1.distalRegion.getGenome(), it1.distalRegion.getChrom(), 
									Math.min(start2, start1), Math.max(end1, end2));
							it1.distalRegion  = r;
							it2.distalRegion  = r;
							Point mid = r.getMidpoint();
							if (mid.distance(it1.distalPoint)<mid.distance(it2.distalPoint))
								it2.distalPoint = it1.distalPoint;
							else
								it1.distalPoint = it2.distalPoint;
						}
					}
				}
			}		
		}

		
		// for each TSS and its distal anchor, count how many read pairs lead to other distal anchors
		// it is like finding the other two edges of the triangle. The count is the sum (over all distal anchors) 
		// of min connecting read count (btw tss-distal2 and distal-distal2).
		System.out.println("\nCount indirect / triangle read pairs, "+CommonUtils.timeElapsed(tic0));
		HashMap<String, ArrayList<Interaction>> gene2it = new HashMap<String, ArrayList<Interaction>>();
		for (Interaction it: interactions){
			if (!gene2it.containsKey(it.geneSymbol))
				gene2it.put(it.geneSymbol, new ArrayList<Interaction>());
			gene2it.get(it.geneSymbol).add(it);
		}
		for (String gene: gene2it.keySet()){
			ArrayList<Interaction> its = gene2it.get(gene);
			for (Interaction it: its){
//				System.out.println(it);
				ArrayList<Integer> indirectCounts = new ArrayList<Integer>();
				Region distal = it.distalRegion.expand(chiapet_radius, chiapet_radius);
				for (Interaction it2: its){
					if (it==it2)
						continue;
					int start = it2.distalRegion.getStart()-chiapet_radius;
					int end = it2.distalRegion.getEnd()+chiapet_radius;
					if (distal.getStart() > start){		// if it1 has higher coord, select by high and then low
						ArrayList<Integer> idx = CommonUtils.getPointsWithinWindow(highEnds, distal);
						ArrayList<ReadPair> rps = new ArrayList<ReadPair> ();
						for (int ii: idx){
							ReadPair rp2 = high.get(ii);
							int pos = rp2.r1.getLocation();
							if(pos>=start && pos<=end)
								rps.add(rp2);
						}
						indirectCounts.add(Math.min(rps.size(), it2.count));		// distal2--distal , it2_count
					}
					else{		// select by low and then high
						ArrayList<Integer> idx = CommonUtils.getPointsWithinWindow(lowEnds, distal);
						ArrayList<ReadPair> rps = new ArrayList<ReadPair> ();
						for (int ii: idx){
							ReadPair rp2 = low.get(ii);
							int pos = rp2.r2.getLocation();
							if(pos>=start && pos<=end)
								rps.add(rp2);
						}
						indirectCounts.add(Math.min(rps.size(), it2.count));		// distal--distal2 , it2_count
					}
				}
				int sum = 0;
				for (int c: indirectCounts)
					sum += c;
				it.indirectCount = sum;
				
//				System.out.println();
//				for (int c: indirectCounts)
//					System.out.print(c);
//				System.out.println();
			}
		}
	
		
		System.out.println("\n\nAnnotate and report, "+CommonUtils.timeElapsed(tic0));
		// report the interactions and annotations
		// annotate the proximal and distal anchors with TF and HM and regions
		StringBuilder sb = new StringBuilder();
		for (Interaction it: interactions){
			ArrayList<Integer> isOverlapped = new ArrayList<Integer>();
			Point mid = it.distalRegion.getMidpoint();
			int radius = it.distalRegion.getWidth()/2+chiapet_radius;
			for (ArrayList<Point> ps : allPeaks){
				ArrayList<Point> p = CommonUtils.getPointsWithinWindow(ps, mid, radius);
				isOverlapped.add(p.size());
				if (!it.isTfAnchord && !p.isEmpty()){	// if the distal point has not been anchored by a TF site
					it.distalPoint = p.get(p.size()/2);
					it.isTfAnchord = true;
				}
			}
			for (List<Region> rs: allRegions){
				isOverlapped.add(CommonUtils.getRegionsOverlapsWindow(rs, it.distalRegion, chiapet_radius).size());
			}
			// proximal
			for (ArrayList<Point> ps : allPeaks){
				ArrayList<Point> p = CommonUtils.getPointsWithinWindow(ps, it.tss, chiapet_radius);
				isOverlapped.add(p.size());
			}
			Region tssRegion = it.tss.expand(chiapet_radius);
			for (List<Region> rs: allRegions){
				isOverlapped.add(CommonUtils.getRegionsOverlapsWindow(rs, tssRegion, chiapet_radius).size());
			}
			// print out TF and region overlaps
			sb.append(it.toString()).append("\t");
			for (int b: isOverlapped)
				sb.append(b).append("\t");
			
			// print ChIA-PET call overlap info
			Point tssPoint = it.tss;
			int tssHalfWidth = it.tssRegion.getWidth()/2+chiapet_radius;
			int distalHalfWidth = it.distalRegion.getWidth()/2+chiapet_radius;
			Point distalPoint = it.distalPoint;
			Point tssLeft = new Point(genome, tssPoint.getChrom(), tssPoint.getLocation()-2000);
			Point tssRight = new Point(genome, tssPoint.getChrom(), tssPoint.getLocation()+2000);
			
			// GERM
			int index = Collections.binarySearch(germTss,  tssLeft);
			if( index < 0 )  							// if key not found
				index = -(index+1); 
			int indexRight = Collections.binarySearch(germTss,  tssRight);
			if( indexRight < 0 )  							// if key not found
				indexRight = -(indexRight+1); 
			// if key match found, continue to search ( binarySearch() give undefined index with multiple matches)
			boolean isGermOverlapped = false;
			indexRange: for (int i=index-1;i<=indexRight+2;i++){
				if (i<0 || i>=germTss.size())
					continue;
				try{
					Point tt = germTss.get(i);
					if (tt.distance(tssPoint)<=tssHalfWidth){ 
						if (!germTss2distals.containsKey(tt))
							continue;
						for (Point d:germTss2distals.get(tt)){
							if (d.distance(distalPoint)<=distalHalfWidth){
								isGermOverlapped = true;
								break indexRange;
							}
						}
					}
				}
				catch (IllegalArgumentException e){	// ignore								
				}								
			}
			sb.append(isGermOverlapped?"1\t":"0\t");
			
			// Mango
			// aPoints and bPoints are the midPoint of the two anchorRegions
			index = Collections.binarySearch(aPoints,  tssLeft);
			if( index < 0 )  							// if key not found
				index = -(index+1); 
			indexRight = Collections.binarySearch(aPoints,  tssRight);
			if( indexRight < 0 )  							// if key not found
				indexRight = -(indexRight+1); 
			// if key match found, continue to search ( binarySearch() give undefined index with multiple matches)
			boolean isMangoOverlapped = false;
			indexA: for (int i=index-1;i<=indexRight+2;i++){
				if (i<0 || i>=aPoints.size())
					continue;
				try{
					Point a = aPoints.get(i);
					if (a.distance(tssPoint)<=tssHalfWidth){ 
						if (!a2bs.containsKey(a))
							continue;
						for (Point b:a2bs.get(a)){
							if (b.distance(distalPoint)<=distalHalfWidth){
								isMangoOverlapped = true;
								break indexA;
							}
						}
					}
				}
				catch (IllegalArgumentException e){	// ignore								
				}
			}
			if (isMangoOverlapped)
				sb.append("1\t");
			else{
				index = Collections.binarySearch(bPoints,  tssLeft);
				if( index < 0 )  							// if key not found
					index = -(index+1); 
				indexRight = Collections.binarySearch(bPoints,  tssRight);
				if( indexRight < 0 )  							// if key not found
					indexRight = -(indexRight+1); 
				// if key match found, continue to search ( binarySearch() give undefined index with multiple matches)
				isMangoOverlapped = false;
				indexB: for (int i=index-1;i<=indexRight+2;i++){
					if (i<0 || i>=bPoints.size())
						continue;
					try{
						Point b = bPoints.get(i);
						if (b.distance(tssPoint)<=tssHalfWidth){ 
							if (!b2as.containsKey(b))
								continue;
							for (Point a:b2as.get(b)){
								if (a.distance(distalPoint)<=distalHalfWidth){
									isMangoOverlapped = true;
									break indexB;
								}
							}
						}
					}
					catch (IllegalArgumentException e){	// ignore								
					}								
				}
				sb.append(isMangoOverlapped?"1\t":"0\t");
			}
			
			CommonUtils.replaceEnd(sb, '\n');
			
		}
		CommonUtils.writeFile(Args.parseString(args, "out", "Result")+".readClusters.txt", sb.toString());
		
		// output BEDPE format
		sb = new StringBuilder();
		for (Interaction it: interactions){
			Region distalLocal = it.distalRegion.expand(read_merge_dist, read_merge_dist);
			Region tssLocal = it.tssRegion.expand(read_merge_dist, read_merge_dist);
			int distalLocalCount, tssLocalCount;
			if (it.distalPoint.offset(it.tss)<0){	// distal is in lower coord
				distalLocalCount = CommonUtils.getPointsWithinWindow(lowEnds, distalLocal).size();
				tssLocalCount = CommonUtils.getPointsWithinWindow(highEnds, tssLocal).size();
				sb.append(String.format("%s\t%s\t%d\t%d\t%d\n", distalLocal.toBED(), tssLocal.toBED(),
						it.count, distalLocalCount, tssLocalCount));
			}
			else{	// distal is in higher coord
				distalLocalCount = CommonUtils.getPointsWithinWindow(highEnds, distalLocal).size();
				tssLocalCount = CommonUtils.getPointsWithinWindow(lowEnds, tssLocal).size();
				sb.append(String.format("%s\t%s\t%d\t%d\t%d\n", tssLocal.toBED(), distalLocal.toBED(),
						it.count, tssLocalCount, distalLocalCount));
			}
		}
		CommonUtils.writeFile(Args.parseString(args, "out", "Result")+".bedpe", sb.toString());
		
		System.out.println("\n\nDone: "+CommonUtils.timeElapsed(tic0));
	}
	private class TSS implements Comparable<TSS>{
		public TSS(String symbol, StrandedPoint coord, int id){
			this.symbol = symbol;
			this.coord = coord;
			this.id = id;
		}
		String symbol;
		int id;
		StrandedPoint coord;
		public int compareTo(TSS t) {
			return coord.compareTo(t.coord);
		}
	}
	
	private class TSSwithReads{
		String symbol;
		int id;
		StrandedPoint coord;
		TreeMap<Integer, ArrayList<Boolean>> reads;		// distal read offset --> binary binding indicator
	}
	
	/** r1 should have lower coordinate than r2 */
	private class ReadPair implements Comparable<ReadPair>{
		Point r1;
		Point r2;
		
		public int compareRead1(ReadPair i) {
			return r1.compareTo(i.r1);
		}
		public int compareRead2(ReadPair i) {
			return r2.compareTo(i.r2);
		}
		@Override
		public int compareTo(ReadPair arg0) {
			// TODO Auto-generated method stub
			return 0;
		}
		public String toString(){
			return r1.toString()+"--"+r2.toString();
		}
	}
	private class ReadPairCluster implements Comparable<ReadPairCluster>{
		int r1min = Integer.MAX_VALUE;
		int r1max = -1;
		int r2min = Integer.MAX_VALUE;
		int r2max = -1;
		private ArrayList<ReadPair> reads = new ArrayList<ReadPair>();
		void addReadPair(ReadPair rp){
			if (r1min>rp.r1.getLocation())
				r1min=rp.r1.getLocation();
			if (r2min>rp.r2.getLocation())
				r2min=rp.r2.getLocation();
			if (r1max<rp.r1.getLocation())
				r1max=rp.r1.getLocation();
			if (r2max<rp.r2.getLocation())
				r2max=rp.r2.getLocation();
			reads.add(rp);
		}
		double getDensity(int padding){
			return reads.size()*10000000.0/((r1max-r1min+padding)*(r2max-r2min+padding));
		}
		@Override
		public int compareTo(ReadPairCluster arg0) {
			// TODO Auto-generated method stub
			return 0;
		}
		void sortByRead1(){
			Collections.sort(reads, new Comparator<ReadPair>(){
	            public int compare(ReadPair o1, ReadPair o2) {
	                return o1.compareRead1(o2);
	            }
	        });
		}
		void sortByRead2(){
			Collections.sort(reads, new Comparator<ReadPair>(){
	            public int compare(ReadPair o1, ReadPair o2) {
	                return o1.compareRead2(o2);
	            }
	        });
		}
		public String toString2(){
			StringBuilder sb = new StringBuilder();
			sb.append(reads.size()).append("=<");
			for (ReadPair rp: reads)
				sb.append(rp.toString()).append(",");
			CommonUtils.replaceEnd(sb, '>');
			return sb.toString();
		}
		public String toString0(){
			StringBuilder sb = new StringBuilder();
			sb.append(reads.size()).append("=<");
			sb.append(reads.get(0).r1.getChrom()).append(":").append(r1min).append("-").append(r1max);
			sb.append("==").append(r2min).append("-").append(r2max).append(">");
			return sb.toString();
		}
		public String toString(){
			StringBuilder sb = new StringBuilder();
			sb.append(reads.size()).append("=<");
			sb.append(r1max-r1min);
			sb.append("==").append(r2max-r2min).append(">");
			return sb.toString();
		}
		
	}
	
	class Interaction{
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
//		double pvalue;
		public String toString(){
//			return String.format("%d %.1f\t< %s %s -- %s >", count, density, geneSymbol, tssRegion, distalRegion);
			int dist = distalPoint.offset(tss);
			int padding = Math.abs(dist/20);
			return String.format("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%.1f", geneSymbol, 
					(tss instanceof StrandedPoint)?(StrandedPoint)tss:tss, tssRegion, distalPoint, distalRegion, 
					tss.getChrom()+":"+(Math.min(Math.min(tssRegion.getStart(), distalRegion.getStart()), tss.getLocation())-padding)+"-"+
					(Math.max(Math.max(tssRegion.getEnd(), distalRegion.getEnd()), tss.getLocation())+padding), 
					tssRegion.getWidth(), distalRegion.getWidth(), dist, count, indirectCount, density);
		}
	}

	class InteractionCall{
		Point tss;
		String tssString;
		String geneSymbol;
		String geneID;
		Region distal;
		double pvalue;
		TreeSet<Point> overlapCoords = new TreeSet<Point>();
		
		public String toString(){
			int start, end;
			if (tss.getLocation() < distal.getMidpoint().getLocation()){	// if TSS is upstream
				start = tss.getLocation();
				end = distal.getEnd();
			}
			else{
				start = distal.getStart();
				end = tss.getLocation();
			}
			Region it = new Region(genome, tss.getChrom(), start, end).expand(2000, 2000);
			return String.format("%s\t%d\t%s\t%s\t%s\t%s\t%.2f", it, distal.getMidpoint().distance(tss), distal.toString(), tssString, geneSymbol, geneID, -Math.log10(pvalue));
		}
	}
}

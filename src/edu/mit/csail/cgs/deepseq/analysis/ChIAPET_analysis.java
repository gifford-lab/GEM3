package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeMap;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.analysis.TFBS_SpaitialAnalysis.Site;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class ChIAPET_analysis {
	Genome genome;
	TreeMap<Region, Interaction> r2it = new TreeMap<Region, Interaction>();
	String fileName = null;
	Set<String> flags;
	String[] args;
	
	public ChIAPET_analysis(String[] args){
		
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

		flags = Args.parseFlags(args);
		this.args = args;
		
		fileName = Args.parseString(args, "bedpe", null);	
	}
	
	public static void main(String args[]){
		ChIAPET_analysis analysis = new ChIAPET_analysis(args);
		int type = Args.parseInteger(args, "type", 0);
		analysis.cleanUpOverlaps();
		analysis.countGenesPerRegion();
		analysis.StatsTAD();
//		switch(type){
//		case 0:
//			analysis.loadFile(fileName);
//			break;
//		case 1:
//			analysis.loadFile(fileName);
//			analysis.stats(fileName);
//			break;
//		}
		
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
			Interaction it = new Interaction();
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
						Interaction it=null;
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
				Interaction it=null;
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
				Interaction it = r2it.get(r);
				sb.append(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%.2e", it.distal.toBED(), it.tss.expand(2000).toBED(), it.tss.toString(), it.geneID, it.geneSymbol, it.distal.toString(), it.distal.getWidth(), it.pvalue)).append("\n");
			}
			CommonUtils.writeFile(fileName.replace(".bedpe", (flags.contains("rm_self")?".rm_self":"") + ".merged.bedpe"), sb.toString());
		}
	}
	
	/**
	 * Count number of genes connected by each region (overlapping regions are merged).
	 * @param fileName
	 */
	void countGenesPerRegion(){
		ArrayList<Region> rs = new ArrayList<Region>();
		rs.addAll(r2it.keySet());		// rs is sorted, r2it is a TreeMap
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
				merged = r;
				HashSet<String> genes = new HashSet<String>();
				
				for (int id:ids){
					sb.append(id).append(",");
					genes.add(r2it.get(rs.get(id)).geneSymbol);
				}
				CommonUtils.replaceEnd(sb, '\t');
				for (String s:genes)
					sb.append(s).append(",");
				CommonUtils.replaceEnd(sb, '\t');
				sb.append(genes.size()).append("\n");
				ids.clear();
				ids.add(i);
			}
		}
		// finish the last merge
		sb.append(merged.toString()).append("\t").append(merged.getWidth()).append("\t");
		HashSet<String> genes = new HashSet<String>();
		for (int id:ids){
			sb.append(id).append(",");
			genes.add(r2it.get(rs.get(id)).geneSymbol);
		}
		CommonUtils.replaceEnd(sb, '\t');
		for (String s:genes)
			sb.append(s).append(",");
		CommonUtils.replaceEnd(sb, '\t');
		sb.append(genes.size()).append("\n");
		
//		System.out.println(sb.toString());
		CommonUtils.writeFile(fileName.replace(".bedpe", (flags.contains("rm_self")?".rm_self":"") + ".per_region_count.txt"), sb.toString());
	}
	
	void StatsTAD(){
		String tad_file = Args.parseString(args, "tad", null);
		ArrayList<Region> tads = CommonUtils.load_BED_regions(genome, tad_file).car();
		Collections.sort(tads);
		ArrayList<Interaction> itNonTAD = new ArrayList<Interaction>();
		ArrayList<Interaction> itSameTAD = new ArrayList<Interaction>();
		ArrayList<Interaction> itCrossTAD = new ArrayList<Interaction>();
		
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
		for (Interaction it: itSameTAD)
			sb.append(it.toString()).append("\t1\n");
		for (Interaction it: itCrossTAD)
			sb.append(it.toString()).append("\t2\n");
		for (Interaction it: itNonTAD)
			sb.append(it.toString()).append("\t0\n");
		CommonUtils.writeFile(fileName.replace(".bedpe", (flags.contains("rm_self")?".rm_self":"") + ".StatsTAD.txt"), sb.toString());
		
	}
	
	
	class Interaction{
		Point tss;
		String tssString;
		String geneSymbol;
		String geneID;
		Region distal;
		double pvalue;
		
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

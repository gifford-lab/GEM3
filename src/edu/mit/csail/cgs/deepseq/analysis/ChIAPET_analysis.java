package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.TreeMap;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class ChIAPET_analysis {
	Genome genome;
	TreeMap<Region, Interaction> r2it = new TreeMap<Region, Interaction>();
	
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
	}
	
	public static void main(String args[]){
		ChIAPET_analysis analysis = new ChIAPET_analysis(args);
		String fileName = Args.parseString(args, "bedpe", null);		
		analysis.loadFile(fileName);
	}
	
	void loadFile(String fileName){
		ArrayList<String> texts = CommonUtils.readTextFile(fileName);
		HashMap<Point, ArrayList<Region>> tmp = new HashMap<Point, ArrayList<Region>>();
		for (String line:texts){
			String f[] = line.trim().split("\t");
			Interaction it = new Interaction();
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
			
			// TODO: filter interactions that have distal regions containing the TSS?

			while (r2it.containsKey(it.distal))
				it.distal = it.distal.expand(-1, -1);
			r2it.put(it.distal, it);
			if (!tmp.containsKey(it.tss))
				tmp.put(it.tss, new ArrayList<Region>() );
			tmp.get(it.tss).add(it.distal);
		}
		
		// for each tss, merge overlapping distal regions (<1kb)
		for (Point tss: tmp.keySet()){
			ArrayList<Region> regions = tmp.get(tss);
			ArrayList<Region> mergedRegions = new ArrayList<Region>();
			Collections.sort(regions);
			Region previous = regions.get(0);
			ArrayList<Region> previousRegions = new ArrayList<Region>();
			previousRegions.add(previous);
			
			for (int i = 1; i < regions.size(); i++) {
	          Region region = regions.get(i);
				// if overlaps with previous region, combine the regions
				if (previous.overlaps(region)){
					previous = previous.combine(region);
					previousRegions.add(region);
				} 
				else{	// not overlap any more, update, then move to next one
					mergedRegions.add(previous);
					previous = region;
					// merge overlapping regions, update interactions
					if (previousRegions.size()>1){	// merged
						Interaction it=null;
						double bestPvalue = 1;
						for (Region r:previousRegions){
							it = r2it.get(r);
							r2it.remove(r);			// remove old one
							bestPvalue = Math.min(bestPvalue, it.pvalue);
						}
						it.distal = previous;
						it.pvalue = bestPvalue;
						r2it.put(previous, it);		// add merged region
					}
					previousRegions.clear();
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
		for (Region r:r2it.keySet()){
			Interaction it = r2it.get(r);
			System.out.println(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%.2e", it.distal.toBED(), it.tss.expand(2000).toBED(), it.tss.toString(), it.geneID, it.geneSymbol, it.distal.toString(), it.distal.getWidth(), it.pvalue));;
		}
	}
	class Interaction{
		Point tss;
		String geneSymbol;
		String geneID;
		Region distal;
		double pvalue;
	}
}

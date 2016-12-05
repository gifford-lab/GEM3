package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeMap;
import java.util.TreeSet;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils.Gene;
import edu.mit.csail.cgs.tools.utils.Args;

public class GeneAnnotationAnalysis {
	public static void main(String args[]){
		int type = Args.parseInteger(args, "type", 0);
		switch(type){
			case 1: calcInterGenicDistance(args);break;
			case 2: getTADStats(args);break;
		}
	}
	/**
	 * Assign TSS to the TADs, and get stats for the TADs
	 * @param args
	 */
	private static void getTADStats(String args[]){
		Genome genome = CommonUtils.parseGenome(args);
		System.out.println("Loading gene annotations ... ");
		String gene_anno = Args.parseString(args, "gene_anno", "");
		ArrayList<Gene> genes = CommonUtils.loadGeneAnnotations(gene_anno);
		
		mergeMultipleTSSs(genes);
		
		String tad_file = Args.parseString(args, "tad", null);
		ArrayList<Region> tads = CommonUtils.load_BED_regions(genome, tad_file).car();
		Collections.sort(tads);
		ArrayList<Region> gene2tad = new ArrayList<Region>();
		for (int i=0;i<genes.size();i++){
			gene2tad.add(null);
		}
		HashMap<Region, ArrayList<Gene>> tad2genes = new HashMap<Region, ArrayList<Gene>> ();
		for (int i=0;i<genes.size();i++){
			Gene g = genes.get(i);
			Point p = g.getTSS(genome);
			for (int j=0;j<tads.size();j++){
				Region tad = tads.get(j);
				if (tad.contains(p)){
					gene2tad.set(i, tad);
					if(!tad2genes.containsKey(tad))
						tad2genes.put(tad, new ArrayList<Gene>());
					tad2genes.get(tad).add(g);
					break;
				}
			}
		}
		// output
		StringBuilder sb = new StringBuilder();
		for (int i=0;i<genes.size();i++){
			Gene g = genes.get(i);
			Region tad = gene2tad.get(i);
			if (tad!=null)
				sb.append(g.name).append("\t").append(tad.toString()).append("\t")
				.append(tad.getWidth()).append("\t")
				.append(tad2genes.get(tad).size()).append("\t")
				.append(tad.getWidth()/tad2genes.get(tad).size()).append("\n");
			else
				sb.append(g.name).append("\t").append(0).append("\t")
				.append(0).append("\t").append(0).append("\n");
		}
		CommonUtils.writeFile(new File(gene_anno).getName().replace(".txt", ".TAD_anno.txt"), sb.toString());
		System.out.println("Done!");

	}
	

	/**
	 * The nearby TSSs are merged if they are within 1000bp of each other. 
	 * The median coordinate is used as the TssMerged, the TESs are recorded. 
	 * Then for each Tss, compute the shortest distance to TSSs with lower and higher coordinates, 
	 * same for TES (except also excluding the TESs of merged genes). 
	 * The chrom, strand and start position are recorded for each gene. 
	 * The strand information can then be used to convert lower/higher to up/downstream TSS/TES using the awk tool.
	 */
	private static void calcInterGenicDistance(String args[]){
		System.out.println("Loading gene annotations ... ");
		String gene_anno = Args.parseString(args, "gene_anno", "");
		ArrayList<Gene> genes = CommonUtils.loadGeneAnnotations(gene_anno);
		int dist = Args.parseInteger(args, "merge_dist", 1000);
		
		mergeNearbyTSS(genes, dist);
		
		ArrayList<Gene> genesByEnds = new ArrayList<Gene>();
		genesByEnds.addAll(genes);
		Collections.sort(genesByEnds, new Comparator<Gene>(){
            public int compare(Gene o1, Gene o2) {
                return o1.compareToByEnd(o2);
            }
		});
		
		// compute distance to nearest TSS and TES in up and downstream
		// do not separate plus and minus strands, use the TSS-strand to determine up or down
		
		// separate the genes by chroms
		TreeMap<String, ArrayList<Gene>> geneByChrom = new TreeMap<String, ArrayList<Gene>>();
		for (Gene g: genesByEnds){
			if (!geneByChrom.containsKey(g.chr))
				geneByChrom.put(g.chr, new ArrayList<Gene>());
			geneByChrom.get(g.chr).add(g);
		}
		TreeMap<String, int[]> endsByChrom = new TreeMap<String, int[]>();
		for (String chr: geneByChrom.keySet()){
			ArrayList<Gene> g_c = geneByChrom.get(chr);
			int[] ends = new int [g_c.size()];
			for (int i=0;i<g_c.size();i++)
				ends[i]=g_c.get(i).end;
			endsByChrom.put(chr, ends);
		}
		
		StringBuilder sb = new StringBuilder();
		sb.append("Gene\tChr\tStart\tisPlus\tlTSS\thTSS\tlTES\thTES\n");	//l: low coord, h: high coord
		for (Gene g: genes){
			sb.append(String.format("%s\t%s\t%d\t%d\t", g.name, g.chr, g.start, g.strand=='+'?1:0));
			
			// nearest TSS
			int tssId = Collections.binarySearch(genes, g, new Comparator<Gene>(){
	            public int compare(Gene o1, Gene o2) {
	                return o1.compareToByStart(o2);
	            }
			});
			
			// no need to check tssID, it is guaranteed to match a TSS
			// lower coord side
			if (tssId==0 || !g.chr.equalsIgnoreCase(genes.get(tssId-1).chr))
				sb.append("-1\t");
			else
				sb.append(g.start-genes.get(tssId-1).start).append("\t");
			// higher coord side
			if (tssId==genes.size()-1 || !g.chr.equalsIgnoreCase(genes.get(tssId+1).chr))
				sb.append("-1\t");
			else
				sb.append(genes.get(tssId+1).start-g.start).append("\t");	
			
			// nearest TES		
			int start = g.start;
			int[] ends = endsByChrom.get(g.chr);
			int tesId = Arrays.binarySearch(ends, start);
			if (tesId<0){
				tesId = -(tesId+1)-1; 		// if not match, set the lowTES to the position before the inserted position
				int highId = tesId+1;
				HashSet<Integer> mergedEnds = g.mergedEnds;
				if (mergedEnds.isEmpty())
					mergedEnds.add(g.end);
				while (tesId>=0 && mergedEnds.contains(ends[tesId]))	// exclude mergedEnds
					tesId--;
				while (highId<ends.length && mergedEnds.contains(ends[highId]))	// exclude mergedEnds
					highId++;
				if (tesId==-1)
					sb.append("-1\t");
				else
					sb.append(start - ends[tesId]).append("\t");	
				if (highId==ends.length)
					sb.append("-1\t");
				else
					sb.append(ends[highId] - start).append("\t");	
			}
			else
				sb.append("0\t0\t");	// TSS overlap with a TES

			sb.append("\n");		
		}
		CommonUtils.writeFile(Args.parseString(args, "out", "GeneAnno")+"_"+dist+"_byCoord.txt", sb.toString());
		
	}
	
	/**
	 * Merge multiple TSS of a gene into the median TSS
	 * @param genes the gene list will be modified
	 */
	private static void mergeMultipleTSSs (ArrayList<Gene> genes){

		System.out.print("Merge "+genes.size()+" TSSs to ");
		TreeMap<String, ArrayList<Gene>> name2genes = new TreeMap<String, ArrayList<Gene>>();
		for (Gene g:genes){
			if (!name2genes.containsKey(g.name))
				name2genes.put(g.name, new ArrayList<Gene>());
			name2genes.get(g.name).add(g);
		}
		ArrayList<Gene> toRemove = new ArrayList<Gene>();
		for (String n: name2genes.keySet()){
			ArrayList<Gene> gs  = name2genes.get(n);
			if (gs.size()>=2){
				Collections.sort(gs);
				gs.remove(gs.size()/2);			// keep the median TSS gene from being removed
				toRemove.addAll(gs);
			}
		}
		genes.removeAll(toRemove);
		genes.trimToSize();
		System.out.println(genes.size());

	}
	
	/**
	 * The nearby TSSs are merged if they are within dist of each other. 
	 * @param genes the gene list will be modified
	 * @param dist
	 */
	private static void mergeNearbyTSS(ArrayList<Gene> genes, int dist){

		System.out.println("TSS merge distance: "+dist);
		Collections.sort(genes);
		
		// merge nearby TSS within "dist", merged genes format: A|B|C
		System.out.println("Total TSSs: "+genes.size());
		ArrayList<HashSet<Integer>> mergedIdx = new ArrayList<HashSet<Integer>>();
		HashSet<Integer> idx = new HashSet<Integer>();
		for (int i=1;i<genes.size();i++){
			Gene g1 = genes.get(i-1);
			Gene g2 = genes.get(i);
			if (g1.chr.equals(g2.chr) && g2.start-g1.start<=dist){ // to merge
				idx.add(i-1);
				idx.add(i);
			}
			else if (!idx.isEmpty()){	// gap
				mergedIdx.add((HashSet<Integer>)idx.clone());
				idx.clear();
			}
		}
		System.out.println("Groups of nearby TSSs: "+mergedIdx.size());
		ArrayList<Integer> toRemove = new ArrayList<Integer>();
		ArrayList<Gene> toAddback = new ArrayList<Gene>();
		
		for (HashSet<Integer> idxs:mergedIdx){	// each group of nearby TSSs
			if (idxs.isEmpty())
				continue;
			TreeSet<String> ids = new TreeSet<String>();
			TreeSet<String> names = new TreeSet<String>();		// same gene names will show once only
			TreeSet<Integer> starts = new TreeSet<Integer>();
			HashSet<Integer> mergedEnds = new HashSet<Integer>();
			Gene g=null;
			for (int ii: idxs){
				g = genes.get(ii);
				ids.add(g.id);
				names.add(g.name);
				starts.add(g.start);
				mergedEnds.add(g.end);
				toRemove.add(ii);		// remove individual TSSs
			}
			// update g to be the representative gene TSS
			// names, could have duplicates
			if (names.size()>1){
				StringBuilder sbNames = new StringBuilder();
				String last = names.pollLast();
				for (String name: names)
					sbNames.append(name).append("|");
				sbNames.append(last);
				g.name = sbNames.toString().trim();
			}
			
			// ids
			if (ids.size()>1){
				StringBuilder sb = new StringBuilder();
				String last = ids.pollLast();
				for (String id: ids)
					sb.append(id).append("|");
				sb.append(last);
				g.id = sb.toString().trim();
			}

			// set the start to be the median TSS start
			int count = starts.size()/2;
			for (int s:starts){
				count--;
				if (count==0){
					g.start = s;
					break;
				}
			}
			g.mergedEnds = mergedEnds;
			toAddback.add(g);				// add genes back to the list
		}
		Collections.sort(toRemove);
		Collections.reverse(toRemove);
		for (int i: toRemove)
			genes.remove(i);
		genes.addAll(toAddback);
		
		genes.trimToSize();
		Collections.sort(genes, new Comparator<Gene>(){
            public int compare(Gene o1, Gene o2) {
                return o1.compareToByStart(o2);
            }
		});
		
		System.out.println("Total TSSs after merging: "+genes.size());
		
	}

	
}

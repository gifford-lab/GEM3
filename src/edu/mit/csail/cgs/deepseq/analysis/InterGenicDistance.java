package edu.mit.csail.cgs.deepseq.analysis;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.TreeMap;
import java.util.TreeSet;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.deepseq.discovery.kmer.Kmer;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.tools.utils.Args;

public class InterGenicDistance {
	public static void main(String args[]){
		InterGenicDistance analysis = new InterGenicDistance();
		System.out.println("Loading GENCODE annotations ... ");
		String gene_anno = Args.parseString(args, "gene_anno", "");
		ArrayList<String> texts = CommonUtils.readTextFile(gene_anno);
		int dist = Args.parseInteger(args, "merge_dist", 1000);
		System.out.println("TSS merge distance: "+dist);
		ArrayList<Gene> genes = new ArrayList<Gene>();
		
		// UCSC Table format : 
		// #bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	name2
		// #bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	geneSymbol	refseq	description
		for (int i=0;i<texts.size();i++){
			String t = texts.get(i);
			if (t.startsWith("#"))
				continue;
			Gene g = analysis.new Gene();
			String f[] = t.split("\t");
			g.id = f[1].trim();
			g.chr = f[2].replace("chr", "");
			g.strand = f[3].charAt(0);
			g.start = Integer.parseInt(f[4]);
			g.end = Integer.parseInt(f[5]);
			g.name = f[12].trim();
			genes.add(g);
		}
		
		Collections.sort(genes, new Comparator<Gene>(){
            public int compare(Gene o1, Gene o2) {
                return o1.compareToByStart(o2);
            }
		});
		
		ArrayList<Gene> genesByEnds = new ArrayList<Gene>();
		genesByEnds.addAll(genes);
		Collections.sort(genesByEnds, new Comparator<Gene>(){
            public int compare(Gene o1, Gene o2) {
                return o1.compareToByEnd(o2);
            }
		});
		
		// merge nearby TSS
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
			
			// no need to check tssID, it is guaranteed to match
			if (tssId==0 || !g.chr.equalsIgnoreCase(genes.get(tssId-1).chr))
				sb.append("-1\t");
			else
				sb.append(g.start-genes.get(tssId-1).start).append("\t");
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
	
	private class Gene implements Comparable<Gene> {
		String id;
		String name;
		String chr;
		int start;
		int end;
		char strand;
		HashSet<Integer> mergedEnds = new HashSet<Integer>();
		
		public int compareToByStart(Gene p) {
		    if (!chr.equals(p.chr)) {
		      return chr.compareTo(p.chr);
		    }
		    if (start < p.start) {
		      return -1;
		    }
		    if (start > p.start) {
		      return 1;
		    }
		    return 0;
		}
		
		public int compareToByEnd(Gene p) {
		    if (!chr.equals(p.chr)) {
		      return chr.compareTo(p.chr);
		    }
		    if (end < p.end) {
		      return -1;
		    }
		    if (end > p.end) {
		      return 1;
		    }
		    return 0;
		}

		@Override
		public int compareTo(Gene arg0) {
			// TODO Auto-generated method stub
			return 0;
		}
		
		public String toString(){
			return name+"-->"+chr+":"+start+"-"+end;
		}
	}
	
}

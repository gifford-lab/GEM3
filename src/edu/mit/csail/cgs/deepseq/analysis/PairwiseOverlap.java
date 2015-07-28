package edu.mit.csail.cgs.deepseq.analysis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeMap;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KMAC;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public class PairwiseOverlap{

	public static void main(String[] args) {
	
		int type = Args.parseInteger(args, "type", 1);
		switch(type){
			case 1: diffBinding(args, parseGenome(args)); break;
			case 2: diffKmer(args); break;
			case 3: diffRegions(args, parseGenome(args)); break;
		}
		
	}
	
	public static Genome parseGenome(String[] args) {
		Genome genome = null;
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
		return genome;
	}
	
	/**
	 * Report the relationship between a pair of binding calls: 
	 * <ol>
	 * <li> The distance between binding calls
	 * <li> The percentages of intersection relative the whole sets
	 * </ol>
	 * example: --species "Mus musculus;mm9" --TF1 "C:\Data\workspace\gse\CTCF_outputs\Ctcf_2_GEM_events.txt" --TF2 "C:\Data\workspace\gse\CTCF_outputs\Ctcf_1_GEM_events.txt" --win 100 --out TF1_vs_TF2.txt
	 */
	public static void diffBinding(String[] args, Genome genome){
		String name1 = Args.parseString(args, "name1", "TF1");
		String name2 = Args.parseString(args, "name2", "TF2");
	    	
		// load binding sites (TF1)
	    ArrayList<String> texts = CommonUtils.readTextFile(Args.parseString(args, "TF1", null));
		// group sites by chrom
		TreeMap<String, ArrayList<Point>> chrom2sites = new TreeMap<String, ArrayList<Point>>();
		ArrayList<Point> tf1_pts = new ArrayList<Point>();
		for (String line:texts){
			if (line.length()==0)
				continue;
			if (line.startsWith("#")||line.startsWith("Position"))
				continue;
			String[] f = line.split("\\s+");		// match one or more white space
			Point p = Point.fromString(genome, f[0]);
			tf1_pts.add(p);
			String chr = p.getChrom();
			if (!chrom2sites.containsKey(chr))
				chrom2sites.put(chr, new ArrayList<Point>());
			chrom2sites.get(chr).add(p);
		}		
		// sort sites in each chrom
		for (String chr: chrom2sites.keySet()){
			ArrayList<Point> sites = chrom2sites.get(chr);
			Collections.sort(sites);
		}
		System.out.println(name1+": "+tf1_pts.size());
			
		int win = Args.parseInteger(args, "win", 100);
		StringBuilder shared = new StringBuilder("# TF1: "+Args.parseString(args, "TF1", null)+"\n");
		shared.append("# TF2: "+Args.parseString(args, "TF2", null)+"\n");
		shared.append("#TF1\tTF2\tOffset\n");
		
		// load binding sites (TF2)
	    texts = CommonUtils.readTextFile(Args.parseString(args, "TF2", null));
	    ArrayList<Point> tf2_only_pts = new ArrayList<Point> ();
	    HashSet<Point> tf1_shared_pts = new HashSet<Point> ();
	    int tf2_total = 0;
	    int shared_count = 0;
		
		for (String line:texts){
			if (line.length()==0)
				continue;
			if (line.startsWith("#")||line.startsWith("Position"))
				continue;
			tf2_total++;
			String[] f = line.split("\\s+");		// match one or more white space
			Point tf2 = Point.fromString(genome, f[0]);
			if (chrom2sites.containsKey(tf2.getChrom())){
				ArrayList<Point> sites = chrom2sites.get(tf2.getChrom());
				ArrayList<Point> results = CommonUtils.getPointsWithinWindow(sites, tf2, win);

				if (results.isEmpty())
					tf2_only_pts.add(tf2);
				else{
					tf1_shared_pts.addAll(results);
					for (Point tf1:results){
						shared.append(tf1.toString()).append("\t").append(tf2.toString()).append("\t").append(tf2.offset(tf1)).append("\n");
						shared_count++;
					}
				}
			}
			else
				tf2_only_pts.add(tf2);
		}
		tf1_pts.removeAll(tf1_shared_pts);
		System.out.println(name2+": "+tf2_total);
		System.out.println(name1+" only: "+tf1_pts.size());
		System.out.println(name2+" only: "+tf2_only_pts.size());
		System.out.println("Shared: "+shared_count);
		System.out.println("(Note: If an event of TF1 overlaps with multiple events of TF2, it will be counted multiple times, vice versa. Therefore the total number may not add up exactly.)");
		
		StringBuilder tf1_only = new StringBuilder("# TF1: "+Args.parseString(args, "TF1", null)+"\n");
		StringBuilder tf2_only = new StringBuilder("# TF2: "+Args.parseString(args, "TF2", null)+"\n");
		for (Point p: tf1_pts)
			tf1_only.append(p.toString()).append("\n");		
		for (Point p: tf2_only_pts)
			tf2_only.append(p.toString()).append("\n");
		
		CommonUtils.writeFile(name1+"_"+name2+"_shared.txt", shared.toString());
		CommonUtils.writeFile(name1+"_only_vs_"+name2+".txt", tf1_only.toString());
		CommonUtils.writeFile(name2+"_only_vs_"+name1+".txt", tf2_only.toString());
	}

	/**
	 * Report the relationship between 2 list of regions: 
	 * <ol>
	 * <li> The distance between binding calls
	 * <li> The percentages of intersection relative the whole sets
	 * </ol>
	 * example: --type 3 --species "Mus musculus;mm9" --RS1 "C:\Data\workspace\gse\CTCF_outputs\Ctcf_2_GEM_events.txt" --RS2 "C:\Data\workspace\gse\CTCF_outputs\Ctcf_1_GEM_events.txt" --win 100 
	 */
	public static void diffRegions(String[] args, Genome genome){
	    	
		// load binding sites (TF1)
	    ArrayList<String> texts = CommonUtils.readTextFile(Args.parseString(args, "RS1", null));
		// group sites by chrom
		TreeMap<String, ArrayList<Region>> chrom2sites = new TreeMap<String, ArrayList<Region>>();
		ArrayList<Region> tf1_pts = new ArrayList<Region>();
		for (String line:texts){
			if (line.length()==0)
				continue;
			if (line.startsWith("#")||line.startsWith("Position"))
				continue;
			String[] f = line.split("\\s+");		// match one or more white space
			Region p = Region.fromString(genome, f[0]);
			tf1_pts.add(p);
			String chr = p.getChrom();
			if (!chrom2sites.containsKey(chr))
				chrom2sites.put(chr, new ArrayList<Region>());
			chrom2sites.get(chr).add(p);
		}		
		// sort sites in each chrom
		for (String chr: chrom2sites.keySet()){
			ArrayList<Region> sites = chrom2sites.get(chr);
			Collections.sort(sites);
		}
			
		int win = Args.parseInteger(args, "win", 50);
		StringBuilder shared = new StringBuilder("# RS1: "+Args.parseString(args, "RS1", null)+"\n");
		shared.append("# RS2: "+Args.parseString(args, "RS2", null)+"\n");
		shared.append("#RS1\tRS2\tRS1_ID\tRS2_ID\n");
		
		// load binding sites (TF2)
	    texts = CommonUtils.readTextFile(Args.parseString(args, "RS2", null));
	    HashSet<Region> tf1_shared_pts = new HashSet<Region> ();
	    int rs2_id=1;
		StringBuilder tf2_only = new StringBuilder("# RS2: "+Args.parseString(args, "RS2", null)+"\n");
		for (String line:texts){
			if (line.length()==0)
				continue;
			if (line.startsWith("#")||line.startsWith("Position"))
				continue;
			String[] f = line.split("\\s+");		// match one or more white space
			Region tf2 = Region.fromString(genome, f[0]);
			if (chrom2sites.containsKey(tf2.getChrom())){
				ArrayList<Region> sites = chrom2sites.get(tf2.getChrom());
				ArrayList<Region> results = CommonUtils.getRegionsOverlapsWindow(sites, tf2, win);

				if (results.isEmpty())
					tf2_only.append(tf2.toString()).append("\t").append(rs2_id).append("\n");
				else{
					tf1_shared_pts.addAll(results);
					for (Region tf1:results){
						shared.append(tf1.toString()).append("\t").append(tf2.toString()).append("\t")
						.append(tf1_pts.indexOf(tf1)+1).append("\t").append(rs2_id).append("\n");
					}
				}
			}
			else
				tf2_only.append(tf2.toString()).append("\t").append(rs2_id).append("\n");
			rs2_id++;
		}
		
		StringBuilder tf1_only = new StringBuilder("# RS1: "+Args.parseString(args, "RS1", null)+"\n");
		for (int i=0;i<tf1_pts.size();i++){
			Region p = tf1_pts.get(i);
			if (!tf1_shared_pts.contains(p))
				tf1_only.append(p.toString()).append("\t").append(i+1).append("\n");		
		}
		
		String name1 = Args.parseString(args, "name1", "RS1");
		String name2 = Args.parseString(args, "name2", "RS2");
		CommonUtils.writeFile(name1+"_"+name2+"_shared.txt", shared.toString());
		CommonUtils.writeFile(name1+"_only_vs_"+name2+".txt", tf1_only.toString());
		CommonUtils.writeFile(name2+"_only_vs_"+name1+".txt", tf2_only.toString());
	}
	
	public static void diffKmer(String[] args){
	    int k_min = Args.parseInteger(args, "k_min", 5);
	    int k_max = Args.parseInteger(args, "k_max", 15);
	    double p = Args.parseDouble(args, "p", -3);
	    
	    ArrayList<String> texts = CommonUtils.readFastaFile(Args.parseString(args, "fasta1", null));
	    String[] seqs = new String[texts.size()];
	    texts.toArray(seqs);	    
	    int t1 = seqs.length;
	    
	    texts = CommonUtils.readFastaFile(Args.parseString(args, "fasta2", null));
	    String[] seqs2 = new String[texts.size()];
	    texts.toArray(seqs2);	  
	    texts = null;
	    int t2 = seqs2.length;
	    
	    StringBuilder sb = new StringBuilder();
	    sb.append("#Fasta1:\t"+Args.parseString(args, "fasta1", null)+"\n");
	    sb.append("#Fasta2:\t"+Args.parseString(args, "fasta2", null)+"\n");
	    sb.append("#Total\t"+t1+"\t"+t2+"\n");
	    for (int k=k_min;k<=k_max;k++){	    	
		    HashMap<String, Integer> m1 = CommonUtils.countKmers(k, seqs);
		    HashMap<String, Integer> m2 = CommonUtils.countKmers(k, seqs2);
		    for (String key:m1.keySet()){
		    	int c1 = m1.get(key);
		    	int c2 = 0;
		    	if (m2.containsKey(key)){
		    		c2 = m2.get(key);
		    		m2.remove(key);
		    	}
		    	else{
		    		key = SequenceUtils.reverseComplement(key);
		    		if (m2.containsKey(key)){
			    		c2 = m2.get(key);
			    		m2.remove(key);
			    	}
		    	}
		    	double hgp = KMAC.computeHGP(t1, t2, c1, c2);
		    	if (hgp<p)
		    		sb.append(String.format("%s\t%d\t%d\t%.1f\t1\n", key, c1, c2, hgp));	  
		    	else{
		    		hgp = KMAC.computeHGP(t2, t1, c2, c1);
			    	if (hgp<p)
			    		sb.append(String.format("%s\t%d\t%d\t%.1f\t2\n", key, c1, c2, hgp));	 
		    	}
		    		
		    }
		    for (String key:m2.keySet()){
		    	int c2 = m2.get(key);
		    	double hgp = KMAC.computeHGP(t2, t1, c2, 0);
		    	if (hgp<p)
		    		sb.append(String.format("%s\t%d\t%d\t%.1f\t2\n", key, 0, c2, hgp));	
		    }
	    }
//		    System.out.print(sb.toString());
	    CommonUtils.writeFile("diff_kmer.txt", sb.toString());
	}
}

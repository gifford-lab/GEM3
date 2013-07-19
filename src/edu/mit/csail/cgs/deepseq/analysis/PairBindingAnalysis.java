package edu.mit.csail.cgs.deepseq.analysis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeMap;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KMAC;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public class PairBindingAnalysis {

	public static void main(String[] args) {
	
		int type = Args.parseInteger(args, "type", 1);
		switch(type){
			case 1: diffBinding(args, parseGenome(args)); break;
			case 2: diffKmer(args); break;
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
			
		int win = Args.parseInteger(args, "win", 100);
		StringBuilder out = new StringBuilder("# TF1: "+Args.parseString(args, "TF1", null)+"\n");
		out.append("# TF2: "+Args.parseString(args, "TF2", null)+"\n");
		out.append("#TF1\tTF2\tOffset\n");
		
		// load binding sites (TF2)
	    texts = CommonUtils.readTextFile(Args.parseString(args, "TF2", null));
	    ArrayList<Point> tf2_only_pts = new ArrayList<Point> ();
	    HashSet<Point> tf1_shared_pts = new HashSet<Point> ();
		for (String line:texts){
			if (line.length()==0)
				continue;
			if (line.startsWith("#")||line.startsWith("Position"))
				continue;
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
						out.append(tf1.toString()).append("\t").append(tf2.toString()).append("\t").append(tf2.offset(tf1)).append("\n");
					}
				}
			}
			else
				tf2_only_pts.add(tf2);
		}
		tf1_pts.removeAll(tf1_shared_pts);
		StringBuilder tf1_only = new StringBuilder("# TF1: "+Args.parseString(args, "TF1", null)+"\n");
		StringBuilder tf2_only = new StringBuilder("# TF2: "+Args.parseString(args, "TF2", null)+"\n");
		for (Point p: tf1_pts)
			tf1_only.append(p.toString()).append("\n");		
		for (Point p: tf2_only_pts)
			tf2_only.append(p.toString()).append("\n");
		
		String name1 = Args.parseString(args, "name1", "TF1");
		String name2 = Args.parseString(args, "name2", "TF2");
		CommonUtils.writeFile(name1+"_"+name2+"_shared.txt", out.toString());
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

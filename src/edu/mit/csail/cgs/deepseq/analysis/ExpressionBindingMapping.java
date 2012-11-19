package edu.mit.csail.cgs.deepseq.analysis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.TreeMap;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.analysis.MultiTF_Binding.Site;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public class ExpressionBindingMapping {

	/**
	 * @param produce the mapping between TSS of genes in the expression data and the binding sites
	 * example: --species "Mus musculus;mm9" --expression iHox.exoncounts.raw.txt.disp0.05.iHoxc6-vs-iHoxc10.edgeR_GLM_DE.txt --tss igenomes_ens_gene.110620.tss --binding iHox.all.events.table --out iHox_ExpressionBindingMap_win20k.txt
	 */
	public static void main(String[] args) {
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
	    
		// load gene names from expression data
		ArrayList<String> geneNames = new ArrayList<String>();
		ArrayList<String> texts = CommonUtils.readTextFile(Args.parseString(args, "expression", null));
		for (String line:texts){
			if (line.length()==0)
				continue;
			if (line.startsWith("#"))
				continue;
			String[] f = line.split("\\s+");		// match one or more white space
			geneNames.add(f[0]);
		}
		geneNames.trimToSize();
		
		// load TSS annotation file
		HashMap<String, String> gene2TSS = new HashMap<String, String>();
		texts = CommonUtils.readTextFile(Args.parseString(args, "tss", null));
		for (String line:texts){
			if (line.length()==0)
				continue;
			if (line.startsWith("#"))
				continue;
			String[] f = line.split("\\s+");		// match one or more white space
			String tss = f[0]+"\t"+f[5];
			String[] f1 = f[1].replace("ID=","").split(";");
			gene2TSS.put(f1[0], tss);
		}
		
		// load binding sites
		texts = CommonUtils.readTextFile(Args.parseString(args, "binding", null));
		// group sites by chrom
		TreeMap<String, ArrayList<Point>> chrom2sites = new TreeMap<String, ArrayList<Point>>();
		for (String line:texts){
			if (line.length()==0)
				continue;
			if (line.startsWith("#"))
				continue;
			String[] f = line.split("\\s+");		// match one or more white space
			Point p = Point.fromString(genome, f[0]);
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
			
		int win = Args.parseInteger(args, "win", 20000);
		StringBuilder out = new StringBuilder("#Gene\tCoord\tStrand\tSite\tOffset\n");
		for (String g:geneNames){
			String tss = gene2TSS.get(g);
			String[] f = tss.split("\t");
			Point t = Point.fromString(genome, f[0]);
			boolean isPlusString = f[1].equals("+");
			if (chrom2sites.containsKey(t.getChrom())){
				ArrayList<Point> sites = chrom2sites.get(t.getChrom());
				ArrayList<Point> results = getPointsWithinWindow(sites, t, win);
				for (Point p:results){
					int offset = isPlusString?p.offset(t):p.offset(t)*-1;			// binding - tss
					out.append(g).append("\t").append(t.toString()).append("\t").append(isPlusString?"+":"-").append("\t");
					out.append(p.toString()).append("\t").append(offset).append("\n");
				}
			}
		}
		
		CommonUtils.writeFile(Args.parseString(args, "out", "out.txt"), out.toString());

	}

	static private ArrayList<Point> getPointsWithinWindow(ArrayList<Point> sites, Point anchor, int win){
		ArrayList<Point> results = new ArrayList<Point>();
		Region r = anchor.expand(win);
		Point start = r.startPoint();
		Point end = r.endPoint();
		int startIndex = -1;
		int endIndex = -1;
		int i = Collections.binarySearch(sites, start);
		if (i<0)
			startIndex=-i-1;		// -index-1, the insertion point
		else
			startIndex = i;
		i = Collections.binarySearch(sites, end);
		if (i<0)
			endIndex=-i-2;			// -index-1-1, the point before the insertion point
		else
			endIndex = i;
		if (startIndex<=endIndex){
			for (int j=startIndex;j<=endIndex;j++){
				results.add(sites.get(j));
			}
		}
		return results;
	}
}

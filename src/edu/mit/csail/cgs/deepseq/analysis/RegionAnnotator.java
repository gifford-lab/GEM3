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
		RegionAnnotator mtb = new RegionAnnotator(args);
		int round = Args.parseInteger(args, "r", 2);
		int type = Args.parseInteger(args, "type", 0);
		switch(type){
		case 0:  	// print out the overlap percentage of each query region for specified annotated regions
			mtb.regionAnnotationOverlaps();
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
	
	private void regionAnnotationOverlaps(){
		ArrayList<Region> queryRegions = CommonUtils.loadRegionFile(Args.parseString(args, "region_file", region_file), genome);
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
}

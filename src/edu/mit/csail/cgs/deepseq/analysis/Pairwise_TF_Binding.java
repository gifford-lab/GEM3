package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrixImport;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSPeak;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class Pairwise_TF_Binding {


	public static void main(String[] args) {
		Genome genome=null;
		ArrayList<String> names = new ArrayList<String>();
		ArrayList<WeightMatrix> pwms = new ArrayList<WeightMatrix>();
		ArrayList<ArrayList<Site>> all_sites = new ArrayList<ArrayList<Site>>();
		double gc = 0.42;//mouse		gc=0.41 for human
		double wm_factor = 0.6;	// PWM threshold, as fraction of max score
		File dir;
		boolean oldFormat =  false;
		SequenceGenerator<Region> seqgen;

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

		Set<String> flags = Args.parseFlags(args);
		wm_factor = Args.parseDouble(args, "pwm_factor", wm_factor);
		gc = Args.parseDouble(args, "gc", gc);
		seqgen = new SequenceGenerator<Region>();
		seqgen.useCache(!flags.contains("no_cache"));
		
		// load Event And Motifs
		String sox2_gem = Args.parseString(args, "sox2_gem", null);
		String sox2_gps = Args.parseString(args, "sox2_gps", null);
		String sox2_pfm = Args.parseString(args, "sox2_pfm", null);
		String oct4_gem = Args.parseString(args, "oct4_gem", null);
		String oct4_gps = Args.parseString(args, "oct4_gps", null);
		String oct4_pfm = Args.parseString(args, "oct4_pfm", null);
		
		WeightMatrix sox2_wm = CommonUtils.loadPWM_PFM_file(sox2_pfm, gc);
		WeightMatrix oct4_wm = CommonUtils.loadPWM_PFM_file(oct4_pfm, gc);
		
		all_sites.add(loadSites(genome, seqgen, sox2_gem, null, 0));
		all_sites.add(loadSites(genome, seqgen, oct4_gem, null, 1));
		all_sites.add(loadSites(genome, seqgen, oct4_gps, null, 2));
		all_sites.add(loadSites(genome, seqgen, oct4_gps, oct4_wm, 3));
		all_sites.add(loadSites(genome, seqgen, sox2_gem, sox2_wm, 4));
		names.add("sox2_gem");
		names.add("oct4_gem");
		names.add("oct4_gps");
		names.add("oct4_motif");
		names.add("sox2_motif");
		printBindingOffsets(genome, seqgen, all_sites, names, sox2_wm, wm_factor);
	}
	
	private static ArrayList<Site> loadSites(Genome genome, SequenceGenerator<Region> seqgen, String filePath, WeightMatrix wm, int tf_id){
		Pairwise_TF_Binding ptb = new Pairwise_TF_Binding();
		int round = 50;
		int posShift=0, negShift=0;
		WeightMatrixScorer scorer=null;
		if (wm!=null){
			scorer = new WeightMatrixScorer(wm);
			posShift = wm.length()/2;				// shift to the middle of motif
			negShift = wm.length()-1-wm.length()/2;
		}
		try{
			List<GPSPeak> gpsPeaks = GPSParser.parseGPSOutput(filePath, genome);
			ArrayList<Site> sites = new ArrayList<Site>();
			for (GPSPeak p:gpsPeaks){
				Site site = ptb.new Site();
				site.tf_id = tf_id;
				
				if (wm!=null){		// use nearest motif as binding site
					Region region = p.expand(round);
					String seq = seqgen.execute(region).toUpperCase();	// here round is the range of window 
					int hit = CommonUtils.scanPWMoutwards(seq, wm, scorer, round, wm.getMaxScore()*0.6).car();
					if (hit!=-999){
						if (hit>=0)
							site.bs = new Point(p.getGenome(), p.getChrom(), region.getStart()+hit+posShift);
						else
							site.bs = new Point(p.getGenome(), p.getChrom(), region.getStart()-hit+negShift);
					}
					else
						site.bs = (Point)p;		// no motif found, still use original GPS call
				}
				else
					site.bs = (Point)p;
				
				sites.add(site);
			}
			return sites;
		}
		catch (IOException e){
			System.out.println(tf_id+" does not have valid GPS/GEM event call file.");
			return null;
		}
	}
	
	
	public Pairwise_TF_Binding(){
	}

	class Site implements Comparable<Site>{
		int tf_id;
		Point bs;
		int id;
		public int compareTo(Site s) {					// descending score
			return(bs.compareTo(s.bs));
		}
		public String toString(){
			return String.format("%s, tf=%d, id=%d", bs.toString(), tf_id, id);
		}
	}

	static private void printBindingOffsets(Genome genome, SequenceGenerator<Region> seqgen, ArrayList<ArrayList<Site>> all_sites, 
			ArrayList<String> names, WeightMatrix wm, double wm_factor){
		// classify sites by chrom
		TreeMap<String, ArrayList<Site>> chrom2sites = new TreeMap<String, ArrayList<Site>>();
		for (ArrayList<Site> sites:all_sites){
			for (Site s:sites){
				String chr = s.bs.getChrom();
				if (!chrom2sites.containsKey(chr))
					chrom2sites.put(chr, new ArrayList<Site>());
				chrom2sites.get(chr).add(s);
			}
		}
		// sort and index sites in each chrom
		for (String chr: chrom2sites.keySet()){
			ArrayList<Site> sites = chrom2sites.get(chr);
			Collections.sort(sites);
			for (int i=0;i<sites.size();i++)
				sites.get(i).id = i;
		}
		int range = 400;
		int seqRange = 20;
		
		int i=0;	// for Sox_GEM as anchor point
		ArrayList<float[]> profiles = new ArrayList<float[]>();
		for (int n=0;n<all_sites.size();n++){
			profiles.add(new float[range*2+1]);
		}
		System.out.println(names.get(i));
		ArrayList<Site> sites_TF = all_sites.get(i);		// all sites of this TF
		WeightMatrixScorer scorer = new WeightMatrixScorer(wm);
		StringBuilder site_sb = new StringBuilder("Site    \torient\t");
		for (int n=0;n<names.size();n++){
			site_sb.append(names.get(n)+"\t");
		}
		site_sb.deleteCharAt(site_sb.length()-1).append("\n");
		for (Site s:sites_TF){
			int id = s.id;			// id of this TF binding in all binding sites in the chrom
			int b = s.bs.getLocation();
			// figure out the direction of TF binding based on sequence motif match
			int direction = 0;		// not sure, because no PWM, or no PWM match
			if (wm!=null){
				String seq = seqgen.execute(s.bs.expand(seqRange)).toUpperCase();
				Pair<Integer, Double> hit = CommonUtils.scanPWMoutwards(seq, wm, scorer, seqRange, wm.getMaxScore()*wm_factor);
				if (hit.car()!=-999){
					if (hit.car()>=0)
						direction = 1;
					else
						direction = -1;
				}
			}
			site_sb.append(s.bs.toString()).append("\t").append(direction).append("\t");
			
			int[] offsets = new int[names.size()];
			for (int n=0;n<offsets.length;n++){
				offsets[n]=999;
			}
			
			// count the nearby binding calls in upstream and downstream direction
			ArrayList<Site> chromSites = chrom2sites.get(s.bs.getChrom());
			for(int p=id;p<chromSites.size();p++){		// downstream (including same BS)
				Site s2 = chromSites.get(p);
				int offset = s2.bs.getLocation()-b;
				if (offset>range)
					break;
				float[] profile = profiles.get(s2.tf_id);
				switch(direction){
					case 0: profile[offset+range]+=0.5;profile[-offset+range]+=0.5;break;
					case 1: profile[offset+range]++;break;
					case -1: profile[-offset+range]++;break;
				}
				if (offsets[s2.tf_id]==999){
					if (s2.tf_id==i && offset==0)	// skip this binding site itself
						continue;
					if (direction==0){
						offsets[s2.tf_id]=offset;
					}
					else
						offsets[s2.tf_id]=offset*direction;
				}
			}
			// upstream is separated to search outwards, starting from BS position
			for(int p=id-1;p>=0;p--){					// upstream
				Site s2 = chromSites.get(p);
				int offset = s2.bs.getLocation()-b;
				if (offset<-range)
					break;
				float[] profile = profiles.get(s2.tf_id);
				switch(direction){
					case 0: profile[offset+range]+=0.5;profile[-offset+range]+=0.5;break;
					case 1: profile[offset+range]++;break;
					case -1: profile[-offset+range]++;break;
				}
				if (offsets[s2.tf_id]==999){
					if (s2.tf_id==i && offset==0)	// skip this binding site itself
						continue;
					if (direction==0){
						offsets[s2.tf_id]=offset;
					}
					else
						offsets[s2.tf_id]=offset*direction;
				}
			}
			for (int n=0;n<offsets.length;n++){
				site_sb.append(offsets[n]).append("\t");
			}
			site_sb.deleteCharAt(site_sb.length()-1).append("\n");
		} // for each site
		
		// output
		String filename1 = names.get(i)+"_site_offsets.txt";
		CommonUtils.writeFile(filename1, site_sb.toString());			
		
		StringBuilder sb = new StringBuilder(names.get(i)+"\t");
		for (int n=0;n<names.size();n++){
			sb.append(names.get(n)+"\t");
		}
		sb.deleteCharAt(sb.length()-1).append("\n");
		
		for (int p=-range;p<=range;p++){
			sb.append(p).append("\t");
			for (int n=0;n<names.size();n++){
				float[] profile = profiles.get(n);
					sb.append(String.format("%.0f\t", profile[p+range]));
			}
			sb.deleteCharAt(sb.length()-1).append("\n");
		}
		String filename = names.get(i)+"_profiles.txt";
		CommonUtils.writeFile(filename, sb.toString());
	}
	

}

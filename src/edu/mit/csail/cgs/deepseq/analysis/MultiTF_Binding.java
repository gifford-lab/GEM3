package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedPoint;
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

public class MultiTF_Binding {

	Genome genome=null;
	ArrayList<String> names = new ArrayList<String>();
	ArrayList<WeightMatrix> pwms = new ArrayList<WeightMatrix>();
	ArrayList<ArrayList<Site>> all_sites = new ArrayList<ArrayList<Site>>();
	double gc = 0.42;//mouse		gc=0.41 for human
	double wm_factor = 0.6;	// PWM threshold, as fraction of max score
	File dir;
	boolean oldFormat =  false;
	boolean useKmerStrand = false;
	private SequenceGenerator<Region> seqgen;

	// command line option:  (the folder contains GEM result folders) 
	// --dir Y:\Tools\GPS\Multi_TFs\oct18_GEM --species "Mus musculus;mm9" --r 2 --pwm_factor 0.6 --expts expt_list.txt [--no_cache --old_format] 
	public static void main(String[] args) {
		MultiTF_Binding mtb = new MultiTF_Binding(args);
		int round = Args.parseInteger(args, "r", 2);
		int prefix = Args.parseInteger(args, "prefix", 0);
		switch(round){
		case 2:
		case 3:
		case 4:mtb.loadEventAndMotifs(2);		// GEM
				mtb.printBindingOffsets(2, prefix);
				break;
		case 1:	mtb.loadEventAndMotifs(1);		// GPS
				mtb.printBindingOffsets(1, prefix);
				break;
		case 9:	mtb.loadEventAndMotifs(1);		// Motif
				mtb.printMotifOffsets();
		default:mtb.loadEventAndMotifs(round);		// Fake GPS calls, but snap to nearest PSSM within 50bp 
				mtb.printBindingOffsets(round, prefix);
				break;
		}
		
	}
	
	public MultiTF_Binding(String[] args){
				
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
		oldFormat = flags.contains("old_format");
		useKmerStrand = flags.contains("kmer_strand");
		dir = new File(Args.parseString(args, "dir", "."));
		names = new ArrayList<String>();
		ArrayList<String> info = CommonUtils.readTextFile(Args.parseString(args, "expts", null));
		for (String txt: info){
			if (!txt.equals("")){
				String[] f = txt.split("\t");
				names.add(f[0]);
			}
		}
		wm_factor = Args.parseDouble(args, "pwm_factor", wm_factor);
		gc = Args.parseDouble(args, "gc", gc);
		seqgen = new SequenceGenerator<Region>();
		seqgen.useCache(!flags.contains("no_cache"));
	}
	
	private void loadEventAndMotifs(int round){

		for (int tf=0;tf<names.size();tf++){
			String name = names.get(tf);

			System.out.println(String.format("TF#%d: loading %s", tf, name));
			
			// load motif files
			WeightMatrix wm = null;
			File dir2= new File(dir, name);
			if (!oldFormat)
				dir2= new File(dir2, name+"_outputs");
			final String suffix = name+"_"+ (round>=2?round:1) +"_PFM";
			File[] files = dir2.listFiles(new FilenameFilter(){
				public boolean accept(File arg0, String arg1) {
					if (arg1.startsWith(suffix))
						return true;
					else
						return false;
				}
			});
			if (files.length==0){
				System.out.println(name+" does not have a motif PFM file.");
				pwms.add(null);
			}
			else{				// if we have valid PFM file
				wm = CommonUtils.loadPWM_PFM_file(files[0].getAbsolutePath(), gc);
				pwms.add( wm );
			}
			
			// load binding event files 
			File gpsFile = new File(dir2, name+"_"+ (round>=2?round:1) +
					(oldFormat?"_GPS_significant.txt":"_GEM_events.txt"));
			String filePath = gpsFile.getAbsolutePath();
			WeightMatrixScorer scorer = null;
			int posShift=0, negShift=0;
			if (round!=1&&round!=2&&round!=9&&wm!=null){	
				scorer = new WeightMatrixScorer(wm);
				posShift = wm.length()/2;				// shift to the middle of motif
				negShift = wm.length()-1-wm.length()/2;
			}
			try{
				List<GPSPeak> gpsPeaks = GPSParser.parseGPSOutput(filePath, genome);
				ArrayList<Site> sites = new ArrayList<Site>();
				for (GPSPeak p:gpsPeaks){
					Site site = new Site();
					site.tf_id = tf;
					site.strand = p.getKmerStrand();
					if (round!=1&&round!=2&&round!=9&&wm!=null){		// use nearest motif as binding site
						Region region = p.expand(round);
						String seq = seqgen.execute(region).toUpperCase();	// here round is the range of window 
						int hit = CommonUtils.scanPWMoutwards(seq, wm, scorer, round, wm.getMaxScore()*wm_factor).car();
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
				all_sites.add(sites);
			}
			catch (IOException e){
				System.out.println(name+" does not have valid GPS/GEM event call file.");
				System.exit(1);
			}
		}
	}
	
	class Site implements Comparable<Site>{
		int tf_id;
		Point bs;
		int id;
		char strand;
		public int compareTo(Site s) {					// descending score
			return(bs.compareTo(s.bs));
		}
	}

	private void printBindingOffsets(int round, int prefix){
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
		// for each TF as anchor point
		for (int i=0;i<names.size();i++){
			ArrayList<float[]> profiles = new ArrayList<float[]>();
			for (int n=0;n<names.size();n++){
				profiles.add(new float[range*2+1]);
			}
			System.out.println(names.get(i));
			ArrayList<Site> sites_TF = all_sites.get(i);		// all sites of this TF
			WeightMatrix wm = pwms.get(i);
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
				if (useKmerStrand && round!=1)
					direction = s.strand=='*'?0:(s.strand=='+'?1:-1);
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
			String filename1 = names.get(i)+(round==2?"":("_"+round))+"_site_offsets.txt";
			CommonUtils.writeFile(new File(dir, filename1).getAbsolutePath(), site_sb.toString());			
			
			StringBuilder sb = new StringBuilder(names.get(i).substring(prefix)+"\t");
			for (int n=0;n<names.size();n++){
				sb.append(names.get(n).substring(prefix)+"\t");
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
			String filename = names.get(i)+(round==2?"":("_"+round))+"_profiles.txt";
			CommonUtils.writeFile(new File(dir, filename).getAbsolutePath(), sb.toString());
		}
	}
	
	private void printAllTFs2SingleSiteOffsets(String siteListFile, String inputTableName){
		// input table has 3 columns: Name:Events[:Motifs]
		// if motif is present, snap to nearest motif
		// otherwise, use event positions directly, either GEM or GPS or list motif instances
		
		// load the list of sites
		ArrayList<StrandedPoint> points = StrandedPoint.fromFile(genome, siteListFile);
		
		// load the input table
		ArrayList<String> rows = CommonUtils.readTextFile(inputTableName);
		HashMap<String, String> name2event = new HashMap<String, String>();
		HashMap<String, String> name2motif = new HashMap<String, String>();
		for (String row: rows){
			String[] f = row.split("\t");
			name2event.put(f[0], f[1]);
			if (f.length==3)
				name2motif.put(f[0], f[2]);
		}
		
		
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
		// for each TF as anchor point
		for (int i=0;i<names.size();i++){
			ArrayList<float[]> profiles = new ArrayList<float[]>();
			for (int n=0;n<names.size();n++){
				profiles.add(new float[range*2+1]);
			}
			System.out.println(names.get(i));
			ArrayList<Site> sites_TF = all_sites.get(i);		// all sites of this TF
			WeightMatrix wm = pwms.get(i);
			WeightMatrixScorer scorer = new WeightMatrixScorer(wm);
			StringBuilder site_sb = new StringBuilder("Site    \torient\t");
			for (int n=0;n<names.size();n++){
				site_sb.append(names.get(n).substring(3)+"\t");
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

		}

	}
	
	/** For each binding event of A, scan with all motif, compute the offset between A and all others */
	private void printMotifOffsets(){
		// classify sites by chrom
		int seqRange = 250;
		for (int i=0;i<names.size();i++){
			System.out.println(names.get(i));
			ArrayList<Site> sites = all_sites.get(i);
			WeightMatrix wm = pwms.get(i);
			if (wm ==null)
				continue;
			String[] seqs = new String[sites.size()];
			for (int j=0;j<sites.size();j++){
				seqs[j] = seqgen.execute(sites.get(j).bs.expand(seqRange)).toUpperCase();
			}
			
			ArrayList[][] hits = new ArrayList[seqs.length][pwms.size()];
			
			for (int j=0;j<pwms.size();j++){
				WeightMatrix pwm = pwms.get(j);
				if (pwm==null)
					continue;
				WeightMatrixScorer scorer = new WeightMatrixScorer(pwm);
				for (int s=0;s<seqs.length;s++){
					hits[s][j]=CommonUtils.getAllPWMHit(seqs[s], pwm.length(), scorer, pwm.getMaxScore()*0.6);
				}
			}
			int seqLen = seqs[0].length();
			for (int j=0;j<pwms.size();j++){
				if (pwms.get(j)==null)
					continue;
				int range = seqLen - pwms.get(i).length()/2 - pwms.get(j).length()/2;
				int[] same = new int[range*2+1];
				int[] diff = new int[range*2+1];
				for (int s=0;s<seqs.length;s++){
					ArrayList<Integer> hitm = hits[s][i];
					ArrayList<Integer> hitj = hits[s][j];
					if (hitm.isEmpty()||hitj.isEmpty())
						continue;
					if (i==j){		//self comparison
						for (int a=0;a<hitm.size();a++){
							int pm = hitm.get(a);
							for (int b=a;b<hitm.size();b++){
								int pj = hitm.get(b);
								if ((pm>=0&&pj>=0) || (pm<0&&pj<0))
									same[pj-pm+range]++;
								else
									diff[-pj-pm+range]++;			// -pj to get the coord on the same strand as pm
							}
						}
					}
					else{
						for (int pm:hitm){
							for (int pj:hitj){
								if ((pm>=0&&pj>=0) || (pm<0&&pj<0))
									same[pj-pm+range]++;
								else
									diff[-pj-pm+range]++;			// -pj to get the coord on the same strand as pm
							}
						}
					}
				}
				StringBuilder sb = new StringBuilder();
				int x[]=new int[range*2+1];
				for (int p=0;p<same.length;p++){
					x[p]=p-range;
					sb.append(String.format("%d\t%d\t%d\n", x[p], same[p], diff[p]));
				}
				String filename = names.get(i)+"_motif_"+names.get(j)+".txt";
				CommonUtils.writeFile(new File(dir, filename).getAbsolutePath(), sb.toString());
			}			
		}
	}
}

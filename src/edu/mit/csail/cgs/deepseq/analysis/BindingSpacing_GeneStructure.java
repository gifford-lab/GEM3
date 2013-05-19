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

/**
 * Compute the spacing between TF binding sites to the gene structures (such as TSS, splice sites, etc.)
 * @author Y
 *
 */
public class BindingSpacing_GeneStructure {

	Genome genome=null;
	ArrayList<String> tf_names = new ArrayList<String>();
	ArrayList<WeightMatrix> pwms = new ArrayList<WeightMatrix>();
	ArrayList<ArrayList<Site>> all_TF_sites = new ArrayList<ArrayList<Site>>();
	ArrayList<ArrayList<Site>> all_ANNO_sites = new ArrayList<ArrayList<Site>>(); // all gene annotations, such as tss, exon start/end, each as a list
	ArrayList<String> annotation_names = new ArrayList<String>();
	double gc = 0.42;//mouse		gc=0.41 for human
	double wm_factor = 0.6;	// PWM threshold, as fraction of max score
	File dir;
	boolean oldFormat =  false;
	private SequenceGenerator<Region> seqgen;

	// command line option:  (the folder contains GEM result folders) 
	// --dir Y:\Tools\GPS\Multi_TFs\oct18_GEM --species "Mus musculus;mm9" --r 2 --pwm_factor 0.6 --expts expt_list.txt [--no_cache --old_format] 
	public static void main(String[] args) {
		BindingSpacing_GeneStructure bsgs = new BindingSpacing_GeneStructure(args);
		int round = Args.parseInteger(args, "r", 2);
		int prefix = Args.parseInteger(args, "prefix", 0);		// # of letters to remove from the expt name, such as 3 for "ES_".
		switch(round){
		case 2:
		case 3:bsgs.loadEvents(round);		// GEM
				bsgs.printBindingOffsets(round, prefix);
				break;
		case 1:	bsgs.loadEvents(1);		// GPS
				bsgs.printBindingOffsets(1, prefix);
				break;
		default:bsgs.loadEvents(round);		// Fake GPS calls, but snap to nearest PSSM within 50bp 
				bsgs.printBindingOffsets(round, prefix);
				break;
		}
		
	}
	
	public BindingSpacing_GeneStructure(String[] args){
				
	    try {
	        Pair<Organism, Genome> pair = Args.parseGenome(args);
	        if(pair != null) {
	            genome = pair.cdr();
	        } else {
	            String genomeString = Args.parseString(args,"g",null);		// text file with chrom lengths
	            if(genomeString != null){
	                genome = new Genome("Genome", new File(genomeString));
	            } else{
	                genome=null;
	            }
	        }
	    } catch (NotFoundException e) {
	      e.printStackTrace();
	    }

		Set<String> flags = Args.parseFlags(args);
		oldFormat = flags.contains("old_format");
		dir = new File(Args.parseString(args, "dir", "."));
		wm_factor = Args.parseDouble(args, "pwm_factor", wm_factor);
		gc = Args.parseDouble(args, "gc", gc);
		seqgen = new SequenceGenerator<Region>();
		seqgen.useCache(!flags.contains("no_cache"));

		tf_names = new ArrayList<String>();
		ArrayList<String> info = CommonUtils.readTextFile(Args.parseString(args, "expts", null));
		for (String txt: info){
			if (!txt.equals("")){
				String[] f = txt.split("\t");
				tf_names.add(f[0]);
			}
		}
		// load gene annotation file
		System.out.println("Loading GENCODE annotations ... ");
		String anno_file = Args.parseString(args, "gene_anno", null);
		ArrayList<String> texts = CommonUtils.readTextFile(anno_file);
		char strand = '*';
		ArrayList<Site> tss = new ArrayList<Site>();					// -1
		ArrayList<Site> first_exon_starts = new ArrayList<Site>();		// -2
		ArrayList<Site> internal_exon_starts = new ArrayList<Site>();		// -3
		ArrayList<Site> internal_exon_ends = new ArrayList<Site>();			// -4
		ArrayList<Site> last_exon_ends = new ArrayList<Site>();			// -5
		for (int i=0;i<texts.size();i++){
			String t = texts.get(i);
			if (t.startsWith("#"))
				continue;
			String f[] = t.split("\t");
			String type = f[2];
			if (type.equals("transcript")){		// only look for transcript
				String chr = f[0].replace("chr", "");
				strand = f[6].charAt(0);
				tss.add(new Site(-1, new Point(genome, chr, Integer.parseInt(f[strand=='+'?3:4])), strand, i));
				ArrayList<Integer> exon_starts = new ArrayList<Integer>();
				ArrayList<Integer> exon_ends = new ArrayList<Integer>();
				// for lines after transcript, should be exons until the next transcript
				for (int j=i+1;j<texts.size();j++){
					String t_e = texts.get(j);
					String f_e[] = t_e.split("\t");
					String type_e = f_e[2];
					if (type_e.equals("exon")){
						if(strand=='+'){
							exon_starts.add(Integer.parseInt(f_e[3]));
							exon_ends.add(Integer.parseInt(f_e[4]));
						}else{
							exon_starts.add(Integer.parseInt(f_e[4]));
							exon_ends.add(Integer.parseInt(f_e[3]));
						}
					}
					if (type_e.equals("transcript")||j==(texts.size()-1)){	// next transcript, or end of file
						first_exon_starts.add(new Site(-2, new Point(genome, chr, exon_starts.get(0)), strand, i+1));
						if (exon_starts.size()>2){
							for (int k=1;k<exon_starts.size();k++)
								internal_exon_starts.add(new Site(-3, new Point(genome, chr, exon_starts.get(k)), strand, i+k+1));
							for (int k=0;k<exon_ends.size()-1;k++)
								internal_exon_ends.add(new Site(-4, new Point(genome, chr, exon_ends.get(k)), strand, i+k+1));
						}
						last_exon_ends.add(new Site(-5, new Point(genome, chr, exon_ends.get(exon_ends.size()-1)), strand, i+exon_ends.size()));							
						break;
					}
				}
			}
		}
		removeDuplicateSites(tss);
		all_ANNO_sites.add(tss);
		annotation_names.add("GENCODE_TSS");
		removeDuplicateSites(first_exon_starts);
		all_ANNO_sites.add(first_exon_starts);
		annotation_names.add("GENCODE_FIRST_EXON_STARTS");
		removeDuplicateSites(internal_exon_starts);
		all_ANNO_sites.add(internal_exon_starts);
		annotation_names.add("GENCODE_INTERNAL_EXON_STARTS");
		removeDuplicateSites(internal_exon_ends);
		all_ANNO_sites.add(internal_exon_ends);
		annotation_names.add("GENCODE_INTERNAL_EXON_ENDS");
		removeDuplicateSites(last_exon_ends);
		all_ANNO_sites.add(last_exon_ends);
		annotation_names.add("GENCODE_LAST_EXON_ENDS");
		
		System.out.println("# of TSS is " + tss.size());
		System.out.println("# of FIRST_EXON_STARTS is " + first_exon_starts.size());
		System.out.println("# of INTERNAL_EXON_STARTS is " + internal_exon_starts.size());
		System.out.println("# of INTERNAL_EXON_ENDS is " + internal_exon_ends.size());
		System.out.println("# of LAST_EXON_ENDS is " + last_exon_ends.size());
	}
	
	// remove duplicate sites (same coord and strand) and sort
	private void removeDuplicateSites(ArrayList<Site> sites){
		HashMap<String, Site> maps = new HashMap<String, Site>();
		for (Site s:sites)
			maps.put(s.coord.toString()+s.strand, s);
		sites.clear();
		sites.addAll(maps.values());
		Collections.sort(sites);
	}
	
	private void loadEvents(int round){

		for (int tf=0;tf<tf_names.size();tf++){
			String name = tf_names.get(tf);

			System.out.println(String.format("TF#%d: loading %s ...", tf, name));
			
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
			if (round!=1&&round!=2&&wm!=null){	// if use motif position as bs
				scorer = new WeightMatrixScorer(wm);
				posShift = wm.length()/2;				// shift to the middle of motif
				negShift = wm.length()-1-wm.length()/2;
			}
			try{
				List<GPSPeak> gpsPeaks = GPSParser.parseGPSOutput(filePath, genome);
				ArrayList<Site> sites = new ArrayList<Site>();
				for (GPSPeak p:gpsPeaks){
					Site site = new Site();
					site.type_id = tf;
					//TODO: if use motif, set the strand 
					if (round>3 && wm!=null){		// use nearest motif as binding site
						Region region = p.expand(round);
						String seq = seqgen.execute(region).toUpperCase();	// here round is the range of window 
						int hit = CommonUtils.scanPWMoutwards(seq, wm, scorer, round, wm.getMaxScore()*wm_factor).car();
						if (hit!=-999){
							if (hit>=0)
								site.coord = new Point(p.getGenome(), p.getChrom(), region.getStart()+hit+posShift);
							else
								site.coord = new Point(p.getGenome(), p.getChrom(), region.getStart()-hit+negShift);
						}
						else{
							site.coord = (Point)p;		// no motif found, still use original GPS call
						}
					}
					else{								// do not use motif, use GEM/GPS call
						site.coord = (Point)p;
						site.strand = p.getKmerStrand();
					}
					
					sites.add(site);
				}
				all_TF_sites.add(sites);
			}
			catch (IOException e){
				System.out.println(name+" does not have valid GPS/GEM event call file.");
				System.exit(1);
			}
		}
	}
	
	class Site implements Comparable<Site>{
		public Site(int type_id, Point coord, char strand, int id){
			this.type_id = type_id;
			this.coord = coord;
			this.strand = strand;
			this.id = id;
		}
		int type_id;			// TF id (positive) or gene annotation id (negative)
		Point coord;
		char strand;
		int id;					// site id in the same dataset (type)
		public int compareTo(Site s) {					// descending coordinates
			return(coord.compareTo(s.coord));
		}
		public Site(){
		}
		public String toString(){
			return String.format("type=%d, coord=%s:%s, id=%d", type_id, coord.toString(), strand, id);
		}
	}

	private void printBindingOffsets(int round, int prefix){
		// classify sites by chrom (including all TF and ANNO sites, sorted, for easy search
		// TODO: this may not be scalable when number of TF is large, for a few hundred expts, it works OK.
		TreeMap<String, ArrayList<Site>> chrom2sites = new TreeMap<String, ArrayList<Site>>();
		for (ArrayList<Site> sites:all_TF_sites){
			for (Site s:sites){
				String chr = s.coord.getChrom();
				if (!chrom2sites.containsKey(chr))
					chrom2sites.put(chr, new ArrayList<Site>());
				chrom2sites.get(chr).add(s);
			}
		}
		for (ArrayList<Site> sites:all_ANNO_sites){
			for (Site s:sites){
				String chr = s.coord.getChrom();
				if (!chrom2sites.containsKey(chr))
					chrom2sites.put(chr, new ArrayList<Site>());
				chrom2sites.get(chr).add(s);
			}
		}
		// sort and index sites in each chrom
		System.out.println("Sorting sites ... ");
		for (String chr: chrom2sites.keySet()){
			ArrayList<Site> sites = chrom2sites.get(chr);
			Collections.sort(sites);
			for (int i=0;i<sites.size();i++)
				sites.get(i).id = i;
		}
		int range = 400;
		int seqRange = 20;
		// for TSS as anchor point
		for (int i=0;i<annotation_names.size();i++){
			ArrayList<float[]> profiles = new ArrayList<float[]>();
			for (int n=0;n<tf_names.size();n++){
				profiles.add(new float[range*2+1]);
			}
			System.out.println("Analyzing "+annotation_names.get(i));
			ArrayList<Site> sites_ANNO = all_ANNO_sites.get(i);		// all sites of this TF
			StringBuilder site_sb = new StringBuilder("Site    \torient\t");
			for (int n=0;n<tf_names.size();n++){
				site_sb.append(tf_names.get(n)+"\t");
			}
			site_sb.deleteCharAt(site_sb.length()-1).append("\n");
			for (Site s:sites_ANNO){
				int id = s.id;			// id of this TF binding in all binding sites in the chrom
				int b = s.coord.getLocation();				
				site_sb.append(s.coord.toString()).append("\t").append(s.strand).append("\t");
				
				int[] offsets = new int[tf_names.size()];
				for (int n=0;n<offsets.length;n++){
					offsets[n]=999;
				}
				
				// count the nearby binding calls in upstream and downstream direction
				ArrayList<Site> chromSites = chrom2sites.get(s.coord.getChrom());
				for(int p=id;p<chromSites.size();p++){		// downstream (including same BS)
					Site s2 = chromSites.get(p);
					if (s2.type_id<0)
						continue;							// do not count the annotations sites
					int offset = s2.coord.getLocation()-b;
					if (offset>range)
						break;
					float[] profile = profiles.get(s2.type_id);
					switch(s.strand){
						case '+': profile[offset+range]++;break;
						case '-': profile[-offset+range]++;break;
					}
					if (offsets[s2.type_id]==999){		// only listed the first site of this type
						offsets[s2.type_id]=offset * (s.strand=='+'?1:-1);
					}
				}
				// upstream is separated to search outwards, starting from BS position
				for(int p=id-1;p>=0;p--){					// upstream
					Site s2 = chromSites.get(p);
					if (s2.type_id<0)
						continue;							// do not count the annotations sites
					int offset = s2.coord.getLocation()-b;
					if (offset<-range)
						break;
					float[] profile = profiles.get(s2.type_id);
					switch(s.strand){
						case '+': profile[offset+range]++;break;
						case '-': profile[-offset+range]++;break;
					}
					if (offsets[s2.type_id]==999){
						offsets[s2.type_id]=offset * (s.strand=='+'?1:-1);
					}
				}
				for (int n=0;n<offsets.length;n++){
					site_sb.append(offsets[n]).append("\t");
				}
				site_sb.deleteCharAt(site_sb.length()-1).append("\n");
			} // for each site
			
			// output
			String filename1 = annotation_names.get(i)+(round==2?"":("_"+round))+"_site_offsets.txt";
			CommonUtils.writeFile(new File(dir, filename1).getAbsolutePath(), site_sb.toString());			
			
			StringBuilder sb = new StringBuilder(annotation_names.get(i).substring(prefix)+"\t");
			for (int n=0;n<tf_names.size();n++){
				sb.append(tf_names.get(n).substring(prefix)+"\t");
			}
			sb.deleteCharAt(sb.length()-1).append("\n");
			sb.append("Total\t");
			for (int n=0;n<tf_names.size();n++){
				sb.append(all_TF_sites.get(n).size()+"\t");
			}
			sb.deleteCharAt(sb.length()-1).append("\n");
			for (int p=-range;p<=range;p++){
				sb.append(p).append("\t");
				for (int n=0;n<tf_names.size();n++){
					float[] profile = profiles.get(n);
						sb.append(String.format("%.0f\t", profile[p+range]));
				}
				sb.deleteCharAt(sb.length()-1).append("\n");
			}
			String filename = annotation_names.get(i)+(round==2?"":("_"+round))+"_profiles.txt";
			CommonUtils.writeFile(new File(dir, filename).getAbsolutePath(), sb.toString());
		}
	}		
}

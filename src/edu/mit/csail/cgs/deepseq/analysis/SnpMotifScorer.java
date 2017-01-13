package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.Set;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.discovery.kmer.GappedKmer;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KsmMotif;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KmerGroup;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;

public class SnpMotifScorer {

	public static void main(String[] args){
		// genome info and binding events
		Genome genome=null;
		Organism org=null;
		ArgParser ap = new ArgParser(args);
		Set<String> flags = Args.parseFlags(args);		
	    try {
	      Pair<Organism, Genome> pair = Args.parseGenome(args);
	      if(pair==null){
	        System.err.println("No genome provided; provide a Gifford lab DB genome name.");;System.exit(1);
	      }else{
	        genome = pair.cdr();
	        org = pair.car();
	      }
	    } catch (NotFoundException e) {
	      e.printStackTrace();
	    }
		SequenceGenerator<Region> seqgen = new SequenceGenerator<Region>();
		seqgen.useCache(!flags.contains("no_cache"));		
		seqgen.useLocalFiles(!flags.contains("use_db_genome"));
		if (!flags.contains("use_db_genome"))
			seqgen.setGenomePath(Args.parseString(args, "genome", ""));
		
		String vcf_flie = Args.parseString(args, "vcf", null);
		ArrayList<String> lines = new ArrayList<String>();
		if (vcf_flie!=null)
			lines = CommonUtils.readTextFile(vcf_flie);
		else
			System.exit(-1);
		
		int windowSize = Args.parseInteger(args, "win", 15);
		double gc = Args.parseDouble(args, "gc", 0.41);   //0.41 human, 0.42 mouse
		double pwm_cutoff = Args.parseDouble(args, "pwm_cutoff", 0.7);   // cutoff as the ratio of max PWM score
//		int width = windowSize*2+1;
//		Random randObj = new Random(Args.parseInteger(args, "rand_seed", 0));

		String motif_path = Args.parseString(args, "motif_dir", null);
		File dir = new File(motif_path);
		String pfm_file = Args.parseString(args, "pfms", null);
		String ksm_file = Args.parseString(args, "ksms", null);
		
		ArrayList<String> ksms = new ArrayList<String>();
		if (ksm_file!=null)
			ksms = CommonUtils.readTextFile(new File(dir, ksm_file).getAbsolutePath());
		ArrayList<KsmPwmScanner> kss = new ArrayList<KsmPwmScanner>();
		int count=0;
		for (String ksmPath : ksms){
			File file = new File(dir, ksmPath);
			KsmMotif ksm = GappedKmer.loadKSM(file);
			KsmPwmScanner scanner = new KsmPwmScanner(args, ksm);

			kss.add(scanner);
			System.out.println(String.format("K%d: %s", count, ksmPath)); 
			count++;
		}
		System.out.println();
		
		count=0;
		ArrayList<String> pfms = new ArrayList<String>();
		ArrayList<WeightMatrix> wms = new ArrayList<WeightMatrix>();
		if (pfm_file!=null)
			pfms = CommonUtils.readTextFile(new File(dir, pfm_file).getAbsolutePath());		
		for (String pfm : pfms){
			File file = new File(dir, pfm);
			WeightMatrix wm = CommonUtils.loadPWM_PFM_file(file.getAbsolutePath(), gc); 
			wms.add(wm);
			System.out.println(String.format("P%d: %s", count, pfm)); count++;
		}
		System.out.println();
		
		count = 0;
		for (String line: lines){
			count++;
			if (line.startsWith("#"))
				continue;
			String f[] = line.split("\t");	
			Point p = new Point(genome, f[0], Integer.parseInt(f[1])-1);		// -1, convert to 0-based coordinate for Point and sequence
			String seq = seqgen.execute(p.expand(windowSize)).toUpperCase();
			if (seq.charAt(windowSize)!=f[3].charAt(0)){
				System.err.println("Ref allele "+f[3]+" not match "+seq.charAt(windowSize)+" in sequence "+seq+", line #"+count);
				continue;
			}
			String alt = seq.substring(windowSize)+f[4]+seq.substring(windowSize+1, seq.length());
			String seqRC = SequenceUtils.reverseComplement(seq);
			String altRC = SequenceUtils.reverseComplement(alt);
			String[] info=f[7].split("=");
			String[] info1=info[1].split(";");
			int refc = Integer.parseInt(info1[0]);
			int altc = Integer.parseInt(info[2]);
			
			for (int i=0;i<kss.size();i++){
				KsmPwmScanner scanner = kss.get(i);
				KmerGroup kg = scanner.getBestKG(seq, seqRC);
				KmerGroup kgN = scanner.getBestKG(alt, altRC);

				double ksm = kg==null?0:kg.getScore();
				double ksm_alt = kgN==null?0:kgN.getScore();
				if (ksm!=ksm_alt){
					String match = "ZZ";
					if (kg!=null){
						Pair<Integer,Integer> ends = kg.getMatchEndIndices();
						int start=ends.car(), end=ends.cdr();
						if (start<0)
							start = 0;
						if (end>seq.length())
							end=seq.length();
						if (start<end){
							match = seq.substring(start,end)+"|"+SequenceUtils.reverseComplement(seq).substring(start,end);
						}
					}
					String match_alt = "ZZ";
					if (kgN!=null){
						Pair<Integer,Integer> ends = kgN.getMatchEndIndices();
						int start=ends.car(), end=ends.cdr();
						if (start<0)
							start = 0;
						if (end>alt.length())
							end=alt.length();
						if (start<end){
							match_alt = alt.substring(start,end)+"|"+SequenceUtils.reverseComplement(alt).substring(start,end);
						}
					}	
					boolean isSameDirection = (refc>altc && ksm>ksm_alt) || (refc<altc && ksm<ksm_alt); 
					System.out.println(String.format("%d\t%s\t%s\t%d\t%d\t%d\tKSM\tK_motif_%d\t%.1f\t%.1f\t%s\t%s", count, seq, f[7], refc, altc, isSameDirection?1:0, i, ksm, ksm_alt, match, match_alt));
				}
			}
			
			for (int i=0;i<wms.size();i++){
				WeightMatrix wm = wms.get(i);
				double pwm = WeightMatrixScorer.getMaxSeqScore(wm, seq, false);
				double pwm_alt = WeightMatrixScorer.getMaxSeqScore(wm, alt, false);
				if (pwm!=pwm_alt && (pwm>wm.getMaxScore()*pwm_cutoff || pwm_alt>wm.getMaxScore()*pwm_cutoff)){
					String match=WeightMatrixScorer.getMaxScoreSequence(wm, seq, -1000, 0);
					String match_alt=WeightMatrixScorer.getMaxScoreSequence(wm, alt, -1000, 0);
					boolean isSameDirection = (refc>altc && pwm>pwm_alt) || (refc<altc && pwm<pwm_alt); 
					System.out.println(String.format("%d\t%s\t%s\t%d\t%d\t%d\tPWM\tP_motif_%d\t%.1f\t%.1f\t%s\t%s", count, seq, f[7], refc, altc, isSameDirection?1:0, i, pwm, pwm_alt, match, match_alt));
				}
			}

		} // each expt
	}
	

}

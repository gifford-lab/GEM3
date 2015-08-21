package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrixImport;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.discovery.kmer.GappedKmer;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KMAC1;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KMAC1.KmerGroup;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KMAC1.MotifThreshold;
import edu.mit.csail.cgs.deepseq.discovery.kmer.Kmer;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSPeak;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.utils.stats.ROC;

public class KmerScanner {
	public static char[] letters = {'A','C','T','G'};
	private KMAC1 kEngine;
	// each element in the list is for one ChIP-Seq method
	
	public KmerScanner(ArrayList<Kmer> kmers, int posSeqCount, int negSeqCount, boolean use_base_kmer){
		kEngine = new KMAC1(kmers, null, use_base_kmer);
		kEngine.setTotalSeqCount(posSeqCount, negSeqCount);
	}
	
	public KmerGroup[] query (String seq){
//		KmerGroup[] kg = kEngine.findUnstrandedKmerHits_old(seq);
//		KmerGroup[] kg2 = kEngine.findUnstrandedKmerHits(seq);
		return kEngine.findUnstrandedKmerHits(seq);
	}
	
	public KmerGroup getBestKG (String seq){
		double bestScore = 0;
		KmerGroup best = null;
		KmerGroup[] kgs = query(seq);
		for (KmerGroup kg:kgs)
			if (bestScore < kg.getScore()){
				bestScore = kg.getScore();
				best = kg;
			}
		return best;
	}
	/**
	 * Find the file in path that match the type string
	 */
	public static String getFileName(String path, String type){
		File pathAll = new File(path);
		String name = pathAll.getName();
		File dir = pathAll.getParentFile();
		
		final String suffix = name + type;
		File[] files = dir.listFiles(new FilenameFilter(){
			public boolean accept(File arg0, String arg1) {
				if (arg1.startsWith(suffix))
					return true;
				else
					return false;
			}
		});
		if (files.length==0){
			System.out.println(name+" does not have a "+type+" file.");
			return null;
		}
		else{				// if we have valid file
			return files[0].getAbsolutePath();
		}
	}
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
		
		// load experiment list
		String path = Args.parseString(args, "path", "./");
		String fasta_path = Args.parseString(args, "fasta_path", "./");
		String fasta_suffix = Args.parseString(args, "fasta_suffix", ".fasta");
		String other_pfm_path = Args.parseString(args, "pfm_path", "./");
		String other_pfm_suffix = Args.parseString(args, "pfm_suffix", "");
		int windowSize = Args.parseInteger(args, "win", 50);
		double fpr = Args.parseDouble(args, "fpr", 0.1);
		double gc = Args.parseDouble(args, "gc", 0.41);   //0.41 human, 0.42 mouse
		int width = windowSize*2+1;
		int top = Args.parseInteger(args, "top", 5000);
		if (top==-1)
			top = Integer.MAX_VALUE;
		Random randObj = new Random(Args.parseInteger(args, "rand_seed", 0));

		ArrayList<String> lines = CommonUtils.readTextFile(Args.parseString(args, "expts", null));
		
		String[] pfm_suffixs = new String[0];
		if (!other_pfm_suffix.equals(""))
			pfm_suffixs = other_pfm_suffix.split(";");
		
		for (String line: lines){
			String f[] = line.split("\t");			
			scanSeqs(f[0], path, fasta_path, fasta_suffix, other_pfm_path, pfm_suffixs,
					flags.contains("use_base_kmer"), gc, top, randObj, width, fpr);
		    
		} // each expt
	}
	
	private static void scanSeqs(String expt, String path, String fasta_path, String fasta_suffix,  
			String other_pfm_path, String[] pfm_suffixs, boolean use_base_kmer, double gc,
			int top, Random randObj, int width, double fpr){
		
		System.out.println("Running "+expt);
		long tic = System.currentTimeMillis();
		String kmer=null, pfm=null, fasta_file=null;
		if (expt!=null){
			kmer = getFileName(path+expt, ".m0.KSM");			// old file name format
			if (kmer==null)
				kmer = getFileName(path+expt, "_KSM");		// new file name format, since May 2012
			pfm = getFileName(path+expt, ".all.PFM");
			fasta_file = fasta_path+expt+fasta_suffix;
		}
		
		long t1 = System.currentTimeMillis();
		File file = new File(kmer);
    	System.err.println(kmer);
		ArrayList<Kmer> kmers = GappedKmer.loadKmers(file);
		Pair<Integer, Integer> c = Kmer.getTotalCounts(file);
		
		KmerScanner scanner = new KmerScanner(kmers, c.car(), c.cdr(), use_base_kmer);
		System.out.println("KSM loading:\t"+CommonUtils.timeElapsed(t1));
	        	    
	    long t = System.currentTimeMillis();
	    WeightMatrix motif = CommonUtils.loadPWM_PFM_file(pfm, gc); 
    	System.err.println(pfm);
	    System.out.println("PWM loading:\t"+CommonUtils.timeElapsed(t));
		
	    // additional pfms
	    WeightMatrix[] otherPwms = new WeightMatrix[pfm_suffixs.length];
		ArrayList<ArrayList<Double>> other_scores = new ArrayList<ArrayList<Double>>();
		ArrayList<ArrayList<Double>> otherN_scores = new ArrayList<ArrayList<Double>>();
	    for (int i=0;i<pfm_suffixs.length;i++){
	    	pfm = other_pfm_path+expt+pfm_suffixs[i];
	    	System.err.println(pfm);
	    	otherPwms[i] = CommonUtils.loadPWM_PFM_file(pfm, gc); 
	    	other_scores.add(new ArrayList<Double>());
	    	otherN_scores.add(new ArrayList<Double>());
	    }
	    
		StringBuilder sb = new StringBuilder();

    	System.err.println(fasta_file);
		ArrayList<String> posSeqs = CommonUtils.loadSeqFromFasta(fasta_file);
		int numSeqToRun = Math.min(top, posSeqs.size());
		System.out.println("Scanning "+numSeqToRun+" regions ...");
				
		ArrayList<Double> pwm_scores = new ArrayList<Double>();
		ArrayList<Double> pwmN_scores = new ArrayList<Double>();
		ArrayList<Double> ksm_scores = new ArrayList<Double>();
		ArrayList<Double> ksmN_scores = new ArrayList<Double>();
		
		int PWM_time = 0;
		int KSM_time = 0;
		for (int i=0;i<numSeqToRun;i++){
			String seq = posSeqs.get(i).toUpperCase();
			int startSeq = seq.length()/2 - width/2; int endSeq =startSeq+width;
			seq = seq.substring(startSeq,endSeq);
			String seqN = SequenceUtils.dinu_shuffle(seq, randObj);

			// PWM
			long pwm_t = System.currentTimeMillis();
			double pwm = WeightMatrixScorer.getMaxSeqScore(motif, seq, false);
			String match=WeightMatrixScorer.getMaxScoreSequence(motif, seq, -1000, 0);
			pwm_scores.add(pwm);
			double pwmN = WeightMatrixScorer.getMaxSeqScore(motif, seqN, false);
			String matchN=WeightMatrixScorer.getMaxScoreSequence(motif, seqN, -1000, 0);
			pwmN_scores.add(pwmN);
			PWM_time += System.currentTimeMillis() - pwm_t;
			
			
			// KSM
			long ksm_t = System.currentTimeMillis();
			KmerGroup kg = scanner.getBestKG(seq);
			KmerGroup kgN = scanner.getBestKG(seqN);
			String matchKSM = "ZZ";
			if (i==31)
				i+=0;
			if (kg!=null){
				Pair<Integer,Integer> ends = kg.getMatchEndIndices();
				int start=ends.car(), end=ends.cdr();
				if (start<0)
					start = 0;
				if (end>seq.length())
					end=seq.length();
				if (start<end){
					matchKSM = seq.substring(start,end)+"|"+SequenceUtils.reverseComplement(seq).substring(start,end);
				}
			}
			String matchNKSM = "ZZ";
			if (kgN!=null){
				Pair<Integer,Integer> ends = kgN.getMatchEndIndices();
				int start=ends.car(), end=ends.cdr();
				if (start<0)
					start = 0;
				if (end>seqN.length())
					end=seqN.length();
				if (start<end){
					matchNKSM = seqN.substring(start,end)+"|"+SequenceUtils.reverseComplement(seqN).substring(start,end);
				}
			}					
						
			ksm_scores.add(kg==null?0:kg.getScore());
			ksmN_scores.add(kgN==null?0:kgN.getScore());
			KSM_time += System.currentTimeMillis() - ksm_t;

			// other PWMs
			String[] other_matches = new String[otherPwms.length];
			String[] otherN_matches = new String[otherPwms.length];
			StringBuilder sb_other_matchString = new StringBuilder();
			StringBuilder sb_other_score = new StringBuilder();
			for (int j=0;j<otherPwms.length;j++){
				other_matches[j]=WeightMatrixScorer.getMaxScoreSequence(otherPwms[j], seq, -1000, 0);
				double score = WeightMatrixScorer.getMaxSeqScore(otherPwms[j], seq, false);
				other_scores.get(j).add(score);
				otherN_matches[j]=WeightMatrixScorer.getMaxScoreSequence(otherPwms[j], seqN, -1000, 0);
				double scoreN = WeightMatrixScorer.getMaxSeqScore(otherPwms[j], seqN, false);
				otherN_scores.get(j).add(scoreN);
				sb_other_matchString.append("\t"+other_matches[j]+"\t"+otherN_matches[j]);
				sb_other_score.append(String.format("\t%.2f\t%.2f", score, scoreN));
			}

			sb.append(String.format("%d\t%s\t%s\t%s\t%s%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d%s\n",  
					i, match, matchN, matchKSM, matchNKSM, sb_other_matchString.toString(), pwm, pwmN, 
					kg==null?0:kg.getScore(), kgN==null?0:kgN.getScore(), 
					kg==null?0:-kg.getBestKmer().getHgp(), kgN==null?0:-kgN.getBestKmer().getHgp(), 
					kg==null?0:kg.getBestKmer().getPosHitCount(), kgN==null?0:kgN.getBestKmer().getPosHitCount(),
					sb_other_score.toString()));
		}
		
		System.out.println("Total PWM scanning time:" + PWM_time);
		System.out.println("Total KSM scanning time:" + KSM_time);
		
		CommonUtils.writeFile(expt+"_w"+width+"_scores.txt", sb.toString());
		System.out.println(expt+"_w"+width+"_scores.txt");
		
		System.out.print(String.format("%s\tPWM_KSM_FPR\t%.2f\t%.2f\t%.2f", 
				expt,evaluateScoreROC(pwm_scores, pwmN_scores, fpr),
				evaluateScoreROC(ksm_scores, ksmN_scores, fpr), fpr));
		for (int j=0;j<otherPwms.length;j++){
			System.out.print(String.format("\t%.2f", 
					evaluateScoreROC(other_scores.get(j), otherN_scores.get(j), fpr)));
		}
		System.out.println();
	}
	private static double evaluateScoreROC(ArrayList<Double> posScores, ArrayList<Double> negScores, double falsePositiveRate){
		double[] p = new double[posScores.size()];
		for (int i=0;i<p.length;i++)
			p[i]=posScores.get(i);
		double[] n = new double[negScores.size()];
		for (int i=0;i<n.length;i++)
			n[i]=negScores.get(i);
		
		ROC roc = new ROC(p, n);
		return roc.partialAUC(falsePositiveRate)/falsePositiveRate*100;
	}
	
	public static void main_old(String[] args){
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
		
		// load experiment list
		String path = Args.parseString(args, "path", "./");
		String fasta_path = Args.parseString(args, "fasta_path", "./");
		ArrayList<String> lines = CommonUtils.readTextFile(Args.parseString(args, "expts", null));
		for (String line: lines){
			String f[] = line.split("\t");
			String expt = f[0];
			System.out.println("Running "+expt);
			long tic = System.currentTimeMillis();
			// k-mer group info
//			String expt= Args.parseString(args, "out", "out");
//			System.out.println(path);
			String kmer=null, pfm=null, event=null;
			if (expt!=null){
				kmer = getFileName(path+expt, ".m0.KSM");			// old file name format
				if (kmer==null)
					kmer = getFileName(path+expt, "_KSM");		// new file name format, since May 2012
				pfm = getFileName(path+expt, ".all.PFM");
				event = path+expt+"_GPS_events.txt";
			}
			long t1 = System.currentTimeMillis();
			File file = new File(Args.parseString(args, "kmer", kmer));
			ArrayList<Kmer> kmers = GappedKmer.loadKmers(file);
	//		int clusterId = Args.parseInteger(args, "class", 0);
	//		ArrayList<Kmer> toRemove = new ArrayList<Kmer>();
	//		for (Kmer km:kmers)
	//			if (km.getClusterId()!=clusterId)
	//				toRemove.add(km);
	//		kmers.removeAll(toRemove);
	//		kmers.trimToSize();
			Pair<Integer, Integer> c = Kmer.getTotalCounts(file);
			KmerScanner scanner = new KmerScanner(kmers, c.car(), c.cdr(), !flags.contains("use_gapped_kmer"));
			System.out.println("KSM loading:\t"+CommonUtils.timeElapsed(t1));
		    
		    // PWM
	//		Pair<WeightMatrix, Double> wm = CommonUtils.loadPWM(args, org.getDBID());
	//		WeightMatrix motif = wm.car();
	//		double motifThreshold = wm.cdr();
		    	    
		    long t = System.currentTimeMillis();
		    WeightMatrix motif = CommonUtils.loadPWM_PFM_file(Args.parseString(args, "pfm", pfm), Args.parseDouble(args, "gc", 0.41)); //0.41 human, 0.42 mouse
		    System.out.println("PWM loading:\t"+CommonUtils.timeElapsed(t));
		    
			// event locations
			int windowSize = Args.parseInteger(args, "win", 50);
			String eventFile = Args.parseString(args, "event", event);
			List<GPSPeak> gpsPeaks = null;
			try{
				gpsPeaks = GPSParser.parseGPSOutput(eventFile, genome);
			}
			catch(IOException e){
				e.printStackTrace();
			}
	
			StringBuilder sb = new StringBuilder();
			ArrayList<Region> posRegions = new ArrayList<Region>();
			for(int i=0;i<gpsPeaks.size();i++){
				GPSPeak p = gpsPeaks.get(i);
				posRegions.add(p.expand(windowSize));
			}
			int width = windowSize*2+1;
			int top = Args.parseInteger(args, "top", 5000);
			TreeMap<Region, Region> reg2reg = new TreeMap<Region, Region>();		
			if (!flags.contains("di-shuffle")){
				int count = 0;
				
				each_event: for(int i=0;i<gpsPeaks.size();i++){
					GPSPeak p = gpsPeaks.get(i);
					if (p.getLocation()<1000+windowSize)
						continue each_event;
					Region rNeg = new Point(genome, p.getChrom(), p.getLocation()-1000).expand(windowSize);
					if (posRegions.get(i).getWidth()==width && rNeg.getWidth()==width){
						for(Region r:posRegions){
							if (rNeg.overlaps(r))
								continue each_event;
						}
						reg2reg.put(posRegions.get(i), rNeg);
						count++;
						if (count==5000)
							break;
					}
				}
				posRegions.clear();
				posRegions.addAll(reg2reg.keySet());
			}
			else{
				ArrayList<Region> tmp = new ArrayList<Region>();
				int num = top>posRegions.size()?posRegions.size():top;
				for (int i=0;i<num;i++){
					tmp.add(posRegions.get(i));
				}
				posRegions = tmp;
				tmp = null;
			}
			posRegions.trimToSize();
			Collections.sort(posRegions);
			System.out.println("Got regions:\t"+CommonUtils.timeElapsed(t));
			
			ArrayList<Double> pwm_scores = new ArrayList<Double>();
			ArrayList<Double> pwmN_scores = new ArrayList<Double>();
			ArrayList<Double> ksm_scores = new ArrayList<Double>();
			ArrayList<Double> ksmN_scores = new ArrayList<Double>();
			
			Random randObj = new Random(Args.parseInteger(args, "rand_seed", 0));
			System.out.println("Scanning "+posRegions.size()+" regions ...");
			int PWM_time = 0;
			int KSM_time = 0;
			for (Region r:posRegions){
				String seq = seqgen.execute(r).toUpperCase();
				
				String name_N = "----";
				String seqN;
	//			if (flags.contains("shuffle"))
	//				seqN = SequenceUtils.shuffle(seq, randObj);
				if (flags.contains("di-shuffle"))
					seqN = SequenceUtils.dinu_shuffle(seq, randObj);
				else{
					Region rN = reg2reg.get(r);
					seqN = seqgen.execute(rN).toUpperCase();
					name_N = rN.toString();
				}
				long pwm_t = System.currentTimeMillis();
				double pwm = WeightMatrixScorer.getMaxSeqScore(motif, seq, false);
				pwm_scores.add(pwm);
				double pwmN = WeightMatrixScorer.getMaxSeqScore(motif, seqN, false);
				pwmN_scores.add(pwmN);
				PWM_time += System.currentTimeMillis() - pwm_t;
				
				long ksm_t = System.currentTimeMillis();
				KmerGroup kg = scanner.getBestKG(seq);
				KmerGroup kgN = scanner.getBestKG(seqN);
				ksm_scores.add(kg==null?0:kg.getScore());
				ksmN_scores.add(kgN==null?0:kgN.getScore());
				KSM_time += System.currentTimeMillis() - ksm_t;
				sb.append(String.format("%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d\n", r.toString(), name_N, pwm, pwmN, 
						kg==null?0:kg.getScore(), kgN==null?0:kgN.getScore(), 
						kg==null?0:-kg.getBestKmer().getHgp(), kgN==null?0:-kgN.getBestKmer().getHgp(), 
						kg==null?0:kg.getBestKmer().getPosHitCount(), kgN==null?0:kgN.getBestKmer().getPosHitCount()));
			}
			
			System.out.println("Total PWM scanning time:" + PWM_time);
			System.out.println("Total KSM scanning time:" + KSM_time);
			
			CommonUtils.writeFile(expt+"_w"+width+"_scores.txt", sb.toString());
			System.out.println(expt+"_w"+width+"_scores.txt");
			
			if (flags.contains("compute_enrichment")){
				ArrayList<ScoreEnrichment> pwm_se = computeScoreEnrichments(pwm_scores, pwmN_scores);
				sb = new StringBuilder();
				for (ScoreEnrichment se: pwm_se){
					sb.append(String.format("%.2f\t%d\t%d\t%.2f\n", se.score, se.posHit, se.negHit, se.hgp));
				}
				CommonUtils.writeFile(expt+"_w"+width+"_pwm_enrichment.txt", sb.toString());
				
				ArrayList<ScoreEnrichment> ksm_se = computeScoreEnrichments(ksm_scores, ksmN_scores);
				sb = new StringBuilder();
				for (ScoreEnrichment se: ksm_se){
					sb.append(String.format("%.2f\t%d\t%d\t%.2f\n", se.score, se.posHit, se.negHit, se.hgp));
				}
				CommonUtils.writeFile(expt+"_w"+width+"_ksm_enrichment.txt", sb.toString());
			}
		    System.out.println(expt + " is done:\t"+CommonUtils.timeElapsed(tic));
		}
//		KmerGroup[] kgs = scanner.query("TATTTACATGCAGTGTCCGGAAGACGCCAGAAGAGGGCAGTAGATGCCCTAGTAGTGGAGC");
//		for (KmerGroup kg:kgs){
//			System.out.println(kg.toString());
//		}
	}

	private static ArrayList<ScoreEnrichment> computeScoreEnrichments(ArrayList<Double> posScores, ArrayList<Double> negScores){
		int total  = posScores.size();		
		double posSeqScores[] = new double[total];
		double negSeqScores[] = new double[total];
		for (int i=0;i<total;i++){
			posSeqScores[i]=posScores.get(i);
			negSeqScores[i]=negScores.get(i);			
		}
		ArrayList<ScoreEnrichment> ses = new ArrayList<ScoreEnrichment> ();
		Arrays.sort(posSeqScores);		
		Arrays.sort(negSeqScores);
		
		// find the threshold motif score
		TreeSet<Double> posScoreUnique = new TreeSet<Double>();
		for (double s:posSeqScores)
			posScoreUnique.add(s);
		Double[] posScores_u = new Double[posScoreUnique.size()];
		posScoreUnique.toArray(posScores_u);
		for (int i=0;i<posScores_u.length;i++){
			double score = posScores_u[i];
			if (score<=0)
				continue;
			ScoreEnrichment se = new ScoreEnrichment();
			se.score = score;
			int index = CommonUtils.findKey(posSeqScores, score);
			se.posHit = posSeqScores.length-index;
			index = CommonUtils.findKey(negSeqScores, score);
			se.negHit = negSeqScores.length-index;
			se.hgp = KMAC1.computeHGP(total, total, se.posHit, se.negHit);
			ses.add(se);
		}
		return ses;
	}
	
	private static class ScoreEnrichment{
		double score;
		int posHit;
		int negHit;
		double hgp;
	}
}

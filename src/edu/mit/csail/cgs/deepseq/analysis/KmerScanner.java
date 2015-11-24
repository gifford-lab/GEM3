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
import edu.mit.csail.cgs.deepseq.analysis.TFBS_SpaitialAnalysis.Site;
import edu.mit.csail.cgs.deepseq.discovery.Config;
import edu.mit.csail.cgs.deepseq.discovery.kmer.GappedKmer;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KMAC1;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KMAC1.KmerGroup;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KMAC1.MotifThreshold;
import edu.mit.csail.cgs.deepseq.discovery.kmer.Kmer;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KsmMotif;
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
	
	public KmerScanner(String[] args, ArrayList<Kmer> kmers, int posSeqCount, int negSeqCount, double[] seq_weights, boolean use_odds_ratio){
		kEngine = new KMAC1(kmers);
        Config config = new Config();
        try{
			config.parseArgs(args);   
		}
		catch (Exception e){
			e.printStackTrace();
    		System.exit(-1);
		}  
		kEngine.setConfig(config, "", "");
		kEngine.setTotalSeqCount(posSeqCount, negSeqCount);
		kEngine.setSequenceWeights(seq_weights);
		kEngine.setUseOddsRatio(use_odds_ratio);
	}
	
	public KmerGroup[] query (String seq){
//		return kEngine.findUnstrandedKsmGroupHits(seq);
		return kEngine.findKsmGroupHits(seq);
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
//		int round = Args.parseInteger(args, "r", 2); 	//GEM output round (1 for GPS, 2 for GEM)
//		int type = Args.parseInteger(args, "type", 999);
//		ArrayList<ArrayList<Site>> clusters=null;
//		switch(type){
//		case 999:	// default: simple file loading for RPD public code
//			analysis.loadBindingEvents();
//			clusters = analysis.mergeTfbsClusters();
//			analysis.outputTFBSclusters(clusters);
//			break;
//		case 0:
//			analysis.loadBindingEvents_old();
//			clusters = analysis.mergeTfbsClusters();
//			analysis.outputTFBSclusters(clusters);
//			break;
//		}
//	}
//		
//	private void scan_KSM_PWM(String[] args){
		Set<String> flags = Args.parseFlags(args);		
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
		double fpr = Args.parseDouble(args, "fpr", 0.1);
		double gc = Args.parseDouble(args, "gc", 0.41);   //0.41 human, 0.42 mouse
		int windowSize = Args.parseInteger(args, "win_size", 101);
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
			if (line.startsWith("#"))
				continue;
			scanSeqs(args, f[0], path, fasta_path, fasta_suffix, other_pfm_path, pfm_suffixs,
					flags.contains("or"), !flags.contains("use_seq_weights"), gc, top, randObj, windowSize, fpr);
		    
		} // each expt
	}
	
	private static void scanSeqs(String[] args, String expt, String path, String fasta_path, String fasta_suffix,  
			String other_pfm_path, String[] pfm_suffixs, boolean use_odds_ratio, boolean ignoreWeights, double gc,
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
		KsmMotif ksm = GappedKmer.loadKSM(file, ignoreWeights);
		KmerScanner scanner = new KmerScanner(args, ksm.kmers, ksm.posSeqCount, ksm.negSeqCount, ksm.seq_weights, use_odds_ratio);
		System.out.println("KSM loading:\t"+CommonUtils.timeElapsed(t1));
	        	    
	    long t = System.currentTimeMillis();
	    WeightMatrix motif = CommonUtils.loadPWM_PFM_file(pfm, gc); 
	    double maxPwmScore = motif.getMaxScore();
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
			if (kg!=null)
				matchKSM = kg.getCoveredSequence();
			String matchNKSM = "ZZ";
			if (kgN!=null)
				matchNKSM = kgN.getCoveredSequence();
						
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
		
		double pwm05 = maxPwmScore * 0.5; double count05=0;	// count of neg scores that pass 0.5*maxScore
		double pwm06 = maxPwmScore * 0.6; double count06=0;
		double pwm07 = maxPwmScore * 0.7; double count07=0;
		for (double s:pwmN_scores){
			if (s>=pwm05)
				count05++;
			if (s>=pwm06)
				count06++;
			if (s>=pwm07)
				count07++;
		}
		System.out.print(String.format("%s\tPWM_KSM_FPR_PWM05_PWM06_PWM07\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f", 
				expt,evaluateScoreROC(pwm_scores, pwmN_scores, fpr),
				evaluateScoreROC(ksm_scores, ksmN_scores, fpr), fpr,
				count05/pwmN_scores.size(), count06/pwmN_scores.size(), count07/pwmN_scores.size()));
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

}

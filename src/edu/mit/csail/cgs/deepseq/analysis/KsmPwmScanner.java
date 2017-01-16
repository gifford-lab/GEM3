package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.Random;
import java.util.Set;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.deepseq.discovery.Config;
import edu.mit.csail.cgs.deepseq.discovery.kmer.GappedKmer;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KMAC;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KmerGroup;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KsmMotif;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.utils.stats.ROC;
import net.sf.samtools.util.SequenceUtil;

public class KsmPwmScanner {
	public static char[] letters = {'A','C','T','G'};
	private KMAC kEngine;
	// each element in the list is for one ChIP-Seq method
	
	public KsmPwmScanner(String[] args, KsmMotif ksm){
        Config config = new Config();
        try{
			config.parseArgs(args);   
		}
		catch (Exception e){
			e.printStackTrace();
    		System.exit(-1);
		}  
		kEngine = new KMAC(ksm.kmers, config);
		kEngine.setTotalSeqCount(ksm.posSeqCount, ksm.negSeqCount);
		if (config.use_weighted_kmer)
			kEngine.setSequenceWeights(ksm.seq_weights);
	}
	
	public KmerGroup getBestKG (String seq, String seq_rc){
		KmerGroup[] kgs = kEngine.findKsmGroupHits(seq, seq_rc);
		if (kgs==null)
			return null;
		else
			return kgs[0];
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
		int type = Args.parseInteger(args, "type", 1);
		switch(type){
		case 1:	// default: simple file loading for RPD public code
			scan_KSM_PWM(args);
			break;
		case 2:
			shuffleFasta(args);
			break;
		case 3:
			rocFromScores(args);
			break;
		}
	}
	
	private static void rocFromScores (String[] args){
		double fpr = Args.parseDouble(args, "fpr", 0.1);
		ArrayList<String> lines = CommonUtils.readTextFile(Args.parseString(args, "pos", null));
		ArrayList<Double> posScores = new ArrayList<Double>();
		for (String s:lines){
			String f[] = s.split("\t");
			posScores.add(Double.parseDouble(f[7]));
		}
		posScores.trimToSize();
		lines = CommonUtils.readTextFile(Args.parseString(args, "neg", null));
		ArrayList<Double> negScores = new ArrayList<Double>();
		for (String s:lines){
			String f[] = s.split("\t");
			negScores.add(Double.parseDouble(f[7]));
		}
		negScores.trimToSize();
		double roc = evaluateScoreROC(posScores, negScores, fpr);
		System.out.println(String.format("%s\tTFFM\t%.2f", Args.parseString(args, "expt", null), roc));
	}
	
	private static void shuffleFasta(String[] args){
		Random rand = new Random(Args.parseInteger(args, "seed", 0));
		String fasta = Args.parseString(args, "fasta", null);
		if (fasta==null)
			System.err.println("File not found: "+fasta);
		ArrayList<String> posSeqs = CommonUtils.loadSeqFromFasta(fasta);
		StringBuilder sb = new StringBuilder();
		int count=1;
		for (String s: posSeqs){
			String seqN = SequenceUtils.dinu_shuffle(s, rand);
			sb.append(">Shuffled_").append(count).append("\n");
			sb.append(seqN).append("\n");
			count++;
		}
		CommonUtils.writeFile(fasta+".shuffled", sb.toString());
	}
	
	private static void scan_KSM_PWM(String[] args){
		// load experiment list
		String motif_path = Args.parseString(args, "path", "./");
		String fasta_path = Args.parseString(args, "fasta_path", "./");
		String fasta_suffix = Args.parseString(args, "fasta_suffix", ".fasta");
		String neg_fasta_suffix = Args.parseString(args, "neg_fasta_suffix", null);
		String other_pfm_path = Args.parseString(args, "pfm_path", "./");
		String other_pfm_suffix = Args.parseString(args, "pfm_suffix", "");
		double fpr = Args.parseDouble(args, "fpr", 0.1);
		double gc = Args.parseDouble(args, "gc", 0.41);   //0.41 human, 0.42 mouse
		int windowSize = Args.parseInteger(args, "win_size", 101);
		int top = Args.parseInteger(args, "top", 5000);
		if (top==-1)
			top = Integer.MAX_VALUE;
		int randSeed = Args.parseInteger(args, "rand_seed", 0);
		int negPosRatio = Args.parseInteger(args, "neg_ratio", 1);

		ArrayList<String> lines = CommonUtils.readTextFile(Args.parseString(args, "expts", null));
		
		String[] pfm_suffixs = new String[0];
		if (!other_pfm_suffix.equals(""))
			pfm_suffixs = other_pfm_suffix.split(";");
		
		for (String line: lines){
			String f[] = line.split("\t");	
			if (line.startsWith("#"))
				continue;
			if (f.length==1){
				scanSeqs(args, f[0], motif_path, f[0], fasta_path, fasta_suffix, neg_fasta_suffix, other_pfm_path, pfm_suffixs,
					gc, top, randSeed, negPosRatio, windowSize, fpr);
			}
			else{
				scanSeqs(args, f[0], motif_path, f[1], fasta_path, fasta_suffix, neg_fasta_suffix, other_pfm_path, pfm_suffixs,
					gc, top, randSeed, negPosRatio, windowSize, fpr);
			}
		    
		} // each expt
	}
	
	private static void scanSeqs(String[] args, String expt, String motif_path, String expt2, String fasta_path, String fasta_suffix,  
			String neg_fasta_suffix, String other_pfm_path, String[] pfm_suffixs, double gc,
			int top, int randSeed, int negPosRatio, int width, double fpr){
		
		System.out.println("Running "+expt);
		Random[] randObjs = new Random[negPosRatio];
		for (int i=0;i<negPosRatio;i++)
			randObjs[i] = new Random(randSeed+i);

		String kmer=null, pfm=null, fasta_file=null, fasta_neg_file=null;
		if (expt!=null){
			kmer = getFileName(motif_path+expt, ".m0.KSM");			// old file name format
			if (kmer==null)
				kmer = getFileName(motif_path+expt, "_KSM");		// new file name format, since May 2012
			pfm = getFileName(motif_path+expt, ".all.PFM");
			
			fasta_file = fasta_path+expt2+fasta_suffix;
			if (neg_fasta_suffix!=null)
				fasta_neg_file = fasta_path+expt2+neg_fasta_suffix;
		}
		
		long t1 = System.currentTimeMillis();
		File file = new File(kmer);
    	System.err.println(kmer);
		KsmMotif ksm = GappedKmer.loadKSM(file);
		KsmPwmScanner scanner = new KsmPwmScanner(args, ksm);
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
		
		ArrayList<String> negSeqs = null;
		if (fasta_neg_file != null )
			negSeqs = CommonUtils.loadSeqFromFasta(fasta_neg_file);
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
			//PWM
			long pwm_t = System.currentTimeMillis();
			double pwm = WeightMatrixScorer.getMaxSeqScore(motif, seq, false);
			String match=WeightMatrixScorer.getMaxScoreSequence(motif, seq, -1000, 0);
			pwm_scores.add(pwm);
			PWM_time += System.currentTimeMillis() - pwm_t;
			// KSM
			long ksm_t = System.currentTimeMillis();
			KmerGroup kg = scanner.getBestKG(seq, SequenceUtil.reverseComplement(seq));
			String matchKSM = "ZZ";
			if (kg!=null)
				matchKSM = kg.getCoveredSequence();
			ksm_scores.add(kg==null?0:kg.getScore());
			KSM_time += System.currentTimeMillis() - ksm_t;
			
			// scoring negative shuffled sequences
			for (int j=0;j<negPosRatio;j++){
				String seqN = fasta_neg_file!=null?
						negSeqs.get(i).toUpperCase().substring(startSeq,endSeq) :
						SequenceUtils.dinu_shuffle(seq, randObjs[j]);
				// PWM
				double pwmN = WeightMatrixScorer.getMaxSeqScore(motif, seqN, false);
				String matchN=WeightMatrixScorer.getMaxScoreSequence(motif, seqN, -1000, 0);
				pwmN_scores.add(pwmN);	
				// KSM
				KmerGroup kgN = scanner.getBestKG(seqN, SequenceUtil.reverseComplement(seqN));
				String matchNKSM = "ZZ";
				if (kgN!=null)
					matchNKSM = kgN.getCoveredSequence();
				ksmN_scores.add(kgN==null?0:kgN.getScore());
				
				sb.append(String.format("%d\t%s\t%s\t%s\t%s%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%d\n",  
						i, match, matchN, matchKSM, matchNKSM, "\tNA", pwm, pwmN, 
						kg==null?0:kg.getScore(), kgN==null?0:kgN.getScore(), 
						kg==null?0:-kg.getBestKmer().getHgp(), kgN==null?0:-kgN.getBestKmer().getHgp(), 
						kg==null?0:kg.getBestKmer().getPosHitCount(), kgN==null?0:kgN.getBestKmer().getPosHitCount()));
			}
			
//			TODO: print out all pos scores and then all neg scores
//			// other PWMs
//			String[] other_matches = new String[otherPwms.length];
//			String[] otherN_matches = new String[otherPwms.length];
//			StringBuilder sb_other_matchString = new StringBuilder();
//			StringBuilder sb_other_score = new StringBuilder();
//			for (int j=0;j<otherPwms.length;j++){
//				other_matches[j]=WeightMatrixScorer.getMaxScoreSequence(otherPwms[j], seq, -1000, 0);
//				double score = WeightMatrixScorer.getMaxSeqScore(otherPwms[j], seq, false);
//				other_scores.get(j).add(score);
//				otherN_matches[j]=WeightMatrixScorer.getMaxScoreSequence(otherPwms[j], seqN, -1000, 0);
//				double scoreN = WeightMatrixScorer.getMaxSeqScore(otherPwms[j], seqN, false);
//				otherN_scores.get(j).add(scoreN);
//				sb_other_matchString.append("\t"+other_matches[j]+"\t"+otherN_matches[j]);
//				sb_other_score.append(String.format("\t%.2f\t%.2f", score, scoreN));
//			}

//			sb.append(String.format("%d\t%s\t%s\t%s\t%s%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d%s\n",  
//					i, match, matchN, matchKSM, matchNKSM, sb_other_matchString.toString(), pwm, pwmN, 
//					kg==null?0:kg.getScore(), kgN==null?0:kgN.getScore(), 
//					kg==null?0:-kg.getBestKmer().getHgp(), kgN==null?0:-kgN.getBestKmer().getHgp(), 
//					kg==null?0:kg.getBestKmer().getPosHitCount(), kgN==null?0:kgN.getBestKmer().getPosHitCount(),
//					sb_other_score.toString()));
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
		System.out.print(String.format("%s\tPWM_KSM_FPR_PWM05_PWM06_PWM07\t%.2f\t%.2f\t%.3f\t%.3f\t%.3f\t%.3f", 
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

package edu.mit.csail.cgs.deepseq.analysis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedPoint;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KMAC;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KMAC1;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KMAC1.KmerGroup;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.utils.stats.StatUtil;
import net.sf.samtools.util.SequenceUtil;

public class MotifScan {

	public static void main(String[] args) {
	
		int type = Args.parseInteger(args, "type", 1);
		switch(type){
			case 1: findMotifInstances(args); break;
			case 9: scanWholeGenome(args, parseGenome(args)); break;
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
	
	/** TODO
	 */
	public static void scanWholeGenome(String[] args, Genome genome){
	    
	}

	public static void findMotifInstances(String[] args){
		Set<String> flags = Args.parseFlags(args);
	    
		String fasta = Args.parseString(args, "fasta", null);
		if (fasta==null)
			return;
	    ArrayList<String> texts = CommonUtils.readTextFile(fasta.trim());
	    int lineNum = texts.size();
	    String[] seqs = new String[lineNum/2];
	    String[] names = new String[lineNum/2];		// name is the first field of fasta header
	    for (int i=0;i<texts.size();i=i+2){
	    	seqs[i/2] = texts.get(i+1).toUpperCase();
	    	String s = texts.get(i);
	    	String f[] = s.substring(1, s.length()).split("\t");
	    	names[i/2] = f[0];
	    }
	    
	    // search for PWM motif matches
		ArrayList<MotifInstance> instances = null;
		ArrayList<Integer> motifLengths = new ArrayList<Integer>();
		ArrayList<Double> motifThresholds = new ArrayList<Double>();
		
		String pfm_fle = Args.parseString(args, "pfm", null);
		if (pfm_fle!=null){		// Load multiple PWMs
			StringBuilder sb_header = new StringBuilder();
		    List<WeightMatrix> pwms = CommonUtils.loadPWMs_PFM_file(pfm_fle, Args.parseDouble(args, "gc", 0.41));
		    if (pwms.isEmpty()){
		    	System.out.println("No motif is loaded from \n"+pfm_fle);
		    	System.exit(-1);
		    }
		    double scoreRatio = Args.parseDouble(args, "pwm_cutoff", 0.6);
		    
		    sb_header.append("# Motif Information\n");
		    sb_header.append("#ID\tName\tLetters\tWidth\tMax\tThresh\n");
		    
		    for (int m=0; m<pwms.size(); m++){
		    	WeightMatrix pwm = pwms.get(m);
		    	motifLengths.add(pwm.length());
		    	double threshold = pwm.getMaxScore()*scoreRatio;
		    	motifThresholds.add(threshold);
		    	sb_header.append("#").append(m).append("\t")
		    	.append(pwm.getName()).append("\t")
		    	.append(WeightMatrix.getMaxLetters(pwm)).append("\t")
		    	.append(pwm.length()).append("\t")
		    	.append(String.format("%.2f", pwm.getMaxScore())).append("\t")
		    	.append(String.format("%.2f", threshold)).append("\n");
		    }
		    sb_header.append("#");
		    System.out.println(sb_header.toString());

			instances = getPWMInstances(seqs, pwms, motifThresholds);
		}
		
	    // search for KSM motif matches
		String ksm_fle = Args.parseString(args, "ksm", null);
		if (ksm_fle!=null){		// Load multiple KSMs
			ArrayList<String> lines = CommonUtils.readTextFile(ksm_fle);
			ArrayList<KMAC1> kmacs = new ArrayList<KMAC1>();
			ArrayList<String> knames = new ArrayList<String>();
			for (String l:lines){
				if (l.startsWith("#"))
					continue;
				String[] f = l.split("\t");
				knames.add(f[0].trim());
				kmacs.add(CommonUtils.loadKsmFile(f[1].trim(), flags.contains("use_seq_weights"), flags.contains("or")));
			}
			instances = getKSMInstances(seqs, kmacs, knames);
		}
		
		// search for exact k-mer match
		String kmer = Args.parseString(args, "kmer", null);
		if (kmer!=null){
			StringBuilder sb_header = new StringBuilder();
		    
		    motifLengths.add(kmer.length());
		    
		    sb_header.append("# Motif Information\n");
		    sb_header.append("#ID\tLetters\tWidth\n");
		    sb_header.append("#").append(0).append("\t").append(kmer).append("\t").append(kmer.length()).append("\n");
		    sb_header.append("#");
		    System.out.println(sb_header.toString());
		    
			instances = getKmerInstances(seqs, kmer);
		}
		
	    // output
		StringBuilder sb = new StringBuilder("#Motif\tSeqID\tMotif_Name\tSeqName\tMatch\tSeqPos\tCoord\tStrand\tScore\n");
		Genome g = CommonUtils.parseGenome(args);
	    for (int i=0;i<instances.size();i++){
	    	MotifInstance mi = instances.get(i);
	    	String f[] = names[mi.seqID].split(" ");
	    	String regionString = null;
	    	if (f.length>1)
	    		regionString = f[1];
	    	else
	    		regionString = names[mi.seqID];
	    	
	    	String coor_string = null;
	    	Region region = Region.fromString(g, regionString);
	    	if (region!=null){
		    	Point startPoint = region.startPoint();
		    	int instanceLength = mi.matchSeq.length();
		    	int pos = mi.position + instanceLength/2;		// relative position from the left coord, adjust to motif midpoint
		    	if (mi.strand=='-')
		    		pos = mi.position + (instanceLength-1) - instanceLength/2;
		    	coor_string = new StrandedPoint(g, startPoint.getChrom(), startPoint.getLocation()+pos, mi.strand).toString();
	    	}
	    	else
	    		coor_string = "N.A.";
	    	sb.append(mi.motifID).append("\t").append(mi.seqID).append("\t").append(mi.motifName).append("\t").append(names[mi.seqID]).append("\t").append(mi.matchSeq).append("\t")
	    	.append(mi.position).append("\t").append(coor_string).append("\t").append(mi.strand).append("\t").append(String.format("%.2f", mi.score)).append("\n");
	    }	    
	    String out = Args.parseString(args, "out", fasta.substring(0, fasta.length()-6));
    	CommonUtils.writeFile(out.concat(".motifInstances.txt"), sb.toString());
	    
    	
	    // skip the matched sequences
	    if (flags.contains("skip")){
	    	HashSet<Integer> seqID_matched = new HashSet<Integer>();
	    	for (int i=0;i<instances.size();i++){
	    		seqID_matched.add(instances.get(i).seqID);
	    	}
	    	StringBuilder sb_skip = new StringBuilder();
	    	for (int s=0;s<seqs.length;s++){
	    		if (!seqID_matched.contains(s))
	    			sb_skip.append(">"+names[s]+"\n"+seqs[s]+"\n");
	    	}
	    	CommonUtils.writeFile(out.concat(".motifHitSkipped.fasta"), sb_skip.toString());
	    }
	    
	    // mask the matched instances
	    if (flags.contains("mask")){
	    	StringBuilder sb_mask = new StringBuilder();
	    	for (int i=0;i<instances.size();i++){
	    		MotifInstance mi = instances.get(i);
	    		StringBuilder masked = new StringBuilder(seqs[mi.seqID]);
	    		int endPos = mi.position + motifLengths.get(mi.motifID);		// exclusive end
	    		for (int idx=mi.position;idx<endPos;idx++)
	    			masked.setCharAt(idx, 'N');
	    		seqs[mi.seqID] = masked.toString();
	    	}
	    	for (int s=0;s<seqs.length;s++){
	    		sb_mask.append(">"+names[s]+"\n"+seqs[s]+"\n");
	    	}
	    	CommonUtils.writeFile(out.concat(".motifHitMasked.fasta"), sb_mask.toString());
	    }
	}
	
	public static ArrayList<MotifInstance> getKSMInstances(String[] seqs, ArrayList<KMAC1> kmacs, ArrayList<String> knames) {
		System.out.println("Scanning KSM motifs ...");
	    ArrayList<MotifInstance> instances = new ArrayList<MotifInstance>();
	    String[] seqs_rc = new String[seqs.length];
	    for (int i=0;i<seqs.length;i++)
	    	seqs_rc[i]=SequenceUtil.reverseComplement(seqs_rc[i]);
	    
	    for (int m=0; m<kmacs.size(); m++){
	    	System.out.println("  ... "+knames.get(m)+" ...");
	    	KMAC1 kmac = kmacs.get(m);
	    	for (int s=0; s<seqs.length;s++){
//	    		System.out.print(s+" ");
//	    		if (s==874) {
//	    			kmac.setIsDebugging(); // debug
//	    			System.out.println();
//	    		}
	    		KmerGroup[] kgs = kmac.findKsmGroupHits(seqs[s], seqs_rc[s]);
	    		for (int i=0;i<kgs.length;i++){
	    			KmerGroup kg = kgs[i];
		    		MotifInstance mi = new MotifInstance();
		    		mi.motifID = m;
		    		mi.motifName = knames.get(m);
		    		mi.score = kg.getScore();
		    		int pos = kg.getPosBS();
		    		if (pos > KMAC1.RC-seqs[s].length()*2){	// RC strand match
		    			mi.position = pos-KMAC1.RC;	
		    			mi.strand = '-';
		    		}
		    		else{
		    			mi.position = pos;	
		    			mi.strand = '+';
		    		}
		    		if (kg.getKmers().isEmpty())
		    			System.err.println("Empty Kmer group match: "+kg.toString());
		    		mi.matchSeq = kg.getCoveredSequence();
		    		mi.seqID = s;
		    		instances.add(mi);
	    		}
	    	}
	    }
		return instances;
	}

	/**
	 * Return a list of motif PWM matches from a list of sequences<br>
	 * If the sequence is shorter than the PWM, report a partial score
	 */
	public static ArrayList<MotifInstance> getPWMInstances(String[] seqs, List<WeightMatrix> pwms, List<Double> motifThresholds){
	    ArrayList<MotifInstance> instances = new ArrayList<MotifInstance>();
	    for (int m=0; m<pwms.size(); m++){
	    	WeightMatrix pwm = pwms.get(m);
	    	WeightMatrixScorer scorer = new WeightMatrixScorer(pwm, true);
	    	double threshold = motifThresholds.get(m);
	    	int width = pwm.length();
		    for (int s=0; s<seqs.length;s++){
		    	String str = seqs[s];
		    	if (str.length()>=width){
		    		WeightMatrixScoreProfile profiler = scorer.execute(str);
					for (int i=0;i<profiler.length();i++){
						double score = profiler.getHigherScore(i);
						if (score >= threshold){
							char strand = profiler.getHigherScoreStrand(i);
							String instance = null;
							if (strand=='+')
								instance = str.substring(i, i+width);
							else
								instance = SequenceUtils.reverseComplement(str.substring(i, i+width));
							MotifInstance mi = new MotifInstance();
							mi.motifID = m;
							mi.motifName = pwm.getName();
							mi.seqID = s;
							mi.matchSeq = instance;
							mi.position = i;
							mi.strand = strand;
							mi.score = score;
							instances.add(mi);
						}
					}
		    	}
		    	else{
		    		Pair<Float,Integer> partialScore = WeightMatrixScorer.scorePartialMatrix(pwm, str.toCharArray(), true);
		    		int p = partialScore.cdr()<WeightMatrixScorer.RC/2?partialScore.cdr():(partialScore.cdr()-WeightMatrixScorer.RC);		    		
		    		double partialThreshold = pwm.getPartialMaxScore(Math.abs(p),Math.abs(p)+str.length())*threshold/pwm.getMaxScore();
		    		if (partialScore.car()>=partialThreshold || partialScore.car() >= threshold/2){
			    		MotifInstance mi = new MotifInstance();
						mi.motifID = m;
						mi.strand = partialScore.cdr()<WeightMatrixScorer.RC/2?'+':'-';
						mi.motifName = pwm.getName()+"_p"+(mi.strand=='+'?partialScore.cdr():(partialScore.cdr()-WeightMatrixScorer.RC));
						mi.seqID = s;
						mi.matchSeq = mi.strand=='+'?str:SequenceUtils.reverseComplement(str);
						mi.position = 0;
						mi.score = partialScore.car();
						instances.add(mi);
		    		}
		    	}
		    }
	    }
		return instances;
	}
	/**
	 * Scan the sequences (forward and RC) for EXACT matches of one single k-mer 
	 * @param seqs
	 * @param kmer
	 * @return
	 */
	public static ArrayList<MotifInstance> getKmerInstances(String[] seqs, String kmer){
	    
	    ArrayList<MotifInstance> instances = new ArrayList<MotifInstance>();	    
	    for (int s=0; s<seqs.length;s++){
	    	String seq = seqs[s];
	    	char strand = '+';
	    	int pos = seq.indexOf(kmer);
	    	if (pos>=0){
	    		while(pos>=0){
		    		MotifInstance mi = new MotifInstance();
					mi.motifID = 0;
					mi.seqID = s;
					mi.matchSeq = kmer;
					mi.position = pos;
					mi.strand = strand;
					mi.score = 1;
					instances.add(mi);
					pos = seq.indexOf(kmer, pos+1);
	    		}
	    	}
	    	String kmer_rc = SequenceUtils.reverseComplement(kmer);
	    	strand = '-';
	    	pos = seq.indexOf(kmer_rc);
	    	if (pos>=0){
	    		while(pos>=0){
		    		MotifInstance mi = new MotifInstance();
					mi.motifID = 0;
					mi.seqID = s;
					mi.matchSeq = kmer;
					mi.position = pos;
					mi.strand = strand;
					mi.score = 1;
					instances.add(mi);
					pos = seq.indexOf(kmer, pos+1);
	    		}
	    	}
	    }
		return instances;
	}
}



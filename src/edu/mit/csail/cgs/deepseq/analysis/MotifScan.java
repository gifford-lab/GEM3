package edu.mit.csail.cgs.deepseq.analysis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KMAC;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public class MotifScan {

	public static void main(String[] args) {
	
		int type = Args.parseInteger(args, "type", 1);
		switch(type){
			case 1: reportMotifMatches(args); break;
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

	public static void reportMotifMatches(String[] args){
		StringBuilder sb_header = new StringBuilder();
		Set<String> flags = Args.parseFlags(args);
			    
		String fasta = Args.parseString(args, "fasta", null);
		if (fasta==null)
			return;
	    ArrayList<String> texts = CommonUtils.readTextFile(fasta);
	    int lineNum = texts.size();
	    String[] seqs = new String[lineNum/2];
	    String[] names = new String[lineNum/2];
	    for (int i=0;i<texts.size();i=i+2){
	    	seqs[i/2] = texts.get(i+1).toUpperCase();
	    	String s = texts.get(i);
	    	String f[] = s.substring(1, s.length()).split("\t");
	    	names[i/2] = f[0];
	    }
	    
	    List<WeightMatrix> pwms = CommonUtils.loadPWMs_PFM_file(Args.parseString(args, "pfm", null), Args.parseDouble(args, "gc", 0.41));
	    double scoreRatio = Args.parseDouble(args, "score_ratio", 0.6);
	    
	    sb_header.append("# Motif Information\n");
	    sb_header.append("#ID\tLetters\tWidth\tMax\tThresh\n");
	    for (int m=0; m<pwms.size(); m++){
	    	WeightMatrix pwm = pwms.get(m);
	    	double threshold = pwm.getMaxScore()*scoreRatio;
	    	sb_header.append("#").append(m).append("\t").append(WeightMatrix.getMaxLetters(pwm)).append("\t").append(pwm.length()).append("\t")
	    	.append(String.format("%.2f", pwm.getMaxScore())).append("\t").append(String.format("%.2f", threshold)).append("\n");
	    }
	    sb_header.append("#");
	    StringBuilder sb = new StringBuilder("#Motif\tSeqID\tSeqName\tMatch\tPos\tStrand\tScore\n");
	    ArrayList<MotifInstance> instances = new ArrayList<MotifInstance>();
	    for (int m=0; m<pwms.size(); m++){
	    	WeightMatrix pwm = pwms.get(m);
	    	WeightMatrixScorer scorer = new WeightMatrixScorer(pwm, true);
	    	double threshold = pwm.getMaxScore()*scoreRatio;
	    	int width = pwm.length();
		    for (int s=0; s<seqs.length;s++){		    	
				WeightMatrixScoreProfile profiler = scorer.execute(seqs[s]);
				for (int i=0;i<profiler.length();i++){
					double score = profiler.getMaxScore(i);
					if (score >= threshold){
						char strand = profiler.getMaxStrand(i);
						String instance = null;
						if (strand=='+')
							instance = seqs[s].substring(i, i+width);
						else
							instance = SequenceUtils.reverseComplement(seqs[s].substring(i, i+width));
						MotifInstance mi = new MotifInstance();
						mi.motifID = m;
						mi.seqID = s;
						mi.matchSeq = instance;
						mi.posisition = i;
						mi.strand = strand;
						mi.score = score;
						instances.add(mi);
					}
				}
		    }
	    }
	    
	    // output
	    for (int i=0;i<instances.size();i++){
	    	MotifInstance mi = instances.get(i);
	    	sb.append(mi.motifID).append("\t").append(mi.seqID).append("\t").append(names[mi.seqID]).append("\t").append(mi.matchSeq).append("\t")
	    	.append(mi.posisition).append("\t").append(mi.strand).append("\t").append(String.format("%.2f", mi.score)).append("\n");
	    }	    
	    System.out.println(sb_header.toString());
	    System.out.println(sb.toString());
	    
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
	    	String outName = fasta.substring(0, fasta.length()-4).concat(".motifHitSkipped.fasta");
	    	CommonUtils.writeFile(outName, sb_skip.toString());
	    }
	    
	    // mask the matched instances
	    if (flags.contains("mask")){
	    	StringBuilder sb_mask = new StringBuilder();
	    	for (int i=0;i<instances.size();i++){
	    		MotifInstance mi = instances.get(i);
	    		StringBuilder masked = new StringBuilder(seqs[mi.seqID]);
	    		int endPos = mi.posisition + pwms.get(mi.motifID).length();		// exclusive
	    		for (int idx=mi.posisition;idx<endPos;idx++)
	    			masked.setCharAt(idx, 'N');	    	
	    		seqs[mi.seqID] = masked.toString();
	    	}
	    	for (int s=0;s<seqs.length;s++){
	    		sb_mask.append(">"+names[s]+"\n"+seqs[s]+"\n");
	    	}
	    	String outName = fasta.substring(0, fasta.length()-6).concat(".motifHitMasked.fasta");
	    	CommonUtils.writeFile(outName, sb_mask.toString());
	    }
	}
	
}
class MotifInstance{
	int motifID;
	int seqID;
	String matchSeq;
	int posisition;
	char strand;
	double score;
}

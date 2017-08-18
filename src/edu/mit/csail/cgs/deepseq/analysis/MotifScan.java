package edu.mit.csail.cgs.deepseq.analysis;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedPoint;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.discovery.Config;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KMAC;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KmerGroup;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import net.sf.samtools.util.SequenceUtil;

public class MotifScan {

	public static void main(String[] args) {
	
		int type = Args.parseInteger(args, "type", 1);
		switch(type){
			case 1: findMotifInstances(args); break;
			case 9: scanWholeGenome(args, CommonUtils.parseGenome(args)); break;
		}
		
	}
	
	/** TODO
	 */
	public static void scanWholeGenome(String[] args, Genome genome){
	    
	}

	public static void findMotifInstances(String[] args){
		Set<String> flags = Args.parseFlags(args);
	    boolean toAddFasta = flags.contains("add_fasta");
	    boolean toMakeMatrix = flags.contains("matrix");		// make seq-motif feature matrix

		String fasta = Args.parseString(args, "fasta", null);
	    String out = Args.parseString(args, "out", fasta.replace(".fasta", "").replace(".fa", ""));

	    ArrayList<String> texts = CommonUtils.loadFasta(fasta.trim());
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
		String header = null;
		String pfm_fle = Args.parseString(args, "pfm", null);
		if (pfm_fle!=null){		// Load multiple PWMs in a single file
			StringBuilder sb_header = new StringBuilder();
		    List<WeightMatrix> pwms = CommonUtils.loadPWMs_PFM_file(pfm_fle, Args.parseDouble(args, "gc", 0.41));
		    if (pwms.isEmpty()){
		    	System.out.println("No motif PFM is loaded from \n"+pfm_fle);
		    	System.exit(-1);
		    }
		    
		    // display motif information
		    double scoreRatio = Args.parseDouble(args, "pwm_cutoff", 0.6);
		    sb_header.append("# Motif Information\n");
		    sb_header.append("#ID\tName\tLetters\tWidth\tMax\tThresh\n");
		    for (int m=0; m<pwms.size(); m++){
		    	WeightMatrix pwm = pwms.get(m);
		    	motifLengths.add(pwm.length());
		    	double threshold = pwm.getMaxScore() * scoreRatio;
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
   
		    if (toMakeMatrix){	// make a score matrix (seqs x motifs)
		    	double[][] matrix = makePwmScoreMatrix(seqs, pwms);
		    	StringBuilder msb = new StringBuilder();
		    	msb.append("Seq").append("\t");
		    	for (int i=0; i<pwms.size(); i++)
		    		msb.append(pwms.get(i).name).append("\t");
		    	CommonUtils.replaceEnd(msb, '\n');
		    	for (int j=0; j<seqs.length; j++){
//		    		msb.append(names[j]).append("\t");
			    	for (int i=0; i<pwms.size(); i++){
			    		if (matrix[j][i]<motifThresholds.get(i) && scoreRatio>=0)		// cutoff
			    			msb.append(0).append("\t");
			    		else
			    			msb.append(String.format("%.4f\t", matrix[j][i]));
			    	}
			    	CommonUtils.replaceEnd(msb, '\n');
		    	}
		    	CommonUtils.writeFile(out.concat(".scoreMatrix.txt"), msb.toString());
		    	System.exit(0);
		    }
		    else{	// only report strong hits
				instances = getPWMInstances(seqs, pwms, motifThresholds);
				header = "# numSequence:"+seqs.length+"\n# numMotif:"+pwms.size();
		    }
		}
		
	    // search for KSM motif matches
		String ksm_string = Args.parseString(args, "ksm", null);
		if (ksm_string!=null){		// Load multiple KSMs
			
	        Config config = new Config();
	        try{
				config.parseArgs(args);   
			}
			catch (Exception e){
				e.printStackTrace();
	    		System.exit(-1);
			}  
	        
			ArrayList<KMAC> kmacs = new ArrayList<KMAC>();
			ArrayList<String> knames = new ArrayList<String>();
			if (ksm_string.contains(",")){	// it is a comma-separated string with name,path format
				String f[]=ksm_string.split(",");
				KMAC kmac = CommonUtils.loadKsmFile(f[1].trim(), config);
				if (kmac==null){
					System.err.println("Error in loading "+f[1]+", exit here!");
					System.exit(-1);
				}
				knames.add(f[0].trim());
				kmacs.add(kmac);
				instances = scan_multiKSMs(seqs, kmacs, knames);
			}
			else{		// it is KSM list file path
				ArrayList<String> lines = CommonUtils.readTextFile(ksm_string);
				if (lines.size()>500 && !toMakeMatrix){
					// scan each KSM one by one, avoid loading all KSMs (using too much memory)
					instances = new ArrayList<MotifInstance>();
					int n = 0;
					for (String l:lines){
						if (l.startsWith("#"))
							continue;
						String[] f = l.split("\t");
						String kname = f[0].trim();
						KMAC kmac = CommonUtils.loadKsmFile(f[1].trim(), config);
						if (kmac==null){
							System.err.println("Error in loading "+f[1]+", skipping it ...");
							continue;
						}
						instances.addAll(scan_singleKSM(seqs, kmac, kname, n));
					}
				}
				else{ 
					for (String l:lines){
						if (l.startsWith("#"))
							continue;
						String[] f = l.split("\t");
						KMAC kmac = CommonUtils.loadKsmFile(f[1].trim(), config);
						if (kmac==null){
							System.err.println("Error in loading "+f[1]+", skipping it ...\n");
							continue;
						}
						knames.add(f[0].trim());
						kmacs.add(kmac);
					}
				    if (toMakeMatrix){
				    	CommonUtils.writeFile(out.concat(".scoreMatrix.txt"), makeScoreMatrix_multiKSMs(seqs, kmacs, knames));
				    	System.exit(0);
				    }
					instances = scan_multiKSMs(seqs, kmacs, knames);
				}
			}
			
			header = "# numSequence:"+seqs.length+"\n# numMotif:"+kmacs.size();
		}
		
		// search for exact k-mer match
		String kmer = Args.parseString(args, "kmer", null);
		if (kmer!=null){
			StringBuilder sb_header = new StringBuilder();
		    
		    motifLengths.add(kmer.length());
		    
		    sb_header.append("# K-mer Information\n");
		    sb_header.append("#ID\tLetters\tWidth\n");
		    sb_header.append("#").append(0).append("\t").append(kmer).append("\t").append(kmer.length());
		    sb_header.append("\n");
		    sb_header.append("#");
		    System.out.println(sb_header.toString());
		    
			instances = getKmerInstances(seqs, kmer);
			header = "# numSequence:"+seqs.length+"\n# numMotif:1";
		}
		
	    // output
		System.out.println("Note: for motif instances on the minus strand, the SeqPos is the position on the reverse compliment of the input sequence.");
	    if (toAddFasta){
	    	CommonUtils.writeFile(out.concat(".motifInstances.txt"), 
    			header+"\nMotif\tSeqID\tMotif_Name\tSeqName\tMatch\tSeqPos\tCoord\tStrand\tScore\tFasta\n"); 	// write first, overwrite if the file exists
	    }
	    else{
	    	CommonUtils.writeFile(out.concat(".motifInstances.txt"), 
	    			header+"\nMotif\tSeqID\tMotif_Name\tSeqName\tMatch\tSeqPos\tCoord\tStrand\tScore\n"); 	// write first, overwrite if the file exists	    	
	    }
		StringBuilder sb = new StringBuilder();
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
	    	sb.append(mi.motifID).append("\t").append(mi.seqID).append("\t").append(mi.motifName).append("\t")
	    	.append(names[mi.seqID]).append("\t").append(mi.matchSeq).append("\t")
	    	.append(mi.position).append("\t").append(coor_string).append("\t")
	    	.append(mi.strand).append("\t").append(String.format("%.2f", mi.score));
	    	if (toAddFasta){
	    		sb.append("\t").append(seqs[mi.seqID]);
	    	}
	    	sb.append("\n");
	    	if (sb.length()>1e7){	// write sb in smaller trunks
		    	CommonUtils.appendFile(out.concat(".motifInstances.txt"), sb.toString());
		    	sb = new StringBuilder();
	    	}	    	
	    }	    
    	CommonUtils.appendFile(out.concat(".motifInstances.txt"), sb.toString());	// write out the remaining texts
	    
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
	
	/**
	 * Return all KSM motif matches (with 1+ component k-mer match) from all sequences for all KSMs
	 * @param seqs
	 * @param kmacs
	 * @param knames
	 * @return
	 */
	public static ArrayList<MotifInstance> scan_multiKSMs(String[] seqs, ArrayList<KMAC> kmacs, ArrayList<String> knames) {
		System.out.println("Scanning KSM motifs ...");
	    ArrayList<MotifInstance> instances = new ArrayList<MotifInstance>();
	    String[] seqs_rc = new String[seqs.length];
	    for (int i=0;i<seqs.length;i++)
	    	seqs_rc[i]=SequenceUtil.reverseComplement(seqs[i]);
	    
	    for (int m=0; m<kmacs.size(); m++){
	    	System.out.println("  ... "+knames.get(m)+" ...");
	    	KMAC kmac = kmacs.get(m);
	    	for (int s=0; s<seqs.length;s++){
//	    		System.out.print(s+" ");
	    		if (s==60) {
	    			kmac.setIsDebugging(); // debug
	    			System.out.println();
	    		}
	    		KmerGroup[] kgs = kmac.findKsmGroupHits(seqs[s], seqs_rc[s]);
	    		if (kgs==null)
	    			continue;
	    		for (int i=0;i<kgs.length;i++){
	    			KmerGroup kg = kgs[i];
		    		MotifInstance mi = new MotifInstance();
		    		mi.motifID = m;
		    		mi.motifName = knames.get(m);
		    		mi.score = kg.getPredictedValue();
		    		int pos = kg.getPosBS();
		    		if (pos > KMAC.RC-seqs[s].length()*2){	// RC strand match
		    			mi.position = pos-KMAC.RC;	
		    			mi.strand = '-';
		    		}
		    		else{
		    			mi.position = pos;	
		    			mi.strand = '+';
		    		}
		    		if (kg.getKmers().isEmpty())
		    			System.err.println("Empty Kmer group match: "+kg.toString());
		    		mi.matchSeq = kg.getCoveredSequence()+":"+kg.getAllKmerString();
		    		mi.seqID = s;
		    		instances.add(mi);
	    		}
	    	}
	    }
		return instances;
	}

	/**
	 * Return a matrix of all KSM motif scores from all sequences for all KSMs <br>
	 * Max score is used if a KSM has multiple matches 
	 * @param seqs
	 * @param kmacs
	 * @param knames
	 * @return
	 */
	public static String makeScoreMatrix_multiKSMs(String[] seqs, ArrayList<KMAC> kmacs, ArrayList<String> knames) {
		System.out.println("Making KSM motif score matrix ...");
	    StringBuilder sb = new StringBuilder();
	    for (int m=0; m<kmacs.size(); m++){
	    	sb.append(knames.get(m)).append("\t");
	    }
	    CommonUtils.replaceEnd(sb, '\n');
	    System.out.println(sb.toString());
	    
	    for (int s=0; s<seqs.length;s++){
	    	String seqs_rc = SequenceUtil.reverseComplement(seqs[s]);
		    for (int m=0; m<kmacs.size(); m++){
	    		KmerGroup[] kgs = kmacs.get(m).findKsmGroupHits(seqs[s], seqs_rc);
	    		if (kgs==null)
	    			sb.append(0).append("\t");
	    		else
	    			sb.append(String.format("%.4f\t", kgs[0].getScore()));
	    	}
		    CommonUtils.replaceEnd(sb, '\n');
	    }
		return sb.toString();
	}
	
	// scan single KSM
	public static ArrayList<MotifInstance> scan_singleKSM(String[] seqs, KMAC kmac, String kname, int motifId) {
		System.out.println("Scanning KSM motifs ...");
	    ArrayList<MotifInstance> instances = new ArrayList<MotifInstance>();
	    String[] seqs_rc = new String[seqs.length];
	    for (int i=0;i<seqs.length;i++)
	    	seqs_rc[i]=SequenceUtil.reverseComplement(seqs[i]);
	    
    	System.out.println("  ... "+kname+" ...");
    	for (int s=0; s<seqs.length;s++){
    		if (s==60) {
    			kmac.setIsDebugging(); // debug
    			System.out.println();
    		}
    		KmerGroup[] kgs = kmac.findKsmGroupHits(seqs[s], seqs_rc[s]);
    		if (kgs==null)
    			continue;
    		for (int i=0;i<kgs.length;i++){
    			KmerGroup kg = kgs[i];
	    		MotifInstance mi = new MotifInstance();
	    		mi.motifID = motifId;
	    		mi.motifName = kname;
	    		mi.score = kg.getPredictedValue();
	    		int pos = kg.getPosBS();
	    		if (pos > KMAC.RC-seqs[s].length()*2){	// RC strand match
	    			mi.position = pos-KMAC.RC;	
	    			mi.strand = '-';
	    		}
	    		else{
	    			mi.position = pos;	
	    			mi.strand = '+';
	    		}
	    		if (kg.getKmers().isEmpty())
	    			System.err.println("Empty Kmer group match: "+kg.toString());
	    		mi.matchSeq = kg.getCoveredSequence()+":"+kg.getAllKmerString();
	    		mi.seqID = s;
	    		instances.add(mi);
    		}
    	}
		return instances;
	}
	
	public static double[][] makePwmScoreMatrix(String[] seqs, List<WeightMatrix> pwms){
		System.out.println("Making PWM motif score matrix ...");
		double[][] matrix = new double[seqs.length][pwms.size()];
		for (int i=0;i<seqs.length;i++){
			for (int j=0; j<pwms.size();j++)
				matrix[i][j] = WeightMatrixScorer.getMaxSeqScore(pwms.get(j), seqs[i], false);
		}
		return matrix;
	}

	/**
	 * Return a list of motif PWM matches (above threshold) from a list of sequences<br>
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
							mi.strand = strand;
							mi.score = score;
							if (strand=='+')
								mi.position = i+pwm.length()/2;
							else
								// minus strand, the SeqPos is the position on the reverse compliment of the input sequence.
								mi.position = str.length() - (i+(pwm.length()-pwm.length()/2-1)) -1 ; 
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



package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.discovery.kmer.Kmer;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KmerEngine;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KmerEngine.KmerGroup;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KmerMotifFinder;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSPeak;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class KmerScanner {

	private ArrayList<Kmer> kmers;
	private KmerEngine kEngine;
	// each element in the list is for one ChIP-Seq method
	private ArrayList<String> methodNames = new ArrayList<String>();
	private ArrayList<ArrayList<Point>> events = new ArrayList<ArrayList<Point>>();
	
	public KmerScanner(ArrayList<Kmer> kmers, int posSeqCount, int negSeqCount){
		this.kmers = kmers;
		kEngine = new KmerEngine(kmers, null, false);
		kEngine.setTotalSeqCount(posSeqCount, negSeqCount);
	}
	
	public KmerGroup[] query (String seq){
		return kEngine.query(seq);
	}
	public KmerGroup getBestKG (String seq){
		double minScore = 0;
		KmerGroup best = null;
		KmerGroup[] kgs = query(seq);
		for (KmerGroup kg:kgs)
			if (minScore >= kg.getHgp()){
				minScore = kg.getHgp();
				best = kg;
			}
		return best;
	}	
	public static void main(String[] args){
		// k-mer group info
		String outName= Args.parseString(args, "out", "out");
		File file = new File(Args.parseString(args, "kmer", null));
		ArrayList<Kmer> kmers = Kmer.loadKmers(file);
		int clusterId = Args.parseInteger(args, "cluster", 0);
		ArrayList<Kmer> toRemove = new ArrayList<Kmer>();
		for (Kmer km:kmers)
			if (km.getClusterId()!=clusterId)
				toRemove.add(km);
		kmers.removeAll(toRemove);
		kmers.trimToSize();
		Pair<Integer, Integer> c = Kmer.getTotalCounts(file);
		KmerScanner scanner = new KmerScanner(kmers, c.car(), c.cdr());
		
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
	    
	    // PWM
		Pair<WeightMatrix, Double> wm = CommonUtils.loadPWM(args, org.getDBID());
		WeightMatrix motif = wm.car();
		double motifThreshold = wm.cdr();
	    	    
		// event locations
		int windowSize = Args.parseInteger(args, "win", 50);
		String eventFile = Args.parseString(args, "event", null);
		List<GPSPeak> gpsPeaks = null;
		try{
			gpsPeaks = GPSParser.parseGPSOutput(eventFile, genome);
		}
		catch(IOException e){
			e.printStackTrace();
		}
		SequenceGenerator<Region> seqgen = new SequenceGenerator<Region>();
		seqgen.useCache(true);		
		
		StringBuilder sb = new StringBuilder();
		ArrayList<Region> posRegions = new ArrayList<Region>();
		for(int i=0;i<gpsPeaks.size();i++){
			GPSPeak p = gpsPeaks.get(i);
			posRegions.add(p.expand(windowSize));
		}
		int width = windowSize*2+1;
		TreeMap<Region, Region> reg2reg = new TreeMap<Region, Region>();
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
			}
		}
		
		ArrayList<Double> pwm_scores = new ArrayList<Double>();
		ArrayList<Double> pwmN_scores = new ArrayList<Double>();
		ArrayList<Double> ksm_scores = new ArrayList<Double>();
		ArrayList<Double> ksmN_scores = new ArrayList<Double>();
		for (Region r:reg2reg.keySet()){
			String seq = seqgen.execute(r).toUpperCase();
			Region rN = reg2reg.get(r);
			String seqN = seqgen.execute(rN).toUpperCase();
			double pwm = WeightMatrixScorer.getMaxSeqScore(motif, seq);
			pwm_scores.add(pwm);
			double pwmN = WeightMatrixScorer.getMaxSeqScore(motif, seqN);
			pwmN_scores.add(pwmN);
			KmerGroup kg = scanner.getBestKG(seq);
			KmerGroup kgN = scanner.getBestKG(seqN);
			ksm_scores.add(kg==null?0:-kg.getHgp());
			ksmN_scores.add(kgN==null?0:-kgN.getHgp());
			sb.append(String.format("%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\n", r.toString(), rN.toString(), pwm, pwmN, kg==null?0:-kg.getHgp(), kgN==null?0:-kgN.getHgp()));
		}
		
		CommonUtils.writeFile(outName+"_w"+width+"_scores.txt", sb.toString());
		System.out.println(outName+"_w"+width+"_scores.txt");
		

		ArrayList<ScoreEnrichment> pwm_se = computeScoreEnrichments(pwm_scores, pwmN_scores);
		sb = new StringBuilder();
		for (ScoreEnrichment se: pwm_se){
			sb.append(String.format("%.2f\t%d\t%d\t%.2f\n", se.score, se.posHit, se.negHit, se.hgp));
		}
		CommonUtils.writeFile(outName+"_w"+width+"_pwm_enrichment.txt", sb.toString());
		
		ArrayList<ScoreEnrichment> ksm_se = computeScoreEnrichments(ksm_scores, ksmN_scores);
		sb = new StringBuilder();
		for (ScoreEnrichment se: ksm_se){
			sb.append(String.format("%.2f\t%d\t%d\t%.2f\n", se.score, se.posHit, se.negHit, se.hgp));
		}
		CommonUtils.writeFile(outName+"_w"+width+"_ksm_enrichment.txt", sb.toString());
		
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
			se.hgp = KmerMotifFinder.computeHGP(total, total, se.posHit, se.negHit);
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

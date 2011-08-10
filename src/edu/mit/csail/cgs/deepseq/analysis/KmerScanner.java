package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSPeak;
import edu.mit.csail.cgs.ewok.verbs.motifs.Kmer;
import edu.mit.csail.cgs.ewok.verbs.motifs.KmerEngine;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.ewok.verbs.motifs.KmerEngine.KmerGroup;
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
	public double getMaxKGscore (String seq){
		double minScore = 0;
		KmerGroup[] kgs = query(seq);
		for (KmerGroup kg:kgs)
			if (minScore >= kg.getHgp())
				minScore = kg.getHgp();
		return minScore;
	}	
	public static void main(String[] args){
		// k-mer group info
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
		for(int i=0;i<gpsPeaks.size();i++){
			GPSPeak p = gpsPeaks.get(i);
			String seq = seqgen.execute(p.expand(windowSize));
			Point pt = new Point(genome, p.getChrom(), p.getLocation()+500);
			String seqN = seqgen.execute(pt.expand(windowSize));
			double pwm = WeightMatrixScorer.getMaxSeqScore(motif, seq);
			double pwmN = WeightMatrixScorer.getMaxSeqScore(motif, seqN);
			double kgs = scanner.getMaxKGscore(seq);
			double kgsN = scanner.getMaxKGscore(seqN);
			sb.append(String.format("%d\t%.2f\t%.2f\t%.2f\t%.2f\n", i, pwm, pwmN, kgs, kgsN));
		}
		System.out.println(sb.toString());
		CommonUtils.writeFile("Ctcf_scores.txt", sb.toString());
		
//		KmerGroup[] kgs = scanner.query("TATTTACATGCAGTGTCCGGAAGACGCCAGAAGAGGGCAGTAGATGCCCTAGTAGTGGAGC");
//		for (KmerGroup kg:kgs){
//			System.out.println(kg.toString());
//		}
	}

}

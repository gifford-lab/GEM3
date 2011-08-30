package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.discovery.kmer.Kmer;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KmerEngine;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KmerEngine.KmerGroup;
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
		
		for (Region r:reg2reg.keySet()){
			String seq = seqgen.execute(r);
			Region rN = reg2reg.get(r);
			String seqN = seqgen.execute(rN);
			double pwm = WeightMatrixScorer.getMaxSeqScore(motif, seq);
			double pwmN = WeightMatrixScorer.getMaxSeqScore(motif, seqN);
			KmerGroup kg = scanner.getBestKG(seq);
			KmerGroup kgN = scanner.getBestKG(seqN);
			sb.append(String.format("%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\n", r.toString(), rN.toString(), pwm, pwmN, kg==null?0:kg.getHgp(), kgN==null?0:kgN.getHgp()));
		}
		CommonUtils.writeFile(outName+"_w"+width+"_scores.txt", sb.toString());
		System.out.println(outName+"_w"+width+"_scores.txt");
		
//		KmerGroup[] kgs = scanner.query("TATTTACATGCAGTGTCCGGAAGACGCCAGAAGAGGGCAGTAGATGCCCTAGTAGTGGAGC");
//		for (KmerGroup kg:kgs){
//			System.out.println(kg.toString());
//		}
	}

}

package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
import java.util.ArrayList;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.ewok.verbs.motifs.Kmer;
import edu.mit.csail.cgs.ewok.verbs.motifs.KmerEngine;
import edu.mit.csail.cgs.ewok.verbs.motifs.KmerEngine.KmerGroup;

public class KmerScanner {
	private Genome genome;
	private Organism org;
	private String[] args;
	private WeightMatrix motif = null;
	private String outName="out";
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
	
	public static void main(String[] args){
		ArrayList<Kmer> kmers = Kmer.loadKmers(new File(args[0]));
		int clusterId = Integer.parseInt(args[1]);
		ArrayList<Kmer> toRemove = new ArrayList<Kmer>();
		for (Kmer km:kmers)
			if (km.getClusterId()!=clusterId)
				toRemove.add(km);
		kmers.removeAll(toRemove);
		kmers.trimToSize();
		KmerScanner scanner = new KmerScanner(kmers, 39070, 39750);
		KmerGroup[] kgs = scanner.query("CACGAAGCACACGCCCGAAAATCCTGAGCACGTGGCTCTACCGAGGGACTGGAAGCGCTCC");
		for (KmerGroup kg:kgs){
			System.out.println(kg.toString());
		}
	}

}

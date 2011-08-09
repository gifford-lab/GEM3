package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
import java.util.ArrayList;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.ewok.verbs.motifs.Kmer;
import edu.mit.csail.cgs.ewok.verbs.motifs.KmerEngine;

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
	
	public KmerScanner(ArrayList<Kmer> kmers){
		this.kmers = kmers;
		kEngine = new KmerEngine(kmers, null, false);
		
	}
	public static void main(String[] args){
		KmerScanner scanner = new KmerScanner(Kmer.loadKmers(new File(args[0])));
	}

}

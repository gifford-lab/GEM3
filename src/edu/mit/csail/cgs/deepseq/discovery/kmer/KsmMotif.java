package edu.mit.csail.cgs.deepseq.discovery.kmer;

import java.util.ArrayList;

public class KsmMotif {
	public ArrayList<Kmer> kmers;
	public int posSeqCount;
	public int negSeqCount;
	public double score;
	public double[] seq_weights;
}

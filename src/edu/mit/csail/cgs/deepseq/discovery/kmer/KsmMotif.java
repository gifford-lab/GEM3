package edu.mit.csail.cgs.deepseq.discovery.kmer;

import java.util.ArrayList;

public class KsmMotif {
	public ArrayList<Kmer> kmers;
	public int posSeqCount;
	public int negSeqCount;
	public double cutoff;
	public double[] seq_weights;
	public int[] posCoveredWidth;
	public int[] negCoveredWidth;
	public String[][] posHitStrings;		// KSM matched seq
	public String[][] negHitStrings;
}

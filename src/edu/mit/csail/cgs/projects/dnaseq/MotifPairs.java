package edu.mit.csail.cgs.projects.dnaseq;

import java.util.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.motifs.*;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.tools.motifs.WeightMatrixScanner;
import edu.mit.csail.cgs.tools.utils.Args;

/**
 * Scans genomic sequence for a collection of motifs and 
 * reports on the frequency with which one motif follows another.
 *
 * This is important to the DNAseq HMM as a prior on the transition
 * probabilities from one motif state to another
 */


public class MotifPairs {

    private Genome genome;
    private List<WeightMatrix> matrices;
    private short[][] follows;
    private int distanceLimit = 5;
    private List<Region> regionsToScan;

    public MotifPairs() {}
    public void parseArgs(String args[]) throws NotFoundException {
        distanceLimit = Args.parseInteger(args,"distance",distanceLimit);
        matrices = new ArrayList<WeightMatrix>();
        matrices.addAll(Args.parseWeightMatrices(args));
        genome = Args.parseGenome(args).cdr();
        follows = new short[matrices.size()][matrices.size()];
        regionsToScan = Args.parseRegionsOrDefault(args);
    }
    public void run() {
        int window = 10000;
        SequenceGenerator seqgen = new SequenceGenerator();
        List<List<WMHit>> hits = new ArrayList<List<WMHit>>(matrices.size());
        float[] cutoffs = new float[matrices.size()];
        for (int i = 0; i < matrices.size(); i++) {
            cutoffs[i] = (float)(matrices.get(i).getMaxScore() * .6);
        }

        seqgen.useCache(true);
        seqgen.useLocalFiles(true);
        for (Region chrom : regionsToScan) {
            int start = chrom.getStart();
            while (start + window < chrom.getEnd()) {
                Region r = new Region(chrom.getGenome(), chrom.getChrom(), start, start + window);
                char[] seq = seqgen.execute(r).toCharArray();
                for (int i = 0; i < matrices.size(); i++) {
                    hits.set(i, WeightMatrixScanner.scanSequence(matrices.get(i),
                                                                 cutoffs[i],
                                                                 seq));
                }
                for (int i = 0; i < matrices.size(); i++) {
                    for (int j = 0; j < matrices.size(); j++) {
                        if (i == j) { continue;}

                    }
                }


                start += window;
            }
            seqgen.clearCache();
        }

    }



}
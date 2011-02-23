package edu.mit.csail.cgs.projects.dnaseq;

import java.io.*;
import java.util.*;
import java.sql.SQLException;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.motifs.*;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.tools.motifs.WeightMatrixScanner;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.projects.readdb.Client;
import edu.mit.csail.cgs.projects.readdb.ClientException;
import edu.mit.csail.cgs.projects.readdb.Aggregator;

/**
 * Test the HMM by calling states on unlabeled data.
 *
 * java edu.mit.csail.cgs.projects.dnaseq.HMMTest --species "$HS;hg19" [--modelfile hmm.model] --dnaseq "Crawford GM12878 DNaseSeq GM12878 against Input;statistical 1/11/11" --region "7:27145469-27218813"
 *
 * [--modelfile hmm.model] model file to load
 *
 * The analysis specified by --dnaseq isn't actually used; the
 * code just uses the same input alignments.
 */

public class HMMTest {
    private ChipSeqLoader loader;
    private HMMReads reads;
    private Genome genome;
    private SequenceGenerator seqgen;
    private List<Region> testRegions;
    private List<ChipSeqAlignment> alignments;
    private String modelFname;

    private HMM hmm;

    public HMMTest() throws IOException, ClientException, SQLException {
        reads = new HMMReads();
        loader = new ChipSeqLoader();
    }
    public void parseArgs(String args[]) throws NotFoundException, SQLException, IOException {
        modelFname = Args.parseString(args,"modelfile","hmm.model");
        genome = Args.parseGenome(args).cdr();
        testRegions = Args.parseRegions(args);
        ChipSeqAnalysis dnaseq = Args.parseChipSeqAnalysis(args,"dnaseq");
        alignments = new ArrayList<ChipSeqAlignment>();
        alignments.addAll(dnaseq.getForeground());
        seqgen = new SequenceGenerator(genome);
        seqgen.useLocalFiles(true);
        seqgen.useCache(true);
        hmm = new HMM(modelFname);
        hmm.toLogProbabilities();
    }
    private Collection<Region> getTestRegions() throws SQLException {
        return testRegions;
    }

    public void test() throws IOException, ClientException, SQLException {
        int motiflength = (hmm.numStates - 2) / 2;
        for (Region region : getTestRegions()) {
            char[] sequence = seqgen.execute(region).toCharArray();
            ReadCounts counts = reads.getReadCounts(region,
                                                    alignments);
            int readCounts[] = counts.getCounts();

            System.err.println("Running Forward-Backward on " + region);
            byte backPointers[][] = new byte[sequence.length][hmm.numStates];
            double lastProb[][] = new double[sequence.length][hmm.numStates];
            for (int i = 0; i < hmm.numStates; i++) {
                lastProb[0][i] = hmm.initialProbabilities[i] + Math.log(hmm.states[i].getProb(sequence[0], readCounts[0]));
            }
            for (int t = 1; t < sequence.length; t++) {
                for (int s = 0; s < hmm.numStates; s++) {
                    double dataprob = Math.log(hmm.states[s].getProb(sequence[t], readCounts[t]) + .000001);
                    if (Double.isNaN(dataprob)) {
                        throw new RuntimeException(String.format("dp is nan in %d from %c, %d",
                                                                 s, sequence[t], readCounts[t]));
                    }

                    byte bestPrevState = -1;
                    double bestPrevProb = Double.NEGATIVE_INFINITY;
                    for (byte y = 0; y < hmm.numStates; y++) {
                        double trans = hmm.transitions[y][s];
                        if (Double.isInfinite(trans)) { continue;}
                        double p = trans + lastProb[t-1][y];
                        if (p > bestPrevProb) {
                            bestPrevProb = p;
                            bestPrevState = y;
                        }
                    }                    
                    lastProb[t][s] = dataprob + bestPrevProb;
                    if (Double.isNaN(lastProb[t][s])) {
                        throw new RuntimeException(String.format("Setting NaN at %d, %d from %f, %f",
                                                                 t,s,dataprob, bestPrevProb));
                    }

                    backPointers[t][s] = bestPrevState;
                }
            }
            byte state = -1;
            double bestprob = Double.NEGATIVE_INFINITY;
            for (byte s = 0; s < hmm.numStates; s++) {
                if (lastProb[lastProb.length - 1][s] > bestprob) {
                    state = s;
                    bestprob = lastProb[lastProb.length - 1][s];
                }
            }
            if (state == -1) {
                System.err.println("Couldn't find a best ending state.  Going to fill your screen now");
                for (int t = 0; t < sequence.length && t < 100; t++) {
                    System.err.print(String.format("%d %c  ",t, sequence[t]));
                    for (int s = 0; s < hmm.numStates; s++) {
                        System.err.print(String.format("  %.2e",lastProb[t][s]));
                    }
                    System.err.println();
                }
            }

            byte bestStateSequence[] = new byte[sequence.length];
            bestStateSequence[bestStateSequence.length - 1] = state;
            int t = sequence.length - 1;
            state = backPointers[t][state];
            t--;
            while (t >= 0) {
                bestStateSequence[t] = state;
                state = backPointers[t][state];
                t--;
            }
            System.err.println();
            for (int i = 0; i < bestStateSequence.length; i++) {
                if (bestStateSequence[i] >= 2 && (i == 0 || bestStateSequence[i-1] < 2)) {
                    System.out.println(String.format("%s:%d-%d",
                                                     region.getChrom(),
                                                     (i + region.getStart()),
                                                     (i + region.getStart() + motiflength)));
                }
            }
        }
    }
    public static void main(String args[]) throws Exception {
        HMMTest test = new HMMTest();
        test.parseArgs(args);
        test.test();
    }



}
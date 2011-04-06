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
 * Run an HMM with multiple motifs to determine binding sites in DNase 
 * hypersensitive regions.  Needs an HMM model from one of the training programs
 * as input to get the nucleotide and read counts for the generic insensitive
 * and sensitive states.  The motif states here use the nucleotide counts from
 * the database and the read counts from the insensitive state.
 *
 * The HMM model allows transitions between sensitive and insensitive.  Sensitive
 * can also transition into a motif, which then runs straight through the motif
 * and back to sensitive.
 *
 * java edu.mit.csail.cgs.projects.dnaseq.MultiMotifHMMTest --species "$HS;hg19" --modelfile chipseq.model --wm "CTCF;JASPAR 11/09 MA0139.1" --wm "SRF;TRANSFAC 10.4 V$SRF_Q4 M00810" --dnaseq "Crawford GM12878 DNaseSeq GM12878 against Input;statistical 1/11/11" --region "7:27m-40m"
 *
 *
 */

public class MultiMotifHMMTest {
    // filled in by constructo
    private ChipSeqLoader loader;
    private HMMReads reads;

    // filled in by parseArgs
    private Genome genome;
    private SequenceGenerator seqgen;
    private List<Region> testRegions;
    private List<ChipSeqAlignment> alignments;
    private String modelFname;
    private List<WeightMatrix> matrices;
    private List<HMM> trainedSingleMotifHMMs;
    private List<String> trainedNames;
    private boolean debug;

    // filled in by setupHMM
    private short numStates;
    private double transitions[][];
    private HMMState states[];
    private double initialProbabilities[];
    private short possiblePreviousStates[][];
    private short firstMotifState[]; // for weight matrices
    private short firstModelState[]; // for models trained and saved as an HMM
    private String stateNames[];

    public MultiMotifHMMTest() throws IOException, ClientException, SQLException {
        loader = new ChipSeqLoader();
        reads = new HMMReads();
    }
    public void parseArgs(String args[]) throws NotFoundException, SQLException, IOException {
        modelFname = Args.parseString(args,"modelfile","hmm.model");
        genome = Args.parseGenome(args).cdr();
        testRegions = Args.parseRegions(args);
        reads.smooth(Args.parseInteger(args,"smooth",0));
        ChipSeqAnalysis dnaseq = Args.parseChipSeqAnalysis(args,"dnaseq");
        alignments = new ArrayList<ChipSeqAlignment>();
        alignments.addAll(dnaseq.getForeground());
        seqgen = new SequenceGenerator(genome);
        seqgen.useLocalFiles(true);
        seqgen.useCache(true);
        debug = Args.parseFlags(args).contains("debug");
        matrices = new ArrayList<WeightMatrix>();
        Collection<WeightMatrix> twm = Args.parseWeightMatrices(args);
        if (twm.size() > 20) {
            System.err.println("Too many weight matrices returned.  Not using them.");
        } else {
            matrices.addAll(twm);
            for (WeightMatrix m : twm) {
                matrices.add(WeightMatrix.reverseComplement(m));
            }            
        }

        MarkovBackgroundModel bgModel = null;
        String bgmodelname = Args.parseString(args,"bgmodel","whole genome zero order");
        BackgroundModelMetadata md = BackgroundModelLoader.getBackgroundModel(bgmodelname,
                                                                              1,
                                                                              "MARKOV",
                                                                              genome.getDBID());
        for (WeightMatrix motif : matrices) {
            if (bgModel == null) {
                motif.toFrequency();
            } else {
                motif.toFrequency(bgModel);
            }
        }
        trainedSingleMotifHMMs = new ArrayList<HMM>();
        trainedNames = new ArrayList<String>();
        trainedNames.addAll(Args.parseStrings(args,"model"));
        for (String f : trainedNames) {
            trainedSingleMotifHMMs.add(new HMM(f));
        }

        System.err.println("There are " +matrices.size() + " motifs");
        System.err.println(matrices.toString());
        System.err.println("There are " + trainedNames.size() + " trained models");
    }
    private Collection<Region> getTestRegions() throws SQLException {
        return testRegions;
    }
    public void setupHMM() throws IOException {
        HMM baseModel = new HMM(modelFname);
        HMMState insensitive = baseModel.states[0];
        HMMState sensitive = baseModel.states[1];
        
        numStates = 2;
        firstMotifState = new short[matrices.size()];
        firstModelState = new short[trainedSingleMotifHMMs.size() * 2];
        for (int i = 0; i < matrices.size(); i++) {
            firstMotifState[i] = numStates;
            numStates += matrices.get(i).length();
        }
        for (int i = 0; i < trainedSingleMotifHMMs.size(); i++) {
            firstModelState[i*2] = numStates;
            numStates += (trainedSingleMotifHMMs.get(i).numStates - 2) / 2;
            firstModelState[i*2+1] = numStates;
            numStates += (trainedSingleMotifHMMs.get(i).numStates - 2) / 2;
        }

        states = new HMMState[numStates];
        initialProbabilities = new double[numStates];
        transitions = new double[numStates][numStates];
        possiblePreviousStates = new short[numStates][];
        stateNames = new String[numStates];

        initialProbabilities[0] = 1;

        states[0] = insensitive;
        states[1] = sensitive;
        possiblePreviousStates[0] = new short[2];
        possiblePreviousStates[0][0] = 0;
        possiblePreviousStates[0][1] = 1;
        transitions[0][0] = .99999;
        transitions[0][1] = .00001;
        stateNames[0] = "unenriched";

        possiblePreviousStates[1] = new short[2 + matrices.size() + 2*trainedSingleMotifHMMs.size()];
        possiblePreviousStates[1][0] = 0;
        possiblePreviousStates[1][1] = 1;
        transitions[1][0] = .005;
        transitions[1][1] = .990;
        for (int i = 0; i < matrices.size(); i++) {
            possiblePreviousStates[1][2+i] = (short)(firstMotifState[i] + matrices.get(i).length() - 1);
            transitions[1][firstMotifState[i]] = .005 / (double)(numStates - 2);
        }
        for (int i = 0; i < trainedSingleMotifHMMs.size(); i++) {
            possiblePreviousStates[1][2+matrices.size()+i*2] = (short)(firstModelState[i*2] + (trainedSingleMotifHMMs.get(i).numStates-2)/2 - 1);
            possiblePreviousStates[1][2+matrices.size()+i*2+1] = (short)(firstModelState[i*2+1] + (trainedSingleMotifHMMs.get(i).numStates-2)/2 - 1);
            transitions[1][firstModelState[i*2]] = .005 / (double)(numStates - 2);
            transitions[1][firstModelState[i*2+1]] = .005 / (double)(numStates - 2);
        }

        stateNames[1] = "enriched";

        int[] sensReadCounts = sensitive.getCounts();
        int[] insensReadCounts = insensitive.getCounts();
        int footprintPrior = 1000;
        
        for (int i = 0; i < matrices.size(); i++) {
            WeightMatrix matrix = matrices.get(i);
            short fms = firstMotifState[i];
            for (short j = 0; j < matrix.length(); j++) {
                int[] footprintReadCounts = new int[insensReadCounts.length];
                // this forumula comes from ~/psrg/projects/dnaseq/ctcf_test/predict_counts_from_bases.m
                // and was based on CTCF and Pax5 data in GM12878
                double avgCount = .41 - .059 * matrix.matrix[j]['A'] + .18 * matrix.matrix[j]['C'] + .27 * matrix.matrix[j]['G'] + .015 * matrix.matrix[j]['T'];
                double p = 1.0 / (1.0 + avgCount);
                for (int k = 0; k < footprintReadCounts.length; k++) {
                    footprintReadCounts[k] = (int)(footprintPrior * p * Math.pow(1-p, k));
                }

                states[fms+j] = new HMMState((int)(matrix.matrix[j]['A'] * footprintPrior),
                                             (int)(matrix.matrix[j]['C'] * footprintPrior),
                                             (int)(matrix.matrix[j]['G'] * footprintPrior),
                                             (int)(matrix.matrix[j]['T'] * footprintPrior),
                                             footprintReadCounts);                
                possiblePreviousStates[fms+j] = new short[1];
                if (j == 0) {
                    possiblePreviousStates[fms+j][0] = 1;
                    transitions[fms+j][fms+1] = 1;
                } else if (j == matrix.length() - 1) {
                    possiblePreviousStates[fms+j][0] = (short)(fms+j-1);
                    transitions[fms+j][1] = 1;  
                } else {
                    possiblePreviousStates[fms+j][0] = (short)(fms + j - 1);                    
                    transitions[fms+j][fms+j+1] = 1;
                }
                stateNames[fms+j] = matrices.get(i).toString() + " pos " + j;
            }
        }
        for (int i = 0; i < trainedSingleMotifHMMs.size(); i++) {
            HMM hmm = trainedSingleMotifHMMs.get(i);
            short fms = firstModelState[i*2];
            int len = (hmm.numStates-2) / 2;
            for (int j = 0; j < len; j++) {
                states[fms+j] = hmm.states[2+j];
                possiblePreviousStates[fms+j] = new short[1];
                if (j == 0) {
                    possiblePreviousStates[fms+j][0] = 1;
                    transitions[fms+j][fms+j+1] = 1;
                } else if (j == len - 1) {
                    possiblePreviousStates[fms+j][0] = (short)(fms+j-1);
                    transitions[fms+j][1] = 1;

                } else {
                    possiblePreviousStates[fms+j][0] = (short)(fms+j-1);
                    transitions[fms+j][fms+j+1] = 1;
                }
                stateNames[fms+j] = trainedNames.get(i) + " pos " + j;
            }
            fms = firstModelState[i*2+1];
            for (int j = 0; j < len; j++) {
                states[fms+j] = hmm.states[2+len+j];
                possiblePreviousStates[fms+j] = new short[1];
                if (j == 0) {
                    possiblePreviousStates[fms+j][0] = 1;
                    transitions[fms+j][fms+j+1] = 1;
                } else if (j == len - 1) {
                    possiblePreviousStates[fms+j][0] = (short)(fms+j-1);
                    transitions[fms+j][1] = 1;

                } else {
                    possiblePreviousStates[fms+j][0] = (short)(fms+j-1);
                    transitions[fms+j][fms+j+1] = 1;
                }
                stateNames[fms+j] = trainedNames.get(i) + " pos " + j;
            }
        }

        for (int i = 0; i < transitions.length; i++) {
            for (int j = 0; j < transitions.length; j++) {
                System.err.print(transitions[i][j] + "  ");
                transitions[i][j] = Math.log(transitions[i][j]);
            }
            System.err.println();
        }
        for (int i = 0; i < numStates; i++) {
            System.err.println(states[i].toString());
        }
    }
    public void test() throws IOException, ClientException, SQLException {
        for (Region region : getTestRegions()) {
            char[] sequence = seqgen.execute(region).toCharArray();
            ReadCounts counts = reads.getReadCounts(region,
                                                    alignments);
            int readCounts[] = counts.getCounts();

            System.err.println("Running Forward-Backward on " + region);
            short backPointers[][] = new short[sequence.length][numStates];
            double lastProb[][] = new double[sequence.length][numStates];
            System.err.println(String.format("%d  %d %d %d",
                                             numStates,
                                             lastProb[0].length,
                                             initialProbabilities.length,
                                             states.length));
            for (int i = 0; i < numStates; i++) {
                lastProb[0][i] = Math.log(initialProbabilities[i] * states[i].getProb(sequence[0], readCounts[0]));
            }
            for (int t = 1; t < sequence.length; t++) {
                if (debug) {
                    System.err.println(String.format("%d  %c  %d: ",t + region.getStart(), sequence[t],readCounts[t]));
                }
                for (short s = 0; s < states.length; s++) {
                    double dataprob = Math.log(states[s].getProb(sequence[t], readCounts[t]) + .01);
                    if (Double.isNaN(dataprob)) {
                        throw new RuntimeException(String.format("dp is nan in %d from %c, %d",
                                                                 s, sequence[t], readCounts[t]));
                    }

                    short bestPrevState = -1;
                    double bestPrevProb = Double.NEGATIVE_INFINITY;
                    for (int k = 0; k < possiblePreviousStates[s].length; k++) {
                        short y = possiblePreviousStates[s][k];
                        double trans = transitions[y][s];
                        if (Double.isInfinite(trans)) { continue;}
                        double p = trans + lastProb[t-1][y];
                        if (debug) {
                            System.err.println(String.format("   %d -> %d,  %.2e,  %.2e + %.2e = %.2e",
                                                             y,s,dataprob,trans,lastProb[t-1][y], p+dataprob));
                        }
                        if (p > bestPrevProb) {
                            bestPrevProb = p;
                            bestPrevState = y;
                        }
                    }                    
                    if (Double.isNaN(lastProb[t][s])) {
                        throw new RuntimeException(String.format("Setting NaN at %d, %d from %f, %f",
                                                                 t,s,dataprob, bestPrevProb));
                    }
                    lastProb[t][s] = dataprob + bestPrevProb;
                    backPointers[t][s] = bestPrevState;
                }
            }
            int state = -1;
            double bestprob = Double.NEGATIVE_INFINITY;
            for (int s = 0; s < states.length; s++) {
                if (lastProb[lastProb.length - 1][s] > bestprob) {
                    state = s;
                    bestprob = lastProb[lastProb.length - 1][s];
                }
            }
            if (state == -1) {
                System.err.println("Couldn't find a best ending state.  Going to fill your screen now");
                for (int t = 0; t < sequence.length && t < 100; t++) {
                    System.err.print(String.format("%d %c  ",t, sequence[t]));
                    for (int s = 0; s < states.length; s++) {
                        System.err.print(String.format("  %.2e",lastProb[t][s]));
                    }
                    System.err.println();
                }
            }

            int bestStateSequence[] = new int[sequence.length];
            bestStateSequence[bestStateSequence.length - 1] = state;
            int t = sequence.length - 1;
            state = backPointers[t][state];
            t--;
            while (t >= 0) {
                bestStateSequence[t] = state;
                state = backPointers[t][state];
                t--;
            }
            for (int i = 1; i < bestStateSequence.length; i++) {
                if (bestStateSequence[i] >= 2 && (i == 0 || bestStateSequence[i-1] < 2)) {
                    int motif = -1;
                    for (int j = 0; j < firstMotifState.length; j++) {
                        if (firstMotifState[j] == bestStateSequence[i]) {
                            motif = j;
                        }
                    }
                    System.out.println(String.format("%s:%d-%d\t%s",
                                                     region.getChrom(),
                                                     (i + region.getStart()),
                                                     (i + region.getStart()),
                                                     stateNames[bestStateSequence[i]]));
                }
                // if (bestStateSequence[i] == 1 && bestStateSequence[i-1] == 0) {
                //     System.err.println("Transition to sensitive at " +(i + region.getStart()));
                // }
                // if (bestStateSequence[i-1] == 1 && bestStateSequence[i] == 0) {
                //     System.err.println("Transition from sensitive at " +(i + region.getStart()));
                // }

            }
        }

    }
    public static void main(String args[]) throws Exception {
        MultiMotifHMMTest test = new MultiMotifHMMTest();
        test.parseArgs(args);
        test.setupHMM();
        test.test();
    }

}
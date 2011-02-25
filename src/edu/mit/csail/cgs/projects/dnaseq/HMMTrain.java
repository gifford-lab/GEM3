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
import edu.mit.csail.cgs.projects.readdb.ClientException;
import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;

/**
 * Superclass for HMM training.  Right now it just needs an implementation of
 * isBound() to determine whether a motif instance is bound.
 */

public abstract class HMMTrain extends DNASeqEnrichmentCaller {

    protected WeightMatrix motif;
    protected SequenceGenerator seqgen;
    private float motifCutoff;
    private WMHitStartComparator hitcomp;
    private List<Region> trainingRegions;
    private String modelFname;
    
    protected int numStates;
    protected int transitions[][];
    protected HMMState insensState, sensState, motifForwStates[], motifRevStates[];
    protected int lastStateNum = -1;

    public HMMTrain() throws IOException, ClientException, SQLException {
        super();
        hitcomp = new WMHitStartComparator();
    }
    public void parseArgs(String args[]) throws NotFoundException, SQLException, IOException {
        super.parseArgs(args);
        modelFname = Args.parseString(args,"modelfile","hmm.model");
        trainingRegions = Args.parseRegions(args);
        Collection<WeightMatrix> matrices = Args.parseWeightMatrices(args);
        Iterator<WeightMatrix> iter = matrices.iterator();
        if (iter.hasNext()) {
            motif = iter.next();
        } else {
            throw new NotFoundException("Couldn't find any motifs in the args");
        }
        MarkovBackgroundModel bgModel = null;
        String bgmodelname = Args.parseString(args,"bgmodel","whole genome zero order");
        BackgroundModelMetadata md = BackgroundModelLoader.getBackgroundModel(bgmodelname,
                                                                              1,
                                                                              "MARKOV",
                                                                              genome.getDBID());
        if (bgModel == null) {
            motif.toLogOdds();
        } else {
            motif.toLogOdds(bgModel);
        }


        if (iter.hasNext()) {
            System.err.println("More than one motif specified in args.  Using the first; "+ motif.toString());
        }        
        motifCutoff = (float)(motif.getMaxScore() * Args.parseDouble(args,"cutoff",.7));
        System.err.println("Motif Cutoff is " + motifCutoff);
        
        seqgen = new SequenceGenerator(genome);
        seqgen.useLocalFiles(true);
        seqgen.useCache(true);

        numStates = 2 + 2 * motif.length();
        transitions = new int[numStates][numStates];
        for (int i = 0; i < numStates; i++) {
            for (int j = 0; j < numStates; j++) {
                transitions[i][j] = 0;
            }
        }

        insensState = new HMMState();
        sensState = new HMMState();
        motifForwStates = new HMMState[motif.length()];
        motifRevStates = new HMMState[motif.length()];
        for (int i = 0; i < motifForwStates.length; i++) {
            motifForwStates[i] = new HMMState();
            motifRevStates[i] = new HMMState();
        }
    }
    private Collection<Region> getTrainingRegions() throws SQLException {
        return trainingRegions;
    }
    protected void newTrainingRegion(Region r) {
    }
    protected abstract boolean hasBinding(int position);
    public void train() throws IOException, SQLException, ClientException {
        int trainedon = 0;
        for (Region region : getTrainingRegions()) {
            newTrainingRegion(region);
            System.err.println("Trainig on " + region);
            ReadCounts readcounts = reads.getReadCounts(region, getAlignments());
            List<? extends Region> dnaseqevents = getHyperSensitiveRegions(region, readcounts);
            int start = region.getStart();
            for (Region sensitiveRegion : dnaseqevents) {
                trainedon++;
                Region insensitiveRegion = new Region(genome, sensitiveRegion.getChrom(), start, sensitiveRegion.getStart());
                addTraining(false, insensitiveRegion, readcounts);
                addTraining(true, sensitiveRegion, readcounts);
                start = sensitiveRegion.getEnd();
            }
            addTraining(false, new Region(genome, region.getChrom(), start, region.getEnd()), readcounts);
        }
        System.err.println("Trained on " + trainedon + " hypersensitive regions");
    }
    private List<ChipSeqAnalysisResult> filterSensitiveRegions(List<ChipSeqAnalysisResult> input) {
        ArrayList<ChipSeqAnalysisResult> output = new ArrayList<ChipSeqAnalysisResult>();
        for (ChipSeqAnalysisResult r : input) {
            if (r.foregroundReadCount > 200 && Math.log(r.pvalue) < - 20) {
                output.add(r);
            }
        }
        return output;
    }
    private void addTraining(boolean sensitive,
                             Region region,
                             ReadCounts counts) throws SQLException {
        char[] sequence = seqgen.execute(region).toCharArray();
        List<WMHit> hits = WeightMatrixScanner.scanSequence(motif, motifCutoff, sequence);
        Collections.sort(hits, hitcomp);
        int motifHitStarts[] = new int[hits.size()];
        for (int i = 0; i < hits.size(); i++) {
            motifHitStarts[i] = hits.get(i).getStart() + region.getStart();
            if (i > 0 && motifHitStarts[i] - motifHitStarts[i-1] < motif.length()) {
                System.err.println("Overlapping motifs " + motifHitStarts[i-1] + "," + motifHitStarts[i]);
            }

        }       
        int motifIndex = -1;
        boolean plusStrand = true;
        boolean isBound = false;
        for (int pos = region.getStart(); pos < region.getEnd(); pos++) {
            char letter = sequence[pos - region.getStart()];

            int mi = Arrays.binarySearch(motifHitStarts, pos);
            boolean ms = mi >= 0;
            if (ms) {
                plusStrand = hits.get(mi).getStrand().equals("+");
            }

            if (motifIndex >= 0) {
                motifIndex++;
            } 
            if (motifIndex == motif.length()) {
                motifIndex = -1;                
            } 
            if (motifIndex == -1 && ms) {
                motifIndex = 0;
                isBound = hasBinding(pos + motif.length() / 2);
            }
            if (motifIndex == -1) {
                isBound = hasBinding(pos);
            }

            int count = counts.getCount(pos);
            HMMState state = null;
            int statenum = 0;
            if (isBound && sensitive && motifIndex >= 0) {
                state = plusStrand ? motifForwStates[motifIndex] : motifRevStates[motifIndex];
                statenum = 2 + (plusStrand ? 0 : motif.length()) + motifIndex;
            } else if (sensitive) {
                state = sensState;
                statenum = 1;
            } else {
                state = insensState;
                statenum = 0;
            }
            if (lastStateNum != -1) {
                transitions[lastStateNum][statenum]++;
            }
            lastStateNum = statenum;
            state.addData(sequence[pos - region.getStart()],
                          count);
        }
    }
    private void printTransitions() {
        for (int i = 0; i < numStates; i++) {
            for (int j = 0; j < numStates; j++) {
                System.out.print(String.format("%d\t",transitions[i][j]));
            }
            System.out.println();
        }
    }
    private void printState(HMMState state, String label) {
        System.out.println("\t" + label);
        System.out.println(state.toString());
    }
    public void saveModel() throws IOException {
        PrintWriter pw = new PrintWriter(modelFname);
        pw.println(numStates);
        pw.print(".98\t.01");
        for (int i = 0; i < motifForwStates.length * 2; i++) {
            pw.print(String.format("\t%.4f",.01 / (2*motifForwStates.length)));
        }

        pw.println();

        pw.println(insensState.serialize());
        pw.println(sensState.serialize());
        for (int i = 0; i < motifForwStates.length; i++) {
            pw.println(motifForwStates[i].serialize());
        }
        for (int i = 0; i < motifRevStates.length; i++) {
            pw.println(motifRevStates[i].serialize());
        }
        for (int i = 0; i < numStates; i++) {
            for (int j = 0; j < numStates; j++) {
                pw.print(String.format("%d\t",transitions[i][j]));
            }
            pw.println();
        }
        pw.close();
    }
    public void printModel() {
        printTransitions();
        printState(insensState,"non sensitive");
        printState(sensState,"sensitive");
        for (int i = 0; i < motifForwStates.length; i++) {
            printState(motifForwStates[i],String.format("motif %d",i));
        }
        for (int i = 0; i < motifRevStates.length; i++) {
            printState(motifRevStates[i],String.format("motif %d",i));
        }

    }
}


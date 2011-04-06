package edu.mit.csail.cgs.projects.dnaseq;

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

/**
 * Print binding calls for comparison to HMM output.  If a motif is provided,
 * then the positions are the motif positions that are close to binding events.
 * If no motif is specified, then just output the binding calls
 */

public class PrintBindingCalls {
    private WeightMatrix motif;
    private ChipSeqLoader loader;
    private ChipSeqAnalysis binding, dnaseq;
    private Genome genome;
    private SequenceGenerator seqgen;
    private float motifCutoff;
    private WMHitStartComparator hitcomp;
    private List<Region> regions;
    private int bindingDistance;

    public PrintBindingCalls() {
        hitcomp = new WMHitStartComparator();
    }
    public void parseArgs(String args[]) throws SQLException, NotFoundException {
        binding = Args.parseChipSeqAnalysis(args,"chipseq");
        bindingDistance = Args.parseInteger(args,"distance",10);
        dnaseq = null;
        try {
            dnaseq = Args.parseChipSeqAnalysis(args,"dnaseq");
        } catch (RuntimeException e) {
            // don't worry, this just means none was specified
        }
        genome = Args.parseGenome(args).cdr();
        Collection<WeightMatrix> matrices = Args.parseWeightMatrices(args);
        Iterator<WeightMatrix> iter = matrices.iterator();
        if (iter.hasNext() && matrices.size() < 10) {
            motif = iter.next();
        } else {
            motif = null;
        }
        if (motif != null) {
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
            motifCutoff = (float)(motif.getMaxScore() * Args.parseDouble(args,"cutoff",.7));
        }
        seqgen = new SequenceGenerator(genome);
        seqgen.useLocalFiles(true);
        seqgen.useCache(true);
        regions = Args.parseRegionsOrDefault(args);
    }
    public void run() throws SQLException {

        for (Region region : regions) {
            List<ChipSeqAnalysisResult> dnaseqResults = null;
            if (dnaseq != null) {
                dnaseqResults = new ArrayList<ChipSeqAnalysisResult>();
                for (ChipSeqAnalysisResult r : dnaseq.getResults(genome, region)) {
                    dnaseqResults.add(r);
                }
            }
            

            if (motif == null) {
                for (ChipSeqAnalysisResult result : binding.getResults(genome, region)) {
                    boolean print = true;
                    if (dnaseq != null) {
                        print = false;
                        for (ChipSeqAnalysisResult d : dnaseqResults) {
                            if (d.overlaps(result)) {
                                print = true;
                                break;
                            }
                        }
                    }
                    if (print) {
                        System.out.println(result.toString());
                    }
                }
            } else {
                char[] sequence = seqgen.execute(region).toCharArray();
                List<WMHit> hits = WeightMatrixScanner.scanSequence(motif, motifCutoff, sequence);
                Collections.sort(hits, hitcomp);
                int[] motifHitStarts = new int[hits.size()];
                for (int i = 0; i < hits.size(); i++) {
                    motifHitStarts[i] = hits.get(i).getStart() + region.getStart();
                }
                List<Region> bindingEvents = new ArrayList<Region>();
                for (ChipSeqAnalysisResult result : binding.getResults(genome, region)) {
                    bindingEvents.add(result.expand(bindingDistance, bindingDistance));
                }
                for (int i = 0; i < motifHitStarts.length; i++) {
                    int pos = motifHitStarts[i] + motif.length();
                    for (Region b : bindingEvents) {
                        if (pos >= b.getStart() && pos <= b.getEnd()) {
                            boolean print = true;
                            if (dnaseq != null) {
                                print = false;
                                for (ChipSeqAnalysisResult d : dnaseqResults) {
                                    if (d.overlaps(b)) {
                                        print = true;
                                        break;
                                    }
                                }
                            }
                            if (print) {                            
                                System.out.println(region.getChrom() + ":" + motifHitStarts[i] + "-" + (motifHitStarts[i] + motif.length()));
                                break;
                            }
                        }
                    }
                }

            }


            

        }

    }

    public static void main(String args[]) throws NotFoundException, SQLException {
        PrintBindingCalls pbc = new PrintBindingCalls();
        pbc.parseArgs(args);
        pbc.run();
    }


}
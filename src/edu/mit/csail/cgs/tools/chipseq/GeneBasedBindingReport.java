package edu.mit.csail.cgs.tools.chipseq;

import java.util.*;
import java.sql.*;
import java.io.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.projects.readdb.*;
import edu.mit.csail.cgs.datasets.function.*;
import edu.mit.csail.cgs.datasets.motifs.*;
import edu.mit.csail.cgs.tools.motifs.WeightMatrixScanner;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.utils.Enrichment;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.tools.utils.Args;


public class GeneBasedBindingReport {

    /**
     * java edu.mit.csail.cgs.tools.chipseq.GeneBasedBindingReport --species "$MM;mm9" \
     * --analysisname "PPG ES iCdx2 p2A 7-28-10 lane 5 (36bp)"          \
     * --analysisversion "vs PPG Day4 null-antiV5 iTF_iOlig2 1 (default params) run 2 round 3" \
     * --genes refGene [--proxup 5000] [--proxdown200] [--up 10000] [--intronlen 10000] [--thresh .001]
     *
     * Output columns are
     * 0) gene name
     * 1) positions of distal binding events
     * 2) positions of proximal binding events
     * 3) positions of intronic or exonic binding events
     */

    private Genome genome;
    private ChipSeqAnalysis analysis;
    private List<RefGeneGenerator> geneGenerators;
    private int proxup, proxdown, up, intronlen, analysisdbid;
    private double thresh;

    public static void main(String args[]) throws Exception {
        GeneBasedBindingReport report = new GeneBasedBindingReport();
        report.parseArgs(args);
        report.run();
    }
    public void parseArgs(String args[]) throws SQLException, NotFoundException {
        genome = Args.parseGenome(args).cdr();
        analysis = null;
        analysis = Args.parseChipSeqAnalysis(args,"analysis");                                                
        geneGenerators = Args.parseGenes(args);
        proxup = Args.parseInteger(args,"proxup",4000);
        up = Args.parseInteger(args,"up",10000);
        intronlen = Args.parseInteger(args,"intronlen",10000);
        proxdown = Args.parseInteger(args,"proxdown",200);
        thresh = Args.parseDouble(args,"thresh",.01);
        if (intronlen < proxdown) {
            intronlen = proxdown + 1;
        }
    }


    public void run() throws SQLException {
        ArrayList<ChipSeqAnalysisResult> proxevents = new ArrayList<ChipSeqAnalysisResult>(), 
            distalevents= new ArrayList<ChipSeqAnalysisResult>() , intronevents = new ArrayList<ChipSeqAnalysisResult>();
        for (RefGeneGenerator generator : geneGenerators) {
            Iterator<Gene> all = generator.getAll();
            while (all.hasNext()) {
                Gene g = all.next();
                proxevents.clear(); distalevents.clear(); intronevents.clear();
                StrandedRegion wholeRegion = g.expand(up,0);
                Region distalPromoter = g.getStrand() == '+' ? new Region(g.getGenome(), g.getChrom(), g.getStart() - up, g.getStart() - proxup) :
                    new Region(g.getGenome(), g.getChrom(), g.getEnd() + proxup, g.getEnd() + up);

                Region proximalPromoter = g.getStrand() == '+' ? new Region(g.getGenome(), g.getChrom(), g.getStart() - proxup, g.getStart() + proxdown) :
                    new Region(g.getGenome(), g.getChrom(), g.getEnd() - proxdown, g.getEnd() + proxup);
                Region intronicRegion  = g.getStrand() == '+' ? new Region(g.getGenome(), g.getChrom(), g.getStart() + proxdown, g.getStart() + intronlen) :
                    new Region(g.getGenome(), g.getChrom(), g.getEnd() - intronlen, g.getEnd() - proxdown);

                Collection<ChipSeqAnalysisResult> allResults = analysis.getResults(genome, wholeRegion);
                for (ChipSeqAnalysisResult r : allResults) {
                    if (!Double.isInfinite(r.pvalue) && r.pvalue > thresh) {
                        System.err.println("Skipping " + r + " because " + r.pvalue);
                        continue;
                    }
                    if (r.overlaps(distalPromoter)) {
                        distalevents.add(r);
                    } else if (r.overlaps(proximalPromoter)) {
                        proxevents.add(r);
                    } else if (r.overlaps(intronicRegion)) {
                        intronevents.add(r);
                    } 

                }
                System.out.println(String.format("%s\t%s\t%s\t%s",
                                                 g.getName(),
                                                 distalevents.toString(),
                                                 proxevents.toString(),
                                                 intronevents.toString()));
            }
        }
    }
}
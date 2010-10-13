package edu.mit.csail.cgs.tools.chipseq;

import java.util.*;
import java.sql.*;
import java.io.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.tools.utils.Args;



/**
 * java edu.mit.csail.cgs.tools.chipseq.GeneBasedBindingReport --species "$MM;mm9" \
 * --analysisname "PPG ES iCdx2 p2A 7-28-10 lane 5 (36bp)"              \
 * --analysisversion "vs PPG Day4 null-antiV5 iTF_iOlig2 1 (default params) run 2 round 3" \
 * --genes refGene [--proxup 5000] [--proxdown200] [--up 10000] [--intronlen 10000] [--thresh .001]
 *
 * Output columns are
 * 0) gene name
 * 1) positions of distal binding events
 * 2) positions of proximal binding events
 * 3) positions of intronic or exonic binding events
 */


public class GeneBasedBindingReport extends GeneBasedReport {

    private ChipSeqAnalysis analysis;
    private double pvalthresh;

    public static void main(String args[]) throws Exception {
        GeneBasedBindingReport report = new GeneBasedBindingReport();
        report.parseArgs(args);
        report.getRegions(args);
        report.report();
    }

    public void parseArgs(String args[]) throws SQLException, NotFoundException {
        super.parseArgs(args);
        analysis = null;
        analysis = Args.parseChipSeqAnalysis(args,"analysis");                                                        
        pvalthresh = Args.parseDouble(args,"thresh",.01);
    }
    public void getRegions(String args[]) {}
    public Collection<ChipSeqAnalysisResult> getOverlappingRegions (Region wholeRegion) throws SQLException {
        ArrayList<ChipSeqAnalysisResult> output = new ArrayList<ChipSeqAnalysisResult>();
        for (ChipSeqAnalysisResult r : analysis.getResults(getGenome(), wholeRegion)) {
            if (!Double.isInfinite(r.pvalue) && r.pvalue > pvalthresh) {
                continue;
            }
            output.add(r);
        }
        return output;
    }

}
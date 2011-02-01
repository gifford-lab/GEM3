package edu.mit.csail.cgs.tools.chipseq;

import java.sql.SQLException;
import java.util.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.*;


/**
 * Base class for comparing two chipseq analyses.
 *
 * main() should instantiate the object, call parseArgs(),
 * and then call printOutputEvents().
 *
 * Subclasses must override getOutputEvents() to provide the
 * list of output events.
 *
 * parseArgs() here looks for --one 'analysisname;version',
 * --two 'analysisname;version' and --maxdist
 */

public abstract class CompareTwoAnalyses {

    private ChipSeqAnalysis one, two;
    private int maxDistance = 20;
    private Genome genome;
    private List<Region> analysisRegions;

    public CompareTwoAnalyses() {}

    public void parseArgs(String args[]) throws NotFoundException, SQLException {
        one = Args.parseChipSeqAnalysis(args,"one");
        two = Args.parseChipSeqAnalysis(args,"two");
        maxDistance = Args.parseInteger(args,"maxdist",maxDistance);
        genome = Args.parseGenome(args).cdr();
        analysisRegions = Args.parseRegionsOrDefault(args);
    }
    public ChipSeqAnalysis getAnalysisOne() {return one;}
    public ChipSeqAnalysis getAnalysisTwo() {return two;}
    public Genome getGenome() {return genome;}
    public List<ChipSeqAnalysisResult> getResultsOne() throws SQLException {
        List<ChipSeqAnalysisResult> output = new ArrayList<ChipSeqAnalysisResult>();
        for (Region r : getAnalysisRegions()) {
            output.addAll(getResultsOne(r));
        }
        return output;
    }
    public List<ChipSeqAnalysisResult> getResultsTwo() throws SQLException {
        List<ChipSeqAnalysisResult> output = new ArrayList<ChipSeqAnalysisResult>();
        for (Region r : getAnalysisRegions()) {
            output.addAll(getResultsTwo(r));
        }
        return output;
    }
    public List<ChipSeqAnalysisResult> getResultsOne(Region region) throws SQLException {
        return one.getResults(region);
    }
    public List<ChipSeqAnalysisResult> getResultsTwo(Region region) throws SQLException {
        return two.getResults(region);
    }
    public List<Region> getAnalysisRegions() {
        return analysisRegions;
    }
    /** returns true if the sorted list contains r, according
        to whatever overlap criteria we're using
     */
    public boolean containsMatch(List<ChipSeqAnalysisResult> list,
                                 ChipSeqAnalysisResult r) {
        int i = Collections.binarySearch(list,r);
        if (i >= 0) {
            return true;
        }
        int inspt = -1 - i;
        System.err.println("Looking for " + r + " gave " + inspt);
        try {
            System.err.println("  " + list.get(inspt));
            System.err.println("  " + list.get(inspt-1));
        } catch (Exception e) {        }

        if (inspt < list.size() && list.get(inspt).distance(r) <= maxDistance) {
            System.err.println("*");
            return true;
        }
        if (inspt > 0 && list.get(inspt-1).distance(r) <= maxDistance) {
            System.err.println("*");
            return true;
        }
        return false;
    }
    public abstract List<ChipSeqAnalysisResult> getOutputEvents() throws SQLException ;
    public void printOutputEvents() throws SQLException {
        for (ChipSeqAnalysisResult r : getOutputEvents()) {
            System.out.println(String.format("%s:%d-%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2e\t%.1f",
                                             r.getChrom(), r.getStart(), r.getEnd(),
                                             r.getPosition(),
                                             r.getFG(), r.getBG(), r.getStrength(), r.getShape(),
                                             r.getPValue(), r.getFoldEnrichment()));
        }
    }
}
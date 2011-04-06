package edu.mit.csail.cgs.projects.dnaseq;

import java.io.*;
import java.util.*;
import java.sql.SQLException;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.projects.readdb.ClientException;
import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;

/**
 * DNASeqEnrichmentCaller --species "$MM;mm9" 
 * 
 * To analyze the same alignments as an existing analysis:
 * --dnaseq "analysisname;analysisversion"
 *
 * --dnaseqfg "alignname;alignversion" --dnaseqbg "alignname;alignversion"
 *
 * --windowsize 40
 * --windowstep 10
 * --combinegap 40
 * --minwidth 60
 * --minfold 1.5
 * --subsample .2
 *
 * Output columns are
 * 0) chrom
 * 1) start
 * 2) end
 * 3) center
 * 4) foreground read count
 * 5) background read count
 * 6) enrichment ratio
 * 7) pvalue
 */

public class DNASeqEnrichmentCaller {

    protected HMMReads reads;    
    protected Genome genome;
    
    protected ChipSeqLoader loader;
    private Poisson poisson;
    private int sensitiveRegionWindowSize = 40;
    private int sensitiveRegionWindowStep = 10;
    private int combineRegionsGap = 40;
    private int minimumRegionWidth = 60;
    private double totalFgReads, totalBgReads;
    private double channelScalingFactor, minFoldChange, subsample;
    private List<ChipSeqAlignment> alignments, bgAlignments;

    private List<Region> regionsToCall;

    public DNASeqEnrichmentCaller() throws IOException, ClientException, SQLException {
        loader = new ChipSeqLoader();
        DRand re = new DRand();
		poisson = new Poisson(0, re);
        reads = new HMMReads();
    }
    public Genome getGenome() {return genome;}
    public HMMReads getReads() {return reads;}
    public List<Region> getRegionsToCall() {return regionsToCall;}
    public List<ChipSeqAlignment> getAlignments() {return alignments;}
    public List<ChipSeqAlignment> getBGAlignments() {return bgAlignments;}
    public List<ChipSeqAlignment> getFGAlignments() {return alignments;}
    public void parseArgs(String args[]) throws NotFoundException, SQLException, IOException {
        genome = Args.parseGenome(args).cdr();
        alignments = new ArrayList<ChipSeqAlignment>();
        bgAlignments = new ArrayList<ChipSeqAlignment>();
        if (Args.parseString(args,"dnaseq",null) != null) {
            try {
                ChipSeqAnalysis dnaseq = Args.parseChipSeqAnalysis(args,"dnaseq");
                alignments.addAll(dnaseq.getForeground());
                bgAlignments.addAll(dnaseq.getBackground());
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        reads.smooth(Args.parseInteger(args,"smooth",0));
        subsample = Args.parseDouble(args,"subsample",2);
        if (subsample < 1) {
            reads.subSample(subsample);
        }

        List<ChipSeqLocator> fg = Args.parseChipSeq(args,"dnaseqfg");
        List<ChipSeqLocator> bg = Args.parseChipSeq(args,"dnaseqbg");
        for (ChipSeqLocator locator : fg) {
            System.err.println("fg locator " + locator);
            alignments.addAll(loader.loadAlignments(locator,genome));
            System.err.println("Alignments now are " + alignments);
        }
        for (ChipSeqLocator locator : bg) {
            bgAlignments.addAll(loader.loadAlignments(locator,genome));
        }
        if (alignments.size() == 0) {
            throw new NotFoundException("Didn't find any foreground chipseq alignments");
        }
        if (bgAlignments.size() == 0) {
            throw new NotFoundException("Didn't find any background chipseq alignments");
        }

        sensitiveRegionWindowSize = Args.parseInteger(args,"windowsize",sensitiveRegionWindowSize);
        sensitiveRegionWindowStep = Args.parseInteger(args,"windowstep",sensitiveRegionWindowStep);
        combineRegionsGap = Args.parseInteger(args,"combinegap",combineRegionsGap);
        minimumRegionWidth = Args.parseInteger(args,"minwidth",minimumRegionWidth);
        minFoldChange = Args.parseDouble(args,"minfold",1.5);

        totalFgReads = 0;
        for (ChipSeqAlignment a : alignments) {
            totalFgReads += loader.countAllHits(a);
            System.err.println(a.getDBID() + " -> " + totalFgReads);
        }
        if (subsample < 1) {
            totalFgReads = (int)(subsample * totalFgReads);
        }

        totalBgReads = 0;
        for (ChipSeqAlignment a : bgAlignments) {
            totalBgReads += loader.countAllHits(a);
            System.err.println(a.getDBID() + " -> " + totalBgReads);
        }
        channelScalingFactor = totalFgReads / totalBgReads;
        System.err.println(String.format("Channel scaling factor is %.0f/%.0f=%.2f", totalFgReads,totalBgReads,channelScalingFactor));
        regionsToCall = Args.parseRegionsOrDefault(args);
    }
    public List<ChipSeqAnalysisResult> getHyperSensitiveRegions(Region queryRegion) throws ClientException, IOException {
        return getHyperSensitiveRegions(queryRegion,
                                        reads.getReadCounts(queryRegion, alignments));
    }
    public List<ChipSeqAnalysisResult> getHyperSensitiveRegions(Region queryRegion, ReadCounts fgReadCounts) throws ClientException, IOException {

        reads.subSample(2);
        ReadCounts bgReadCounts = reads.getReadCounts(queryRegion,
                                                      bgAlignments);
        reads.subSample(subsample);

        List<Region> enriched = new ArrayList<Region>();
        int start = queryRegion.getStart();
        int lastGoodStart = -1;

        while (start + sensitiveRegionWindowSize < queryRegion.getEnd()) {
            int end = start + sensitiveRegionWindowSize;
            int fgcount = 0;
            int bgcount = 0;
            for (int i = start; i < end; i++) {
                fgcount += fgReadCounts.getCount(i);
                bgcount += bgReadCounts.getCount(i);
            }
            poisson.setMean(bgcount * channelScalingFactor * minFoldChange + 1);
            double bgcdf = poisson.cdf(fgcount - 1);
            poisson.setMean(1 + minFoldChange * totalFgReads * sensitiveRegionWindowSize / 2080000000.0);
            double unicdf = poisson.cdf(fgcount - 1);
            double cdf = Math.min(unicdf, bgcdf);
            //            System.err.println(String.format("\t %d to +%d :: %d / %d -> %.2f %.2f", start, sensitiveRegionWindowSize, fgcount, bgcount, bgcdf, unicdf));
            //            double cdf = bgcdf;
            if (cdf >= .99) {
                if (lastGoodStart == -1) {
                    lastGoodStart = start;
                }
            } else {
                if (lastGoodStart > 0) {
                    enriched.add(new Region(queryRegion.getGenome(), queryRegion.getChrom(), lastGoodStart, start - sensitiveRegionWindowStep + sensitiveRegionWindowSize));
                }
                lastGoodStart = -1;
            }
            start += sensitiveRegionWindowStep;
        }
        if (lastGoodStart != -1) {
            enriched.add(new Region(queryRegion.getGenome(), queryRegion.getChrom(), lastGoodStart, queryRegion.getEnd()));
        }
        return combineRegions(enriched, fgReadCounts, bgReadCounts);
    }
    private List<ChipSeqAnalysisResult> combineRegions(List<Region> input, ReadCounts fgCounts, ReadCounts bgCounts) {
        List<ChipSeqAnalysisResult> output = new ArrayList<ChipSeqAnalysisResult>();
        Collections.sort(input);
        int i = 0;
        Region accumulating = null;
        while (i < input.size()) {
            if (accumulating == null) {
                accumulating = input.get(i);
            } else {
                if (accumulating.getChrom().equals(input.get(i).getChrom()) &&
                    accumulating.distance(input.get(i)) < combineRegionsGap) {
                    accumulating = new Region(accumulating.getGenome(),
                                              accumulating.getChrom(),
                                              accumulating.getStart(),
                                              input.get(i).getEnd());
                } else {
                    if (accumulating.getWidth() > minimumRegionWidth) {
                        output.add(regionToDNASeq(accumulating, fgCounts, bgCounts));
                    }
                    accumulating = input.get(i);
                }
            }
            i++;
        }
        if (accumulating != null) {
            if (accumulating.getWidth() > minimumRegionWidth) {
                output.add(regionToDNASeq(accumulating, fgCounts, bgCounts));
            }
        }
        return output;
    }
    private ChipSeqAnalysisResult regionToDNASeq(Region r, ReadCounts fgCounts, ReadCounts bgCounts) {
        double fc = 0, bc = 0;
        for (int i = r.getStart(); i < r.getEnd(); i++) {
            fc += fgCounts.getCount(i);
            bc += bgCounts.getCount(i);
        }
        poisson.setMean(bc * channelScalingFactor * minFoldChange);
        double bgcdf = poisson.cdf((int)(fc - 1));
        poisson.setMean(minFoldChange * totalFgReads * sensitiveRegionWindowSize / 2080000000.0);
        double unicdf = poisson.cdf((int)(fc - 1));
        double cdf = Math.min(unicdf, bgcdf);

        double fold = Math.min(fc / (bc * channelScalingFactor),
                               fc / (totalFgReads * sensitiveRegionWindowSize / 2080000000.0));

        return new ChipSeqAnalysisResult(r.getGenome(), r.getChrom(), r.getStart(), r.getEnd(), 
                                         (r.getStart() + r.getEnd()) / 2,                                         
                                         fc, bc,
                                         0.0, //strength
                                         0.0, // shape
                                         1 - cdf, // pvalue
                                         fold); //foldchange
    }
    public static void main(String args[]) throws Exception {
        DNASeqEnrichmentCaller caller = new DNASeqEnrichmentCaller();
        caller.parseArgs(args);
        System.out.println("Chrom\tStart\tEnd\tCenter\tFG\tBG\tRatio\tPvalue");
        for (Region region : caller.getRegionsToCall()) {
            try {
                for (ChipSeqAnalysisResult call : caller.getHyperSensitiveRegions(region)) {
                    System.out.println(String.format("%s\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%f",
                                                     call.getChrom(),
                                                     call.getStart(),
                                                     call.getEnd(),
                                                     (call.getStart() + call.getEnd())/2,
                                                     call.getFG(),
                                                     call.getBG(),
                                                     call.getFoldEnrichment(),
                                                     call.getPValue()));
                }
            } catch (Exception e) {
                System.err.println("Error on " + region + ".  Skipping");
                System.err.println(e.toString());
            }
        }
    }
    

}
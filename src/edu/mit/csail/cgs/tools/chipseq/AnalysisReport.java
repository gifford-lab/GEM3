package edu.mit.csail.cgs.tools.chipseq;

import java.util.*;
import java.sql.*;
import java.io.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.datasets.function.*;
import edu.mit.csail.cgs.datasets.motifs.*;
import edu.mit.csail.cgs.tools.motifs.WeightMatrixScanner;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.utils.Enrichment;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.tools.utils.Args;

// need to add motif comparison!

/**
 * Print a report about a chipseq analysis run.
 *
 * java edu.mit.csail.cgs.tools.chipseq.AnalysisReport --species "$MM;mm9" \
 * --analysisname "PPG ES iCdx2 p2A 7-28-10 lane 5 (36bp)" \
 * --analysisversion "vs PPG Day4 null-antiV5 iTF_iOlig2 1 (default params) run 2 round 3" \
 * --genes refGene [--up 5000] [--down 200] [--thresh .001] [--nogenes] [--nogo]
 * [--output base_file_name] [--dumpevents] [--dumpgenes] [--dumpgo]
 * [--topevents 10000]
 *
 * --up and --down specify how far upstream and downstream from a gene's TSS should
 * make a binding event count as binding that gene.
 * --thresh is the p-value threshold for GO enrichment
 * --nogo and --nogenes suppress reporting of GO cats and bound genes
 */

public class AnalysisReport {

    private Genome genome;
    private ChipSeqAnalysis analysis;
    private List<RefGenePromoterGenerator> promoterGenerators;
    private Collection<WeightMatrix> matrices;
    private int up, down, topEvents;
    private double thresh, wmcutoff;
    private boolean dumpEvents, dumpGenes, dumpGO, html, noGO, noGenes;
    private String outputBase;

    // these all get filled in by run() and
    // used by report() 
    private Collection<Region> regions;
    private List<ChipSeqAnalysisResult> events;
    private Set<Gene> boundGenes;
    private List<RefGeneGenerator> geneGenerators;
    private List<Enrichment> enrichments;

    public static void main(String args[]) throws Exception {
        AnalysisReport report = new AnalysisReport();
        report.parseArgs(args);
        report.run();
        report.printReport();
        report.runMotifReport();
    }

    public AnalysisReport() {
    }
    public void parseArgs(String args[]) throws SQLException, NotFoundException, IOException {
        genome = Args.parseGenome(args).cdr();
        analysis = null;
        analysis = Args.parseChipSeqAnalysis(args,"analysis");                                                
        regions = Args.parseRegionsOrDefault(args);
        geneGenerators = Args.parseGenes(args);
        promoterGenerators = new ArrayList<RefGenePromoterGenerator>();
        matrices = Args.parseWeightMatrices(args);
        up = Args.parseInteger(args,"up",4000);
        down = Args.parseInteger(args,"down",200);
        thresh = Math.log(Args.parseDouble(args,"thresh",.01));
        wmcutoff = Args.parseDouble(args,"wmcutoff",.4);
        topEvents = Args.parseInteger(args,"topevents",-1);
        for (RefGeneGenerator rg : geneGenerators) {
            promoterGenerators.add(new RefGenePromoterGenerator(genome,
                                                                rg.getTable(),
                                                                up,
                                                                down));
        }
        outputBase = analysis.getName() + "___" + analysis.getVersion();
        outputBase = outputBase.replaceAll("[^\\w]+","_");
        outputBase = Args.parseString(args,"output",outputBase);
        dumpEvents  = Args.parseFlags(args).contains("dumpevents");
        dumpGenes = Args.parseFlags(args).contains("dumpgenes");
        dumpGO = Args.parseFlags(args).contains("dumpgo");
        html = Args.parseFlags(args).contains("html");
        noGO = Args.parseFlags(args).contains("nogo");
        noGenes = Args.parseFlags(args).contains("nogenes");
        if (noGenes) {
            noGO = true;
        }
        events = new ArrayList<ChipSeqAnalysisResult>();

        MarkovBackgroundModel bgModel = null;
        String bgmodelname = Args.parseString(args,"bgmodel","whole genome zero order");
        BackgroundModelMetadata md = BackgroundModelLoader.getBackgroundModel(bgmodelname,
                                                                              1,
                                                                              "MARKOV",
                                                                              genome.getDBID());
        if (md != null) {
            bgModel = BackgroundModelLoader.getMarkovModel(md);
        } else {
            System.err.println("Couldn't get metadata for " + bgmodelname);
        }
        if (bgModel == null) {
            for (WeightMatrix m : matrices) {
                m.toLogOdds();
            }            
        } else {
            for (WeightMatrix m : matrices) {
                m.toLogOdds(bgModel);
            }            
        }        
        dumpParams(args);
    }
    private void getEvents() throws SQLException {
        for (Region r : regions) {
            events.addAll(analysis.getResults(genome,r));
        }
    }
    private void getBoundGenes() throws SQLException {
        boundGenes = new HashSet<Gene>();
        for (RefGenePromoterGenerator pg : promoterGenerators) {
            for (ChipSeqAnalysisResult event : events) {
                Iterator<Gene> iter = pg.execute(event);
                while (iter.hasNext()) {
                    boundGenes.add(iter.next());
                }
            }
        }
    }
    private void getEnrichedCategories() throws SQLException {
        
        FunctionLoader funcloader = new GOFunctionLoader(GOFunctionLoader.getDefaultDBName());
        FunctionalUtils funcutils = new FunctionalUtils(funcloader, genome.getSpecies());
        
        Set<String> allGeneNames = new HashSet<String>();
        Map<String,Integer> chromlengths = genome.getChromLengthMap();
        for (Region r : regions) {
            for (RefGeneGenerator generator : geneGenerators) {
                Iterator<Gene> genes = generator.execute(r);
                while (genes.hasNext()) {
                    allGeneNames.add(genes.next().getName());
                }
            }
        }
        
        Set<String> foreground = new HashSet<String>();
        for (Gene g : boundGenes) {
            foreground.add(g.getName());
        }

        Map<String,Enrichment> output = funcutils.calculateTotalEnrichments(allGeneNames, foreground);
        enrichments = new ArrayList<Enrichment>();
        for (String c : output.keySet()) {
            Enrichment e = output.get(c);
            if (e.getLogPValue() <= thresh) {
                enrichments.add(e);
            }
        }
    }

    public void run() throws SQLException {
        System.err.println("Getting Events");
        getEvents();
        Collections.sort(events, new ChipSeqAnalysisPValueComparator());
        if (topEvents > 0 && topEvents < events.size()) {
            events = events.subList(events.size()-topEvents,events.size());
        }
        if (noGenes) {
            boundGenes = new HashSet<Gene>();
            enrichments = new ArrayList<Enrichment>();
            return;
        }
        System.err.println("Getting Bound Genes");
        getBoundGenes();
        if (noGO) {
            enrichments = new ArrayList<Enrichment>();
            return;
        }
        System.err.println("Getting enriched categories");
        getEnrichedCategories();
        Collections.sort(enrichments, new EnrichmentPvalueComparator());
    }
    public void dumpParams(String args[]) throws IOException {
        PrintWriter pw = new PrintWriter(outputBase + ".params");
        pw.print("cmdline=java edu.mit.csail.cgs.tools.chipseq.AnalysisReport ");
        for (int i = 0; i < args.length; i++) {
            pw.print(args[i] + " ");
        }
        pw.println();
        pw.println("genome=" + genome.toString());
        pw.println("up=" + up);
        pw.println("down=" + down);
        pw.println("topEvents=" + topEvents);
        pw.println("thresh=" + thresh);
        pw.println("wmcutoff=" + wmcutoff);
        for (WeightMatrix m : matrices) {
            pw.println("matrix="+m.toString());
        }
        pw.println("analysis=" + analysis.toString());


        
    }
    public void dumpEvents() throws IOException {
        if (!dumpEvents) {
            return;
        }
        PrintWriter pw = new PrintWriter(outputBase + ".events");
        for (ChipSeqAnalysisResult a : events) {
            pw.println(String.format("%s\t%.2f\t%.2e\t%.2f\t%.2f\t%.2f",
                                     a.toString(),
                                     a.strength,
                                     a.pvalue,
                                     a.foldEnrichment,
                                     a.foregroundReadCount,
                                     a.backgroundReadCount));
        }
        pw.close();

    }
    public void dumpGenes() throws IOException {
        if (!dumpGenes) {
            return;
        }
        PrintWriter pw = new PrintWriter(outputBase + ".genes");
        for (Gene g : boundGenes) {
            pw.println(g.toString());
        }
        pw.close();

    }
    public void dumpGO() throws IOException {
        if (!dumpGO) {
            return;
        }
        PrintWriter pw = new PrintWriter(outputBase + ".go");
        for (Enrichment e : enrichments) {
            pw.println(e.toString());
        }
        pw.close();

    }
    public void printReport() throws IOException {
        dumpEvents();
        dumpGenes();
        dumpGO();
        if (html) {
            printHTMLReport();
        } else {
            printTextReport();
        }
    }
    public void printHTMLReport() {
        System.err.println("<tr><th>Analysis</th><th>Total Events</th><th>Bound Genes</th><th>GO Cats</th></tr>");
        System.out.println(String.format("<tr><td>%s,%s,%s</td><td>%d</td><td>%d</td><td>%d</td></tr>\n",
                                         analysis.getName(), 
                                         analysis.getVersion(),
                                         analysis.getProgramName(),
                                         events.size(),
                                         boundGenes.size(),
                                         enrichments.size()));
        
    }
    public void printTextReport() {
        System.out.println(String.format("%s\t%s\t%s",analysis.getName(), 
                                         analysis.getVersion(),
                                         analysis.getProgramName()));
        System.out.println("Total Events\t" + events.size());
        System.out.println("Total Bound Genes\t" + boundGenes.size());
        System.out.println("GO enriched categories\t" + enrichments.size());
    }
    public void runMotifReport() throws IOException {
        if (matrices.size() == 0) {
            System.err.println("No weight matrices available.  Not running report");
            return;
        }

        PrintWriter pw = new PrintWriter(outputBase + ".motifs");
        pw.print("BindingLocation\tBindingStrength\tBindingPval\tBindingFoldEnrichment\tBindingFGCount\tBindingBGCount");
        for (WeightMatrix wm : matrices) {
            pw.print("\t" + wm.toString());
        }
        pw.println();
        SequenceGenerator seqgen = new SequenceGenerator();
        for (ChipSeqAnalysisResult a : events) {
            Region toScan = a.expand(100,100);
            char[] seq = seqgen.execute(toScan).toCharArray();
            pw.print(String.format("%s\t%.2f\t%.2e\t%.2f\t%.2f\t%.2f",
                                     a.toString(),
                                     a.strength,
                                     a.pvalue,
                                     a.foldEnrichment,
                                     a.foregroundReadCount,
                                     a.backgroundReadCount));

            for (WeightMatrix wm : matrices) {
                List<WMHit> hits = WeightMatrixScanner.scanSequence(wm,
                                                                    (float)(wm.getMaxScore() * wmcutoff),
                                                                    seq);
                float bestscore = 0;
                for (WMHit h : hits) {
                    if (h.score > bestscore) {
                        bestscore = h.score;
                    }
                }
                pw.print(String.format("\t%.2f",bestscore));                
            }
            pw.println();
        }
        pw.close();
    }
   

}

class ChipSeqAnalysisRatioComparator implements Comparator<ChipSeqAnalysisResult> {
    public int compare(ChipSeqAnalysisResult a, ChipSeqAnalysisResult b) {
        return Double.compare(a.foldEnrichment,b.foldEnrichment);
    }
    public boolean equals(Object o) {
        return o instanceof ChipSeqAnalysisRatioComparator;
    }
}
class ChipSeqAnalysisPValueComparator implements Comparator<ChipSeqAnalysisResult> {
    public int compare(ChipSeqAnalysisResult a, ChipSeqAnalysisResult b) {
        int r = Double.compare(b.pvalue,a.pvalue);
        if (r == 0) {
            r = Double.compare(a.foldEnrichment,b.foldEnrichment);
        }
        return r;
    }
    public boolean equals(Object o) {
        return o instanceof ChipSeqAnalysisRatioComparator;
    }
}

class EnrichmentPvalueComparator implements Comparator<Enrichment> {

    public boolean equals(Object o) {
        return o instanceof EnrichmentPvalueComparator;
    }
    public int compare(Enrichment a, Enrichment b) {
        double diff = b.getLogPValue() - a.getLogPValue();
        if (diff < 0) {
            return -1;
        } else if (diff > 0) {
            return 1;
        } else {
            return 0;
        }

    }
}
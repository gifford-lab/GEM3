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

// need to add motif comparison!

/**
 * Print a report about a chipseq analysis run.
 *
 * java edu.mit.csail.cgs.tools.chipseq.AnalysisReport --species "$MM;mm9" \
 * --analysisname "PPG ES iCdx2 p2A 7-28-10 lane 5 (36bp)" \
 * --analysisversion "vs PPG Day4 null-antiV5 iTF_iOlig2 1 (default params) run 2 round 3" \
 * --genes refGene [--up 5000] [--down 200] [--eventthresh .001] [--gothresh .001] [--nogenes] [--nogo]
 * [--output base_file_name] [--dumpevents] [--dumpgenes] [--dumpgo] [--dumptrack]
 * [--topevents 10000] [--project syscode]
 *
 * --up and --down specify how far upstream and downstream from a gene's TSS should
 * make a binding event count as binding that gene.
 * --thresh is the p-value threshold for GO enrichment
 * --nogo and --nogenes suppress reporting of GO cats and bound genes
 * --project is used to determine the URL on nanog for the bigbed and bigwig files
 */

public class AnalysisReport {

    private Genome genome;
    private ChipSeqAnalysis analysis;
    private Collection<WeightMatrix> matrices;
    private int up, down, topEvents, analysisdbid;
    private double eventthresh, gothresh, wmcutoff;
    private boolean dumpEvents, dumpGenes, dumpGO, dumpMotifs, html, noGO, noGenes, dumpTrack;
    private String outputBase, projectName;

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
        analysisdbid = analysis.getDBID();
        regions = Args.parseRegionsOrDefault(args);
        geneGenerators = Args.parseGenes(args);
        matrices = Args.parseWeightMatrices(args);
        up = Args.parseInteger(args,"up",4000);
        down = Args.parseInteger(args,"down",200);
        gothresh = Math.log(Args.parseDouble(args,"gothresh",.01));
        eventthresh = Args.parseDouble(args,"eventthresh",.01);
        wmcutoff = Args.parseDouble(args,"wmcutoff",.4);
        topEvents = Args.parseInteger(args,"topevents",-1);
        outputBase = analysis.getName() + "___" + analysis.getVersion();
        outputBase = outputBase.replaceAll("[^\\w]+","_");
        outputBase = Args.parseString(args,"output",outputBase);
        dumpEvents  = Args.parseFlags(args).contains("dumpevents");
        dumpGenes = Args.parseFlags(args).contains("dumpgenes");
        dumpTrack = Args.parseFlags(args).contains("dumptrack");
        dumpMotifs = Args.parseFlags(args).contains("dumpmotifs");
        dumpGO = Args.parseFlags(args).contains("dumpgo");
        html = Args.parseFlags(args).contains("html");
        noGO = Args.parseFlags(args).contains("nogo");
        noGenes = Args.parseFlags(args).contains("nogenes");
        projectName = Args.parseString(args,"project","syscode");
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
            for (ChipSeqAnalysisResult result : analysis.getResults(genome,r)) {
                if (result.pvalue <= eventthresh) {
                    events.add(result);
                }
            }
        }
    }
    private void getBoundGenes() throws SQLException {
        boundGenes = new HashSet<Gene>();
        Map<String,ArrayList<ChipSeqAnalysisResult>> map = new HashMap<String, ArrayList<ChipSeqAnalysisResult>>();
        for (ChipSeqAnalysisResult event : events) {
            if (!map.containsKey(event.getChrom())) {
                map.put(event.getChrom(), new ArrayList<ChipSeqAnalysisResult>());
            }
            map.get(event.getChrom()).add(event);
        }        
        for (RefGeneGenerator generator : geneGenerators) {
            Iterator<Gene> all = generator.getAll();
            while (all.hasNext()) {
                Gene g = all.next();
                Region r = g.expand(up, down - g.getWidth());
                if (!map.containsKey(g.getChrom())) {
                    continue;
                }
                for (ChipSeqAnalysisResult event : map.get(g.getChrom())) {
                    if (event.overlaps(r)) {
                        boundGenes.add(g);
                        break;
                    }
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
            if (e.getLogPValue() <= gothresh) {
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
        pw.println("gothresh=" + gothresh);
        pw.println("eventthresh=" + eventthresh);
        pw.println("wmcutoff=" + wmcutoff);
        for (WeightMatrix m : matrices) {
            pw.println("matrix="+m.toString());
        }
        pw.println("analysis=" + analysis.toString());


        
    }
    public void dumpTrack() throws IOException {
        if (!dumpTrack) {
            return;
        }
        PrintWriter bed = new PrintWriter(analysisdbid + ".bedtrack");
        PrintWriter wiggle = new PrintWriter(analysisdbid + ".wigtrack");

        bed.println(String.format("track type=bigBed name=\"%s calls\" description=\"Calls for %s\" visibility=full db=%s bigDataUrl=http://nanog.csail.mit.edu/%s_tracks/%d.bb",
                                  analysis.getName(), analysis.toString(), genome.getVersion(),projectName,analysisdbid));
        wiggle.println(String.format("track type=bigWig name=\"%s data\" description=\"Data for %s\" visibility=full db=%s bigDataUrl=http://nanog.csail.mit.edu/%s_tracks/%d.bw",
                                     analysis.getName(), analysis.toString(), genome.getVersion(),projectName,analysisdbid));
        bed.close();
        wiggle.close();
        bed = new PrintWriter(analysisdbid + ".bed");
        wiggle = new PrintWriter(analysisdbid + ".wig");

        for (ChipSeqAnalysisResult a : events) {
            bed.println(String.format("chr%s\t%d\t%d\tbinding\t%d",
                                      a.getChrom(),
                                      a.getStart(),
                                      a.getEnd(),
                                      (int)Math.round(a.foldEnrichment*10)));
        }
        bed.close();
        Set<ChipSeqAlignment> fg = analysis.getForeground();
        Collection<String> alignids = new ArrayList<String>();
        for (ChipSeqAlignment align : fg) {
            alignids.add(Integer.toString(align.getDBID()));
        }
        Client client = null;
        try {
            client = new Client();
        } catch (Exception e) {
            e.printStackTrace();
            wiggle.close();
            return;
        }
        Map<String,Integer> chroms = genome.getChromIDMap();
        Map<String,Integer> lengths = genome.getChromLengthMap();
        int binsize = 50;
        for (String chromname : chroms.keySet()) {
            try {
                int chromid = chroms.get(chromname);

                Map<Integer,Integer> plus = client.getHistogram(alignids,
                                                                chromid,
                                                                false,
                                                                false,
                                                                binsize,
                                                                0,1000000000,
                                                                null, // minweight
                                                                true); // plus strand
                Map<Integer,Integer> minus = client.getHistogram(alignids,
                                                                 chromid,
                                                                 false,false,
                                                                 binsize,
                                                                 0,1000000000,
                                                                 null, // minweight
                                                                 false); // plus strand
                wiggle.println(String.format("variableStep chrom=chr%s span=%d",
                                             chromname, binsize/2));
                TreeSet<Integer> allPos = new TreeSet<Integer>();
                allPos.addAll(plus.keySet());
                allPos.addAll(minus.keySet());
                for (int pos : allPos) {
                    if (pos + binsize > lengths.get(chromname)) {
                        break;
                    }
                    if (plus.containsKey(pos)) {
                        wiggle.println(String.format("%d\t%d",
                                                     pos,plus.get(pos)));
                    }
                    if (minus.containsKey(pos)) {
                        wiggle.println(String.format("%d\t-%d",
                                                     pos+binsize/2,minus.get(pos)));
                    }                    
                }
            } catch (ClientException e) {
                e.printStackTrace();
            }                
        }
        if (client != null) {
            client.close();
        }
        wiggle.close();

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
            pw.println(e.reportLine());
        }
        pw.close();

    }
    public void printReport() throws IOException {
        dumpEvents();
        dumpGenes();
        dumpGO();
        dumpTrack();
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
        if (!dumpMotifs) {
            return;
        }
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
        seqgen.useCache(true);
        seqgen.useLocalFiles(true);
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


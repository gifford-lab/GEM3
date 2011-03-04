package edu.mit.csail.cgs.projects.dnaseq;

import java.sql.SQLException;
import java.util.*;
import java.io.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.motifs.*;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.tools.motifs.*;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.projects.readdb.*;

/**
 * Generates the input files for CENTIPEDE given
 * - a motif to scan for
 * - a dnaseq experiment
 * - a set of gene annotations to use
 * - a conservation track
 *
 *
 * The cuts file has one row per motif instance.  The number of columns is motif width + 400- it
 * gives the number of dnaseq cuts over and on 200bp either side of the motif
 *
 * The annots file contains the chromosome, start, stop, and strand of the motif as well
 * as the match score, conservation score, and distance to the nearest start site
 */

public class CentipedeFiles {

    private Genome genome;
    private ChipSeqLoader loader;
    private List<ChipSeqAlignment> alignments;
    private List<RefGeneGenerator> refgenes;
    private SequenceGenerator seqgen;
    private PhastConsGenerator phastcons;
    private WeightMatrix motif;
    private float motifCutoff;
    private List<Region> regions;
    private int centipedeWidth = 100;
    private PrintWriter annotsPW, cutsPW;
    private Client client;


    public CentipedeFiles() {}
    public void parseArgs(String args[]) throws NotFoundException, IOException, ClientException, SQLException {
        loader = new ChipSeqLoader();
        client = new Client();
        genome = Args.parseGenome(args).cdr();
        List<ChipSeqLocator> cslocators = Args.parseChipSeq(args);       
        alignments = new ArrayList<ChipSeqAlignment>();
        for (ChipSeqLocator l : cslocators) {
            alignments.addAll(loader.loadAlignments(l, genome));
        }
        refgenes = Args.parseGenes(args);
        seqgen = new SequenceGenerator();
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
        motifCutoff = (float)(motif.getMaxScore() * Args.parseDouble(args,"cutoff",.7));
        regions = splitRegions(Args.parseRegionsOrDefault(args));

        String outputBase = Args.parseString(args,"outputBase",motif.getName());
        annotsPW = new PrintWriter(outputBase + ".annots");
        cutsPW = new PrintWriter(outputBase + ".cuts");
        annotsPW.println("chrom\tStart\tEnd\tStrand\tPWMscore\tConsScore\tTSSdist");
        phastcons = new PhastConsGenerator(genome, genome.getSpecies().equals("Mus musculus") ? "phastCons30way" : "phastCons46way");
    }
    private List<Region> splitRegions(List<Region> in)  {
        int width = 10000000;
        List<Region> out = new ArrayList<Region>();
        for (Region r : in) {
            for (int start = r.getStart(); start < r.getEnd(); start += width) {
                int end = start + width < r.getEnd() ? start + width : r.getEnd();
                out.add(new Region(r.getGenome(),
                                   r.getChrom(),
                                   start,
                                   end));
            }
        }
        return out;
    }
    public void run() throws ClientException, SQLException, IOException {
        int motifCount = 0;
        for (Region region : regions) {
            try {
                int chromID = genome.getChromID(region.getChrom());
                int starts[] = getStartSites(region);
                char[] sequence = seqgen.execute(region).toCharArray();
                List<WMHit> hits = WeightMatrixScanner.scanSequence(motif, motifCutoff, sequence);
                for (WMHit hit : hits) {
                    if (hit.getStart() + region.getStart() < centipedeWidth) { continue;}

                    motifCount++;
                    int distanceToTSS = closestDistance(region.getStart() + hit.getStart(),
                                                        region.getStart() + hit.getEnd(),
                                                        starts);
                    Region trainingRegion = new Region(region.getGenome(),
                                                       region.getChrom(),
                                                       region.getStart() + hit.getStart() - centipedeWidth,
                                                       region.getStart() + hit.getEnd() + centipedeWidth);
                    Region consRegion = new Region(region.getGenome(),
                                                   region.getChrom(),
                                                   region.getStart() + hit.getStart() - 20,
                                                   region.getStart() + hit.getEnd() + 20);

                    double conservation  =0;
                    Iterator<ScoredRegion> iter = phastcons.execute(consRegion);
                    while (iter.hasNext()) {
                        ScoredRegion sr = iter.next();
                        if (sr.getScore() > conservation) {
                            conservation = sr.getScore();
                        }
                    }                

                    int plus[] = getCounts(trainingRegion, alignments, true);
                    int minus[] = getCounts(trainingRegion, alignments, false);
                    StringBuilder builder = new StringBuilder();
                    for (int i = 0; i < plus.length; i++) {
                        if (i != 0) {
                            builder.append("\t");
                        }
                        builder.append(plus[i]);
                    }
                    for (int i = 0; i < minus.length; i++) {
                        builder.append("\t" + minus[i]);
                    }
                    annotsPW.println(String.format("%s\t%d\t%d\t%s\t%.4f\t%.4f\t%d",
                                                   region.getChrom(),
                                                   region.getStart() + hit.getStart(),
                                                   region.getStart() + hit.getEnd(),
                                                   hit.getStrand(),
                                                   hit.getScore(),
                                                   conservation,
                                                   distanceToTSS));
                    cutsPW.println(builder.toString());
                }
            } catch (Exception e) {
                e.printStackTrace();
            }

        }
        System.err.println("Found " + motifCount + " motif examples");
    }
    private int[] getCounts(Region region,
                            Collection<ChipSeqAlignment> alignments,
                            boolean plusStrand) throws IOException, ClientException {
        int[] output = new int[region.getWidth()];
        int regionstart = region.getStart();
        for (int i = 0; i < output.length; i++) {
            output[i] = 0;
        }
        for (ChipSeqAlignment a : alignments) {
            TreeMap<Integer,Integer> m = client.getHistogram(Integer.toString(a.getDBID()),
                                                             region.getGenome().getChromID(region.getChrom()),
                                                             false, // paired
                                                             false, // read extension
                                                             1, // binsize
                                                             10, // dedup
                                                             region.getStart(),
                                                             region.getEnd(),
                                                             null,  // minweight
                                                             plusStrand,
                                                             true); //isLeft
            for (int pos : m.keySet()) {
                output[pos - regionstart] += m.get(pos);
            }
        }
        return output;        
    }
    private int[] getStartSites(Region r) {
        List<Integer> starts = new ArrayList<Integer>();
        r = r.expand(1000000,1000000);
        for (RefGeneGenerator gen : refgenes) {
            Iterator<Gene> iter = gen.execute(r);
            while (iter.hasNext()) {
                starts.add(iter.next().getFivePrime());
            }
        }
        int[] out = new int[starts.size()];
        for (int i = 0; i < starts.size(); i++) {
            out[i] = starts.get(i);
        }
        Arrays.sort(out);
        return out;
    }
    // returns shortest distance from r to an element of l
    private int closestDistance(int start, int end, int[] l) {
        int inspoint = Arrays.binarySearch(l, start);
        if (inspoint >= 0) {
            return 0;
        } else {
            inspoint = -1 - inspoint;
        }
        int mindist = Integer.MAX_VALUE;
        // fudge this.  we need to look at a few more points because they may be closer
        // to the end of the region
        for (int j = inspoint; j < inspoint + 5 && j < l.length; j++) {
            int d = Math.min(Math.abs(l[j] - start), Math.abs(l[j] - end));
            if (d < mindist) {
                mindist = d;
            }
        }
        return mindist;
    }
    public void close() throws IOException {
        annotsPW.close();
        cutsPW.close();
    }
    public static void main(String args[]) throws Exception {
        CentipedeFiles files = new CentipedeFiles();
        files.parseArgs(args);
        files.run();
        files.close();
    }

}
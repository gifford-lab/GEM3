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


public class GeneBasedReport {

    /**
     * java edu.mit.csail.cgs.tools.chipseq.GeneBasedReport --species "$MM;mm9" \
     * --genes refGene [--proxup 5000] [--proxdown200] [--up 10000] [--intronlen 10000]
     * --regions foo.txt
     *
     * --regions - to read from STDIN
     *
     * Output columns are
     * 0) gene name
     * 1) positions of distal binding events
     * 2) positions of proximal binding events
     * 3) positions of intronic or exonic binding events
     */

    private Genome genome;
    private List<Region> allRegions;
    private List<RefGeneGenerator> geneGenerators;
    private int proxup, proxdown, up, intronlen, analysisdbid;
    private boolean firstIntron;

    public static void main(String args[]) throws Exception {
        GeneBasedReport report = new GeneBasedReport();
        report.parseArgs(args);
        report.getRegions(args);
        report.report();
    }
    public Genome getGenome() {return genome;}
    public void parseArgs(String args[]) throws SQLException, NotFoundException {
        genome = Args.parseGenome(args).cdr();
        geneGenerators = Args.parseGenes(args);
        proxup = Args.parseInteger(args,"proxup",4000);
        up = Args.parseInteger(args,"up",10000);
        intronlen = Args.parseInteger(args,"intronlen",10000);
        proxdown = Args.parseInteger(args,"proxdown",200);
        firstIntron = Args.parseFlags(args).contains("firstintron");
        if (intronlen < proxdown) {
            intronlen = proxdown + 1;
        }
    }
    public void getRegions (String args[]) throws IOException, NotFoundException{
        allRegions = Args.readLocations(args,"regions");
    }

    public Collection<? extends Region> getOverlappingRegions(Region target) throws Exception {
        ArrayList<Region> output = new ArrayList<Region>();
        for (Region r : allRegions) {
            if (r.overlaps(target)) {
                output.add(r);
            }
        }
        return output;
    }

    public void report() throws Exception {
        ArrayList<Region> proxevents = new ArrayList<Region>(), 
            distalevents= new ArrayList<Region>() , intronevents = new ArrayList<Region>();
        for (RefGeneGenerator generator : geneGenerators) {
            generator.retrieveExons(firstIntron);
            Iterator<Gene> all = generator.getAll();
            while (all.hasNext()) {
                Gene g = all.next();
                proxevents.clear(); distalevents.clear(); intronevents.clear();
                StrandedRegion wholeRegion = g.expand(up,0);
                int thisintronlen = Math.min(intronlen, g.getWidth());
                Region distalPromoter = g.getStrand() == '+' ? new Region(g.getGenome(), g.getChrom(), g.getStart() - up, g.getStart() - proxup) :
                    new Region(g.getGenome(), g.getChrom(), g.getEnd() + proxup, g.getEnd() + up);

                Region proximalPromoter = g.getStrand() == '+' ? new Region(g.getGenome(), g.getChrom(), g.getStart() - proxup, g.getStart() + proxdown) :
                    new Region(g.getGenome(), g.getChrom(), g.getEnd() - proxdown, g.getEnd() + proxup);
                Region intronicRegion = null;
                if (firstIntron ) {
                    Iterator<Region> iter = ((ExonicGene)g).getExons();
                    while (iter.hasNext()) {
                        Region intron = iter.next();
                        if (intronicRegion == null || (intron.distance(proximalPromoter) < intronicRegion.distance(proximalPromoter))) {
                            intronicRegion = intron;
                        }
                    }
                    if (intronicRegion == null) {
                        intronicRegion = new Region(g.getGenome(), g.getChrom(), 1,1);
                    }

                } else {
                    intronicRegion = g.getStrand() == '+' ? new Region(g.getGenome(), g.getChrom(), g.getStart() + proxdown, g.getStart() + proxdown + thisintronlen) :
                        new Region(g.getGenome(), g.getChrom(), g.getEnd() - proxdown - thisintronlen, g.getEnd() - proxdown);
                }

                for (Region r : getOverlappingRegions(wholeRegion)) {
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
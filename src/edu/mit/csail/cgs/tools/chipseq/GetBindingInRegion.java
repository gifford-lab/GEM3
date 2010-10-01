package edu.mit.csail.cgs.tools.chipseq;

import java.util.*;
import java.sql.*;
import java.io.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.warpdrive.components.Snapshot;


/**
 * Reads regions on STDIN, from a file, or uses a list such as all genes
 * and produces a line on STDOUT for each that lists either the analyses and binding
 * events for each.
 * 
 * java edu.mit.csail.cgs.tools.chipseq.GetBindingInRegion --species "$MM;mm9" --genes refGene
 * java edu.mit.csail.cgs.tools.chipseq.GetBindingInRegion --species "$MM;mm9" --region 10:0-1000000
 * java edu.mit.csail.cgs.tools.chipseq.GetBindingInRegion --species "$MM;mm9" --regionfile foo.txt
 * java edu.mit.csail.cgs.tools.chipseq.GetBindingInRegion --species "$MM;mm9"
 *     this form above reads regions from STDIN
 * 
 * --plots causes plots with the chipseq data to be produced. 
 *
 */


public class GetBindingInRegion {

    private ChipSeqLoader loader;
    private Genome genome;
    private List<RefGeneGenerator> geneGenerators;
    private List<Region> regions;
    private int plotExpand;
    private boolean plots;
    private String plotPrefix = "";

    public GetBindingInRegion() throws SQLException, IOException { 
        geneGenerators = new ArrayList<RefGeneGenerator>();
        regions = new ArrayList<Region>();
        loader = new ChipSeqLoader();
        plots = false;
    }
    public void parseArgs(String[] args) throws NotFoundException, IOException{
        genome = Args.parseGenome(args).cdr();
        geneGenerators = Args.parseGenes(args);
        List<Region> fromfile = Args.readLocations(args,"regionfile");
        if (fromfile != null) {
            regions.addAll(fromfile);
        }
        regions.addAll(Args.parseRegions(args));
        plots = Args.parseFlags(args).contains("plots");
        plotExpand = Args.parseInteger(args,"plotexpand",2500);
        plotPrefix = Args.parseString(args,"plotprefix","");
    }

    public void setGenome(Genome g) {genome = g;}
    public void setGenes(Collection<RefGeneGenerator> g) {geneGenerators.clear(); geneGenerators.addAll(g);}
    public void setRegions(Collection<Region> r) {regions.clear(); regions.addAll(r);}

    public void computeRegions() throws IOException, SQLException  {
        for (RefGeneGenerator generator : geneGenerators) {
            Iterator<Gene> iter = generator.getAll();
            while (iter.hasNext()) {
                regions.add(iter.next());
            }
        }
        if (regions.size() == 0) {
            System.err.println("Reading regions from STDIN");
            BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
            String line = null;
            while ((line = reader.readLine()) != null) {
                regions.add(Region.fromString(genome, line));
            }
        }
        Collections.sort(regions);
    }
    public List<Region> getRegions() {return regions;}
    public Collection<ChipSeqAnalysis> getAnalysesForRegion(Region r) throws SQLException {
        return ChipSeqAnalysis.withResultsIn(loader, r);
    }
    public void print(Region r) throws SQLException {
        System.out.print(r.toString());
        for (ChipSeqAnalysis a : getAnalysesForRegion(r)) {
            System.out.print("\t" + a.getName() + ";" + a.getVersion());
        }
        System.out.println();
    }
    public void plot(Region r) throws SQLException, NotFoundException, IOException {
        Collection<ChipSeqAnalysis> analyses = getAnalysesForRegion(r);
        if (analyses.size() == 0) { return ;}
        r = r.expand(plotExpand,plotExpand);
        System.err.println("Doing plot for " + r);
        ArrayList<String> args = new ArrayList<String>();
        args.add("--noexit");
        args.add("--species");
        args.add(String.format("%s;%s", genome.getSpecies(), genome.getVersion()));
        args.add("--picture");
        args.add(plotPrefix + (r.toString().replaceAll(":","_")) + ".png");
        args.add("--genes");
        args.add("refGene");
        args.add("--chrom");
        args.add(String.format("%s:%d-%d",r.getChrom(), r.getStart(),r.getEnd()));
        Map<Pair<String,String>, List<String>> chipseqtracks = new HashMap<Pair<String,String>, List<String>>();
        for (ChipSeqAnalysis analysis : analyses) {
            args.add("--chipseqanalysis");
            args.add(analysis.getName() + ";" + analysis.getVersion());
            for (ChipSeqAlignment align : analysis.getForeground()) {
                Pair<String,String> namever = new Pair<String,String>(align.getExpt().getName(),
                                                                      align.getName());
                if (!chipseqtracks.containsKey(namever)) {
                    chipseqtracks.put(namever, new ArrayList<String>());
                }
                chipseqtracks.get(namever).add(align.getExpt().getReplicate());
            }
        }
        for (Pair<String,String> p : chipseqtracks.keySet()) {
            args.add("--chipseq");
            String arg = p.car();
            for (String rep : chipseqtracks.get(p)) {
                arg = arg + ";" + rep;
            }
            arg = arg + ";" + p.cdr();
            args.add(arg);
        }

        System.err.println("Args are " + args);
        String[] asarray = new String[args.size()];
        for (int i = 0; i < args.size(); i++) {
            asarray[i] = args.get(i);
        }

        Snapshot.main(asarray);
    }
    public void close() {
        loader.close();
        loader = null;
        for (RefGeneGenerator g : geneGenerators) {
            g.close();
        }
    }

    public static void main(String args[]) {
        try {
            GetBindingInRegion get = new GetBindingInRegion();
            get.parseArgs(args);
            get.computeRegions();
            for (Region r : get.getRegions()) {
                if (get.plots) {
                    get.plot(r);
                } else {
                    get.print(r);
                }
            }
            get.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.exit(0);
    }
}
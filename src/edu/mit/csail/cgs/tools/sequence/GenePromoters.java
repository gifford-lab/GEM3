package edu.mit.csail.cgs.tools.sequence;

import java.util.*;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.IOException;
import java.sql.SQLException;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.ewok.verbs.RefGeneGenerator;
import edu.mit.csail.cgs.ewok.verbs.GeneToPromoter;
import edu.mit.csail.cgs.ewok.verbs.FastaWriter;
import edu.mit.csail.cgs.ewok.verbs.ChromRegionIterator;
import edu.mit.csail.cgs.tools.utils.Args;

/**
 * cat gene_names.txt | java GenePromoters --species "$SC;Sigmav7" --genes s288cMapped --upstream 500 --downstream 100
 * cat gene_names.txt | java GenePromoters --species "$SC;Sigmav7" --genes s288cMapped --upstream 500 --downstream 100 --fasta
 * java GenePromoters --species "$SC;Sigmav7" --genes s288cMapped --upstream 500 --downstream 100 --fasta --allgenes
 *
 *
 */

public class GenePromoters {

    private int upstream, downstream;
    private List<RefGeneGenerator> geneGenerators;
    private Genome genome;
    private boolean allGenes, toFasta, dontOverlapOrfs;
    private GeneToPromoter promgen;
    private FastaWriter<NamedStrandedRegion> fwriter;

    public static void main(String args[]) throws Exception {
        GenePromoters gp = new GenePromoters();
        gp.parseArgs(args);        
        gp.run();
    }
    public GenePromoters() {}
    public void parseArgs(String args[]) throws NotFoundException {
        geneGenerators = Args.parseGenes(args);
        for (RefGeneGenerator r : geneGenerators) {
            r.retrieveExons(false);
        }        
        upstream = Args.parseInteger(args,"upstream",10000);
        downstream = Args.parseInteger(args,"downstream",2000);
        genome = Args.parseGenome(args).getLast();        
        allGenes = Args.parseFlags(args).contains("allgenes");
        toFasta = Args.parseFlags(args).contains("fasta");
        dontOverlapOrfs = Args.parseFlags(args).contains("dontoverlaporfs");
        if (toFasta) {
            fwriter = new FastaWriter<NamedStrandedRegion>(System.out);
        }

        if (dontOverlapOrfs) {
            promgen = new GeneToPromoter(upstream, downstream, geneGenerators);
        } else {
            promgen = new GeneToPromoter(upstream, downstream);
        }
    }
    public void run() throws IOException {
        if (allGenes) {
            ChromRegionIterator chroms = new ChromRegionIterator(genome);
            while (chroms.hasNext()) {
                Region chrom = chroms.next();
                /* we'll use all the gene generators provided but don't want to output duplicate regions.
                   seen is keyed on 5' and contains a list of 3' ends that have already been output.
                */
                Map<Integer,List<Integer>> seen = new HashMap<Integer,List<Integer>>();
                for (RefGeneGenerator refgene : geneGenerators) {
                    Iterator<Gene> iter = refgene.execute(chrom);
                    while (iter.hasNext()) {
                        Gene g = iter.next();
                        if (seen.containsKey(g.getFivePrime()) &&
                            seen.get(g.getFivePrime()).contains(g.getThreePrime())) {
                            continue;
                        } else {
                            if (!seen.containsKey(g.getFivePrime())) {
                                seen.put(g.getFivePrime(), new ArrayList<Integer>());
                            }
                            seen.get(g.getFivePrime()).add(g.getThreePrime());
                        }
                        output(promgen.execute(g));
                    }
                }
            }
        } else {
            BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
            String line = null;
            while ((line = reader.readLine()) != null) {
                Gene g = null;
                String pieces[] = line.split("\\t");
                for (int i = 0; i < pieces.length; i++) {
                    for (RefGeneGenerator refgene : geneGenerators) {
                        Iterator<Gene> iter = refgene.byName(pieces[i]);
                        while (iter.hasNext()) {
                            if (g == null) {
                                g = iter.next();
                            } else {
                                iter.next();
                            }

                        }
                        if (g != null) {
                            break;
                        }
                    }
                    if (g != null) {
                            break;
                    }                    
                }
                if (g != null) {
                    g.setName(line);
                    output(promgen.execute(g));
                }
            }
        }
    }
    public void output(NamedStrandedRegion r) {
        if (toFasta) {
            fwriter.consume(r);
        } else {
            System.out.println(String.format("%s\t%s:%d-%d:%s",
                                             r.toString(),
                                             r.getChrom(), r.getStart(), r.getEnd(), r.getStrand()));
        }
    }

}
package edu.mit.csail.cgs.tools.core;

import java.util.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.tools.utils.Args;


/**
 * Dumps gene annotations in GFF format
 *
 */

public class DumpGeneGFF {

    public static void main(String args[]) throws Exception {
        Genome genome = Args.parseGenome(args).cdr();
        RefGeneGenerator genes = Args.parseGenes(args).get(0);
        ChromRegionIterator chroms = new ChromRegionIterator(genome);
        while (chroms.hasNext()) {
            Region chrom = chroms.next();
            Iterator<Gene> geneiter = genes.execute(chrom);
            while (geneiter.hasNext()) {
                Gene g = geneiter.next();
                System.out.println(String.format("%s\t%s\tCDS\t%d\t%d\t%d\t%s\t0\t",
                                                 g.getChrom(),
                                                 "database",
                                                 g.getStart(),
                                                 g.getEnd(),
                                                 g.getWidth(),
                                                 g.getStrand()));
                                                 

            }
        }            
    }
}
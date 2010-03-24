package edu.mit.csail.cgs.tools.rnaseq;

import java.util.*;
import java.io.IOException;
import java.sql.SQLException;
import cern.jet.random.Gamma;
import cern.jet.random.Binomial;
import cern.jet.random.engine.DRand;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.ewok.verbs.*;

/**
 * Print read counts for genes.
 *
 * Usage:
 * java DumpExpression --species "$SC;Sigmav6" --genes s288cMapped --align "Sc TAP-tagged Khd1 against 10560-6B RNA;trial 10;Novoalign trial 10 uploaded 310109"
 *
 * Output format is <gene name>\t<number of reads>
 * The number of reads reported is only reads on the same strand
 * as the gene's ORF.
 */

public class DumpExpression {

    ChipSeqLoader loader;
    List<ChipSeqAlignment> aligns;
    RefGeneGenerator genes;
    Genome genome;

    public static void main(String args[]) throws Exception {
        DumpExpression d = new DumpExpression();
        d.parseArgs(args);
        d.run();
        d.loader.close();
    }

    public DumpExpression() throws SQLException, IOException {
        loader = new ChipSeqLoader();
        aligns = new ArrayList<ChipSeqAlignment>();
    }
    public void parseArgs(String args[]) throws SQLException, NotFoundException {
        genome = Args.parseGenome(args).cdr();
        List<ChipSeqLocator> locators = Args.parseChipSeq(args,"align");
        for (ChipSeqLocator locator : locators) {
            aligns.addAll(loader.loadAlignments(locator,genome));
        }
        // parseGenes returns a list of genes; just take the first one
        genes = Args.parseGenes(args).get(0);
        if (aligns.size() == 0) {
            throw new NotFoundException("--aligns didn't match any alignments");
        }
    }

    public void run() throws SQLException, IOException {
        ChromRegionIterator chroms = new ChromRegionIterator(genome);

        while (chroms.hasNext()) {
            Region chrom = chroms.next();
            Iterator<Gene> geneiter = genes.execute(chrom);
            while (geneiter.hasNext()) {
                Gene g = geneiter.next();
                int count = loader.countByRegion(aligns,g);
                System.out.println(g.toString() + "\t" + count);
            }

        }

    }
}
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
 * DifferentialExpression --species "$SC;Sigmav6" --one "Sigma polyA RNA, haploid from 2/9/09;3/17/09;bowtie --best -m 100 -k 100" \
 *                        --two "Sigma polyA RNA, tetraploid from 2/9/09;3/17/09;bowtie --best -m 100 -k 100" --genes sgdGene \
 *                        [--bothstrands] [--byweight]
 *
 * Output columns are
 * - gene name
 * - count one
 * - count two
 * - frequency in sample count (count / total count)
 * - frequency two
 * - pval of count two given frequency one
 * - pval of count one given frequency two
 */

public class DifferentialExpression {

    ChipSeqLoader loader;
    List<ChipSeqAlignment> one, two;
    RefGeneGenerator genes;
    Genome genome;
    boolean bothstrands, byweight;

    public static void main(String args[]) throws Exception {
        DifferentialExpression d = new DifferentialExpression();
        d.parseArgs(args);
        d.run();
        d.loader.close();
    }

    public DifferentialExpression() throws SQLException, IOException {
        loader = new ChipSeqLoader();
        one = new ArrayList<ChipSeqAlignment>();
        two = new ArrayList<ChipSeqAlignment>();
    }
    public void parseArgs(String args[]) throws SQLException, NotFoundException {
        genome = Args.parseGenome(args).cdr();
        List<ChipSeqLocator> locators = Args.parseChipSeq(args,"one");
        for (ChipSeqLocator locator : locators) {
            one.addAll(loader.loadAlignments(locator));
        }
        locators = Args.parseChipSeq(args,"two");
        for (ChipSeqLocator locator : locators) {
            two.addAll(loader.loadAlignments(locator));
        }
        // parseGenes returns a list of genes; just take the first one
        genes = Args.parseGenes(args).get(0);
        if (one.size() == 0) {
            throw new NotFoundException("--one didn't match any alignments");
        }
        if (two.size() == 0) {
            throw new NotFoundException("--two didn't match any alignments");
        }
        bothstrands = Args.parseFlags(args).contains("bothstrands");
        byweight = Args.parseFlags(args).contains("byweight");
    }

    public void run() throws SQLException, IOException {
        ChromRegionIterator chroms = new ChromRegionIterator(genome);
        double totalweightone = getWeight(one), totalweighttwo = getWeight(two);
        double totalcountone, totalcounttwo;
        if (byweight) {
            totalcountone = getWeight(one);
            totalcounttwo = getWeight(two);
        } else {
            totalcountone = getCount(one);
            totalcounttwo = getCount(two);
        }

        System.err.println("Total weight one is " + totalweightone + " hits one is " + totalcountone);
        System.err.println("Total weight two is " + totalweighttwo + " hits two is " + totalcounttwo);
        
        Binomial binomial = new Binomial(100,.1,new DRand((int)(System.currentTimeMillis() % 0xFFFFFFFF)));
        
        while (chroms.hasNext()) {
            Region chrom = chroms.next();
            Iterator<Gene> geneiter = genes.execute(chrom);
            while (geneiter.hasNext()) {
                Gene g = geneiter.next();
                double countone, counttwo;
                if (bothstrands) {
                    if (byweight) {
                        countone = loader.weightByRegion(one,(Region)g);
                        counttwo = loader.weightByRegion(two,(Region)g);
                    } else {
                        countone = loader.countByRegion(one,(Region)g);
                        counttwo = loader.countByRegion(two,(Region)g);
                    }
                } else {
                    if (byweight) {
                        countone = loader.weightByRegion(one,g);
                        counttwo = loader.weightByRegion(two,g);
                    } else {
                        countone = loader.countByRegion(one,g);
                        counttwo = loader.countByRegion(two,g);
                    }
                }

                if (countone < 2 && counttwo < 2) { continue; }
                if (countone < 2) {countone = 2;}
                if (counttwo < 2) {counttwo = 2;}
                
                double pone = countone / totalcountone;
                double ptwo = counttwo / totalcounttwo;

                binomial.setNandP((int)totalcountone, ptwo);
                double cdf = binomial.cdf((int)countone);
                double pvalonegiventwo = Math.min(cdf, 1 - cdf);
                binomial.setNandP((int)totalcounttwo, pone);
                cdf = binomial.cdf((int)counttwo);
                double pvaltwogivenone = Math.min(cdf, 1 - cdf);
                System.out.println(String.format("%s\t%.0f\t%.0f\t%.4e\t%.4e\t%.4e\t%.4e", g.toString(),
                                                 countone,
                                                 counttwo,
                                                 pone, ptwo,
                                                 pvaltwogivenone,
                                                 pvalonegiventwo));                
            }

        }

    }
    private double getWeight(Collection<ChipSeqAlignment> alignments) throws SQLException, IOException {
        double weight = 0;
        for (ChipSeqAlignment a : alignments) {
            weight += loader.weighAllHits(a);
        }
        return weight;
    }
    private int getCount(Collection<ChipSeqAlignment> alignments) throws SQLException, IOException {
        int count = 0;
        for (ChipSeqAlignment a : alignments) {
            count += loader.countAllHits(a);
        }
        return count;
    }

}
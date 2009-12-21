package edu.mit.csail.cgs.tools.motifs;

/**
 * Reads a list of n motif names on STDIN (each line contains a set of one or more names for the same gene).
 * Outputs on STDOUT an n x n matrix of motif presence.  The value in each entry (row i column j) is the highest
 * score of motif i in the promoter of gene j.  An entry of zero indicates that the motif or gene couldn't
 * be found or that no motif was found.
 *
 * Input lines should contain tab-separated aliases for the same gene.  Output will be tab separated
 *
 * cat aliases.txt | java MotifOccurrenceMatrix --species "$MM;mm8" --genes refGene --genes ensGene --upstream 10000 --downstream 2000
 */

import java.io.*;
import java.util.*;
import java.sql.SQLException;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.datasets.motifs.*;
import edu.mit.csail.cgs.ewok.verbs.RefGeneGenerator;
import edu.mit.csail.cgs.tools.motifs.WeightMatrixScanner;
import edu.mit.csail.cgs.tools.utils.Args;

public class MotifOccurrenceMatrix {

    private List<List<String>> genes;
    private List<WeightMatrix> matrices;
    private double[][] bestScores;
    private int upstream, downstream;
    private List<RefGeneGenerator> geneGenerators;
    private Genome genome;
    private WeightMatrixLoader loader;

    public static void main(String args[]) throws Exception {
        MotifOccurrenceMatrix mom = new MotifOccurrenceMatrix();
        mom.parseArgs(args);
        mom.readGenes();
        mom.getMatrices();
        mom.computeMatrix();
        mom.printMatrix();
    }

    public MotifOccurrenceMatrix() {
        loader = new WeightMatrixLoader();
    }
    public void parseArgs(String args[]) throws IOException, NotFoundException {
        geneGenerators = Args.parseGenes(args);
        upstream = Args.parseInteger(args,"upstream",10000);
        downstream = Args.parseInteger(args,"downstream",2000);
        genome = Args.parseGenome(args).getLast();        
    }
    /**
     * reads the list of genes from STDIN.  Each line is one elemen of genes
     * and contains the tab-separated gene names
     */
    public void readGenes() throws IOException {
        genes = new ArrayList<List<String>>();
        BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
        String line = null;
        while ((line = reader.readLine()) != null) {
            ArrayList<String> aliases = new ArrayList<String>();
            String[] split = line.split("\\t");
            for (int i = 0; i < split.length; i++) {
                aliases.add(split[i]);
            }
            genes.add(aliases);            
        }
        reader.close();
        bestScores = new double[genes.size()][genes.size()];
    }
    /**
     * Finds a WeightMatrix for each gene by searching through the aliases
     * in order.
     *
     * It might eventually make sense to store a set of matrices and
     * then score them all and take the highest score for each target gene.
     */
    public void getMatrices() {
        List<String> types = new ArrayList<String>();
        types.add("TRANSFAC");
        types.add("MEME");
        types.add("TAMO");
        types.add(null);
        matrices = new ArrayList<WeightMatrix>();
        for (int i = 0; i < genes.size(); i++) {
            List<String> aliases = genes.get(i);
            WeightMatrix matrix = null;
            // first try the aliases listed in the input
            for (String a : aliases) {
                for (String t : types) {
                    Collection<WeightMatrix> collection = loader.query(a, null, t);
                    for (WeightMatrix m : collection) {
                        if (m != null) {
                            matrix = m;
                            break;
                        }
                    }
                    if (matrix != null) {break;}
                }
                if (matrix != null) {break;}

            }            
            if (matrix == null) {
                // now see what aliases we have in our database for the gene and try all of those
                HashSet<String> morealiases = new HashSet<String>();
                for (String a : aliases) {
                    for (RefGeneGenerator g : geneGenerators) {
                        Iterator<Gene> iter = g.byName(a);
                        while (iter.hasNext()) {
                            Gene gene = iter.next();
                            morealiases.addAll(gene.getAliases());
                        }
                    }
                }
                for (String a : morealiases) {
                    for (String t : types) {
                        Collection<WeightMatrix> collection = loader.query(a, null, t);
                        for (WeightMatrix m : collection) {
                            if (m != null) {
                                matrix = m;
                                break;
                            }
                        }
                        if (matrix != null) {break;}
                    }
                    if (matrix != null) {break;}                    
                }            
            }
            if (matrix == null) {
                System.err.println("No Matrix for " + aliases);
            }

            matrices.add(matrix);
        }
    }
    /**
     * Returns the genomic region to search for a given set of aliases.
     * Looks through the aliases in order, trying all elements of geneGenerators
     * on the first before going on to the second.
     *
     * Returns null if no promoter region can be found.
     */
    public Region getPromoterForGene(List<String> aliases) {
        for (String a : aliases) {
            for (RefGeneGenerator g : geneGenerators) {
                Iterator<Gene> iter = g.byName(a);
                if (iter.hasNext()) {
                    Gene gene = iter.next();
                    StrandedRegion sr = new StrandedRegion(gene.getGenome(),
                                                           gene.getChrom(),
                                                           gene.getFivePrime(),
                                                           gene.getFivePrime(),
                                                           gene.getStrand());
                    sr = sr.expand(upstream, downstream);
                    System.err.println(a + " -> " + sr);
                    return sr;
                }
            }
        }
        System.err.println("No promoter for " + aliases);
        return null;
    }
    public void computeMatrix() throws SQLException {
        WMHitScoreComparator comparator = new WMHitScoreComparator();
        List<Region> promoters = new ArrayList<Region>();
        for (int j = 0; j < genes.size(); j++) {
            Region promoter = getPromoterForGene(genes.get(j));
            promoters.add(promoter);
        }
        for (int i = 0; i < genes.size(); i++) {
            WeightMatrix m = matrices.get(i);
            if (m == null) {
                continue;
            }
            for (int j = 0; j < genes.size(); j++) {
                Region promoter = promoters.get(j);
                if (promoter == null) {
                    continue;
                }

                Genome.ChromosomeInfo chrinfo = genome.getChrom(promoter.getChrom());
                String seq = genome.getChromosomeSequence(chrinfo, promoter.getStart(), promoter.getEnd());
                List<WMHit> hits = WeightMatrixScanner.scanSequence(m,
                                                                    (float)(m.getMaxScore() * .75),
                                                                    seq.toCharArray());
                Collections.sort(hits,comparator);
                if (hits.size() > 0) {
                    bestScores[i][j] = hits.get(0).score;
                } else {
                    bestScores[i][j] = 0;
                }
            }
        }
    }
    public void printMatrix () {
        for (int i = 0; i < bestScores.length; i++) {
            for (int j = 0; j < bestScores[i].length - 1; j++) {
                System.out.print(bestScores[i][j] + "\t");
            } 
            System.out.print(bestScores[i][bestScores[i].length - 1] + "\n");
        }
    }

}

class WMHitScoreComparator implements Comparator<WMHit> {
    public int compare(WMHit a, WMHit b) {
        return (int)(100 * (b.score - a.score));
    }
}
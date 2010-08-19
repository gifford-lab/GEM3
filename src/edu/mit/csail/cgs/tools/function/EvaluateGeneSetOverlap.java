package edu.mit.csail.cgs.tools.function;


import java.util.*;
import java.io.*;
import java.sql.SQLException;
import cern.jet.random.Gamma;
import cern.jet.random.Binomial;
import cern.jet.random.engine.DRand;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.ewok.verbs.RefGeneGenerator;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;

/**
 * Determines whether the overlap between two sets of genes
 * is significant given the size of each set relative to
 * the set of all genes.
 *
 * java edu.mit.csail.cgs.tools.function.EvaluateGeneSetOverlap --species "$MM;mm9" --genes refGene --one gene_names_one.txt --two gene_names_two.txt
 *
 */

public class EvaluateGeneSetOverlap {

    private Genome genome;
    private Set<String> allNames, namesOne, namesTwo;
    private RefGeneGenerator refgene;    

    private double freqone, freqtwo;
    private int overlap;

    public EvaluateGeneSetOverlap(Genome g) throws SQLException {
        genome = g;
        setGeneTable("refGene");
    }
    public EvaluateGeneSetOverlap(Genome g, String genes) throws SQLException {
        genome = g;
        setGeneTable(genes);        
    }
    public void setGeneTable(String table) throws SQLException {
        refgene = new RefGeneGenerator(genome, table);
        allNames = new HashSet<String>();
        Iterator<Gene> all = refgene.getAll();
        while (all.hasNext()) {
            allNames.add(all.next().getName());
        }
    }
    public void setSetOne(Set<String> names) {
        namesOne = names;
        namesOne.retainAll(allNames);
    }
    public void setSetTwo(Set<String> names) {
        namesTwo = names;
        namesTwo.retainAll(allNames);
    }
    public void compute() {
        freqone = (double)namesOne.size() / (double)allNames.size();
        freqtwo = (double)namesTwo.size() / (double)allNames.size();
        overlap = 0;
        for (String s : namesOne) {
            if (namesTwo.contains(s)) {
                overlap++;
            }
        }
    }
    /* prob that a gene is in set one */
    public double getFreqOne() {return freqone;}
    /* prob that a gene is in set two */
    public double getFreqTwo() {return freqtwo;}
    /* gene number of genes in one and two */
    public int getOverlap() {return overlap;}

    public double getPvalOneGivenTwo() {
        Binomial binomial = new Binomial(namesOne.size(),freqtwo,new DRand((int)(System.currentTimeMillis() % 0xFFFFFFFF)));
        return 1 - binomial.cdf(overlap);
    }
    public double getPvalTwoGivenOne() {
        Binomial binomial = new Binomial(namesTwo.size(),freqone,new DRand((int)(System.currentTimeMillis() % 0xFFFFFFFF)));
        return 1 - binomial.cdf(overlap);
    }


    public static void main(String args[]) throws IOException, NotFoundException, SQLException {
        String genes = Args.parseString(args,"genes","refGene");
        EvaluateGeneSetOverlap overlap = new EvaluateGeneSetOverlap(Args.parseGenome(args).cdr(), genes);
        overlap.setSetOne(parseFile(Args.parseString(args,"one",null)));
        overlap.setSetTwo(parseFile(Args.parseString(args,"two",null)));
        overlap.compute();
   
        System.out.println(String.format("All %d, One %d, Two %d, Overlap %d",
                                         overlap.allNames.size(), overlap.namesOne.size(), 
                                         overlap.namesTwo.size(), overlap.getOverlap()));
        System.out.println(String.format("E[overlap] %.2f",
                                         overlap.getFreqOne() * overlap.getFreqTwo() * overlap.allNames.size()));
        System.out.println(String.format("pvalues %.4e, %.4e", overlap.getPvalOneGivenTwo(), 
                                         overlap.getPvalTwoGivenOne()));

    }
    public static Set<String> parseFile(String filename) throws IOException {
        HashSet<String> out = new HashSet<String>();
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        String line = null;
        while ((line = reader.readLine()) != null) {
            out.add(line);
        }
        reader.close();
        return out;
    }


}
package edu.mit.csail.cgs.tools.chipchip;
import java.util.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.datasets.chipchip.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.binding.*;

/**
 * Scans for binding (using JBD output) around genes.
 *
 * GeneJBDCalls --species "$SC;SGDv1" --genes sgdGene --jbd "analysis;version" --upstream 500 --downstream 100 
 *
 * --peaks: passes the peaks option to BayesBindingGenerator (tries to get the whole peak region?)
 * --includeORF: by default GeneJBDCalls only looks at the specified region upstream and downstream of the start site. 
 *               This option makes GeneJBDCalls look n bases upstream of the ORF through m bases downstream (ie, it
 *               looks for binding in the ORF too).
 * --noorfoverlap: don't allow the region you get from expanding upstream and downstream to overlap another ORF.
 *
 */

public class GeneJBDCalls {

    private BayesBindingGenerator bayes;
    private RefGeneGenerator genes;
    private Genome genome;
    private int upstream, downstream;
    private boolean includeORF, dontOverlapORFs;

    public static void main(String args[]) throws NotFoundException {
        GeneJBDCalls caller = new GeneJBDCalls();
        caller.parseArgs(args);
        caller.run();
    }
    public GeneJBDCalls() {}
    public void parseArgs(String args[]) throws NotFoundException {
        genome = Args.parseGenome(args).cdr();
        genes = Args.parseGenes(args).get(0);
        String exptname = Args.parseString(args,"jbd",null);
        String exptparts[] = exptname.split(";");
        ChipChipDataset dataset = new ChipChipDataset(genome);
        ChipChipBayes b = dataset.getBayes(exptparts[0], exptparts[1]);
        bayes = new BayesBindingGenerator(b,
                                          Args.parseDouble(args,"probthresh",.2),
                                          Args.parseDouble(args,"sizethresh",1),
                                          Args.parseFlags(args).contains("peaks"));        
        upstream = Args.parseInteger(args,"upstream",500);
        downstream = Args.parseInteger(args,"downstream",100);
        includeORF = Args.parseFlags(args).contains("includeORF"); // include the coding region of the gene in question
        dontOverlapORFs = Args.parseFlags(args).contains("noorfoverlap"); // don't allow the upstream or downstream values to cause overlap with another orf
    }
    public void run() {
        ChromRegionIterator chroms = new ChromRegionIterator(genome);
        while (chroms.hasNext()) {
            Region chrom = chroms.next();
            Iterator<Gene> geneiter = genes.execute(chrom);
            while (geneiter.hasNext()) {
                StringBuffer buffer = new StringBuffer();
                Gene g = geneiter.next();
                Region r = null;
                if (includeORF) {
                    r = g.expand(upstream,downstream);
                } else {
                    r = g.expand(upstream,downstream - g.getWidth());
                }
                if (dontOverlapORFs) {
                    ArrayList<Integer> positions = new ArrayList<Integer>();
                    int tss = g.getStrand() == '+' ? g.getStart() : g.getEnd();
                    positions.add(tss);
                    positions.add(r.getStart());
                    positions.add(r.getEnd());
                    Iterator<Gene> otherGenes = genes.execute(r);
                    while (otherGenes.hasNext()) {
                        Gene o = otherGenes.next();
                        if (o.getID().equals(g.getID()) || o.overlaps(g)) {
                            continue;
                        }
                        positions.add(o.getStart());
                        positions.add(o.getEnd());
                    }
                    Collections.sort(positions);
                    int i = positions.indexOf(tss);
                    if (i == 0) { i = 1;}
                    int start = positions.get(i-1);
                    int end = positions.get(i+1);
                    //                    System.err.println("  " + r + " -> " + start + "," + end);
                    r = new Region(r.getGenome(), r.getChrom(), start,end);

                }

                buffer.append(g.toString());
                Iterator<BindingExtent> bindingiter = bayes.execute(r);
                while (bindingiter.hasNext()) {
                    buffer.append("\t" + bindingiter.next().toString());
                }
                System.out.println(buffer.toString());
            }
        }                              
    }

}
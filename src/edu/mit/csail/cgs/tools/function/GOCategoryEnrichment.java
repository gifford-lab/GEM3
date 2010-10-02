package edu.mit.csail.cgs.tools.function;

import java.io.*;
import java.util.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.function.*;
import edu.mit.csail.cgs.utils.Enrichment;
import edu.mit.csail.cgs.ewok.verbs.RefGeneGenerator;
import edu.mit.csail.cgs.tools.utils.Args;

/**
 * cat genelist.txt | java GOCategoryEnrichment [--thresh .01] [--mincount 0] [--dumpgenes]
 *
 * --thresh is pvalue threshold for enriched categories to be included in output
 * --mincount is minimum number of genes in genelist-category overlap for a category to
 *            be included in output
 * --dumpgenes causes the GO category line to be followed by the list of genes 
 *             in the genelist-category overlap
 */

public class GOCategoryEnrichment {

    public static void main(String args[]) throws Exception {
        Genome genome = Args.parseGenome(args).cdr();
        FunctionLoader funcloader = new GOFunctionLoader(GOFunctionLoader.getDefaultDBName());
        FunctionalUtils funcutils = new FunctionalUtils(funcloader, genome.getSpecies());
        double thresh = Math.log(Args.parseDouble(args,"thresh",.01));
        int mincount = Args.parseInteger(args,"mincount",0);
        boolean dumpgenes = Args.parseFlags(args).contains("dumpgenes");

        Set<String> allGeneNames = new HashSet<String>();
        Map<String,Integer> chromlengths = genome.getChromLengthMap();
        List<RefGeneGenerator> genegens = Args.parseGenes(args);
        for (String cname : chromlengths.keySet()) {
            Region chrom = new Region(genome, cname, 0, chromlengths.get(cname));
            for (RefGeneGenerator generator : genegens) {
                Iterator<Gene> genes = generator.execute(chrom);
                while (genes.hasNext()) {
                    allGeneNames.add(genes.next().getName());
                }
            }
        }
        
        Set<String> foreground = new HashSet<String>();
        BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
        String line = null;
        while ((line = reader.readLine()) != null) {
            foreground.add(line);
        }
        //        System.err.println("FG " + foreground);
        //        System.err.println("BG " + allGeneNames);

        Map<String,Enrichment> output = funcutils.calculateTotalEnrichments(allGeneNames, foreground);
        ArrayList<Enrichment> enrichments = new ArrayList<Enrichment>();
        enrichments.addAll(output.values());
        Collections.sort(enrichments, new EnrichmentPvalueComparator());
        HashSet<String> overlap = new HashSet<String>();
        for (Enrichment e : enrichments) {
            if (e.getLogPValue() < thresh && e.getx() > mincount) {
                System.out.println(e.reportLine());
                if (dumpgenes) {
                    String catname = e.getCategory().replaceAll("\\t.*","");
                    Category cat = funcutils.getCategory(catname);
                    if (cat == null) {
                        System.err.println("got null cat for " + e.getCategory());
                    }

                    overlap.clear();
                    for (Assignment a : funcloader.getAssignments(cat) ) {
                        if (foreground.contains(a.getObject())) {
                            overlap.add(a.getObject().toString());
                        }
                    }
                    System.out.println(overlap);
                }
            }
        }
    }
}
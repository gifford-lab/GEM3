package edu.mit.csail.cgs.tools.function;

import java.io.*;
import java.util.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.function.*;
import edu.mit.csail.cgs.utils.Enrichment;
import edu.mit.csail.cgs.ewok.verbs.RefGeneGenerator;
import edu.mit.csail.cgs.tools.utils.Args;

public class GOCategoryEnrichment {

    public static void main(String args[]) throws Exception {
        Genome genome = Args.parseGenome(args).cdr();
        FunctionLoader funcloader = new GOFunctionLoader(GOFunctionLoader.getDefaultDBName());
        FunctionalUtils funcutils = new FunctionalUtils(funcloader, genome.getSpecies());
        double thresh = Math.log(Args.parseDouble(args,"thresh",.01));
        int mincount = Args.parseInteger(args,"mincount",0);

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
        ArrayList<String> keys = new ArrayList<String>();
        keys.addAll(output.keySet());
        Collections.sort(keys, new EnrichmentComparator(output));
        for (String c : keys) {
            Enrichment e = output.get(c);
            if (e.getLogPValue() < thresh && e.getx() > mincount) {
                System.out.println(String.format("%s\t%d/%d\t%d/%d\t%.5f", c, 
                                                 e.getx(), e.getn(), e.getTheta(), e.getN(), e.getLogPValue()));
            }
        }
    }
}

class EnrichmentComparator implements Comparator<String> {

    private Map<String, Enrichment> map;
    public EnrichmentComparator(Map<String, Enrichment> m) {
        map = m;
    }

    public int compare(String a, String b) {
        return Double.compare(map.get(a).getLogPValue(), map.get(b).getLogPValue());
    }

}
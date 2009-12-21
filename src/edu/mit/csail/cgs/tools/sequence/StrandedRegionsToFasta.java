package edu.mit.csail.cgs.tools.sequence;

import java.io.*;

import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.ewok.verbs.StrandedRegionParser;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;

/* reads a list of stranded regions (one per line) on STDIN.
   Produces on stdout a fasta file containing those regions.
   The genome is specified on the command line as
   --species "Mm;mm8"
*/

public class StrandedRegionsToFasta {

    public static void main(String args[]) {
        try {
            Pair<Organism,Genome> pair = Args.parseGenome(args);
            Organism organism = pair.car();
            Genome genome = pair.cdr();
            StrandedRegionParser parser = new StrandedRegionParser(genome);
            SequenceGenerator seqgen = new SequenceGenerator();
            boolean cache = Args.parseFlags(args).contains("cache");
            seqgen.useCache(cache);
            String line;
            BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                StrandedRegion r = parser.execute(line);
                String seq = seqgen.execute(r);
                if(r.getStrand()=='-'){
                	seq = SequenceUtils.reverseComplement(seq);
                }
                System.out.println(">"+r.toString()+"\n"+seq);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

}

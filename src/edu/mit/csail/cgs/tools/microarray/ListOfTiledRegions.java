package edu.mit.csail.cgs.tools.microarray;

import edu.mit.csail.cgs.ewok.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.Pair;

import java.util.*;
import java.io.*;

/* Prints a list of the regions tiled by the specified array design.

   usage:
   java edu.mit.csail.cgs.tools.microarray.FastaOfTiledRegions \
    --species "$SC;SGDv1"
    --design "Sc 244k"
    [--spacing 300] [--mincount 8]
 */

public class ListOfTiledRegions {
    
    public static void main (String args[]) {
        String designName = Args.parseString(args,"design",null);
        int spacing = Args.parseInteger(args,"spacing",200);
        int mincount = Args.parseInteger(args,"mincount",10);

        try {
            Pair<Organism,Genome> pair = Args.parseGenome(args);
            Organism organism = pair.car();
            Genome genome = pair.cdr();       
            TiledRegionGenerator gen = new TiledRegionGenerator(designName,spacing,mincount);
            List<String> chromnames = genome.getChromList();
            for (int i = 0; i < chromnames.size(); i++) {
                Region region = new Region(genome,chromnames.get(i),1,genome.getChromLength(chromnames.get(i)));
                Iterator<Region> tiled = gen.execute(region);
                while (tiled.hasNext()) {
                    Region r = tiled.next();
                    System.out.println(r.toString());
                }
            }
        
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}

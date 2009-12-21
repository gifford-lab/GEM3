package edu.mit.csail.cgs.tools.alignments;

import java.util.List;
import java.util.Iterator;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.alignments.*;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.tools.utils.Args;

public class PrintIndels {

    public static void main(String args[]) throws NotFoundException {
        int splitPenalty, chromPenalty;
        splitPenalty = Args.parseInteger(args,"splitpenalty",-1);
        chromPenalty = Args.parseInteger(args,"chrompenalty",-1);
        String s = Args.parseString(args,"one",null);
        if (s == null) {
            System.err.println("Must supply two genome names as command line args as --one and --two");
            System.exit(1);
        }
        String g1[] = s.split(";");
        s = Args.parseString(args,"two",null);
        if (s == null) {
            System.err.println("Must supply two genome names as command line args as --one and --two");
            System.exit(1);
        }
        String g2[] = s.split(";");
        Genome one = new Genome(g1[0],g1[1]);
        Genome two = new Genome(g2[0],g2[1]);
        String alignType = Args.parseString(args,"type","blast");
        IndelGenerator generator = new IndelGenerator(two);
        generator.setAlignPrefix(alignType);
        if (splitPenalty > 0) {
            generator.setSplitPenalty(splitPenalty);
        }
        if (chromPenalty > 0) {
            generator.setChromPenalty(chromPenalty);
        }
        
        List<String> chromnames = one.getChromList();
        for (String name : chromnames) {

            try {
                System.err.println(" ======== " + name + " ========\n");
                Region chrom = new Region(one,
                                          name,
                                          0,
                                          one.getChromLength(name));
                Iterator<Indel> indels = generator.execute(chrom);
                while (indels.hasNext()) {
                    Indel i = indels.next();
                    System.out.println(String.format("%s\t%d\t%s",
                                                     i.regionString(),
                                                     i.getSize(),
                                                     i.getType()));
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }
}

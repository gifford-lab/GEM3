package edu.mit.csail.cgs.tools.alignments;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.mit.csail.cgs.datasets.alignments.AlignmentStitcher;
import edu.mit.csail.cgs.datasets.alignments.MultiZAlignRegion;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.parsing.alignment.PSL;
import edu.mit.csail.cgs.utils.parsing.alignment.PSLHit;

/**
 * Reads PSL formatted input on STDIN.  Prints, on STDOUT,
 * the set of input lines corresponding to the best single
 * alignment of the query to the target.  Each base in the 
 * query sequence will be present at most once in the set of
 * output lines.  
 */
public class BestBlatAlignment {

    public static void main(String args[]) throws Exception {
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
        
        PSL psl = new PSL(new BufferedReader(new InputStreamReader(System.in)));
        Map<String,List<MZARS>> regions = new HashMap<String,List<MZARS>>();
        while (psl.hasNext()) {
            PSLHit hit = psl.next();
            if (!regions.containsKey(hit.qname)) {
                regions.put(hit.qname, new ArrayList<MZARS>());
            }

            regions.get(hit.qname).add(new MZARS(hit,
                                                 one, hit.qname, hit.qstart, hit.qend,
                                                 two, hit.tname, hit.tstart, hit.tend,
                                                 hit.match - hit.mismatch, hit.strand));
        }
        for (String q : regions.keySet()) {
            AlignmentStitcher<MZARS> stitcher = new AlignmentStitcher<MZARS>(regions.get(q));
            stitcher.stitch();
            List<MZARS> alignedto = stitcher.getBestAlignment(q);
            for (MZARS m : alignedto) {
                System.out.println(m.getHit().toString());
            }
        }
    }

}

class MZARS extends MultiZAlignRegion {

    private PSLHit hit;
    public MZARS (PSLHit hit,
                  Genome baseGenome, String baseChrom, int baseStart, int baseEnd,
                  Genome otherGenome, String otherChrom, int otherStart, int otherEnd,
                  double score, char strand
                  ) {
        super(baseGenome, baseChrom, baseStart, baseEnd,
              otherGenome, otherChrom, otherStart, otherEnd, score, strand);
        this.hit = hit;       
    }
    public PSLHit getHit() { return hit;}

}
package edu.mit.csail.cgs.tools.chipseq;

import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.chipseq.*;
import java.net.*;
import java.util.*;
import java.io.*;
import java.sql.SQLException;

/**
 * Test program to compare readdb to our Oracle chipseq setup.
 * Reads chrom:start-stop:strand values from stdin.  Repeats them
 * on stdout along with the hit positions
 *
 * Usage:
 * QueryOracle --species "$SC;Sigmav6" --align "Sigma polyA RNA, tetraploid from 2/16/09;3/17/09;bowtie --best -m 100 -k 100" [--quiet]
 */

public class QueryOracle {

    private ChipSeqLoader loader;
    private List<ChipSeqAlignment> alignments;
    private Genome genome;
    private boolean quiet;

    public static void main(String args[]) throws Exception {
        QueryOracle query = new QueryOracle();
        query.parseArgs(args);
        query.run(System.in);
    }

    public void parseArgs(String args[]) throws NotFoundException, SQLException, IOException {        
        quiet = Args.parseFlags(args).contains("quiet");
        genome = Args.parseGenome(args).cdr();
        loader = new ChipSeqLoader();
        List<ChipSeqLocator> locators = Args.parseChipSeq(args,"align");
        alignments = new ArrayList<ChipSeqAlignment>();
        for (ChipSeqLocator l : locators) {
            alignments.addAll(loader.loadAlignments(l,genome));
        }        
    }
    public void run(InputStream instream) throws IOException, SQLException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(instream));
        String line;
        while ((line = reader.readLine()) != null) {
            if (!quiet) { System.out.println(line); }


            try {
                StrandedRegion r = (StrandedRegion)Region.fromString(genome,line);
                Collection<ChipSeqHit> hits = loader.loadByRegion(alignments,r);
                if (!quiet) {
                    for (ChipSeqHit h : hits) {
                        if (h.getStrand() == r.getStrand()) {
                            System.out.println(h.getMidpoint().getLocation());
                        }
                    }
                }
            } catch (Exception e) {
                e.printStackTrace();
            }            
        }
    }

}
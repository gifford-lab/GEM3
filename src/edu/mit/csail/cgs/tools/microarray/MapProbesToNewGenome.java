package edu.mit.csail.cgs.tools.microarray;

import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.utils.io.parsing.alignment.PSL;
import edu.mit.csail.cgs.utils.io.parsing.alignment.PSLHit;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.datasets.chipchip.*;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.Pair;

import java.util.*;
import java.sql.*;
import java.io.*;

/**
 * Maps probes to a new genome based on a PSL file.
 * Reads, on STDIN, a PSL file in which the query sequence names
 * were the DBIDs of the probes.  Adds one position in the specified
 * genome for each line of input (so do your own input filtering).  The
 * input must be sorted such that all hits for a probe are on consecutive lines.
 *
 */

public class MapProbesToNewGenome {

    private static PreparedStatement insert;
    private static Genome genome;
    private static Map<String,Integer> chrommap;
    private static Map<String,Integer> chromids;

    public static void main(String args[]) throws Exception {
        genome = Args.parseGenome(args).cdr();
        chrommap = new HashMap<String,Integer>();
        chromids = genome.getChromIDMap();
		java.sql.Connection c = 
			DatabaseFactory.getConnection("chipchip");
        insert = c.prepareStatement("insert into probelocation(id, chromosome, startpos, stoppos, strand, loccount, bitscore) " +
                                                    " values (?,?,?,?,?,?,?)");
        PSL pslfile = new PSL(new BufferedReader(new InputStreamReader(System.in)));
        ArrayList<PSLHit> hits = new ArrayList<PSLHit>();
        String lastprobe = null;
        while (pslfile.hasNext()) {
            PSLHit hit = pslfile.next();
            if (lastprobe != null &&
                !lastprobe.equals(hit.qname)) {
                addHits(hits);
                hits.clear();
            }
            hits.add(hit);
            lastprobe = hit.qname;
        }
        addHits(hits);
        c.commit();
        c.close();
    }
    private static void addHits(Collection<PSLHit> hits) throws SQLException {
        for (PSLHit hit : hits) {
            if (!chrommap.containsKey(hit.tname)) {
                int id;
                String newchrom = hit.tname.replaceFirst("\\.fa.{0,3}$","").replaceFirst("^chr","");
                if (newchrom.equals("Mito")) {newchrom = "mt";}
                if (chromids.containsKey(newchrom)) {
                    id = chromids.get(newchrom);
                    chrommap.put(hit.tname,id);
                } else {
                    for (String s : chromids.keySet()) {
                        System.err.println("Knew about " + s);
                    }
                    throw new IllegalArgumentException("Can't figure out " + hit.tname + "," + newchrom);
                }
            }
            try {
                int chromid = chrommap.get(hit.tname);
                insert.setInt(1, Integer.parseInt(hit.qname));
                insert.setInt(2,chromid);
                insert.setInt(3,hit.tstart);
                insert.setInt(4,hit.tend);
                insert.setString(5,Character.toString(hit.strand));
                insert.setInt(6,hits.size());
                insert.setDouble(7,2*hit.match - hit.mismatch);
                insert.execute();
            } catch (SQLException e) {
                System.err.println("Hit " + hit.qname);
                e.printStackTrace();
            }
        }

    }

}
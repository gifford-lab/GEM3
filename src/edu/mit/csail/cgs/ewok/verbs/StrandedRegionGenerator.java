package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;
import java.sql.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.database.*;

/**
 *  Generator that returns StrandedRegion objects from the specified
 *  table.  The table must have columns: chrom, chromStart, chromEnd,
 *  name, strand.  The regions are returned sorted by ascending
 *  chromStart.  StrandedRegionGenerator is intended for use with the
 * MySQL UCSC annotation tables.
 *
 * It might make sense to generify the output parameter from StrandedRegion to 
 * Y extends StrandedRegion so this class can be more easily extended.  If you do that,
 * also:
 * 1) write a method that returns the SQL statement (or maybe just the list of columns)
 * 2) write a metho that parses a ResultSet into a Y
*/

public class StrandedRegionGenerator<X extends Region> implements Expander<X,StrandedRegion> {

    private Genome genome;
    private String tablename;

    public StrandedRegionGenerator(Genome g, String t) {
        genome = g;
        tablename = t;
    }

    public Iterator<StrandedRegion> execute(X region) {
        try {
            java.sql.Connection cxn =
                genome.getUcscConnection();
            PreparedStatement ps = cxn.prepareStatement("select strand, chromStart, chromEnd from " + tablename + " where chrom = ? and " +
                                                        "((chromStart <= ? and chromEnd >= ?) or (chromStart >= ? and chromStart <= ?)) order by chromStart");
            String chr = region.getChrom();
            if (!chr.matches("^(chr|scaffold).*")) {
                chr = "chr" + chr;
            }
            ps.setString(1,chr);
            ps.setInt(2,region.getStart());
            ps.setInt(3,region.getStart());
            ps.setInt(4,region.getStart());
            ps.setInt(5,region.getEnd());
            ResultSet rs = ps.executeQuery();
            ArrayList<StrandedRegion> results = new ArrayList<StrandedRegion>();
            while (rs.next()) {
                results.add(new StrandedRegion(genome,
                                               region.getChrom(),
                                               rs.getInt(2),
                                               rs.getInt(3),
                                               rs.getString(1).charAt(0)));
                                     
            }
            rs.close();
            ps.close();
            DatabaseFactory.freeConnection(cxn);
            return results.iterator();
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't get UCSC RefGenes",ex);
        }

    }


}

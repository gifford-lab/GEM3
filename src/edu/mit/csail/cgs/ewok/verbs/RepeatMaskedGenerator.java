package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;
import java.sql.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.RepeatMaskedRegion;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.database.*;

/**
 * Maps a Region to RepeatMaskedRegions from a UCSC annotations table.
 */
public class RepeatMaskedGenerator<X extends Region> implements Expander<X,RepeatMaskedRegion> {

    private Genome genome;
    private boolean onetable;
    
    public RepeatMaskedGenerator(Genome g) {
        genome = g;
        onetable = false;
    }

    public Iterator<RepeatMaskedRegion> execute(X region) {
        try {
            java.sql.Connection cxn =
                genome.getUcscConnection();
            String chr = region.getChrom();
            if (!chr.matches("^(chr|scaffold).*")) {
                chr = "chr" + chr;
            }
            String tablename = onetable ? "rmsk" : (chr + "_rmsk");
            PreparedStatement ps = cxn.prepareStatement("select genoStart, genoEnd, strand, repName, repClass, repFamily, swScore from " + tablename + 
                                                        " where genoName = ? and ((genoStart <= ? and genoEnd >= ?) or (genoStart >= ? and genoStart <= ?)) order by genoStart");
            ps.setString(1,chr);
            ps.setInt(2,region.getStart());
            ps.setInt(3,region.getStart());
            ps.setInt(4,region.getStart());
            ps.setInt(5,region.getEnd());
            ResultSet rs = ps.executeQuery();
            ArrayList<RepeatMaskedRegion> results = new ArrayList<RepeatMaskedRegion>();
            while (rs.next()) {
                RepeatMaskedRegion r = new RepeatMaskedRegion(genome,region.getChrom(),
                                                              rs.getInt(1),
                                                              rs.getInt(2),
                                                              rs.getString(4),
                                                              rs.getString(5),
                                                              rs.getString(6),
                                                              (double)rs.getInt(7),
                                                              rs.getString(3).charAt(0));
                results.add(r);
            }
            rs.close();
            ps.close();
            DatabaseFactory.freeConnection(cxn);
            return results.iterator();
        } catch (SQLException ex) {
            ex.printStackTrace();
            if (!onetable) {
                onetable = true;
                return execute(region);
            } else {
                throw new DatabaseException("Couldn't get UCSC repeat masked regions",ex);
            }
        }
    }
}


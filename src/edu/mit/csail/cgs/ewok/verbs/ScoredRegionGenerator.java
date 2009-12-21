package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;
import java.sql.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.ScoredRegion;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.database.*;

/**
 * ScoreDRegionGenerator parses rows from MySQL UCSC style annotation tables
 * into ScoredRegions.    StrandedRegionGenerator provides a similar function for
 * a slightly different table schema (and has further useful documentation).
 */
public class ScoredRegionGenerator<X extends Region> implements Expander<X,ScoredRegion> {

    private Genome genome;
    private String tablename;    
    
    public ScoredRegionGenerator(Genome g, String tablename) {
        genome = g;
        this.tablename = tablename;
    }

    public void setTablename(String t) {tablename = t;}
    public Genome getGenome() {return genome;}
    public String getTablename() {return tablename;}

    public Iterator<ScoredRegion> execute(X region) {
        try {
            java.sql.Connection cxn =
                getGenome().getUcscConnection();
            PreparedStatement ps = cxn.prepareStatement("select chromStart, chromEnd, score from " + getTablename() + " where chrom = ? and " +
                                                        "((chromStart <= ? and chromEnd >= ?) or (chromStart >= ? and chromStart <= ?)) order by chromStart");
            ps.setFetchSize(1000);
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
            ArrayList<ScoredRegion> results = new ArrayList<ScoredRegion>();
            while (rs.next()) {
                ScoredRegion r = new ScoredRegion(genome,
                                                  region.getChrom(),
                                                  rs.getInt(1),
                                                  rs.getInt(2),
                                                  rs.getFloat(3));
                results.add(r);
            }
            rs.close();
            ps.close();
            DatabaseFactory.freeConnection(cxn);
            return results.iterator();
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't get UCSC MultiZ for " + tablename,ex);
        }
    }
}


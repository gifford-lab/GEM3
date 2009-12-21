package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;
import java.sql.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.general.NamedTypedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.database.*;

/** 
 * Generator that returns NamedTypedRegion objects from the specified table. 
 * The table must have columns: chrom, chromStart, chromEnd, name, strand, and type .  
 * This parses tables like the UCSC sgdOther table from the sacCer1 annotation database.
 */

public class SGDOtherGenerator<X extends Region> implements Expander<X,NamedTypedRegion> {

    private Genome genome;
    private String tablename;

    public SGDOtherGenerator(Genome g, String t) {
        genome = g;
        tablename = t;
    }

    public Iterator<NamedTypedRegion> byName(String name) {
        try {
            java.sql.Connection cxn =
                genome.getUcscConnection();
            PreparedStatement ps = cxn.prepareStatement("select name, chrom, strand, type, chromStart, chromEnd from " + tablename + " where name = ?");
            ps.setString(1,name);
            Iterator<NamedTypedRegion> results = parseResults(ps);
            DatabaseFactory.freeConnection(cxn);
            return results;
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't get UCSC RefGenes",ex);
        }
    }


    public Iterator<NamedTypedRegion> execute(X region) {
        try {
            java.sql.Connection cxn =
                genome.getUcscConnection();
            PreparedStatement ps = cxn.prepareStatement("select name, chrom, strand, type, chromStart, chromEnd from " + tablename + " where chrom = ? and " +
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
            Iterator<NamedTypedRegion> results = parseResults(ps);
            DatabaseFactory.freeConnection(cxn);
            return results;
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't get UCSC SGDGenes",ex);
        }
    }

    private Iterator<NamedTypedRegion> parseResults(PreparedStatement ps) throws SQLException {
        try {
            ResultSet rs = ps.executeQuery();
            ArrayList<NamedTypedRegion> results = new ArrayList<NamedTypedRegion>();        
            while (rs.next()) {
                String chr = rs.getString(2);
                chr = chr.replaceFirst("^(chr|scaffold)","");
                String strand = rs.getString(3);
                results.add(new NamedTypedRegion(genome,
                                                 chr,
                                                 rs.getInt(5),
                                                 rs.getInt(6),
                                                 rs.getString(1),
                                                 rs.getString(4),
                                                 strand.length() > 0 ? strand.charAt(0) : ' '));
            }
            rs.close();
            ps.close();
            return results.iterator();
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't get UCSC RefGenes",ex);
        }
    }       
}

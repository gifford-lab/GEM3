/*
 * Created on Sep 12, 2006
 */
package edu.mit.csail.cgs.datasets.chippet;

import java.util.*;
import java.io.*;
import java.sql.*;

import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.datasets.species.Genome;

/**
 * @author tdanford
 */
public class ChipPetDatum extends StrandedRegion {
    
    private ChipPetExpt expt;
    private int count;

    public ChipPetDatum(ResultSet rs, Genome g, ChipPetExpt e) throws SQLException {
        super(g, g.getChromName(rs.getInt(2)), rs.getInt(3), rs.getInt(4), rs.getString(5).charAt(0));
        int exptID = rs.getInt(1);
        
        if(exptID != e.getDBID()) { throw new IllegalArgumentException(exptID + " != " + e.getDBID()); }
        expt = e;
        count = rs.getInt(6);
    }
    
    // for use, in a slightly-optimized loader in the ChipPetLoader class.
    ChipPetDatum(Genome g, String chrom, int start, int end, char str, ChipPetExpt expt, int ct) { 
        super(g, chrom, start, end, str);
        this.expt = expt;
        count = ct;
    }
    
    public ChipPetExpt getExpt() { return expt; }
    public int getCount() { return count; }
    public String toString() { return getLocationString() + " {" + expt.toString() + "}"; } 

    public int hashCode() { 
        int code = 17;
        code += expt.hashCode(); code *= 37;
        code += count; code *= 37;
        code += super.hashCode(); code *= 37;
        return code; 
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof ChipPetDatum)) { return false; }
        ChipPetDatum d = (ChipPetDatum)o;
        if(!expt.equals(d.expt)) { return false; }
        if(count != d.count) { return false; }
        if(!super.equals(d)) { return false; }
        return true;
    }

    public static PreparedStatement prepareLoadByExptAndGenome(java.sql.Connection c) throws SQLException { 
        return c.prepareStatement("select d.expt, d.chromosome, d.startpos, d.stoppos, d.strand, d.peakOverlap from chippetdata d, " +
                "chromosome c, genome g where d.chromosome = c.id and " +
                "c.genome = g.id and d.expt=? and g.id=?");
    }
    
    public static PreparedStatement createLoadByRegion(java.sql.Connection c) throws SQLException { 
        return c.prepareStatement("select d.expt, d.chromosome, d.startpos, d.stoppos, d.strand, d.peakOverlap from chippetdata d, " +
        "chromosome c, genome g where d.chromosome = c.id and " +
        "c.genome = g.id and d.expt=? and g.id=? and d.chromosome=? " + 
        "and ((d.startpos >= ? and d.startpos <= ?) or (d.startpos <= ? and d.stoppos >= ?))");
    }
    
    public static PreparedStatement createInsertStatement(java.sql.Connection c) throws SQLException { 
        return c.prepareStatement("insert into chippetdata " +
        		"(expt, chromosome, startpos, stoppos, strand, peakOverlap) " +
        		"values (?, ?, ?, ?, ?, ?)");
    }
    
}

package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;
import java.sql.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.database.*;

public class ByteSequenceGenerator<X extends Region> implements Mapper<X,Byte[]> {

    String tablename, keycol, valcol;

    // no longer used, but kept for compatibility 
    public ByteSequenceGenerator (Genome g, String table) {
        tablename = table;
        keycol = "chromosome";
        valcol = "quality";
    }
    public ByteSequenceGenerator (String table) {
        tablename = table;
        keycol = "chromosome";
        valcol = "quality";
    }
    public ByteSequenceGenerator(String table, String key, String val) {
        tablename = table;
        keycol = key;
        valcol = val;
    }

    public Byte[] execute(X region) {
        Byte[] result = null;
        String chromname = region.getChrom();
        if (!chromname.matches("^chr") && !chromname.matches("^scaffold")) {
            chromname = "chr" + chromname;
        }

        try {
            Genome genome = region.getGenome();
            java.sql.Connection cxn = genome.getUcscConnection();
            PreparedStatement ps = cxn.prepareStatement("select substr(" + valcol + ",?,?) from " + tablename + 
                                                        " where " + keycol + " = ?");
            int start = Math.max(region.getStart() - 1,0);
            ps.setInt(1,start);
            ps.setInt(2,region.getEnd() - region.getStart() + 1);
            ps.setString(3,chromname);
            ResultSet rs = ps.executeQuery();
            if (rs.next()) {
                byte[] row = rs.getBytes(1);
                result = new Byte[row.length];
                for (int i = 0; i < row.length; i++) {
                    result[i] = row[i];
                }
            } 
            rs.close();
            ps.close();
            DatabaseFactory.freeConnection(cxn);
        } catch (SQLException ex) {
            ex.printStackTrace();           
        } catch (UnknownRoleException ex) {
            ex.printStackTrace();
            throw new DatabaseException("Couldn't connect to core",ex);
        }
        if (result.length != region.getWidth()) {
            throw new RuntimeException("Wanted " + region + " but only got " + result.length);
        }

        return result;
    }

}

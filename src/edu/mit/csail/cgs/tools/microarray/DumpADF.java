package edu.mit.csail.cgs.tools.microarray;

import java.util.*;
import java.sql.*;
import java.io.*;
import edu.mit.csail.cgs.utils.database.*;

public class DumpADF {

    public static void main(String args[]) throws Exception {
        java.sql.Connection cxn = DatabaseFactory.getConnection("chipchip");
        String designname = args[0];
        PreparedStatement getDesignID = cxn.prepareStatement("select id from arraydesign where name = ?");
        getDesignID.setString(1,designname);
        ResultSet rs = getDesignID.executeQuery();
        rs.next();
        int designid = rs.getInt(1);

        PreparedStatement getProbes = cxn.prepareStatement("select blockno, colno, rowno, probeid, probename, sequence, type, galfiles.name " +
                                                           " from probedesign, galfiles where galfiles.id = probedesign.galfile and " +
                                                           " probedesign.arraydesign = ? order by blockno, colno, rowno");
        getProbes.setInt(1, designid);
        rs = getProbes.executeQuery();

        Map<String,PrintWriter> pws = new HashMap<String,PrintWriter>();

        while (rs.next()) {
            String galfile = rs.getString(8);
            String fname = designname + "_" + galfile + ".adf";
            PrintWriter pw = null;
            if (pws.containsKey(fname)) {
                pw = pws.get(fname);
            } else {
                pw = new PrintWriter(fname);
                pw.println("MetaColumn\tMetaRow\tColumn\tRow\tReporter Identifier\tReporter Name\tReporter BioSequence [Actual Sequence]\tReporter BioSequence Type\tReporter BioSequence Polymer Type\tReporter Group [role]\tReporter Control Type");
                pws.put(fname,pw);
            }
            String group = rs.getString(7);
            String type = "";
            if (group.equals("false")) {
                group = "Experimental";                                    
            } else {
                group = "Control";
                type = "control_biosequence";
            }

            String name = rs.getString(5);
            if (name == null) {
                name = rs.getString(4);
            }
            String sequence = rs.getString(6);
            if (sequence == null) {
                sequence = "";
            }


            pw.println(String.format("%d\t%d\t%d\t%d\t%s\t%s\t%s\tss_oligo\tDNA\t%s\t%s",
                                     1, 
                                     rs.getInt(1),
                                     rs.getInt(2),
                                     rs.getInt(3),
                                     rs.getString(4),
                                     name,
                                     sequence,
                                     group,
                                     type));            
        }
        for (String f : pws.keySet()) {
            pws.get(f).close();
        }


        
        
    }




}
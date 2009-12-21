package edu.mit.csail.cgs.tools.microarray;

import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.datasets.chipchip.*;
import edu.mit.csail.cgs.tools.utils.Args;

import java.util.*;
import java.sql.*;

/**
 * Generates FASTA output on stdout of probe sequences.  The sequence
 * names are the probe DBIDs, suitable for later import when mapped to
 * a new genome.  Specify either a designname or a species name or you'll
 * get *all* the probes in the DB.
 *
 * FastaOfProbeSequences [--design 'Hox Array'] [--species 'Mus Musculus']
 */

public class DumpProbeDesign {

    public static void main(String args[]) throws Exception {
        String designName = Args.parseString(args,"design",null);
        String speciesName = Args.parseString(args,"species",null);
        StringBuffer query = new StringBuffer("select id, probename, probeid, type, sequence from probedesign");
        if (speciesName != null) {
            ChipChipMetadataLoader loader = new ChipChipMetadataLoader();
            ArrayList<Integer> dbids = new ArrayList<Integer>();
            for (ArrayDesign d : loader.loadAllArrayDesigns()) {
                if (d.getGenome().getSpecies().equals(speciesName)) {
                    dbids.add(d.getDBID());
                }
            }
            query.append(" where arraydesign in (");
            for (int i = 0; i < dbids.size(); i++) {
                if (i == 0) {
                    query.append(dbids.get(i));
                } else {
                    query.append("," + dbids.get(i));
                }
            }
            query.append(")");
        }
        if (designName != null) {
            ChipChipMetadataLoader loader = new ChipChipMetadataLoader();
            ArrayDesign design = null;
            for (ArrayDesign d : loader.loadAllArrayDesigns()) {
                if (d.getName().equals(designName)) {
                    design = d;
                    break;
                }
            }
            if (design == null) {
                System.err.println("Couldn't find array design with name " + designName);
            } else {
                query.append(" where arraydesign = " + design.getDBID());
            }
        }
		java.sql.Connection c = 
			DatabaseFactory.getConnection("chipchip");
        Statement stmt = c.createStatement();
        System.err.println(query);
        ResultSet rs = stmt.executeQuery(query.toString());
        while (rs.next()) {
            System.out.println(String.format("%d\t%s\t%s\t%s\t%s",
                                             rs.getInt(1),
                                             rs.getString(2),
                                             rs.getString(3),
                                             rs.getString(4),
                                             rs.getString(5)));
        }
        rs.close();
        stmt.close();

    }

}
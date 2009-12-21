package edu.mit.csail.cgs.tools.chipchip;

import java.sql.*;
import java.util.*;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.chipchip.*;
import edu.mit.csail.cgs.tools.utils.Args;

/**
 * Dumps a tab separated file with the following columns
 * - probe sequence
 * - experiment id
 * - cy5
 * - cy3
 *
 * command line args are --species "$SC;SGDv1" --expt "name;version;replicate"
 */
public class DumpProbesIntensities {

    public static void main(String args[]) throws Exception {
        List<ExptNameVersion> envs = Args.parseENV(args);
        ChipChipMetadataLoader metadata = new ChipChipMetadataLoader();
        Collection<Experiment> expts = new ArrayList<Experiment>();
        for (ExptNameVersion env : envs) {
            if (env.getReplicate() == null) {
                expts.addAll(metadata.loadExperiment(env.getName(), env.getVersion()));
            } else {
                expts.add(metadata.loadExperiment(env.getName(), env.getVersion(), env.getReplicate()));
            }
        }
        ArrayList<Integer> exptids = new ArrayList<Integer>();
        for (Experiment e : expts) {
            exptids.add(e.getDBID());
        }
        java.sql.Connection cxn = DatabaseFactory.getConnection("chipchip");
        Statement stmt = cxn.createStatement();
        String sql = "select probedesign .sequence, data.experiment, data.channelone, data.channeltwo from probedesign, data " +
                      "where data.probe = probedesign.id and data.experiment in (";
        for (int i = 0; i < exptids.size(); i++) {
            sql += (i == 0 ? exptids.get(i) : (", " + exptids.get(i)));
        }
        sql += ") order by data.experiment, data.probe";
        ResultSet rs = stmt.executeQuery(sql);
        while (rs.next()) {
            System.out.println(rs.getString(1) + "\t" +
                               rs.getInt(2) + "\t" +
                               rs.getDouble(3) + "\t" +
                               rs.getDouble(4));
        }
        rs.close();
        stmt.close();
    }                   
}
package edu.mit.csail.cgs.tools.chipchip;

import java.sql.*;
import java.util.*;
import java.io.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.chipchip.*;
import edu.mit.csail.cgs.tools.chipchip.NormalizeExperiment;

/**
 * Control experiment normalization 
 * generates a new experiment in which the ratios
 * are the ratios in the old experiment minus the ratio
 * in the control experiment + 1:
 * new = old - control + 1
 *
 * This is useful if you think the experiment was a mix of signals
 * eg, specific + nonspecific, and the control experiment measures
 * only the nonspecific part.
 */


public class ControlExptNormalization extends NormalizeExperiment {

    private String controlname, controlversion, controlreplicate;
    private int oldexptid, newexptid;
    private Map<Integer,Double> ipvals, wcevals;
    private java.sql.Connection cxn;
    private double ipmean = 0, wcemean = 0, ratiomean = 0;
    private boolean normratio = false, multiplicative = true, usewce = true, usebothchannels = false;

    public ControlExptNormalization (Experiment oldexpt,
                                     Experiment newexpt,
                                     String cn, String cv, String cr,
                                     boolean normratio,
                                     boolean multiplicative,
                                     boolean usewce,
                                     boolean usebothchannels) throws SQLException, NotFoundException, IOException {
        this.normratio = normratio;
        this.multiplicative = multiplicative;
        this.usewce = usewce;
        this.usebothchannels = usebothchannels;
        ipvals = new HashMap<Integer,Double>();
        wcevals = new HashMap<Integer,Double>();
        ChipChipMetadataLoader loader = new ChipChipMetadataLoader();
        oldexptid = oldexpt.getDBID();
        newexptid = newexpt.getDBID();        
        String sql;
        ResultSet rs;
        cxn = DatabaseFactory.getConnection("chipchip");

        if (cr != null) {            
            Experiment controlexpt = loader.loadExperiment(cn,cv,cr);
            sql = "select probe, avg(channelone), avg(channeltwo) from data where experiment = ? group by probe";
            PreparedStatement ps = cxn.prepareStatement(sql);
            ps.setInt(1,controlexpt.getDBID());
            rs = ps.executeQuery();
        } else {
            Collection<Experiment> controlexpts = loader.loadExperiment(cn,cv);
            Iterator<Experiment> iter = controlexpts.iterator();
            String exptstring = ((Integer)(iter.next()).getDBID()).toString();
            while (iter.hasNext()) {
                exptstring += "," + iter.next().getDBID();
            }
            sql = "select probe, avg(channelone), avg(channeltwo) from data where experiment in (" + exptstring + ") group by probe";
            Statement stmt = cxn.createStatement();
            rs = stmt.executeQuery(exptstring);
        }
        double ipsum = 0, wcesum = 0, ratiosum = 0;
        int count = 0;
        while (rs.next()) {
            ipvals.put(rs.getInt(1),
                       rs.getDouble(2));
            wcevals.put(rs.getInt(1),
                        rs.getDouble(3));
            ipsum += rs.getDouble(2);
            wcesum += rs.getDouble(2);
            ratiosum += rs.getDouble(2) / rs.getDouble(3);
            count++;
        }
        ipmean = ipsum / count;
        wcemean = wcesum / count;
        ratiomean = ratiosum / count;
        rs.close();
    }

    public void doNorm() throws SQLException {
        cxn.setAutoCommit(false);
        PreparedStatement getdata = cxn.prepareStatement("select data.probe, data.channelone, data.channeltwo from data where " +
                                                         " data.experiment = " + oldexptid);
        PreparedStatement insertdata = cxn.prepareStatement("insert into data (experiment,probe,channelone,channeltwo,mor,channelratio,ratio) values (" + newexptid + ",?,?,?,?,?,?)");
        ResultSet rs = getdata.executeQuery();
        int inserted = 0;
        System.err.println("Normratio is " + normratio + " and multiplicative is " + multiplicative + " and usewce is " + usewce + " and useboth is " + usebothchannels);
        while (rs.next()) {
            int probe = rs.getInt(1);
            double valone = Math.max(rs.getDouble(2),10);            
            double valtwo = Math.max(rs.getDouble(3),10);
            double mor, channelratio, ratio;
            if (Double.isNaN(valone)) {
                valone = 10;
            }
            if (Double.isNaN(valtwo)) {
                valtwo = 10;
            }
            if (normratio){
                if (multiplicative) {
                    channelratio = valone / valtwo;
                    ratio = channelratio * ratiomean / (ipvals.get(probe) / wcevals.get(probe));
                    mor = ratio;
                } else {
                    channelratio = valone / valtwo;
                    ratio = channelratio + (ratiomean - (ipvals.get(probe) / wcevals.get(probe)));
                    mor = ratio;
                }
            } else {
                if (multiplicative) {
                    if (usebothchannels) {
                        valone *= ipmean / ipvals.get(probe);
                        valtwo *= wcemean / wcevals.get(probe);
                    } else if (usewce ) {
                        valone *= wcemean / wcevals.get(probe);
                        valtwo *= wcemean / wcevals.get(probe);
                    } else {
                        valone *= ipmean / ipvals.get(probe);
                        valtwo *= ipmean / ipvals.get(probe);
                    }
                } else {
                    if (usebothchannels) {
                        valone += ipmean - ipvals.get(probe);
                        valtwo += wcemean - wcevals.get(probe);
                    } else if (usewce ) {
                        valone += wcemean - wcevals.get(probe);
                        valtwo += wcemean - wcevals.get(probe);
                    } else {
                        valone += ipmean - ipvals.get(probe);
                        valtwo += ipmean - ipvals.get(probe);
                    }
                }
                ratio = valone / valtwo;
                channelratio = mor = ratio;
            }          
            insertdata.setInt(1,probe);
            insertdata.setDouble(2,valone);
            insertdata.setDouble(3,valtwo);
            insertdata.setDouble(4,ratio);
            insertdata.setDouble(5,ratio);
            insertdata.setDouble(6,ratio);
            insertdata.execute();
            inserted++;
        }
        cxn.commit();
        System.err.println("Inserted " + inserted + " probe values");
        cxn.commit();
    }
}

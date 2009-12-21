package edu.mit.csail.cgs.tools.chipchip;

import java.sql.*;
import java.util.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.chipchip.*;
import edu.mit.csail.cgs.tools.chipchip.NormalizeExperiment;

/**
 * CrossTalkNormalization removes some component of channel a from 
 * channel b and vice versa.
 *
 * iptowce is the amount of bleedthrough from cy5 to cy3.
 */
public class CrossTalkNormalization extends NormalizeExperiment {



    private int oldexptid, newexptid;
    private Organism org;
    private Genome genome;
    private double ipToWCE, wceToIP;
    
    public CrossTalkNormalization (Experiment oldexpt,
                                   Experiment newexpt,
                                   double iptowce, double wcetoip) throws SQLException, NotFoundException {
        ipToWCE = iptowce;
        wceToIP = wcetoip;
        oldexptid = oldexpt.getDBID();
        newexptid = newexpt.getDBID();
    }

    public void doNorm() throws SQLException {        
        java.sql.Connection cxn = DatabaseFactory.getConnection("chipchip");
        Statement stmt = cxn.createStatement();
        String sql = String.format("insert into data(experiment, probe, channelone, channeltwo, channelratio) " +
                                   "(select %d as experiemnt, probe, greatest(channelone - (%f * channeltwo), 10), greatest(channeltwo - (%f * channelone),10), channelratio " + 
                                   "from data where experiment = %d)",
                                   newexptid,
                                   wceToIP, ipToWCE,
                                   oldexptid);        
        stmt.execute(sql);
        sql = String.format("update data set ratio = channelone / channeltwo, mor = channelone / channeltwo where experiment = %d",newexptid);
        stmt.execute(sql);
    }

}
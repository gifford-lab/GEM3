/*
 * Created on Nov 8, 2006
 */
package edu.mit.csail.cgs.datasets.chipchip;

import java.util.*;
import java.sql.*;

import edu.mit.csail.cgs.datasets.general.Cells;
import edu.mit.csail.cgs.datasets.general.Condition;
import edu.mit.csail.cgs.datasets.general.Factor;
import edu.mit.csail.cgs.datasets.locators.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.Closeable;
import edu.mit.csail.cgs.utils.database.*;


/**
 * @author tdanford
 */
public class LocatorLoader implements Closeable {
    
    private java.sql.Connection cxn;
    private Genome genome;

    public LocatorLoader(Genome g) throws SQLException {
        genome = g;
        cxn = DatabaseFactory.getConnection("chipchip");
    }

    public void close() {
        DatabaseFactory.freeConnection(cxn);
        cxn = null;
    }

    public boolean isClosed() {
        return cxn == null;
    }

    public Collection<ExptLocator> loadLocators(Cells cells, Condition cond, Factor factor) throws SQLException { 
        LinkedList<ExptLocator> locs = new LinkedList<ExptLocator>();
        Statement s = cxn.createStatement();
        String bayesquery = "select a.name, a.version from bayesanalysis a, bayesanalysisinputs inp, bayesToGenome b2g where " +
                "a.id=b2g.analysis and b2g.genome=? and a.id=inp.analysis and inp.experiment=?"; 
        String mspquery = "select a.name, a.version from rosettaanalysis a, rosettaanalysisinputs inp, " +
                "rosettaToGenome b2g where a.id=b2g.analysis and b2g.genome=? and a.id=inp.analysis and inp.experiment=?";
        PreparedStatement bps = cxn.prepareStatement(bayesquery);
        PreparedStatement mps = cxn.prepareStatement(mspquery);
        
        String query = "select e.id, e.name, e.version from experiment e, exptToGenome e2g where " +
                "e.id=e2g.experiment and e2g.genome=" + genome.getDBID() + " and e.active=1";
        if(cells != null) { 
            int dbid = cells.getDBID();
            query += " and (e.cellsone=" + dbid + " or e.cellstwo=" + dbid + ")";
        }
        if(cond != null) { 
            int dbid = cond.getDBID();
            query += " and (e.conditionone=" + dbid + " or e.conditiontwo=" + dbid + ")";
        }
        if(factor != null) { 
            int dbid = factor.getDBID();
            query += " and (e.factorone=" + dbid + " or e.factortwo=" + dbid + ")";
        }
        ResultSet rs = s.executeQuery(query);
        
        while(rs.next()) {
            int exptID = rs.getInt(1);
            String name = rs.getString(2); 
            String version = rs.getString(3);
            ChipChipLocator loc = new ChipChipLocator(genome, name, version);
            locs.addLast(loc);
            
            bps.setInt(1,genome.getDBID());
            bps.setInt(2,exptID);
            ResultSet ars = bps.executeQuery();
            while(ars.next()) { 
                String aname = ars.getString(1); 
                String aversion = ars.getString(2);
                BayesLocator aloc = new BayesLocator(genome, aname, aversion);
                locs.addLast(aloc);
            }
            ars.close();

            mps.setInt(1,genome.getDBID());
            mps.setInt(2,exptID);
            ars = mps.executeQuery();
            while(ars.next()) { 
                String aname = ars.getString(1); 
                String aversion = ars.getString(2);
                MSPLocator aloc = new MSPLocator(genome, aname, aversion);
                locs.addLast(aloc);
            }
            ars.close();
        }
        
        rs.close();
        s.close();
        mps.close();
        bps.close();
        return locs;
    }
}

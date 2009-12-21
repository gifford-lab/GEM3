/*
 * Created on Feb 2, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.datasets.locators;

import java.io.*;
import java.util.*;
import java.sql.*;

import edu.mit.csail.cgs.datasets.chipchip.*;
import edu.mit.csail.cgs.datasets.general.Cells;
import edu.mit.csail.cgs.datasets.general.Condition;
import edu.mit.csail.cgs.datasets.general.Factor;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.utils.*;

public class ExptLocatorFilter implements edu.mit.csail.cgs.utils.Closeable {
    
    private java.sql.Connection cxn;
    private Genome genome;
    
    public ExptLocatorFilter() {
        try { 
            cxn = DatabaseFactory.getConnection(ExptLocator.dbRole);
            genome = null;
        } catch(Exception ue) { 
            ue.printStackTrace(System.err);
            throw new RuntimeException("Unknown Role: " + ExptLocator.dbRole);
        }
    }
    
    public void setGenome(Genome g) { 
        genome = g;
    }
    
    public Genome getGenome() { return genome; }

    public void close() {
        DatabaseFactory.freeConnection(cxn);
        cxn = null;
    }

    public boolean isClosed() {
        return cxn == null;
    }

    public Collection<ExptLocator> findLocators(ExptLocatorFilterOptions opts) throws SQLException {
        LinkedList<ExptLocator> locs = new LinkedList<ExptLocator>();
        System.err.println("findLocators : " + genome);
        if(genome == null) { return locs; }

        System.err.println("opts is " + opts);
        int expt_type = opts.getExptType();
        Cells c = opts.getCells();
        Condition d = opts.getCondition();
        Factor f = opts.getFactor();
        
        String query = null;
        Statement stmt = null;
        ResultSet rs = null;
        System.err.println("TYPE is " + expt_type);
        switch(expt_type) {             
        case ExptLocatorFilterOptions.AGILENT_TYPE:
            
            query = "select e.name, e.version, e.replicate from experiment e, exptToGenome e2g " +
                    "where e.id=e2g.experiment and e.active = 1";
            query += " and e2g.genome=" + genome.getDBID();
            query += c != null ? " and (e.cellsone=" + c.getDBID() + " or e.cellstwo=" + c.getDBID() + ")" : "";
            query += d != null ? " and (e.conditionone=" + d.getDBID() + " or e.conditiontwo=" + d.getDBID() + ")" : "";
            query += f != null ? " and (e.factorone=" + f.getDBID() + " or e.factortwo=" + f.getDBID() + ")" : "";
            query += "order by e.name, e.version";
            stmt = cxn.createStatement();
            rs = stmt.executeQuery(query);
            while(rs.next()) { 
                String name = rs.getString(1), version = rs.getString(2), replicate = rs.getString(3);
                ChipChipLocator loc = new ChipChipLocator(genome, name, version, replicate);
                locs.addLast(loc);
            }
            rs.close();
            stmt.close();
            
            break;
            
        case ExptLocatorFilterOptions.MSP_TYPE:
            
            query = "select a.name, a.version from rosettaanalysis a, rosettaToGenome a2g " +
                    " where a.id=a2g.analysis and a.active = 1";
            query += " and a2g.genome=" + genome.getDBID();
            query += c != null ? " and (a.cellsone=" + c.getDBID() + " or a.cellstwo=" + c.getDBID() + ")" : "";
            query += d != null ? " and (a.conditionone=" + d.getDBID() + " or a.conditiontwo=" + d.getDBID() + ")" : "";
            query += f != null ? " and (a.factorone=" + f.getDBID() + " or a.factortwo=" + f.getDBID() + ")" : "";            
            query += "order by a.name, a.version";

            stmt = cxn.createStatement();
            rs = stmt.executeQuery(query);
            while(rs.next()) { 
                String name = rs.getString(1), version = rs.getString(2);
                MSPLocator loc = new MSPLocator(genome, name, version);
                locs.addLast(loc);
            }
            rs.close();
            stmt.close();
            

            break;
        case ExptLocatorFilterOptions.BAYES_TYPE:

            query = "select a.name, a.version from bayesanalysis a, bayesToGenome a2g, bayesanalysisinputs a2e, " +
            "experiment e where a.id=a2g.analysis and a.id=a2e.analysis and a2e.experiment=e.id and a.active=1";
            query += " and a2g.genome=" + genome.getDBID();
            query += c != null ? " and (e.cellsone=" + c.getDBID() + " or e.cellstwo=" + c.getDBID() + ")" : "";
            query += d != null ? " and (e.conditionone=" + d.getDBID() + " or e.conditiontwo=" + d.getDBID() + ")" : "";
            query += f != null ? " and (e.factorone=" + f.getDBID() + " or e.factortwo=" + f.getDBID() + ")" : "";            
            query += "order by a.name, a.version";

            stmt = cxn.createStatement();
            rs = stmt.executeQuery(query);
            while(rs.next()) { 
                String name = rs.getString(1), version = rs.getString(2);
                BayesLocator loc = new BayesLocator(genome, name, version);
                locs.addLast(loc);
            }
            rs.close();
            stmt.close();
            break;
        }
        
        return locs;
    }
}

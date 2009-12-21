/**
 * 
 */
package edu.mit.csail.cgs.datasets.locators;

import java.util.*;
import java.sql.*;
import java.io.*;
import javax.swing.*;
import java.awt.*;
import javax.swing.tree.*;

import edu.mit.csail.cgs.datasets.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.utils.datastructures.TaxonomyImpl;

/**
 * @author Timothy Danford
 */
public class ExptLocatorTree extends TaxonomyImpl<ExptLocator> {
	
    public ExptLocatorTree(Genome g) {
        super();
        try { 
            java.sql.Connection c = DatabaseFactory.getConnection(ExptLocator.dbRole);
            int species=g.getSpeciesDBID();
            int genome = g.getDBID();

            Statement s = c.createStatement();
            ResultSet rs = null;

            rs = s.executeQuery("select e.name, e.version from experiment e, exptToGenome eg where e.active=1 and " + 
                                "e.id=eg.experiment and eg.genome=" + genome);
            while(rs.next()) { 
                String name = rs.getString(1);
                String version = rs.getString(2);
                ChipChipLocator loc = new ChipChipLocator(g, name, version);
                this.addElement(loc.getTreeAddr(), loc);
            }
            rs.close();
                
            rs = s.executeQuery("select ra.name, ra.version from rosettaanalysis ra, rosettaToGenome rg where " + 
                                "ra.id = rg.analysis and ra.active=1 and rg.genome=" + genome);
            while(rs.next()) { 
                String name = rs.getString(1);
                String version = rs.getString(2);
                MSPLocator msp = new MSPLocator(g, name, version);
                this.addElement(msp.getTreeAddr(), msp);
            }
            rs.close();
			
            rs = s.executeQuery("select ra.name, ra.version from bayesanalysis ra, bayesToGenome rg where " + 
                                "ra.id = rg.analysis and ra.active=1 and rg.genome=" + genome);
            while(rs.next()) { 
                String name2 = rs.getString(1);
                String version2 = rs.getString(2);
                ExptLocator loc2 = new BayesLocator(g, name2, version2);
                addElement(loc2.getTreeAddr(), loc2);
            }
            rs.close();
            s.close();

            DatabaseFactory.freeConnection(c);
        } catch(SQLException se) { 
            se.printStackTrace(System.err);
            throw new RuntimeException(se);
        } catch (UnknownRoleException e) {
            e.printStackTrace();
        }
    }
}

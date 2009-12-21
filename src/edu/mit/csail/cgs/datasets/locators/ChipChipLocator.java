/**
 * 
 */
package edu.mit.csail.cgs.datasets.locators;

import java.io.*;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.sql.PreparedStatement;
import java.util.Collection;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;
import java.util.HashSet;

import edu.mit.csail.cgs.datasets.*;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipData;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipDataset;
import edu.mit.csail.cgs.datasets.chipchip.ExptNameVersion;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.preferences.*;

import edu.mit.csail.cgs.utils.database.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

/**
 * @author Timothy Danford
 *
 */
public class ChipChipLocator 
    extends ExptNameVersion 
    implements ExptLocator {
    
    public static Collection<ChipChipLocator> getLocators(java.sql.Connection cxn, Genome genome) 
        throws SQLException {
        
        LinkedList<ChipChipLocator> locs = new LinkedList<ChipChipLocator>();        
        Statement s = cxn.createStatement();

        int species = genome.getSpeciesDBID();
        ResultSet rs = null;
        
        rs = s.executeQuery("select name, version, id from experiment where active=1 and species=" + species);
        
        while(rs.next()) { 
            String name = rs.getString(1);
            String version = rs.getString(2);
            ChipChipLocator loc = new ChipChipLocator(genome, name, version);
            locs.addLast(loc);
        }
        
        rs.close();
        s.close();
        
        return locs;
    }
	
    private ChipChipDataset ds;

    public ChipChipLocator(ChipChipDataset baseDS, String name, String version, String replicate) {
        super(name, version,replicate);
        ds = baseDS;
    }

    public ChipChipLocator(ChipChipDataset baseDS, String name, String version) {
        super(name, version);
        ds = baseDS;
        replicate = null;
    }

    public ChipChipLocator(Genome g, String name, String version, String replicate) { 
        super(name, version,replicate);
        ds = g.getChipChipDataset();
    }
	
    public ChipChipLocator(Genome g, String name, String version) { 
        super(name, version);
        ds = g.getChipChipDataset();
        replicate = null;
    }
	
    public ChipChipLocator(Genome g, DataInputStream dis) 
        throws IOException { 
        super(dis);
        ds = g.getChipChipDataset();
    }
    
    public LinkedList<String> getTreeAddr() { 
        LinkedList<String> lst = new LinkedList<String>();
        lst.addLast("agilent");
        lst.addLast(version);
        return lst;
    }

    public ExptNameVersion getNameVersion() { return this; }

    public Set<String> getReplicates() {
        try {
            java.sql.Connection chipcxn;
            chipcxn = DatabaseFactory.getConnection("chipchip");
            PreparedStatement stmt = chipcxn.prepareStatement("select unique(replicate) from experiment where name = ? and version = ?");
            ResultSet rs = stmt.executeQuery();
            HashSet<String> results = new HashSet<String>();
            while (rs.next()) {
                results.add(rs.getString(1));
            }
            rs.close();
            stmt.close();
            DatabaseFactory.freeConnection(chipcxn);
            return results;
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't connect with role chipchip",ex);
        } catch (SQLException ex) {
            throw new DatabaseException(ex.toString(),ex);
        } 
    }       

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.preferences.Preferences#getName()
     */
//     public String getName() {
//         if (replicate == null) {
//             return "ChipChip(" + super.toString() + ")";
//         } else {
//             return "ChipChip(" + super.toString() + "," + replicate + ")";
//         }
//     }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.preferences.Preferences#createPanel()
     */
    public PreferencesPanel createPanel() {
        return new ChipChipPanel(this);
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.preferences.Preferences#saveFromPanel(edu.mit.csail.cgs.utils.preferences.PreferencesPanel)
     */
    public void saveFromPanel(PreferencesPanel pp) {
        // do nothing, yet.
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.Factory#createObject()
     */
    public ChipChipData createObject() {
        try { 
            return ds.getData(this);
        } catch(NotFoundException nfe) { 
        	nfe.printStackTrace(System.err);
            throw new IllegalArgumentException(nfe.getMessage());
        }
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.Saveable#save(java.io.DataOutputStream)
     */
    public void save(DataOutputStream dos) throws IOException {
        dos.writeInt(0);
        super.save(dos);
    }

    public static class ChipChipPanel extends PreferencesPanel {
		
        private JLabel name, version;
		
        public ChipChipPanel(ChipChipLocator bl) { 
            super();
            setLayout(new BorderLayout());
            JPanel innerPanel = new JPanel(); 
            innerPanel.setLayout(new GridLayout(2, 2));
            innerPanel.add(new JLabel("Expt Name:")); 
            innerPanel.add((name = new JLabel(bl.name)));
            innerPanel.add(new JLabel("Expt Version:"));
            innerPanel.add((version = new JLabel(bl.version)));
        }
		
        public void saveValues() { 
            super.saveValues();
        }
    }
}

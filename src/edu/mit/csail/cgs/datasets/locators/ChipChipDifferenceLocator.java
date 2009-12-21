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
import edu.mit.csail.cgs.datasets.chipchip.ChipChipDifferenceData;
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
public class ChipChipDifferenceLocator 
    extends ChipChipLocator {
    
    private ChipChipLocator loc1, loc2;

    public ChipChipDifferenceLocator(Genome g, ChipChipLocator loc1, ChipChipLocator loc2) {  
        super(g.getChipChipDataset(), 
                String.format("Difference(%s-%s)", loc1.name, loc2.name), 
                String.format("%s:%s", loc1.version, loc2.version));
        this.loc1 = loc1;
        this.loc2 = loc2;
    }
	
    public LinkedList<String> getTreeAddr() { 
        LinkedList<String> lst = new LinkedList<String>();
        lst.addLast("difference");
        return lst;
    }

    public ExptNameVersion getNameVersion() { return this; }

    public boolean equals(Object o) { 
        if(!(o instanceof ChipChipDifferenceLocator)) { return false; }
        ChipChipDifferenceLocator loc = (ChipChipDifferenceLocator)o;
        return loc1.equals(loc.loc1) && loc2.equals(loc.loc2);
    }
    
    public int hashCode() { 
        int code = replicate == null ? 17 : replicate.hashCode();
        code += loc1.hashCode(); 
        code += loc2.hashCode(); 
        code *= 37;
        return code;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.preferences.Preferences#getName()
     */
    public String getName() {
        return super.toString();
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.Factory#createObject()
     */
    public ChipChipData createObject() {
        ChipChipData d1 = loc1.createObject();
        ChipChipData d2 = loc2.createObject();
        return new ChipChipDifferenceData(d1, d2);
    }

    public PreferencesPanel createPanel() {
        return null;
    }

    public void saveFromPanel(PreferencesPanel pp) {
    }
}

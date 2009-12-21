/**
 * 
 */
package edu.mit.csail.cgs.datasets.locators;

import java.io.*;
import java.util.*;

import edu.mit.csail.cgs.datasets.*;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipBayes;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipDataset;
import edu.mit.csail.cgs.datasets.chipchip.AnalysisNameVersion;
import edu.mit.csail.cgs.datasets.chipchip.NameVersion;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.preferences.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

/**
 * @author Timothy Danford
 *
 */
public class BayesLocator 
	extends AnalysisNameVersion 
	implements ExptLocator {
	
	private ChipChipDataset ds;

	public BayesLocator(ChipChipDataset baseDS, String name, String version) {
		super(name, version);
		ds = baseDS;
	}
	
	public BayesLocator(Genome g, String name, String version) {
		super(name, version);
		ds = g.getChipChipDataset();
	}

	public BayesLocator(Genome g, DataInputStream dis) 
		throws IOException { 
		super(dis);
		ds = g.getChipChipDataset();
	}
    
    public boolean equals(Object o) { 
        if(!(o instanceof BayesLocator)) { return false; }
        BayesLocator loc = (BayesLocator)o;
        if(!super.equals(loc)) { return false; }
        return true;
    }
    
    public int hashCode() { 
        int code = 17;
        code += super.hashCode(); code *= 37;
        return code;
    }
    
	public LinkedList<String> getTreeAddr() { 
		LinkedList<String> lst = new LinkedList<String>();
        lst.addLast("bayes");
		lst.addLast(version);
		return lst;
	}

    public NameVersion getNameVersion() { return this; }

	/* (non-Javadoc)
	 * @see edu.mit.csail.cgs.utils.preferences.Preferences#getName()
	 */
	public String getName() {
		return "Bayes(" + super.toString() + ")";
	}

	/* (non-Javadoc)
	 * @see edu.mit.csail.cgs.utils.preferences.Preferences#createPanel()
	 */
	public PreferencesPanel createPanel() {
		return new BayesPanel(this);
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
	public ChipChipBayes createObject() {
		try { 
			return ds.getBayes(name, version);
		} catch(NotFoundException nfe) { 
            nfe.printStackTrace(System.err);
			throw new IllegalArgumentException(super.toString());
		}
	}

	/* (non-Javadoc)
	 * @see edu.mit.csail.cgs.utils.Saveable#save(java.io.DataOutputStream)
	 */
	public void save(DataOutputStream dos) throws IOException {
		dos.writeInt(1);
		super.save(dos);
	}

	public static class BayesPanel extends PreferencesPanel {
		
		private JLabel name, version, resultsversion;
		
		public BayesPanel(BayesLocator bl) { 
			super();
			setLayout(new BorderLayout());
			JPanel innerPanel = new JPanel(); 
			innerPanel.setLayout(new GridLayout(2, 2));
			innerPanel.add(new JLabel("Bayes Name:")); 
			innerPanel.add((name = new JLabel(bl.name)));
			innerPanel.add(new JLabel("Bayes Version:"));
			innerPanel.add((version = new JLabel(bl.version)));
		}
		
		public void saveValues() { 
			super.saveValues();
		}
	}
}

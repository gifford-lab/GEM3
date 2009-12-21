/**
 * 
 */
package edu.mit.csail.cgs.viz.utils;

import java.util.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;

import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.utils.*;

/**
 * @author Timothy Danford
 *
 */
public class GenomeSelectPanel extends JPanel implements EventSource<ActionEvent> {
	
	public static void main(String[] args) { 
		JFrame f = new JFrame("Genome Select Test");
		Container c = (Container)f.getContentPane();
		c.setLayout(new BorderLayout());
		GenomeSelectPanel p = new GenomeSelectPanel("Homo sapiens", "hg17");
		p.addActionListener(new ActionListener() { 
			public void actionPerformed(ActionEvent e) {
				GenomeSelectPanel gsp = (GenomeSelectPanel)e.getSource();
				String s = gsp.getSpecies();
				String g = gsp.getGenome();
				System.out.println(s + "\t" + g);
			}
		});
		
		c.add(p, BorderLayout.NORTH);
		f.setLocation(100, 100);
		f.setSize(200, 100);
		f.setVisible(true);
		f.pack();
		f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
	
	private String defaultSpecies, defaultGenome;
	private DefaultComboBoxModel speciesModel, genomeModel;
	private JComboBox speciesBox, genomeBox;
    private EventSource.Default<ActionEvent> src;
    
    public GenomeSelectPanel(Genome g) { 
    	super();
    	defaultGenome = g.getVersion();
    	defaultSpecies = g.getSpecies();
    	init();
    }
    
    public GenomeSelectPanel() {
    	super();
    	defaultSpecies = defaultGenome = null;
    	init();
    }

	public GenomeSelectPanel(String species, String genome) { 
		super();
		defaultSpecies = species;
		defaultGenome = genome;
		init();
	}
		
	private void init() {
        src = new EventSource.Default<ActionEvent>(this);
		setLayout(new BorderLayout());
		
		speciesModel = new DefaultComboBoxModel();
		genomeModel = new DefaultComboBoxModel();
		
		speciesBox = new JComboBox(speciesModel);
		genomeBox = new JComboBox(genomeModel);
		speciesBox.setLightWeightPopupEnabled(false);
		genomeBox.setLightWeightPopupEnabled(false);

		JLabel speciesLabel = new JLabel("Species:");
		JLabel genomeLabel = new JLabel("Genome Version:");
		
		JPanel pairPanel = new JPanel();
		pairPanel.setLayout(new GridLayout(2, 2));
		add(pairPanel, BorderLayout.NORTH);
		pairPanel.add(speciesLabel);
		pairPanel.add(speciesBox);
		pairPanel.add(genomeLabel);
		pairPanel.add(genomeBox);
		
		populateSpecies();
		
		speciesBox.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) { 
				populateGenome();
			}
		});
		
		genomeBox.addActionListener(new ActionListener() { 
			public void actionPerformed(ActionEvent e) { 
				JComboBox b = (JComboBox)e.getSource();
				ComboBoxModel m = (ComboBoxModel)b.getModel();
				String s = (String)m.getSelectedItem();
				if(s != null) { fireActionEvent(); }
			}
		});
	}
    
    public void addEventListener(Listener<ActionEvent> e) { 
        src.addEventListener(e);
    }
    
    public void removeEventListener(Listener<ActionEvent> e) { 
        src.removeEventListener(e);
    }
	
	public void addActionListener(ActionListener l) { 
        src.addEventListener(new ActionListenerWrapper(l));
	}
	
	public void removeActionListener(ActionListener l) { 
        src.removeEventListener(new ActionListenerWrapper(l));
	}
	
	private void fireActionEvent() { 
		ActionEvent e = new ActionEvent(this, ActionEvent.ACTION_PERFORMED, "Selection Made");
        src.fireEvent(e);
	}
	
	private void populateSpecies() { 
		TreeSet<String> names = new TreeSet<String>(Organism.getOrganismNames());
		String sel = null;
		speciesModel.removeAllElements();
		for(String n : names) { 
			speciesModel.addElement(n);
			if((defaultSpecies==null && sel == null) || n.equals(defaultSpecies)) { 
				sel = n;
			} 
		}
		if(sel != null) { speciesModel.setSelectedItem(sel); } 
		populateGenome();
	}
	
	private void populateGenome() {
		try { 
			Organism org = Organism.getOrganism(getSpecies());
			TreeSet<String> names = new TreeSet<String>(org.getGenomeNames());
			genomeModel.removeAllElements();
			String sel = null;
			for(String n : names) { 
				genomeModel.addElement(n);
				if((defaultGenome == null && sel == null) || n.equals(defaultGenome)) { 
					sel = n;
				}
			}
			if(sel != null) { genomeModel.setSelectedItem(sel); } 
		} catch(NotFoundException nfe) { 
			nfe.printStackTrace(System.err);
			throw new IllegalStateException(nfe);
		}
	}
	
	public String getSpecies() { return (String)speciesModel.getSelectedItem(); }
	public String getGenome() { return (String)genomeModel.getSelectedItem(); }
	
	public Genome findGenome() { 
		try { 
			Organism org = Organism.getOrganism(getSpecies());
			Genome g = org.getGenome(getGenome());
			return g;
		} catch(NotFoundException e) { 
			return null;
		}
	}

    public boolean hasListeners() {
        return src.hasListeners();
    }
}

class ActionListenerWrapper implements Listener<ActionEvent> { 
    private ActionListener listener;
    public ActionListenerWrapper(ActionListener l) { 
        listener = l;
    }
    
    public void eventRegistered(ActionEvent e) { 
        listener.actionPerformed(e);
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof ActionListenerWrapper)) { return false; }
        ActionListenerWrapper w = (ActionListenerWrapper)o;
        return listener == w.listener;
    }
    
    public int hashCode() { return listener.hashCode(); }
}

package edu.mit.csail.cgs.echo.gui;

import javax.swing.*;

import edu.mit.csail.cgs.datasets.binding.BindingScan;
import edu.mit.csail.cgs.datasets.locators.ExptLocator;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.viz.components.BindingScanSelectPanel;
import edu.mit.csail.cgs.viz.utils.GenomeSelectPanel;

import java.awt.*;
import java.awt.event.*;
import java.sql.SQLException;
import java.util.*;

public class BindingSelectFrame extends JFrame 
	implements EventSource<CreationEvent<BindingScan>> {

	private GenomeSelectPanel genomeSelect;
	private BindingScanSelectPanel selectPanel;
	private JButton ok, cancel;
	private EventSource.Default<CreationEvent<BindingScan>> src;
	
	public BindingSelectFrame(Genome g) throws SQLException { 
		super("Experiment Selection");
		
		src = new EventSource.Default<CreationEvent<BindingScan>>();
        
        String gs = g != null ? g.getSpecies() : null; 
        String gv = g != null ? g.getVersion() : null;
        
		genomeSelect = new GenomeSelectPanel(gs, gv);
		selectPanel = new BindingScanSelectPanel();

		selectPanel.setGenome(g);
		
		ok = new JButton("Ok");
		cancel = new JButton("Cancel");
		
		genomeSelect.addEventListener(new Listener<ActionEvent>() {
			public void eventRegistered(ActionEvent e) { 
				selectPanel.setGenome(genomeSelect.findGenome());
			}
		});
		
		ok.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) { okay(); }
		});
		cancel.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) { cancel(); }
		});
		
		Container c = (Container)getContentPane();
		c.setLayout(new BorderLayout());
		
		c.add(genomeSelect, BorderLayout.NORTH);
		c.add(selectPanel, BorderLayout.CENTER);
		
		JPanel buttons = new JPanel();
		buttons.setLayout(new FlowLayout());
		buttons.add(ok);
		buttons.add(cancel);
		c.add(buttons, BorderLayout.SOUTH);
		
		setVisible(true);
		pack();
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
	}
	
	public void okay() {
		LinkedList<BindingScan> locs = 
			new LinkedList<BindingScan>(selectPanel.getSelected());
		BindingScan scan = locs.isEmpty() ? null : locs.getFirst();
		CreationEvent<BindingScan> evt = 
			new CreationEvent<BindingScan>(this, scan);
		System.out.println("Selected: " + scan);
		src.fireEvent(evt);
		dispose();
	}
	
	public void cancel() { 
		CreationEvent<BindingScan> evt = 
			new CreationEvent<BindingScan>(this, null);
		src.fireEvent(evt);
		dispose();
	}

	public void addEventListener(Listener<CreationEvent<BindingScan>> el) {
		src.addEventListener(el);
	}

	public boolean hasListeners() {
		return src.hasListeners();
	}

	public void removeEventListener(Listener<CreationEvent<BindingScan>> el) {
		src.removeEventListener(el);
	}
}

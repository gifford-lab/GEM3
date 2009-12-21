package edu.mit.csail.cgs.echo.gui;

import javax.swing.*;

import edu.mit.csail.cgs.datasets.locators.ExptLocator;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.viz.utils.GenomeSelectPanel;
import edu.mit.csail.cgs.viz.components.ExptSelectPanel;

import java.awt.*;
import java.awt.event.*;
import java.sql.SQLException;
import java.util.*;

public class ExptSelectFrame extends JFrame implements EventSource<CreationEvent<ExptLocator>> {

	private GenomeSelectPanel genomeSelect;
	private ExptSelectPanel selectPanel;
	private JButton ok, cancel;
	private EventSource.Default<CreationEvent<ExptLocator>> src;
	
	public ExptSelectFrame(Genome g) { 
		super("Experiment Selection");
		
		src = new EventSource.Default<CreationEvent<ExptLocator>>();
        
        String gs = g != null ? g.getSpecies() : null; 
        String gv = g != null ? g.getVersion() : null;
        
		genomeSelect = new GenomeSelectPanel(gs, gv);
		selectPanel = new ExptSelectPanel(genomeSelect.findGenome());
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
		LinkedList<ExptLocator> locs = new LinkedList<ExptLocator>(selectPanel.getSelected());
		ExptLocator exptLoc = locs.isEmpty() ? null : locs.getFirst();
		CreationEvent<ExptLocator> evt = new CreationEvent<ExptLocator>(this, exptLoc);
		System.out.println("Selected: " + exptLoc);
		src.fireEvent(evt);
		dispose();
	}
	
	public void cancel() { 
		CreationEvent<ExptLocator> evt = new CreationEvent<ExptLocator>(this, null);
		src.fireEvent(evt);
		dispose();
	}

	public void addEventListener(Listener<CreationEvent<ExptLocator>> el) {
		src.addEventListener(el);
	}

	public boolean hasListeners() {
		return src.hasListeners();
	}

	public void removeEventListener(Listener<CreationEvent<ExptLocator>> el) {
		src.removeEventListener(el);
	}
}

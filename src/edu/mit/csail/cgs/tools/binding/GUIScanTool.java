package edu.mit.csail.cgs.tools.binding;

import java.util.*;
import java.io.*;
import java.sql.*;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.border.*;

import edu.mit.csail.cgs.utils.Listener;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.viz.utils.GenomeSelectPanel;
import edu.mit.csail.cgs.viz.components.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.verbs.binding.CallerMapper;
import edu.mit.csail.cgs.datasets.locators.*;
import edu.mit.csail.cgs.datasets.binding.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;

public class GUIScanTool extends JFrame implements Listener<ActionEvent> {

	public static void main(String[] args) { 
		try {
			GUIScanTool tool = new GUIScanTool();
			tool.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		} catch (SQLException e) {
			e.printStackTrace(System.err);
		} catch (UnknownRoleException e) {
			e.printStackTrace(System.err);
		}
	}
	
	private BindingScanLoader loader;
	private JButton scanButton;	
	private GenomeSelectPanel genomeSelector;
	private ExptSelectPanel selectPanel;
    private JTextField paramsField;

	public GUIScanTool() throws SQLException, UnknownRoleException { 
		super("Scan Tool");
		
		loader = new BindingScanLoader();
		
		Container c = (Container)getContentPane();
		c.setLayout(new BorderLayout());
		c.add(genomeSelector = new GenomeSelectPanel("Homo sapiens", "hg17"), 
				BorderLayout.NORTH);
		c.add(selectPanel = new ExptSelectPanel(getGenome()), 
				BorderLayout.CENTER);
		
		JPanel scanPanel = new JPanel(); 
		scanPanel.setLayout(new GridLayout(1, 1));
		scanPanel.add(scanButton = new JButton("Scan"));
        
        paramsField = new JTextField();
        JPanel paramsPanel = new JPanel(); paramsPanel.setLayout(new BorderLayout());
        paramsPanel.add(paramsField, BorderLayout.NORTH);
        paramsPanel.setBorder(new TitledBorder("Parameters"));
        
        JPanel bottomPanel = new JPanel();
        bottomPanel.setLayout(new BorderLayout());
        bottomPanel.add(scanPanel, BorderLayout.SOUTH);
        bottomPanel.add(paramsPanel, BorderLayout.NORTH);
		
        c.add(bottomPanel, BorderLayout.SOUTH);
        
        scanButton.addActionListener(new ActionListener() { 
			public void actionPerformed(ActionEvent e) { 
				try {
					scan();
				} catch (SQLException e1) {
					e1.printStackTrace(System.err);
				} catch (UnknownRoleException e1) {
					e1.printStackTrace(System.err);
				}
			}
		});
        
        selectPanel.setGenome(getGenome());
        genomeSelector.addEventListener(this);
		
		setVisible(true);
		pack();
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
	}
	
	public void scan() throws SQLException, UnknownRoleException { 
		while(selectPanel.getNumSelected() > 0) {
            
            BindingParameters params = new BindingParameters(paramsField.getText());
			ExptLocator loc = selectPanel.removeFirstSelected();
			CallerMapper cm = null;

			if(loc instanceof ChipChipLocator) {
			    if(params.containsKey("type") && params.get("type").equals("hmm")) { 
			        cm = new ScanTool.HMMCallerMapper(params);
			    } else { 
			        cm = new ScanTool.ChipChipCallerMapper(params);
			    }
			} else if (loc instanceof MSPLocator) {
			    if(params.containsKey("type") && params.get("type").equals("simple")) { 
			        cm = new ScanTool.SimpleMSPCallerMapper(params);
			    } else { 
			        cm = new ScanTool.MSPCallerMapper(params);
			    }
			} else if (loc instanceof BayesLocator) { 
			    cm = new ScanTool.JBDCallerMapper(params);
			}

			if(cm != null && params.containsKey("type") && params.get("type").equals("domain")) { 
			    cm = new ScanTool.DomainCallerMapper(params, cm);
			}

			if(cm != null) { 
			    ScanTool scanner = new ScanTool(getGenome(), loader, cm);
			    scanner.runScan(loc);
			}
		}
	}
	
	public Genome getGenome() {
		try { 
			return Organism.findGenome(genomeSelector.getGenome());
		} catch(NotFoundException nfe) { 
			throw new RuntimeException("Unknown genome", nfe);
		}
	}

	public void eventRegistered(ActionEvent e) {
	    selectPanel.setGenome(getGenome());
	    System.out.println("New Genome: " + getGenome().getVersion());
	}
}

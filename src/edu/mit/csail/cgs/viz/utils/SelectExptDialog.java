/**
 * 
 */
package edu.mit.csail.cgs.viz.utils;

import java.util.*;
import java.io.*;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.border.*;
import javax.swing.tree.*;

import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.datastructures.TaxonomyTreeModel;
import edu.mit.csail.cgs.datasets.locators.ExptLocator;
import edu.mit.csail.cgs.datasets.locators.ExptLocatorTree;
import edu.mit.csail.cgs.datasets.species.*;

/**
 * @author Timothy Danford
 *
 */
public class SelectExptDialog extends JDialog 
	implements EventSource<SelectionEvent> {
	
	private JTree tree;
	private ExptLocatorTree exptLocTree;
	private String speciesName, genomeName;
	private EventSource.Default<SelectionEvent> src;
	private Genome genome;
	
	private LinkedList<ExptLocator> selectedExpts;
	private JButton okButton, cancelButton;
	
	/**
	 * @throws HeadlessException
	 * @throws NotFoundException 
	 */
	public SelectExptDialog(String s, String g) throws HeadlessException, NotFoundException {
		super();
		speciesName = s; genomeName = g;
		Organism org = Organism.getOrganism(speciesName);
		genome = org.getGenome(genomeName);
		src = new EventSource.Default<SelectionEvent>(this);
		selectedExpts = new LinkedList<ExptLocator>();
		layoutDialog();
	}

	/**
	 * @param arg0
	 * @throws HeadlessException
	 * @throws NotFoundException 
	 */
	public SelectExptDialog(Frame arg0, String s, String g) throws HeadlessException, NotFoundException {
		super(arg0);
		speciesName = s; genomeName = g;
		Organism org = Organism.getOrganism(speciesName);
		genome = org.getGenome(genomeName);
		src = new EventSource.Default<SelectionEvent>(this);
		selectedExpts = new LinkedList<ExptLocator>();
		layoutDialog();
	}

	/**
	 * @param arg0
	 * @throws HeadlessException
	 * @throws NotFoundException 
	 */
	public SelectExptDialog(Dialog arg0, String s, String g) throws HeadlessException, NotFoundException {
		super(arg0);
		speciesName = s; genomeName = g;
		Organism org = Organism.getOrganism(speciesName);
		genome = org.getGenome(genomeName);
		src = new EventSource.Default<SelectionEvent>(this);
		selectedExpts = new LinkedList<ExptLocator>();
		layoutDialog();
	}

	private void layoutDialog() { 
		Container c = (Container)getContentPane();
		c.setLayout(new BorderLayout());
		
		exptLocTree = new ExptLocatorTree(genome);
		tree = new JTree(new TaxonomyTreeModel<ExptLocator>(exptLocTree));
		add(new JScrollPane(tree), BorderLayout.CENTER);
		
		okButton = new JButton("Ok");
		cancelButton = new JButton("Cancel");
		JPanel buttonPanel = new JPanel();
		buttonPanel.setLayout(new GridLayout(1, 2));
		buttonPanel.add(okButton);
		buttonPanel.add(cancelButton);
		
		add(buttonPanel, BorderLayout.SOUTH);
		
		okButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				updateSelectedExpts();
				SelectionEvent se = new SelectionEvent(this, SelectionEvent.OK, selectedExpts);
				dispose();
				src.fireEvent(se);
			} 
		});

		cancelButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SelectionEvent se = new SelectionEvent(this, SelectionEvent.CANCEL, null);
				dispose();
				src.fireEvent(se);
			} 
		});
		
		tree.addMouseListener(new MouseAdapter() { 
			public void mouseClicked(MouseEvent e) { 
				if(e.getClickCount() == 2 && e.getButton() == MouseEvent.BUTTON1) {
					TreePath path = tree.getPathForLocation(e.getX(), e.getY());
					Object val = path.getLastPathComponent();
					if(val instanceof ExptLocator) { 
						selectedExpts.clear();
						ExptLocator loc = (ExptLocator)val;
						selectedExpts.addLast(loc);
						dispose();
						SelectionEvent se = new SelectionEvent(this, SelectionEvent.OK, selectedExpts);
						src.fireEvent(se);
					}
				}
			}
		});
	}
	
	public Collection<ExptLocator> getSelectedExpts() { return selectedExpts; }
	
	private void updateSelectedExpts() { 
		TreePath[] paths = tree.getSelectionPaths();
		selectedExpts.clear();
		for(int i = 0; i < paths.length; i++) {
			Object val = paths[i].getLastPathComponent();
			if(val instanceof ExptLocator) { 
				ExptLocator loc = (ExptLocator)val;
				selectedExpts.addLast(loc);
			}			
		}
	}

	public void addEventListener(Listener<SelectionEvent> el) {
		src.addEventListener(el);
	}

	public void removeEventListener(Listener<SelectionEvent> el) {
		src.removeEventListener(el);
	}

    public boolean hasListeners() {
        return src.hasListeners();
    }
}

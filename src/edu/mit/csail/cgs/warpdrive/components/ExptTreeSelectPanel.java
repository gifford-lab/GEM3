/*
 * Created on Aug 24, 2005
 */
package edu.mit.csail.cgs.warpdrive.components;

import java.util.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.tree.*;
import java.awt.*;
import java.awt.event.*;

import edu.mit.csail.cgs.datasets.locators.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.viz.utils.SimpleListPanel;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.datastructures.TaxonomyTreeModel;
import edu.mit.csail.cgs.utils.preferences.*;

/**
 * @author tdanford
 */
public class ExptTreeSelectPanel extends JSplitPane {
    
    private SimpleListPanel<ExptLocator> selected;
    private ExptLocatorTree totalTree;
    private JTree eTree;
    private JButton addButton, deleteButton;
    
    public ExptTreeSelectPanel(Genome genome) {
        super(JSplitPane.VERTICAL_SPLIT);
        setDividerLocation(200);
        selected = new SimpleListPanel<ExptLocator>();
        totalTree = genome != null ? new ExptLocatorTree(genome) : null;
        eTree = genome != null ? new JTree(new TaxonomyTreeModel<ExptLocator>(totalTree)) : new JTree(new NullTreeModel());
        
        addButton = new JButton("Add");
        deleteButton = new JButton("Remove");

        JPanel selPanel = new JPanel();
        selPanel.setLayout(new BorderLayout());
        selPanel.add(selected, BorderLayout.CENTER);
        
        selPanel.add(new JLabel("Selected Expts:"), BorderLayout.NORTH);

        JPanel totalPanel = new JPanel();
        totalPanel.setLayout(new BorderLayout());
        totalPanel.add(new JScrollPane(eTree), BorderLayout.CENTER);
        totalPanel.add(new JLabel("Total Expts:"), BorderLayout.NORTH);
        JPanel buttonPanel = new JPanel();
        buttonPanel.setLayout(new GridBagLayout());
        buttonPanel.add(addButton);
        buttonPanel.add(deleteButton);
        totalPanel.add(buttonPanel, BorderLayout.SOUTH);

        add(selPanel);
        add(totalPanel);
        
        Action addAction = new AbstractAction("Add Experiment") { 
        	public void actionPerformed(ActionEvent e) { 
                TreePath[] paths = eTree.getSelectionModel().getSelectionPaths();
                for(TreePath path : paths) {
                    ExptLocator v = (ExptLocator)path.getLastPathComponent();
                    selected.addValue(v);
                }                        		
        	}
        };
        
        addButton.addActionListener(addAction);
        
        eTree.addMouseListener(new MouseAdapter() {
        	public void mouseClicked(MouseEvent e) { 
        		if(e.getButton() == MouseEvent.BUTTON1 && 
        				e.getClickCount() == 2) {
        			int x = e.getX(), y = e.getY();
        			TreePath selPath = eTree.getPathForLocation(x, y);
        			if(selPath != null && selPath.getLastPathComponent() instanceof ExptLocator) { 
                        ExptLocator v = (ExptLocator)selPath.getLastPathComponent();
                        selected.addValue(v);
        			}
        		}
        	}
        });

        deleteButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                int[] inds = selected.getSelectedIndices();
                Arrays.sort(inds);
                Collection<ExptLocator> selValues = selected.getSelectedValues();
                selected.removeSelectedValues();
            }
        });
    }
    
    public void setGenome(Genome g) { 
        selected.clear(); 

        if(g != null) { 
            totalTree = new ExptLocatorTree(g);
            eTree.setModel(new TaxonomyTreeModel<ExptLocator>(totalTree));
        } else { 
            totalTree = null;
            eTree.setModel(new NullTreeModel());
        }
    }
    
    public Collection<ExptLocator> getSelected() { return selected.getAllValues(); }
    public void addSelected(ExptLocator loc) { selected.addValue(loc); }
    public void addSelected(Collection<ExptLocator> locs) { 
        for(ExptLocator loc : locs) { addSelected(loc); }
    }
}

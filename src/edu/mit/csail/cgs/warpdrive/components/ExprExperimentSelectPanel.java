/*
 * Created on Nov 6, 2006
 */
package edu.mit.csail.cgs.warpdrive.components;

import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.sql.SQLException;
import java.util.Collection;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.DefaultComboBoxModel;
import javax.swing.JComboBox;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.border.TitledBorder;

import edu.mit.csail.cgs.datasets.expression.Experiment;
import edu.mit.csail.cgs.datasets.expression.ExperimentFilter;
import edu.mit.csail.cgs.datasets.expression.ExpressionLoader;
import edu.mit.csail.cgs.datasets.expression.ExpressionMetadataLoader;
import edu.mit.csail.cgs.datasets.general.Cells;
import edu.mit.csail.cgs.datasets.general.Condition;
import edu.mit.csail.cgs.datasets.general.MetadataLoader;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.viz.components.GenericSelectPanel;

/**
 * @author tdanford
 */
public class ExprExperimentSelectPanel extends GenericSelectPanel<Experiment> {
	
    private ExpressionLoader loader;
    private MetadataLoader metaloader;
    private ExpressionMetadataLoader exprMetaLoader;
    
    private ExprExperimentTableModel selectedModel, filteredModel;
    private DefaultComboBoxModel cellsModel, condModel;
    private JComboBox cellsBox, condBox;
    private JTextField regex;
    private ExperimentFilter filter;
    private Wrapper<Cells> noCells;
    private Wrapper<Condition> noCond;
    private Collection<Cells> allCells;
    private Collection<Condition> allConds;

    public ExprExperimentSelectPanel() {
        super();
        try {
        	loader = new ExpressionLoader();
        	metaloader = new MetadataLoader();
        	exprMetaLoader = new ExpressionMetadataLoader(metaloader);

        	selectedModel = new ExprExperimentTableModel();        
        	filteredModel = new ExprExperimentTableModel();
        	noCells = new Wrapper<Cells>("<NONE>", null);
        	noCond = new Wrapper<Condition>("<NONE>", null);

        	init(filteredModel,selectedModel);
        	
        } catch (SQLException e) {
        	e.printStackTrace();
        }
    }
    
    public JPanel getInputsPanel() {
        JPanel inputPanel = new JPanel(); inputPanel.setLayout(new BorderLayout());		
        JPanel boxPanel = new JPanel(); boxPanel.setLayout(new GridLayout(2, 1));
        JPanel cellsPanel = new JPanel(); cellsPanel.setLayout(new BorderLayout());
        JPanel condPanel = new JPanel(); condPanel.setLayout(new BorderLayout());
        JPanel regexPanel = new JPanel(); regexPanel.setLayout(new BorderLayout());
		
        
        JPanel comboPanel = new JPanel(); comboPanel.setLayout(new GridLayout(1, 2));
        comboPanel.add(cellsPanel);
        comboPanel.add(condPanel);
        
        boxPanel.add(comboPanel);
        boxPanel.add(regexPanel);
		
        cellsModel = new DefaultComboBoxModel();
        condModel = new DefaultComboBoxModel();
        regex = new JTextField();
		
        cellsBox = new JComboBox(cellsModel);
        condBox = new JComboBox(condModel);
		
        cellsPanel.add(cellsBox, BorderLayout.CENTER);
        condPanel.add(condBox, BorderLayout.CENTER);
        regexPanel.add(regex, BorderLayout.CENTER);
		
        cellsPanel.setBorder(new TitledBorder("Cells"));
        condPanel.setBorder(new TitledBorder("Condition"));
        regexPanel.setBorder(new TitledBorder("Match"));
		
        inputPanel.add(boxPanel, BorderLayout.CENTER);
        setBorder(new TitledBorder("Binding Scans:"));
        return inputPanel;
    }

    public ExpressionLoader getLoader() { return loader; }
    
    public void clearSelected() { selectedModel.clear(); }

    public void filter() {
        Cells cells = ((Wrapper<Cells>)(cellsModel.getSelectedItem())).value;
        Condition cond = ((Wrapper<Condition>)(condModel.getSelectedItem())).value;

        String reg = regex.getText().trim();
        Pattern patt = null;
        if(reg.length() > 0) {
            patt = Pattern.compile(reg);
        }
        filteredModel.clear();
        try { 
            filter = new ExperimentFilter(getGenome(), loader);
            Collection<Experiment> expts = filter.findExpts(cells, cond);
            for(Experiment expt : expts) {
                String str = expt.toString();
                Matcher m = patt != null ? patt.matcher(str) : null;
                if(m == null || m.find()) { 
                    filteredModel.addObject(expt);
                }
            }
        } catch(Exception se) { 
            se.printStackTrace(System.err);
        }        
    }
    public void retrieveData() {
        try {
            allCells = new TreeSet<Cells>(exprMetaLoader.getAllCells());
            allConds = new TreeSet<Condition>(exprMetaLoader.getAllConditions());
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    public void updateComponents() {
        selectedModel.clear();
        if(cellsModel.getSize() > 0) { cellsModel.removeAllElements(); } 
        if(condModel.getSize() > 0) { condModel.removeAllElements(); } 

        cellsModel.addElement(noCells);
        condModel.addElement(noCond);
        for(Cells cells : allCells) { 
            Wrapper<Cells> wrapper = new Wrapper<Cells>(cells.getName(), cells);
            cellsModel.addElement(wrapper);
        }
            
        for(Condition cond : allConds) { 
            Wrapper<Condition> wrapper = 
                new Wrapper<Condition>(cond.getName(), cond);
            condModel.addElement(wrapper);
        }
            
        cellsModel.setSelectedItem(noCells);
        condModel.setSelectedItem(noCond);
		
        filter();
    }
    
    public void close() {
        if (!exprMetaLoader.isClosed()) {
        	exprMetaLoader.close();
        }
        if (!metaloader.isClosed()) {
            metaloader.close();
        }
    }
    
    public boolean isClosed() { 
    	return exprMetaLoader.isClosed() && metaloader.isClosed();
    }
    
    private static class Wrapper<X> { 
        public X value;
        public String name;
        public Wrapper(String n, X v) { name = n; value = v; }
        public String toString() { return name; }
        public int hashCode() { return name.hashCode(); }
        public boolean equals(Object o) { 
            if(!(o instanceof Wrapper)) { return false; }
            Wrapper w = (Wrapper)o;
            return w.name.equals(name);
        }
    }

}

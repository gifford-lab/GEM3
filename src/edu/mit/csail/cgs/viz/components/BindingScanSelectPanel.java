/*
 * Created on Nov 6, 2006
 */
package edu.mit.csail.cgs.viz.components;

import java.util.*;
import java.sql.SQLException;
import java.util.regex.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.*;

import edu.mit.csail.cgs.datasets.binding.*;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.chipchip.*;
import edu.mit.csail.cgs.datasets.general.Cells;
import edu.mit.csail.cgs.datasets.general.MetadataLoader;
import edu.mit.csail.cgs.datasets.general.Condition;
import edu.mit.csail.cgs.datasets.general.Factor;
import edu.mit.csail.cgs.datasets.species.Genome;

/**
 * @author tdanford
 */
public class BindingScanSelectPanel extends GenericSelectPanel<BindingScan> {
	
    private BindingScanLoader bindingloader;
    private MetadataLoader metaloader;
    private ChipChipMetadataLoader chipchiploader;
    
    private BindingScanTableModel selectedModel, filteredModel;
    private DefaultComboBoxModel cellsModel, condModel, factorModel;
    private JComboBox cellsBox, condBox, factorBox;
    private JTextField regex;
    private BindingScanFilter filter;
    private Wrapper<Cells> noCells;
    private Wrapper<Condition> noCond;
    private Wrapper<Factor> noFactor;
    private Collection<Cells> allCells;
    private Collection<Condition> allConds;
    private Collection<Factor> allFactors;
    private Collection<BindingScan> scans;

    public BindingScanSelectPanel() 
    	throws SQLException, UnknownRoleException {
        super();
        bindingloader = new BindingScanLoader();
        metaloader = new MetadataLoader();
        chipchiploader = new ChipChipMetadataLoader(metaloader);

        selectedModel = new BindingScanTableModel();        
        filteredModel = new BindingScanTableModel();
        noCells = new Wrapper<Cells>("<NONE>", null);
        noCond = new Wrapper<Condition>("<NONE>", null);
        noFactor = new Wrapper<Factor>("<NONE>", null);	
        scans = new ArrayList<BindingScan>();
        init(filteredModel,selectedModel);
    }
    
    public JPanel getInputsPanel() {
        JPanel inputPanel = new JPanel(); inputPanel.setLayout(new BorderLayout());		
        JPanel boxPanel = new JPanel(); boxPanel.setLayout(new GridLayout(2, 1));
        JPanel cellsPanel = new JPanel(); cellsPanel.setLayout(new BorderLayout());
        JPanel condPanel = new JPanel(); condPanel.setLayout(new BorderLayout());
        JPanel regexPanel = new JPanel(); regexPanel.setLayout(new BorderLayout());
        JPanel factorPanel = new JPanel(); factorPanel.setLayout(new BorderLayout());
		
        
        JPanel comboPanel = new JPanel(); comboPanel.setLayout(new GridLayout(1, 3));
        comboPanel.add(cellsPanel);
        comboPanel.add(condPanel);
        comboPanel.add(factorPanel);
        
        boxPanel.add(comboPanel);
        boxPanel.add(regexPanel);
		
        cellsModel = new DefaultComboBoxModel();
        condModel = new DefaultComboBoxModel();
        factorModel = new DefaultComboBoxModel();
        regex = new JTextField();
		
        cellsBox = new JComboBox(cellsModel);
        condBox = new JComboBox(condModel);
        factorBox = new JComboBox(factorModel);
		
        cellsPanel.add(cellsBox, BorderLayout.CENTER);
        condPanel.add(condBox, BorderLayout.CENTER);
        factorPanel.add(factorBox, BorderLayout.CENTER);
        regexPanel.add(regex, BorderLayout.CENTER);
		
        cellsPanel.setBorder(new TitledBorder("Cells"));
        condPanel.setBorder(new TitledBorder("Condition"));
        factorPanel.setBorder(new TitledBorder("Factor"));
        regexPanel.setBorder(new TitledBorder("Match"));
		
        inputPanel.add(boxPanel, BorderLayout.CENTER);
        setBorder(new TitledBorder("Binding Scans:"));
        return inputPanel;
    }

    public BindingScanLoader getBindingLoader() { return bindingloader; }
    
    public void clearSelected() { selectedModel.clear(); }

    public void filter() {
        Cells cells = ((Wrapper<Cells>)(cellsModel.getSelectedItem())).value;
        Condition cond = ((Wrapper<Condition>)(condModel.getSelectedItem())).value;
        Factor factor = ((Wrapper<Factor>)(factorModel.getSelectedItem())).value;

        filter(cells,cond,factor);
        filteredModel.clear();
        for (BindingScan bs : scans) {
            filteredModel.addObject(bs);
        }
    }
    public void filter(Cells cells, Condition cond, Factor factor) {
        String reg = regex.getText().trim();
        Pattern patt = null;
        if(reg.length() > 0) {
            patt = Pattern.compile(reg);
        }        
        try { 
            synchronized(scans) {
                filter = new BindingScanFilter(getGenome(),bindingloader);
                Collection<BindingScan> tempscans = filter.findScans(cells, cond, factor);        
                scans.clear();
                for(BindingScan bs : tempscans) {
                    String str = bs.toString();
                    Matcher m = patt != null ? patt.matcher(str) : null;
                    if(m == null || m.find()) {             
                        scans.add(bs);
                    }
                }
            }
        } catch(Exception se) { 
            se.printStackTrace(System.err);
        }        
    }

    public void retrieveData() {
        try {
            allCells = new TreeSet<Cells>(chipchiploader.loadAllCells(getGenome()));
            allConds = new TreeSet<Condition>(chipchiploader.loadAllConditions(getGenome()));
            allFactors = new TreeSet<Factor>(chipchiploader.loadAllFactors(getGenome()));
            filter(null,null,null);
        } catch (SQLException e) {
            e.printStackTrace();
        }
        
    }
    public void updateComponents() {
        synchronized(scans) {
            filteredModel.clear();
            for (BindingScan bs : scans) {
                filteredModel.addObject(bs);
            }
        }
        selectedModel.clear();
        if(cellsModel.getSize() > 0) { cellsModel.removeAllElements(); } 
        if(condModel.getSize() > 0) { condModel.removeAllElements(); } 
        if(factorModel.getSize() > 0) { factorModel.removeAllElements(); }		
        cellsModel.addElement(noCells);
        condModel.addElement(noCond);
        factorModel.addElement(noFactor);
        for(Cells cells : allCells) { 
            Wrapper<Cells> wrapper = new Wrapper<Cells>(cells.getName(), cells);
            cellsModel.addElement(wrapper);
        }
            
        for(Condition cond : allConds) { 
            Wrapper<Condition> wrapper = 
                new Wrapper<Condition>(cond.getName(), cond);
            condModel.addElement(wrapper);
        }
            
        for(Factor factor : allFactors){ 
            Wrapper<Factor> wrapper = 
                new Wrapper<Factor>(factor.getName(), factor);
            factorModel.addElement(wrapper);
        }
		
        cellsModel.setSelectedItem(noCells);
        condModel.setSelectedItem(noCond);
        factorModel.setSelectedItem(noFactor);
    }
    
    public void close() {
        if (!bindingloader.isClosed()) {
            bindingloader.close();
        }
        if (!metaloader.isClosed()) {
            metaloader.close();
        }
    }
    
    public boolean isClosed() { 
    	return bindingloader.isClosed() && metaloader.isClosed();
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

/*
 * Created on Nov 6, 2006
 */
package edu.mit.csail.cgs.viz.components;

import java.util.*;
import java.util.regex.*;
import java.sql.SQLException;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.*;

import edu.mit.csail.cgs.datasets.locators.*;
import edu.mit.csail.cgs.datasets.chipchip.*;
import edu.mit.csail.cgs.datasets.general.Cells;
import edu.mit.csail.cgs.datasets.general.MetadataLoader;
import edu.mit.csail.cgs.datasets.general.Condition;
import edu.mit.csail.cgs.datasets.general.Factor;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.species.Genome;


/**
 * @author tdanford
 */
public class ExptSelectPanel extends GenericSelectPanel<ExptLocator> {
    
    private MetadataLoader loader;
    private ChipChipMetadataLoader chipLoader;
    
    private ExptLocatorFilter filter;
    private DefaultComboBoxModel cellsModel, condModel, factorModel;
    private JComboBox cellsBox, condBox, factorBox;
    private JTextField regex;
    private boolean hasAgilent, hasReps, hasMSP, hasBayes;
    private JRadioButton agilentButton, agilentReplicateButton, mspButton, bayesButton;
    private ButtonGroup typeGroup;
    private Wrapper<Cells> noCells;
    private Wrapper<Condition> noCond;
    private Wrapper<Factor> noFactor;
    private ExptTableModel selectedModel, filteredModel;
    private Collection<Cells> allCells;
    private Collection<Condition> allConds;
    private Collection<Factor> allFactors;

    private Collection<ExptLocator> scans;    
    private JButton addDiffButton;

    public ExptSelectPanel(Genome g) {
        super();
        hasAgilent = true;
        hasReps = true;
        hasMSP = true;
        hasBayes = true;
        init(g);       
    }

    /* this constructor lets you determine which experiment type selection buttons are shown.  The choices
       are experiment, experiment w/ replicate, MSP, and Bayes
    */
    public ExptSelectPanel(Genome g,
                           boolean hasExperimentButton,
                           boolean hasReplicateButton,
                           boolean hasMSPButton,
                           boolean hasBayesButton) {
        super();
        hasAgilent = hasExperimentButton;
        hasReps = hasReplicateButton;
        hasMSP = hasMSPButton;
        hasBayes = hasBayesButton;
        init(g);       
    }

    public void init (Genome g) {
        scans = new TreeSet<ExptLocator>();
        try { 
            loader = new MetadataLoader();
            chipLoader = new ChipChipMetadataLoader(loader);
        } catch(Exception e) { 
            e.printStackTrace(System.err);
            throw new RuntimeException(e.getMessage(), e);
        }
        filter = new ExptLocatorFilter();
        if (g != null) {filter.setGenome(g);}
        noCells = new Wrapper<Cells>("<NONE>", null);
        noCond = new Wrapper<Condition>("<NONE>", null);
        noFactor = new Wrapper<Factor>("<NONE>", null);		
        setBorder(new TitledBorder("Selected Experiments:"));
        selectedModel = new ExptTableModel();
        filteredModel = new ExptTableModel();
        
        super.init(filteredModel, selectedModel);

        buttonPanel.add(addDiffButton = new JButton("AddDiff"));
        addDiffButton.addActionListener(new ActionListener() { 
            public void actionPerformed(ActionEvent e) { 
                addDiff();
            }
        });
        setGenome(g);

    }
    
    public void addDiff() { 
        int[] inds = filteredList.getSelectedRows();
        if(inds.length==2) { 
            ExptLocator e1 = filteredModel.getObject(inds[0]);
            ExptLocator e2 = filteredModel.getObject(inds[1]);
            
            if(e1 instanceof ChipChipLocator && e2 instanceof ChipChipLocator) {
                ChipChipLocator al1 = (ChipChipLocator)e1;
                ChipChipLocator al2 = (ChipChipLocator)e2;

                ChipChipLocator dl = new ChipChipDifferenceLocator(getGenome(), al1, al2);
                if(!selectedModel.contains(dl)) { 
                    selectedModel.addObject(dl);
                }
            }
        }
    }
    
    public JPanel getInputsPanel() {
        JPanel inputPanel = new JPanel(); inputPanel.setLayout(new BorderLayout());

        JPanel boxPanel = new JPanel(); boxPanel.setLayout(new GridLayout(3, 1));
        
        JPanel buttonPanel = new JPanel(); buttonPanel.setLayout(new GridLayout(1, 3));
        JPanel cellsPanel = new JPanel(); cellsPanel.setLayout(new BorderLayout());
        JPanel condPanel = new JPanel(); condPanel.setLayout(new BorderLayout());
        JPanel regexPanel = new JPanel(); regexPanel.setLayout(new BorderLayout());
        JPanel factorPanel = new JPanel(); factorPanel.setLayout(new BorderLayout());
        
        JPanel comboPanel = new JPanel(); comboPanel.setLayout(new GridLayout(1, 3));
        comboPanel.add(cellsPanel);
        comboPanel.add(condPanel);
        comboPanel.add(factorPanel);
		
        boxPanel.add(buttonPanel);
        boxPanel.add(comboPanel);
        boxPanel.add(regexPanel);
        
        agilentButton = new JRadioButton("ChipChip");
        agilentReplicateButton = new JRadioButton("ChipChip w/ Replicate");
        mspButton = new JRadioButton("Rosetta");
        bayesButton = new JRadioButton("Bayes");
        typeGroup = new ButtonGroup();
        boolean selected = false;
        if (hasAgilent) {
            typeGroup.add(agilentButton);
            agilentButton.setSelected(true);
            selected = true;
            buttonPanel.add(agilentButton);
        }
        if (hasReps) {
            typeGroup.add(agilentReplicateButton);
            buttonPanel.add(agilentReplicateButton);
            if (!selected) {
                agilentReplicateButton.setSelected(true);
                selected = true;
            }
        }
        if (hasMSP) {
            typeGroup.add(mspButton);
            buttonPanel.add(mspButton);
            if (!selected) {
                mspButton.setSelected(true);
                selected = true;
            }
        }
        if (hasBayes) {
            typeGroup.add(bayesButton);
            buttonPanel.add(bayesButton);
            if (!selected) {
                bayesButton.setSelected(true);
            }

        }
		
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
		
        buttonPanel.setBorder(new TitledBorder("Experiment Type"));
        cellsPanel.setBorder(new TitledBorder("Cells"));
        condPanel.setBorder(new TitledBorder("Condition"));
        factorPanel.setBorder(new TitledBorder("Factor"));
        regexPanel.setBorder(new TitledBorder("Match"));
		
        inputPanel.add(boxPanel, BorderLayout.CENTER);
        return inputPanel;
    }
    /* this filter() only gets called from the swing thread */
    public void filter() {
        Cells cells = ((Wrapper<Cells>)(cellsModel.getSelectedItem())).value;
        Condition cond = ((Wrapper<Condition>)(condModel.getSelectedItem())).value;
        Factor factor = ((Wrapper<Factor>)(factorModel.getSelectedItem())).value;
        filter(cells,cond,factor);

        filteredModel.clear();
        synchronized(scans) {
            for(ExptLocator bs : scans) {
                filteredModel.addObject(bs);
            }
        }
    }
    /* this gets called from anywhere and just updates scans */
    public void filter(Cells cells, Condition cond, Factor factor) {
        synchronized(scans) {
            scans.clear();
            try { 
                ExptLocatorFilterOptions opts = new ExptLocatorFilterOptions();
                // use this as the default in case nothing is set yet
                opts.setChipChipType();
                if(mspButton.isSelected()) { opts.setMSPType(); }
                if(bayesButton.isSelected()) { opts.setBayesType(); }
                if(cells != null) { opts.setCells(cells); } 
                if(cond != null) { opts.setCondition(cond); }
                if(factor != null) { opts.setFactor(factor); }
            
                Collection<ExptLocator> scanstemp = filter.findLocators(opts);
                String reg = regex.getText().trim();
                Pattern patt = null;
                if(reg != null && reg.length() > 0) {
                    patt = Pattern.compile(reg);
                }
                for(ExptLocator bs : scanstemp) {
                    if (agilentButton.isSelected() && bs instanceof ChipChipLocator) {
                        ((ChipChipLocator)bs).clearReplicate();
                    }
                    String str = bs.toString();
                    Matcher m = patt != null ? patt.matcher(str) : null;
                    if(m == null || m.find()) {                  
                        scans.add(bs);
                    }
                }
            } catch(SQLException se) { 
                se.printStackTrace(System.err);
            }        
        }
    }
    
    public void retrieveData() {
        try {
            filter.setGenome(getGenome());
            allCells = new TreeSet<Cells>();
            allConds = new TreeSet<Condition>();
            allFactors = new TreeSet<Factor>();
            allCells.addAll(chipLoader.loadAllCells(getGenome()));
            allConds.addAll(chipLoader.loadAllConditions(getGenome()));
            allFactors.addAll(chipLoader.loadAllFactors(getGenome()));
            filter(null,null,null);
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }
    public void updateComponents() {
        synchronized(scans) {
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
		
            for(ExptLocator bs : scans) {
                filteredModel.addObject(bs);
            }
            if (agilentButton.isSelected()) {
                for(ExptLocator bs : scans) {
                    if (bs instanceof ChipChipLocator) {
                        ((ChipChipLocator)bs).clearReplicate();
                    }
                }
            }
        }
    }
        
    public Collection<ExptLocator> getSelected() { 
        LinkedList<ExptLocator> scans = new LinkedList<ExptLocator>();
        for(int i = 0; i < selectedModel.getSize(); i++) { 
            scans.addLast(selectedModel.getObject(i));
        }
        return scans;
    }
    
    public void close() {
        if (!loader.isClosed()) {
            loader.close();
        }
        if (!filter.isClosed()) {
            filter.close();
        }
    }
    
    public boolean isClosed() { 
        return loader.isClosed() && filter.isClosed();
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

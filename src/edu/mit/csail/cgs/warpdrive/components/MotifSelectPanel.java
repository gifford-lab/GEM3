package edu.mit.csail.cgs.warpdrive.components;

import java.util.*;
import java.util.regex.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.*;
import javax.swing.table.*;

import edu.mit.csail.cgs.datasets.motifs.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.viz.components.GenericSelectPanel;
import edu.mit.csail.cgs.viz.components.SelectionEvent;

public class MotifSelectPanel extends GenericSelectPanel<WeightMatrix> {

    private WeightMatrixLoader loader;
    private MotifTableModel filteredModel, selectedModel;
    private JComboBox nameBox, versionBox, typeBox;
    private JTextField regex;
    private TreeSet<String> names, versions, types;
    private DefaultComboBoxModel nameModel, versionModel, typeModel;

    public MotifSelectPanel() {
        super();
        loader = new WeightMatrixLoader();
        filteredModel = new MotifTableModel();
        selectedModel = new MotifTableModel();
        init(filteredModel,selectedModel);
        names = new TreeSet<String>();
        versions = new TreeSet<String>();
        types = new TreeSet<String>();
    }
    public JPanel getInputsPanel() {
        JPanel inputsPanel = new JPanel(); inputsPanel.setLayout(new BorderLayout());
        JPanel namePanel = new JPanel(); namePanel.setLayout(new BorderLayout());
        JPanel versionPanel = new JPanel(); versionPanel.setLayout(new BorderLayout());
        JPanel typePanel = new JPanel(); typePanel.setLayout(new BorderLayout());
        JPanel regexPanel = new JPanel(); regexPanel.setLayout(new BorderLayout());
        JPanel comboPanel = new JPanel(); comboPanel.setLayout(new GridLayout(2,2));
        comboPanel.add(namePanel);
        comboPanel.add(versionPanel);
        comboPanel.add(typePanel);
        comboPanel.add(regexPanel);        
        

        nameModel = new DefaultComboBoxModel();
        versionModel = new DefaultComboBoxModel();
        typeModel = new DefaultComboBoxModel();
        regex = new JTextField();
        nameModel.addElement("<None>");
        versionModel.addElement("<None>");
        typeModel.addElement("<None>");
        nameBox = new JComboBox(nameModel);
        versionBox = new JComboBox(versionModel);
        typeBox = new JComboBox(typeModel);
        Dimension d = new Dimension(100,30);
        nameBox.setPreferredSize(d);
        versionBox.setPreferredSize(d);
        typeBox.setPreferredSize(d);
        namePanel.add(nameBox,BorderLayout.CENTER);
        versionPanel.add(versionBox,BorderLayout.CENTER);
        typePanel.add(typeBox,BorderLayout.CENTER);
        regexPanel.add(regex,BorderLayout.CENTER);
        namePanel.setBorder(new TitledBorder("Name"));
        versionPanel.setBorder(new TitledBorder("Version"));
        typePanel.setBorder(new TitledBorder("Type"));        
        regexPanel.setBorder(new TitledBorder("Match Motif"));
        inputsPanel.add(comboPanel,BorderLayout.CENTER);
        setBorder(new TitledBorder("Weight Matrices:"));
        return inputsPanel;
    }

    public void filter() {
        filteredModel.clear();
        String n, v, t;
        n = (String)nameBox.getSelectedItem();
        v = (String)versionBox.getSelectedItem();
        t = (String)typeBox.getSelectedItem();
        if (n.equals("<None>")) {n = null;}
        if (v.equals("<None>")) {v = null;}
        if (t.equals("<None>")) {t = null;}
        Pattern patt = null;
        String reg = regex.getText().trim();
        if(reg != null && reg.length() > 0) {
            patt = Pattern.compile(reg);
        }
        for (WeightMatrix wm : loader.query(n,v,t)) {
            if (patt == null || patt.matcher(wm.toString()).find()) {
                filteredModel.addObject(wm);
            }
        }
    }
    public void close() {super.close(); loader.close();}
    public boolean isClosed() {return super.isClosed() && loader.isClosed();}
    public void retrieveData() { 
        names.clear();
        names.addAll(loader.getNames());
        versions.clear();
        versions.addAll(loader.getVersions());
        types.clear();
        types.addAll(loader.getTypes());
    }
    public void updateComponents() { 
        nameModel.removeAllElements();
        nameModel.addElement("<None>");
        for (String s : names) {
            nameModel.addElement(s);
        }

        versionModel.removeAllElements();
        versionModel.addElement("<None>");
        for (String s : versions) {
            versionModel.addElement(s);
        }

        typeModel.removeAllElements();
        typeModel.addElement("<None>");
        for (String s : types) {
            typeModel.addElement(s);
        }

    }
}

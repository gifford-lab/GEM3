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

public class MotifScanSelectPanel extends GenericSelectPanel<WeightMatrixScan> {

    private MotifScanTableModel filteredModel, selectedModel;
    private DefaultComboBoxModel nameModel, versionModel, typeModel, scannameModel;
    private JComboBox nameBox, versionBox, typeBox, scannameBox;
    private JTextField regex;
    private TreeSet<String> names, versions, types, scannames;
    private Collection<WeightMatrixScan> scans;

    public MotifScanSelectPanel() {
        super();
        filteredModel = new MotifScanTableModel();
        selectedModel = new MotifScanTableModel();
        scans = new ArrayList<WeightMatrixScan>();
        super.init(filteredModel,selectedModel);
    }
    public JPanel getInputsPanel() {
        JPanel namePanel = new JPanel(); namePanel.setLayout(new BorderLayout());
        JPanel versionPanel = new JPanel(); versionPanel.setLayout(new BorderLayout());
        JPanel typePanel = new JPanel(); typePanel.setLayout(new BorderLayout());
        JPanel scannamePanel = new JPanel(); scannamePanel.setLayout(new BorderLayout());
        JPanel regexPanel = new JPanel(); regexPanel.setLayout(new BorderLayout());
        GridBagLayout gridbag = new GridBagLayout();
        JPanel comboPanel = new JPanel(); comboPanel.setLayout(gridbag);

        nameModel = new DefaultComboBoxModel();
        versionModel = new DefaultComboBoxModel();
        typeModel = new DefaultComboBoxModel();
        scannameModel = new DefaultComboBoxModel();
        regex = new JTextField();
        nameModel.addElement("<None>");
        versionModel.addElement("<None>");
        typeModel.addElement("<None>");
        scannameModel.addElement("<None>");
        nameBox = new JComboBox(nameModel);
        versionBox = new JComboBox(versionModel);
        typeBox = new JComboBox(typeModel);
        scannameBox = new JComboBox(scannameModel);
        Dimension d = new Dimension(100,30);
        nameBox.setPreferredSize(d);
        versionBox.setPreferredSize(d);
        Dimension d2 = new Dimension(60,30);
        typeBox.setPreferredSize(d2);
        scannameBox.setPreferredSize(d2);
        namePanel.add(nameBox,BorderLayout.CENTER);
        versionPanel.add(versionBox,BorderLayout.CENTER);
        typePanel.add(typeBox,BorderLayout.CENTER);
        scannamePanel.add(scannameBox,BorderLayout.CENTER);
        regexPanel.add(regex, BorderLayout.CENTER);
        namePanel.setBorder(new TitledBorder("Name"));
        versionPanel.setBorder(new TitledBorder("Version"));
        typePanel.setBorder(new TitledBorder("Type"));     
        scannamePanel.setBorder(new TitledBorder("Scan Name"));
        setBorder(new TitledBorder("Weight Matrices:"));
        regexPanel.setBorder(new TitledBorder("Match Motif Scan"));


        GridBagConstraints c = new GridBagConstraints();
        c.weightx = 1.0;
        c.weighty = 1.0;
        c.fill = GridBagConstraints.BOTH;

        c.gridwidth = 1;
        gridbag.setConstraints(namePanel,c);
        comboPanel.add(namePanel);        

        c.gridwidth = GridBagConstraints.REMAINDER;
        gridbag.setConstraints(versionPanel,c);
        comboPanel.add(versionPanel);

        c.gridwidth = 1;
        gridbag.setConstraints(typePanel,c);
        comboPanel.add(typePanel);

        c.gridwidth = 1;
        gridbag.setConstraints(scannamePanel,c);
        comboPanel.add(scannamePanel);

        c.gridwidth = GridBagConstraints.REMAINDER;
        gridbag.setConstraints(regexPanel,c);
        comboPanel.add(regexPanel);        

        return comboPanel;
    }
    public void filter() {
        if (getGenome() == null) {return;}
        filteredModel.clear();
        String n, v, t, s;
        n = (String)nameBox.getSelectedItem();
        v = (String)versionBox.getSelectedItem();
        t = (String)typeBox.getSelectedItem();
        s = (String)scannameBox.getSelectedItem();
        if (n.equals("<None>")) {n = null;}
        if (v.equals("<None>")) {v = null;}
        if (t.equals("<None>")) {t = null;}
        if (s.equals("<None>")) {s = null;}
        filter(n,v,t,s);
        for (WeightMatrixScan wm : scans) {
            filteredModel.addObject(wm);
        }
    }
    public void filter(String n, String v, String t, String s) {
        scans.clear();
        String reg = regex.getText().trim();
        if(reg != null && reg.length() > 0) {
            Pattern patt = null;
            patt = Pattern.compile(reg);
            for (WeightMatrixScan scan : WeightMatrixScan.getScansForGenome(getGenome().getDBID(),n,v,t,s)) {
                if (patt.matcher(scan.toString()).find()) {
                    scans.add(scan);
                }
            }
        } else {
            scans.addAll(WeightMatrixScan.getScansForGenome(getGenome().getDBID(),n,v,t,s));
        }
    }

    public void retrieveData() {
        names = new TreeSet<String>();
        versions = new TreeSet<String>();
        types = new TreeSet<String>();
        scannames = new TreeSet<String>();
        names.addAll(WeightMatrixScan.getNames(getGenome().getDBID()));
        versions.addAll(WeightMatrixScan.getVersions(getGenome().getDBID()));
        types.addAll(WeightMatrixScan.getTypes(getGenome().getDBID()));
        scannames.addAll(WeightMatrixScan.getScanNames(getGenome().getDBID()));
    }
    public void updateComponents() {
        nameModel.removeAllElements();
        nameModel.addElement("<None>");
        for (String name : names) {
            nameModel.addElement(name);
        }
        versionModel.removeAllElements();
        versionModel.addElement("<None>");
        for (String version : versions) {
            versionModel.addElement(version);
        }
        typeModel.removeAllElements();
        typeModel.addElement("<None>");
        for (String type : types) {
            typeModel.addElement(type);
        }
        scannameModel.removeAllElements();
        scannameModel.addElement("<None>");
        for (String name: scannames) {
            scannameModel.addElement(name);
        }
        for (WeightMatrixScan wm : scans) {
            filteredModel.addObject(wm);
        }
    }
}

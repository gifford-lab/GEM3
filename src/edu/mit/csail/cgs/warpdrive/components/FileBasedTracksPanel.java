package edu.mit.csail.cgs.warpdrive.components;

import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.io.File;


public class FileBasedTracksPanel extends JPanel implements ActionListener {

    private java.util.List<JTextField> fnames, labels;
    private java.util.List<JButton> browseButtons;
    private GridBagLayout layout;

    public FileBasedTracksPanel() {
        super();
        fnames = new ArrayList<JTextField>();
        labels = new ArrayList<JTextField>();
        browseButtons = new ArrayList<JButton>();
        layout = new GridBagLayout();
        setLayout(layout);

        GridBagConstraints constraints = new GridBagConstraints();
        constraints.weightx = 1.0;        
        constraints.fill = GridBagConstraints.BOTH;
        constraints.gridwidth = 1;
        JLabel jl = new JLabel("File Name");
        layout.setConstraints(jl,constraints);
        add(jl);
        jl = new JLabel("");
        layout.setConstraints(jl,constraints);
        add(jl);
        jl = new JLabel("Track Label");
        constraints.fill = GridBagConstraints.BOTH;
        constraints.gridwidth = GridBagConstraints.REMAINDER;
        layout.setConstraints(jl,constraints);
        add(jl);
        for (int i = 0; i < 6; i++) {
            addRow();
        }        
    }

    public void addRow() {
        GridBagConstraints constraints = new GridBagConstraints();
        constraints.weightx = 1.0;        
        JTextField namefield = new JTextField("");
        JTextField labelfield = new JTextField("");
        JButton button = new JButton("Browse");
        fnames.add(namefield);
        labels.add(labelfield);
        browseButtons.add(button);
        button.addActionListener(this);
        namefield.setColumns(20);
        labelfield.setColumns(10);

        constraints.fill = GridBagConstraints.BOTH;
        constraints.gridwidth = 1;
        layout.setConstraints(namefield,constraints);
        add(namefield);
        constraints.fill = GridBagConstraints.BOTH;
        constraints.gridwidth = 1;
        layout.setConstraints(button,constraints);
        add(button);
        constraints.fill = GridBagConstraints.BOTH;
        constraints.gridwidth = GridBagConstraints.REMAINDER;
        layout.setConstraints(labelfield,constraints);
        add(labelfield);            
        
        invalidate();
    }
    
    public void actionPerformed(ActionEvent e) {
        int index = browseButtons.indexOf(e.getSource());
        if (index == -1) {
            return;
        }       
        JFileChooser chooser;
        chooser = new JFileChooser(new File(System.getProperty("user.dir")));
        int v = chooser.showOpenDialog(null);
        if(v == JFileChooser.APPROVE_OPTION) { 
            File f = chooser.getSelectedFile();
            fnames.get(index).setText(f.getAbsolutePath());
        }
    }

    /* parses the frame into the specified map and returns it.  If
       the input is null, a new map is created */
    public Map<String,String> parse(Map<String,String> map) {
        if (map == null) {
            map = new HashMap<String,String>();
        }
        for (int i = 0; i < fnames.size(); i++) {
            String name = fnames.get(i).getText();
            String label = labels.get(i).getText();
            if (name.length() > 0) {
                if (label.length() == 0) {
                    label = name;
                }
                map.put(name, label);
            }
        }
        return map;
    }

    public void fill(Map<String,String> map) {
        int index = 0;
        for (String k : map.keySet()) {
            if (index >= fnames.size()) {
                addRow();
            }
            fnames.get(index).setText(k);
            labels.get(index).setText(map.get(k));
            index++;
        }

    }

    

    
}
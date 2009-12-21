/*
 * Created on Oct 20, 2005
 */
package edu.mit.csail.cgs.datasets.orthology;

import java.util.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.border.*;
import javax.swing.table.*;

import edu.mit.csail.cgs.utils.*;

/**
 * @author tdanford
 */
public class OrthologyPanel extends JPanel {
    
    private TotalOrthologyMapping totalMapping;
    private JList list;
    private DefaultListModel model;
    private JLabel mappingLabel;

    public OrthologyPanel(TotalOrthologyMapping tom) {
        super();
        totalMapping = tom;
        setLayout(new BorderLayout());
        
        model = new DefaultListModel();
        list = new JList(model);
        
        add(new JScrollPane(list), BorderLayout.CENTER);
        mappingLabel = new JLabel("");
        add(mappingLabel, BorderLayout.NORTH);
        
        update();
    }

    public void update() {
        OrthologyMapping base = totalMapping.getMapping();
        int size = totalMapping.size();
        mappingLabel.setText(base.getName() + "/" + base.getVersion() + ": " + size + " pairs");
        model.clear();
        for(OrthologyPair p : totalMapping.getPairs()) { 
            model.addElement(p);
        }
    }
}

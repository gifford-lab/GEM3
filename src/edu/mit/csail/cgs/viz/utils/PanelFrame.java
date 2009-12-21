/*
 * Created on Oct 20, 2005
 */
package edu.mit.csail.cgs.viz.utils;

import java.util.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;

/**
 * @author tdanford
 */
public class PanelFrame extends JFrame {

    private JPanel panel;
    
    public PanelFrame(String title, JPanel p) {
        super(title);
        panel = p;
        Container c = (Container)getContentPane();
        c.setLayout(new BorderLayout());
        
        c.add(panel, BorderLayout.CENTER);
        
        setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        setVisible(true);
        pack();
        setLocation(getX() + 100, getY() + 100);
    }

    public JPanel getPanel() { return panel; }
}

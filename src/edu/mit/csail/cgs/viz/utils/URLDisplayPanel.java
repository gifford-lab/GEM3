/*
 * Created on Oct 21, 2005
 */
package edu.mit.csail.cgs.viz.utils;

import java.util.*;
import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

/**
 * @author tdanford
 */
public class URLDisplayPanel extends JPanel {
    
    private JTextField text;
    private Action copyAction;
    private JButton copyButton;

    public URLDisplayPanel(String urlString) {
        super();
        setLayout(new BorderLayout());
        text = new JTextField(urlString);
        text.setEditable(false);
        add(text, BorderLayout.CENTER);
        
        copyAction = new AbstractAction("Copy URL") {
            public void actionPerformed(ActionEvent arg0) {
                copyURL();
            } 
        };
        copyButton = new JButton(copyAction);
        add(copyButton, BorderLayout.SOUTH);
    }

    public void copyURL() { 
        text.selectAll();
        text.copy();
    }

    public String getURLString() { return text.getText(); }
}

package edu.mit.csail.cgs.warpdrive.components;

import java.util.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;


public class HashtableConfigurationFrame extends JFrame implements ActionListener {
    JButton ok, cancel;
    HashtableConfigurationPanel panel;
    Component comp;
    public HashtableConfigurationFrame(String name, HashtableConfigurationPanel p, Component comp) {
        panel = p;
        Container cf = getContentPane();
        this.comp = comp;
        JScrollPane scrollpane = new JScrollPane(p);
        scrollpane.setPreferredSize(p.getPreferredSize());
        cf.setLayout(new BorderLayout());
        cf.add(scrollpane,BorderLayout.CENTER);
        cf.add(new JLabel("Configuring " + name),BorderLayout.NORTH);
        cf.add(scrollpane);
        JPanel buttonPanel = new JPanel();
        buttonPanel.setLayout(new GridBagLayout());
        Dimension buttonSize = new Dimension(40,25);
        ok = new JButton("OK");
        cancel = new JButton("Cancel");
        ok.setMaximumSize(buttonSize);
        cancel.setMaximumSize(buttonSize);
        ok.addActionListener(this);
        cancel.addActionListener(this);
        buttonPanel.add(ok);
        buttonPanel.add(cancel);
        cf.add(buttonPanel,BorderLayout.SOUTH);
    }

    public void actionPerformed (ActionEvent e) {
        comp.repaint();
        if (e.getSource() == ok) {
            panel.parse();
            this.dispose();
        } else if (e.getSource() == cancel) {
            panel.clearConfiguring();
            this.dispose();
        }
    }
}

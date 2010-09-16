package edu.mit.csail.cgs.warpdrive.components;

import java.util.*;
import java.util.regex.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;
import java.sql.*;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.viz.components.GenericSelectPanel;
import edu.mit.csail.cgs.tools.utils.Args;

/**
 * This is a GUI app for managing ChipSeqAnalysis objects that lets you set their
 * active/inactive flag
 */


public class ChipSeqAnalysisFrame extends JFrame implements ActionListener {

    private ChipSeqAnalysisSelectPanel panel;
    private JPanel buttonPanel;
    private JButton activeButton, inactiveButton;
    private JMenuBar menuBar;
    private JMenu fileMenu;
    private JMenuItem exitItem;

    public ChipSeqAnalysisFrame(Genome g) {
        super();
        panel = new ChipSeqAnalysisSelectPanel(g);

        JPanel buttonPanel = new JPanel();
        buttonPanel.setLayout(new GridBagLayout());
        Dimension buttonSize = new Dimension(30,20);
        activeButton = new JButton("Set as Active");
        inactiveButton = new JButton("Set as Inactive");
        activeButton.setMaximumSize(buttonSize);
        inactiveButton.setMaximumSize(buttonSize);
        buttonPanel.add(activeButton);
        buttonPanel.add(inactiveButton);
        activeButton.addActionListener(this);
        inactiveButton.addActionListener(this);
        Container content = getContentPane();
        content.setLayout(new BorderLayout());
        content.add(buttonPanel,BorderLayout.SOUTH);
        content.add(panel,BorderLayout.CENTER);

        menuBar = new JMenuBar();
    	fileMenu = new JMenu("File");
    	fileMenu.setMnemonic(KeyEvent.VK_F);
    	menuBar.add(fileMenu);
    	exitItem = new JMenuItem("Exit", KeyEvent.VK_X);
    	exitItem.addActionListener(new ActionListener() {
    		public void actionPerformed(ActionEvent e) {
    			exit_actionPerformed(e);
    		}
    	});
    	fileMenu.add(exitItem);
    	this.setJMenuBar(menuBar);


        setSize(600,500);
        setVisible(true);
    }
    public void actionPerformed (ActionEvent e) {
        if (e.getSource() == activeButton) {
            for (ChipSeqAnalysis a : panel.getSelected()) {
                a.setActive(true);
                try {
                    a.storeActiveDB();   
                } catch (SQLException ex) {
                    ex.printStackTrace();
                }
            }
        } else if (e.getSource() == inactiveButton) {
            for (ChipSeqAnalysis a : panel.getSelected()) {
                a.setActive(false);
                try {
                    a.storeActiveDB();   
                } catch (SQLException ex) {
                    ex.printStackTrace();
                }
            }
        }
    }
    void exit_actionPerformed(ActionEvent e) {
        panel.close();
        try {
            Thread.sleep(400);
        } catch (Exception ex) {

        }
        System.exit(0);
    }
    public static void main(String args[]) throws Exception {
        ChipSeqAnalysisFrame f = new ChipSeqAnalysisFrame(Args.parseGenome(args).cdr());
    }

}
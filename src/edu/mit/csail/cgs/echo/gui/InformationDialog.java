package edu.mit.csail.cgs.echo.gui;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.util.*;

public class InformationDialog<X> extends JDialog {

	private JComponent info;
	private JButton ok;

	public InformationDialog(JComponent comp) { 
		super();
		
		Container c = (Container)getContentPane();
		c.setLayout(new BorderLayout());
		
		JPanel buttonPanel = new JPanel(); buttonPanel.setLayout(new FlowLayout());
		buttonPanel.add(ok = new JButton("Ok"));
		
		c.add(buttonPanel, BorderLayout.SOUTH);
		c.add(comp, BorderLayout.CENTER);
		
		ok.addActionListener(new ActionListener() { 
			public void actionPerformed(ActionEvent e) { 
				setVisible(false);
				dispose();
			}
		});
        
		setVisible(true);
		pack();
	}
}

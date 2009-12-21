package edu.mit.csail.cgs.echo.gui;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import java.util.*;

public class StringSelectionPanel extends JPanel implements SwingSelectionComponent<String> {
	
	private JTextField entryField;
	
	public StringSelectionPanel(String a) {
		super();
		setLayout(new BorderLayout());
		add(entryField = new JTextField(a != null ? a : ""), BorderLayout.NORTH);
	}

	public JComponent asJComponent() {
		return this;
	}

	public String getSelectedValue() {
		String txt = entryField.getText();
		return txt;
	}
}

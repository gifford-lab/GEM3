package edu.mit.csail.cgs.echo.gui;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import java.util.*;

public class StringArraySelectionPanel extends JPanel implements SwingSelectionComponent<String> {
	
	private ButtonGroup group;
	private JRadioButton[] buttons;
	private String[] array;
	
	public StringArraySelectionPanel(String[] a) {
		super();
		array = (String[])a.clone();
		setLayout(new GridLayout(array.length, 1));
		group = new ButtonGroup();
		buttons = new JRadioButton[array.length];
		for(int i = 0; i < buttons.length; i++) { 
			buttons[i] = new JRadioButton(array[i]);
			group.add(buttons[i]);
			if(i == 0) { buttons[i].setSelected(true); }
			add(buttons[i]);
		}
	}

	public JComponent asJComponent() {
		return this;
	}

	public String getSelectedValue() {
		for(int i = 0; i < array.length; i++) { 
			if(buttons[i].isSelected()) { return array[i]; }
		}
		return null;
	}
}

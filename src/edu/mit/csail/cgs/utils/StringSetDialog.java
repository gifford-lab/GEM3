package edu.mit.csail.cgs.utils;

import java.lang.*;
import java.io.*;
import java.util.*;

import javax.swing.*;
import javax.swing.event.*;

import java.awt.*;
import java.awt.event.*;

public class StringSetDialog extends JDialog { 

    private String[] fValues;
    
    public StringSetDialog(String[] array) { 
	super();
	fValues = (String[])array.clone();
	Arrays.sort(fValues); 

	Container c = getContentPane();
	c.setLayout(new BorderLayout());

	Vector<String> vec = new Vector<String>();
	for(int i = 0; i < fValues.length; i++) { 
	    vec.add(fValues[i]); 
	}
	JList lst = new JList(vec);
	c.add(new JScrollPane(lst), BorderLayout.CENTER);

	setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
	setJMenuBar(createMenuBar());

	pack();
    }

    public StringSetDialog(String[] array, boolean sortValues) { 
	super();
	fValues = (String[])array.clone();
	if(sortValues) { Arrays.sort(fValues); } 

	Container c = getContentPane();
	c.setLayout(new BorderLayout());

	Vector<String> vec = new Vector<String>();
	for(int i = 0; i < fValues.length; i++) { 
	    vec.add(fValues[i]); 
	}
	JList lst = new JList(vec);
	c.add(new JScrollPane(lst), BorderLayout.CENTER);

	setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
	setJMenuBar(createMenuBar());

	pack();
    }

    private JMenuBar createMenuBar() { 
	JMenuBar jmb = new JMenuBar();
	JMenu menu; JMenuItem item;

	jmb.add((menu = new JMenu("File")));
	menu.add((item = new JMenuItem("Close")));
	item.addActionListener(new ActionListener() { 
		public void actionPerformed(ActionEvent e) { 
		    dispose();
		}
	    });
	
	return jmb;
    }
}


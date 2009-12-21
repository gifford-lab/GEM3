
package edu.mit.csail.cgs.echo.components;

import java.util.*;
import java.io.*;
import java.awt.*;

import javax.swing.*;

import edu.mit.csail.cgs.echo.gui.InformationDialog;
import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.Sink;

public class ListingSink implements SelfDescribingVerb, Sink { 

	private static final EchoType[] paramClassArray = null;
	private static final String[] paramNameArray = null;
    
    private static final EchoType[] inputClasses = { EchoType.OBJECT_CLASS };
    private static final String[] inputNames = { "Objects" };

	private Map<String,Object> params;
	private int count;
	private DefaultListModel model;
    
	public ListingSink() { 
		params = new HashMap<String,Object>();
		init();
	}
	
	public void init(Map pms) { 
		params = new HashMap<String,Object>(pms); 
		init();
	}

	public void init() {
		count = 0;
		model = new DefaultListModel(); 
	}
	
	public EchoType[] getParameterClasses() { return paramClassArray; }
	public String[] getParameterNames() { return paramNameArray; }
    
    public EchoType[] getInputClasses() { return inputClasses; }
    public String[] getInputNames() { return inputNames; }
    
    public EchoType getOutputClass() { return null; }
    
	public void consume(Object val) {
		model.addElement(val);
		count++;
	}
    
    public void consume(Iterator itr) { 
		init();
        while(itr.hasNext()) { 
			consume(itr.next());
        }
		finish();
    }

	public void finish() {
        JList list = new JList(model);
        JPanel p = new JPanel(); p.setLayout(new BorderLayout());
        JLabel lbl = new JLabel("# Items: " + count);
        p.add(new JScrollPane(list), BorderLayout.CENTER);
        p.add(lbl, BorderLayout.NORTH);
        
        new InformationDialog(p);
	}
}



package edu.mit.csail.cgs.echo.components;

import java.util.*;
import java.io.*;

import javax.swing.JLabel;

import edu.mit.csail.cgs.echo.gui.InformationDialog;
import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.Sink;

public class CountingSink implements SelfDescribingVerb, Sink { 

	private static final EchoType[] paramClassArray = null;
	private static final String[] paramNameArray = null;
    
    private static final String[] inputNames = { "Objects" };
    private static final EchoType[] inputClasses = { EchoType.OBJECT_CLASS };

	private Map<String,Object> params;
	private int count;
	
	public CountingSink() { 
		params = new HashMap<String,Object>();
		init();
	}
	
	public void init(Map<String,Object> pms) { 
		params = new HashMap<String,Object>(pms); 
		init();
	}
	
	public EchoType[] getParameterClasses() { return paramClassArray; }
	public String[] getParameterNames() { return paramNameArray; }

    public EchoType[] getInputClasses() { return inputClasses; }
    public String[] getInputNames() { return inputNames; }
    
    public EchoType getOutputClass() { return null; }

	public void init() { count = 0; }

	public void consume(Iterator itr) {
		init();
		while(itr.hasNext()) {
			consume(itr.next());
		}
		finish();
	}

	public void consume(Object v) { 
		count += 1;
	}

	public void finish() {
		JLabel lbl = new JLabel("# Items: " + count);
		new InformationDialog(lbl);
	}
}


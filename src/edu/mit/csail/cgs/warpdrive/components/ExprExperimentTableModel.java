/*
 * Created on Nov 20, 2006
 */
package edu.mit.csail.cgs.warpdrive.components;

import javax.swing.event.*;
import javax.swing.table.*;
import javax.swing.*;
import java.util.*;
import edu.mit.csail.cgs.datasets.expression.*;
import edu.mit.csail.cgs.viz.components.ObjectTableModel;

public class ExprExperimentTableModel extends ObjectTableModel<Experiment> {
	
	private boolean sortByName;
	
	public ExprExperimentTableModel() { 
		sortByName = true;
	}
    
    public int getColumnCount() {
        return 1;
    }
    
    public Class getColumnClass(int i) {
        if(i==0) { return String.class; }
        if(i==1) { return String.class; }
        return null;
    }
    
    public String getColumnName(int i) {
        if(i==0) { return "Name"; }
        return null;
    }

    public Object getValueAt(int rowIndex, int c) {
        if(c==0) { return getObject(rowIndex).getName(); }
        return null;
    }

    public void sortByName() {
        sort(new ExperimentNameComparator());
        sortByName = true;
    }
    
}

class ExperimentNameComparator implements Comparator<Experiment> {
    public int compare(Experiment a, Experiment b) {
        return a.getName().compareTo(b.getName());
    }
}


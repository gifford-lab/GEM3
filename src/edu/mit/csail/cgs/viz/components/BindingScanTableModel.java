/*
 * Created on Nov 20, 2006
 */
package edu.mit.csail.cgs.viz.components;

import javax.swing.event.*;
import javax.swing.table.*;
import javax.swing.*;
import java.util.*;
import edu.mit.csail.cgs.datasets.binding.*;

public class BindingScanTableModel extends ObjectTableModel<BindingScan> {
    
    private boolean sortByVersion = true;
    
    public int getColumnCount() {
        return 2;
    }
    public Class getColumnClass(int i) {
        if(i==0) { return String.class; }
        if(i==1) { return String.class; }
        return null;
    }
    public String getColumnName(int i) {
        if(i==0) { return "Version"; }
        if(i==1) { return "Type"; }
        return null;
    }

    public Object getValueAt(int rowIndex, int c) {
        if(c==0) { return getObject(rowIndex).getVersion(); }
        if(c==1) { return getObject(rowIndex).getType(); }
        return null;
    }

    public void sortByVersion() {
        sort(new BindingScanVersionComparator());
        sortByVersion = true;
    }
    
    public void sortByType() {
        sort(new BindingScanTypeComparator());
        sortByVersion = false;
    }    
}

class BindingScanVersionComparator implements Comparator<BindingScan> {
    public int compare(BindingScan a, BindingScan b) {
        return a.getVersion().compareTo(b.getVersion());
    }
}
class BindingScanTypeComparator implements Comparator<BindingScan> {
    public int compare(BindingScan a, BindingScan b) {
        return a.getType().compareTo(b.getType());
    }
}


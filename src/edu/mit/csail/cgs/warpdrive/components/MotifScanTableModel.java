package edu.mit.csail.cgs.warpdrive.components;

import javax.swing.event.*;
import javax.swing.table.*;
import javax.swing.*;
import java.util.*;

import edu.mit.csail.cgs.datasets.motifs.WeightMatrixScan;
import edu.mit.csail.cgs.viz.components.ObjectTableModel;

public class MotifScanTableModel extends ObjectTableModel<WeightMatrixScan> {

    private boolean sortByVersion = false;

    public int getColumnCount() {
        return 4;
    }
    public Class getColumnClass(int i) {
        if(i==0) { return String.class; }
        if(i==1) { return String.class; }
        if(i==2) { return String.class; }
        if(i==3) { return String.class; }
        return null;
    }
    public String getColumnName(int i) {
        if(i==0) { return "Name";}
        if(i==1) { return "Version"; }
        if(i==2) { return "Type"; }
        if (i==3) {return "Scan Name";}
        return null;
    }

    public Object getValueAt(int rowIndex, int c) {
        if(c==0) { return getObject(rowIndex).matrix.name; }
        if(c==1) { return getObject(rowIndex).matrix.version; }
        if(c==2) { return getObject(rowIndex).matrix.type; }
        if(c==3) { return getObject(rowIndex).scanname; }
        return null;
    }

    public void sortByVersion() {
        sort(new WeightMatrixScanVersionComparator());
        sortByVersion = true;
    }
    
    public void sortByName() {
        sort(new WeightMatrixScanNameComparator());
        sortByVersion = false;
    }
    
}

class WeightMatrixScanVersionComparator implements Comparator<WeightMatrixScan> {
    public int compare(WeightMatrixScan a, WeightMatrixScan b) {
        return a.matrix.version.compareTo(b.matrix.version);
    }
}
class WeightMatrixScanNameComparator implements Comparator<WeightMatrixScan> {
    public int compare(WeightMatrixScan a, WeightMatrixScan b) {
        return a.matrix.name.compareTo(b.matrix.name);
    }
}



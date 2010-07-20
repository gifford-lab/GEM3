package edu.mit.csail.cgs.warpdrive.components;

import javax.swing.event.*;
import javax.swing.table.*;
import javax.swing.*;
import java.util.*;
import edu.mit.csail.cgs.datasets.locators.*;
import edu.mit.csail.cgs.datasets.chipseq.ChipSeqAnalysis;
import edu.mit.csail.cgs.viz.components.ObjectTableModel;

public class ChipSeqAnalysisTableModel extends ObjectTableModel<ChipSeqAnalysis> {
    
    private int compareNameVersions(String n1, String v1, String n2, String v2) { 
        int c = v1.compareTo(v2);
        if(c != 0) { return c; }
        return n1.compareTo(n2);
    }
    
    private int findNewIndex(ChipSeqAnalysis bs) {
        String n = bs.getName(), v = bs.getVersion();
        for(int i = 0; i < getSize(); i++) { 
            ChipSeqAnalysis os = getObject(i);
            String on = os.getName(), ov = os.getVersion();
            int c = compareNameVersions(n, v, on, ov);
            if(c <= 0) { return i; }
        }
        return getSize();
    }

    public int getColumnCount() {
        return 3;
    }
    
    public Class getColumnClass(int i) {
        if(i==0) { return String.class; }
        if(i==1) { return String.class; }
        if(i==2) { return String.class; }
        return null;
    }


    public String getColumnName(int i) {
        if(i==0) { return "Analysis Name"; }
        if(i==1) { return "Analysis Version"; }
        if(i==2) { return "Program"; }
        return null;
    }

    public Object getValueAt(int rowIndex, int c) {
        if(c==0) { return getObject(rowIndex).getName(); }
        if(c==1) { return getObject(rowIndex).getVersion(); }
        if(c==2) { return getObject(rowIndex).getProgramName(); }
        return null;
    }

}



/*
 * Created on Nov 20, 2006
 */
package edu.mit.csail.cgs.viz.components;

import javax.swing.event.*;
import javax.swing.table.*;
import javax.swing.*;
import java.util.*;
import edu.mit.csail.cgs.datasets.locators.*;

public class ExptTableModel extends ObjectTableModel<ExptLocator> {
    
    private int compareNameVersions(String n1, String v1, String n2, String v2) { 
        int c = v1.compareTo(v2);
        if(c != 0) { return c; }
        return n1.compareTo(n2);
    }
    
    private int findNewIndex(ExptLocator bs) {
        String n = bs.getNameVersion().name, v = bs.getNameVersion().version;
        for(int i = 0; i < getSize(); i++) { 
            ExptLocator os = getObject(i);
            String on = os.getNameVersion().name, ov = os.getNameVersion().version;
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
        if(i==0) { return "Version"; }
        if(i==1) { return "Type"; }
        if(i==2) { return "Replicate"; }
        return null;
    }

    public Object getValueAt(int rowIndex, int c) {
        if(c==0) { return getObject(rowIndex).getNameVersion().name; }
        if(c==1) { return getObject(rowIndex).getNameVersion().version; }
        if(c==2) { 
            ExptLocator l = getObject(rowIndex);
            if (l instanceof ChipChipLocator) {
                return ((ChipChipLocator)l).getReplicate();
            } else {
                return null;
            }
        }
        return null;
    }

}



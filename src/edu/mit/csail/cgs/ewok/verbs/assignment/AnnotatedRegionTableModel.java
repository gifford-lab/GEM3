/*
 * Created on Dec 5, 2006
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.ewok.verbs.assignment;

import java.awt.BorderLayout;
import java.awt.Container;
import java.sql.SQLException;
import java.util.*;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.binding.BindingScan;
import edu.mit.csail.cgs.datasets.binding.BindingScanLoader;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.ewok.verbs.binding.BindingScanGenerator;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

public class AnnotatedRegionTableModel implements TableModel {

    // columns:
    // 0: chrom
    // 1: start
    // 2: stop
    // 3: bit-vector

    private Annotations<Region, BindingEvent> annotations;
    private LinkedList<TableModelListener> listeners;
    
    public static final int CHROMCOL = 0;
    public static final int STARTCOL = 1;
    public static final int ENDCOL = 2;
    public static final int BITSCOL = 3;
    
    public AnnotatedRegionTableModel(Annotations<Region, BindingEvent> annots) { 
        annotations = annots;
        listeners = new LinkedList<TableModelListener>();
    }

    public void addTableModelListener(TableModelListener tml) {
        listeners.add(tml);
    }

    public void removeTableModelListener(TableModelListener tml) {
        listeners.remove(tml);
    }

    public Class< ? > getColumnClass(int c) {
        switch(c) { 
        case CHROMCOL: return String.class;
        case STARTCOL: return Integer.class;
        case ENDCOL: return Integer.class;
        case BITSCOL: return String.class;
        }
        return null;
    }

    public int getColumnCount() {
        return 4;
    }

    public String getColumnName(int c) {
        switch(c) { 
        case CHROMCOL: return "Chrom:";
        case STARTCOL: return "Start:";
        case ENDCOL: return "End:";
        case BITSCOL: return "Binding:";
        }
        return null;
    }

    public int getRowCount() {
        return annotations.getNumItems();
    }

    public Object getValueAt(int r, int c) {
        switch(c) { 
        case CHROMCOL: return annotations.getItem(r).getChrom();
        case STARTCOL: return annotations.getItem(r).getStart();
        case ENDCOL: return annotations.getItem(r).getEnd();
        case BITSCOL: return getBits(r);
        }
        return null;
    }
    
    private String getBits(int r) {
        return annotations.getAnnotationBitVector(annotations.getItem(r));
    }

    public boolean isCellEditable(int arg0, int arg1) {
        return false;
    }

    public void setValueAt(Object arg0, int arg1, int arg2) {
        throw new UnsupportedOperationException("Can't edit cell.");
    }
}

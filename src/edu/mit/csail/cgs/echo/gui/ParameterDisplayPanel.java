/*
 * Created on Feb 21, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.echo.gui;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;

import edu.mit.csail.cgs.ewok.types.EchoType;

import java.util.*;

public class ParameterDisplayPanel extends JPanel {
    
    public ParameterDisplayPanel(String[] names, EchoType[] classes) {
        super();
        setLayout(new BorderLayout());
        
        TableModel model = new ParameterTableModel(names, classes);
        JTable table = new JTable(model);
        add(new JScrollPane(table), BorderLayout.CENTER);
    }

    public static class ParameterTableModel implements TableModel {
        
        private String[] pnames;
        private EchoType[] pclasses;
        
        private LinkedList<TableModelListener> listeners;
        
        public ParameterTableModel(String[] n, EchoType[] c) { 
            pnames = n;
            pclasses = c;
            listeners = new LinkedList<TableModelListener>();
        }

        public void addTableModelListener(TableModelListener arg0) {
            listeners.addLast(arg0);
        }

        public Class< ? > getColumnClass(int columnIndex) {
            return String.class;
        }

        public int getColumnCount() {
            return 2;
        }

        public String getColumnName(int columnIndex) {
            switch(columnIndex) { 
            case 0: return "Name";
            case 1: return "Type";
            }
            return null;
        }

        public int getRowCount() {
            return pnames != null ? pnames.length : 0;
        }

        public Object getValueAt(int rowIndex, int columnIndex) {
            switch(columnIndex) { 
            case 0: return pnames[rowIndex];
            case 1: 
            	return pclasses[rowIndex] != null ? pclasses[rowIndex].getName() : "null";
            }
            return null;
        }

        public boolean isCellEditable(int rowIndex, int columnIndex) {
            return false;
        }

        public void removeTableModelListener(TableModelListener l) {
            listeners.remove(l);
        }

        public void setValueAt(Object aValue, int rowIndex, int columnIndex) {
            throw new UnsupportedOperationException();
        } 
    }
}

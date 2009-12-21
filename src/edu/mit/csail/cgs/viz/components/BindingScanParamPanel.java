/*
 * Created on Feb 15, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.viz.components;

import java.util.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.TableModelListener;
import javax.swing.table.*;

public class BindingScanParamPanel extends JPanel {
    
    public static class Frame extends JFrame { 
        
        public Frame(Map<String,String> p) { 
            super("BindingScan Parameters");
            BindingScanParamPanel pp = new BindingScanParamPanel(p);
            Container c = (Container)getContentPane();
            c.setLayout(new BorderLayout());
            c.add(pp, BorderLayout.CENTER);
            
            setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
            setVisible(true);
            pack();
        }
    }
    
    private Map<String,String> params;
    private JTable table;
    
    public BindingScanParamPanel(Map<String,String> p) { 
        super();
        params = new HashMap<String,String>(p);

        setLayout(new BorderLayout());
        TableModel model = new ParamTableModel(params);
        table = new JTable(model);
        
        add(new JScrollPane(table), BorderLayout.CENTER);
    }
    
    private static class ParamTableModel implements TableModel {
        
        private LinkedList<TableModelListener> listeners;
        private Map<String,String> params;
        private String[] keys;
        
        public ParamTableModel(Map<String,String> p) { 
            listeners = new LinkedList<TableModelListener>();
            params = p;
            Vector<String> keyv = new Vector<String>(params.keySet());
            keys = keyv.toArray(new String[keyv.size()]);
            Arrays.sort(keys);
        }

        public void addTableModelListener(TableModelListener tml) {
            if(!listeners.contains(tml)) { listeners.addLast(tml); } 
        }

        public Class<?> getColumnClass(int ind) {
            return String.class;
        }

        public int getColumnCount() {
            return 2;
        }

        public String getColumnName(int columnIndex) {
            switch(columnIndex) { 
            case 0: return "Parameter";
            case 1: return "Value";
            default: return null;
            }
        }

        public int getRowCount() {
            return keys.length;
        }

        public Object getValueAt(int rowIndex, int columnIndex) {
            switch(columnIndex) { 
            case 0: return keys[rowIndex];
            case 1: return params.get(keys[rowIndex]);
            default: return null;
            }
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

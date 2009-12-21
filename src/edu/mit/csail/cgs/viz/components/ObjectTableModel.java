package edu.mit.csail.cgs.viz.components;

import javax.swing.event.*;
import javax.swing.table.*;
import javax.swing.*;
import java.util.*;

public class ObjectTableModel<X> implements TableModel {

    private Vector<X> objects;
    private LinkedList<TableModelListener> listeners;

    public ObjectTableModel() {
        objects = new Vector<X>();
        listeners = new LinkedList<TableModelListener>();
    }
    public void clear() {
        if (objects.size() > 0) {
            TableModelEvent evt = 
                new TableModelEvent(this, 0, objects.size()-1, TableModelEvent.ALL_COLUMNS, TableModelEvent.DELETE);
            objects.clear();
            notifyListeners(evt);
        }
    }
    public boolean contains(X o) {
        return objects.contains(o);
    }
    public int indexOf(X o) { return objects.indexOf(o); }    
    public int getSize() { return objects.size(); }
    
    public void addObject(X o) { 
        if (!contains(o)) {
            objects.add(o);
            sort();
            TableModelEvent evt = new TableModelEvent(this, objects.size()-1, objects.size()-1, TableModelEvent.ALL_COLUMNS, TableModelEvent.INSERT);
            notifyListeners(evt);
        }
    }
    
    public X getObject(int i) { 
        return objects.get(i);
    }
    public Collection<X> getObjects() {
        return (Collection<X>)objects.clone();
    }
    
    public void deleteObject(int i) { 
        objects.remove(i);
        TableModelEvent evt = new TableModelEvent(this, i, i, TableModelEvent.ALL_COLUMNS, TableModelEvent.DELETE);
        notifyListeners(evt);
    }
    
    private void notifyListeners(TableModelEvent evt) { 
        for(TableModelListener tml : listeners) { 
            tml.tableChanged(evt);
        }
    }

    public void addTableModelListener(TableModelListener arg0) {
        listeners.addLast(arg0);
    }

    public int getColumnCount() {
        return 1;
    }

    public Class getColumnClass(int i) {
        if(i==0) { return String.class; }
        return null;
    }


    public String getColumnName(int i) {
        if(i==0) { return "Name";}
        return null;
    }

    public int getRowCount() {
        return objects.size();
    }
    
    public Object getValueAt(int rowIndex, int c) {
        if(c==0) { return objects.get(rowIndex).toString(); }
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

    public void sortByColumn(int index) {
        Collections.sort(objects,new OTMColumnComparator(this,objects,index));
    }
    
    public void sort() {
    }

    public void sort(Comparator<X> comparator) {
        Collections.sort(objects,comparator);
    }
}

class OTMColumnComparator implements Comparator {
    private int col;
    private TableModel model;
    private Vector vector;
    public OTMColumnComparator(TableModel tm, Vector v, int c) {
        col = c;
        model = tm;
        vector = v;
    }
    public int compare(Object first, Object second) {
        int listposfirst = vector.indexOf(first);
        int listpossecond = vector.indexOf(second);
        if (listposfirst == -1 || listpossecond == -1) {
            return 0;
        }
        first = model.getValueAt(listposfirst,col);
        second = model.getValueAt(listpossecond,col);
        if (first instanceof Comparable) {
            return ((Comparable)first).compareTo(second);
        } else {
            return 0;
        }


    }
    
    
}



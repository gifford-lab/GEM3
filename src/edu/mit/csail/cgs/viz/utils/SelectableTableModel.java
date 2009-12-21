/**
 * 
 */
package edu.mit.csail.cgs.viz.utils;

import java.util.Collection;
import java.util.LinkedList;
import java.util.Vector;

import javax.swing.table.AbstractTableModel;

/**
 * @author Timothy Danford
 *
 */
public class SelectableTableModel<X> extends AbstractTableModel { 

	private boolean isMultiSelectable;
	private Vector<Boolean> selections;
	private Vector<X> values;
	private Class[] columnTypes;
	private String[] columnNames;
	
	public SelectableTableModel(String colName, Class colClass) { 
		selections = new Vector<Boolean>();
		values = new Vector<X>();
		columnTypes = new Class[2]; 
		columnTypes[0] = Boolean.class;
		columnTypes[1] = colClass;
		columnNames = new String[2];
		columnNames[0] = "?";
		columnNames[1] = colName;
		isMultiSelectable = true;
	}
	
	public int getRowCount() { return values.size(); }
	public boolean isMultiSelectable() { return isMultiSelectable; }
	public void setMultiSelectable(boolean v) { isMultiSelectable = v; }
	
	public int getColumnCount() { return 2; }
	public Class getColumnClass(int c) { return columnTypes[c]; }
	public String getColumnName(int c) { return columnNames[c]; }
	
	public boolean isCellEditable(int r, int c) { 
		return r >= 0 && r < values.size() && c == 0;
	}
	
	public void addValue(X v) { 
		values.add(v); selections.add(false);
		fireTableRowsInserted(values.size()-1, values.size()-1);
	}
	
	public X getValue(int r) { return values.get(r); }
    
    public void setValue(int r, X v) { 
        values.set(r, v);
        fireTableRowsUpdated(r, r);
    }
    
	public boolean isSelected(int r) { return selections.get(r); }

	public void setAllSelections(boolean v) { 
		for(int i = 0; i < selections.size(); i++) { 
			selections.set(i, v);
		}
		fireTableRowsUpdated(0, selections.size()-1);
	}
    
    public X getFirstSelected() {
        for(int i = 0; i < selections.size(); i++) { 
            if(selections.get(i)) { 
                return values.get(i);
            }
        }
        return null;
    }

	public Collection<X> getSelected() {
		LinkedList<X> lst = new LinkedList<X>();
		for(int i = 0; i< selections.size(); i++) { 
			if(selections.get(i)) { 
				lst.addLast(values.get(i));
			}
		}
		return lst;
	}
	
	public void clear() { 
		int end = values.size()-1;
		selections.clear();
		values.clear();
		if(end >= 0) { fireTableRowsDeleted(0, end); } 
	}
	
	public void setValueAt(Object v, int r, int c) { 
		if(c == 0) { 
			if(!isMultiSelectable && !selections.get(r)) { 
				for(int i = 0; i < selections.size(); i++) { 
					if(selections.get(i)) { 
						selections.set(i, new Boolean(false));
						fireTableRowsUpdated(i, i);
					}
				}
			}
			
			selections.set(r, (Boolean)v);
			fireTableRowsUpdated(r, r);
		}
	}
	
	public Object getValueAt(int r, int c) { 
		if(c == 0) { return selections.get(r); }
		if(c == 1) { return values.get(r); }
		return null;
	}
    
    public int[] getSelectionIndices() { 
        int count = 0; 
        for(int i = 0; i < selections.size(); i++) { 
            if(selections.get(i)) { count += 1; }
        }
        int[] array = new int[count];
        int j = 0;
        for(int i = 0; i < selections.size(); i++) { 
            if(selections.get(i)) { array[j++] = i; }
        }
        return array;
    }
}

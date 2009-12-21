/**
 * @author tdanford
 */
package edu.mit.csail.cgs.viz.utils;

import javax.swing.*;
import javax.swing.table.*;

/**
 * @author Timothy Danford
 *
 */
public class SelectableTable<X> extends JTable {

	private SelectableTableModel<X> myModel;
	
	public SelectableTable(String colName, Class colClass) { 
		super();
		myModel = new SelectableTableModel<X>(colName, colClass);
		setModel(myModel);
		if(getModel() == null) { throw new IllegalArgumentException(); }
		setDefaultEditor(Boolean.class, new DefaultCellEditor(new JCheckBox()));
		TableColumnModel cModel = getColumnModel();
		TableColumn col = cModel.getColumn(0);
		col.setMaxWidth(20);
	}

	public SelectableTableModel<X> getModel() { return myModel; }
	public void setMultiSelectable(boolean v) { myModel.setMultiSelectable(v); }
	public boolean isMultiSelectable() { return myModel.isMultiSelectable(); }
}

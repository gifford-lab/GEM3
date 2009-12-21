package edu.mit.csail.cgs.datasets.expression.gui;

import java.util.*;

import javax.swing.*;
import javax.swing.event.TableModelListener;
import javax.swing.table.TableModel;

import edu.mit.csail.cgs.datasets.expression.Experiment;
import edu.mit.csail.cgs.utils.Pair;

import java.awt.*;
import java.awt.event.*;

public class ExperimentPanel {
	
	private Experiment expt;
	
	private JTextField nameField, typeField, cellsField, conditionField;
	private JTable paramsTable;
	private ParamsTableModel paramsModel;
	
	public ExperimentPanel(Experiment expt) { 
		super();
		this.expt = expt;
		
		nameField = new JTextField();
		typeField = new JTextField();
		cellsField = new JTextField();
		conditionField = new JTextField();
		
		nameField.setEditable(false);
		typeField.setEditable(false);
		cellsField.setEditable(false);
		conditionField.setEditable(false);
		
		nameField.setText(expt.getName());
		typeField.setText(String.valueOf(expt.getValueType()));
		cellsField.setText(expt.getCells().getName());
		conditionField.setText(expt.getCondition().getName());
		
		JLabel nameLabel = new JLabel("Name:"),
			typeLabel = new JLabel("Value Type:"), 
			cellsLabel = new JLabel("Cells:"), 
			conditionLabel = new JLabel("Condition:");
		
		
		paramsModel = new ParamsTableModel(expt);
		paramsTable = new JTable(paramsModel);
	}
	
	private static class ParamsTableModel implements TableModel {
		
		private LinkedList<TableModelListener> listeners;
		private Vector<Pair<String,String>> values;
		
		public ParamsTableModel(Experiment expt) { 
			values = new Vector<Pair<String,String>>();
			listeners = new LinkedList<TableModelListener>();
			for(String p : expt.getParams()) { 
				values.add(new Pair<String,String>(p, expt.getParamValue(p)));
			}
		}

		public void addTableModelListener(TableModelListener tml) {
			listeners.addLast(tml);
		}

		public Class getColumnClass(int index) {
			switch(index) { 
			case 0: 
			case 1:
				return String.class;
			default:
				return null;
			}
		}

		public int getColumnCount() { return 2; }

		public String getColumnName(int index) {
			switch(index) { 
			case 0: return "Key:";
			case 1: return "Value:";
			default: return null;
			}
		}

		public int getRowCount() { return values.size(); }

		public Object getValueAt(int r, int c) {
			switch(c) { 
			case 0: return values.get(r).getFirst();
			case 1: return values.get(r).getLast();
			default: return null;
			}
		}

		public boolean isCellEditable(int arg0, int arg1) { return false; }

		public void removeTableModelListener(TableModelListener tml) {
			listeners.remove(tml);
		}

		public void setValueAt(Object arg0, int arg1, int arg2) {
			throw new UnsupportedOperationException();
		} 
		
	}
}

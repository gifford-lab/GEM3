/*
 * Created on Nov 9, 2006
 */
package edu.mit.csail.cgs.warpdrive.components;

import java.sql.SQLException;
import java.util.*;
import java.io.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;
import javax.swing.border.*;
import java.awt.*;
import java.awt.event.*;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.function.*;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.tools.utils.Args;

/**
 * @author tdanford
 */
public class GOAnnotationPanel extends JPanel {
    
    public static void main(String[] args) { 
        try {
            DatabaseFunctionLoader loader = new DatabaseFunctionLoader();
            GOAnnotationPanel panel = new GOAnnotationPanel(loader, Args.parseGenome(args).cdr());
            Frame f = new Frame(panel);
            f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        } catch (SQLException e) {
            e.printStackTrace();
        } catch (UnknownRoleException e) {
            e.printStackTrace();
        } catch (NotFoundException e) {
            e.printStackTrace();
        }

    }
    
    private FunctionLoader loader;
    private Map<String,FunctionVersion> versionMap;
    
    private DefaultComboBoxModel versionModel;
    private JComboBox versionBox;
    
    private AnnotationTableModel annotationModel;
    private JTable annotationTable;
    
    private JTextField goField;
    private JButton annotateButton;
    
    public GOAnnotationPanel(FunctionLoader l, Genome g) throws SQLException {
        super();
        
        loader = l;
        FunctionVersion currentVersion = l.getVersion(g.getVersion());
        
        versionModel = new DefaultComboBoxModel();
        
        for(FunctionVersion v : loader.getAllVersions()) { 
            versionMap.put(v.getName(), v);
            currentVersion = v;
            versionModel.addElement(v.getName());
        }
        
        if(currentVersion != null) { versionModel.setSelectedItem(currentVersion.getName()); }
        
        annotationModel = new AnnotationTableModel();

        annotationTable = new JTable(annotationModel);
        versionBox = new JComboBox(versionModel);
        
        JPanel tablePanel = new JPanel(); tablePanel.setLayout(new BorderLayout());
        tablePanel.add(new JScrollPane(annotationTable), BorderLayout.CENTER);
        tablePanel.setBorder(new TitledBorder("Annotations:"));
        
        JPanel entryPanel = new JPanel(); entryPanel.setLayout(new GridLayout(2, 2));
        entryPanel.add(new JLabel("Version:"));
        entryPanel.add(versionBox);
        entryPanel.add(new JLabel("ID:"));
        entryPanel.add(goField = new JTextField());
        
        annotateButton = new JButton("Annotate");
        annotateButton.addActionListener(new ActionListener() { 
            public void actionPerformed(ActionEvent e) { 
                annotate();
            }
        });
        
        setLayout(new BorderLayout());
        add(entryPanel, BorderLayout.NORTH);
        add(tablePanel, BorderLayout.CENTER);
        add(annotateButton, BorderLayout.SOUTH);
    }
    
    public void setVersion(String v) { 
        if(versionMap.containsKey(v)) { 
            versionModel.setSelectedItem(v);
        }
    }
    
    public void setID(String id) { 
        goField.setText(id);
    }
    
    public void annotate() { 
        annotationModel.clear();
        String versionName = (String)versionModel.getSelectedItem();
        if(versionName != null) { 
            FunctionVersion version = versionMap.get(versionName);
            String id = goField.getText();
            if(id != null) { 
                try {
                    Collection<Assignment> assigns = loader.getAssignments(id, version);
                    annotationModel.addValues(assigns);
                } catch (SQLException e) {
                    e.printStackTrace();
                }
            }
        }
    }
    
    public static class Frame extends JFrame {
        
        private GOAnnotationPanel panel;
        
        public Frame(GOAnnotationPanel p) { 
            super("GO Annotation");
            panel = p;
            
            Container c = (Container)getContentPane();
            c.setLayout(new BorderLayout());
            c.add(p, BorderLayout.CENTER);
            
            setVisible(true);
            pack();
            setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
        }
        
        public GOAnnotationPanel getPanel() { return panel; }
    }
    
    private static class AnnotationTableModel implements TableModel {
        
        private Vector<Assignment> assigns;
        private LinkedList<TableModelListener> listeners;
        
        public AnnotationTableModel() { 
            assigns = new Vector<Assignment>();
            listeners = new LinkedList<TableModelListener>();
        }
        
        public void addValues(Collection<Assignment> as) { 
            assigns.addAll(as);
            TableModelEvent evt = new TableModelEvent(this);
            for(TableModelListener list : listeners) { 
                list.tableChanged(evt);
            }            
        }
        
        public void clear() { 
            assigns.clear();
            TableModelEvent evt = new TableModelEvent(this);
            for(TableModelListener list : listeners) { 
                list.tableChanged(evt);
            }
        }

        public int getRowCount() {
            return assigns.size();
        }

        public int getColumnCount() {
            return 2;
        }

        public String getColumnName(int i) {
            switch(i) {
            case 0:
                return "GO ID:";
            case 1:
                return "Description:";
            default:
                return "";
            }
        }

        public Class getColumnClass(int index) {
            return String.class;
        }

        public boolean isCellEditable(int arg0, int arg1) {
            return false;
        }

        public Object getValueAt(int row, int col) {
            switch(col) {
            case 0: return assigns.get(row).getCategory().getName();
            case 1: return assigns.get(row).getCategory().getDescription();
            default:
                return null;
            }
        }

        public void setValueAt(Object arg0, int arg1, int arg2) {
            throw new UnsupportedOperationException();
        }

        public void addTableModelListener(TableModelListener list) {
            listeners.addLast(list);
        }

        public void removeTableModelListener(TableModelListener list) {
            listeners.remove(list);
        } 
        
    }
    
    
}

/*
 * Created on Dec 5, 2006
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.ewok.verbs.assignment;

import java.sql.SQLException;
import java.util.*;
import javax.swing.*;
import java.awt.*;
import javax.swing.event.*;
import javax.swing.table.*;

import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.ewok.verbs.binding.BindingScanGenerator;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

import edu.mit.csail.cgs.datasets.binding.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;

public class AnnotatedBindingTableModel implements TableModel {
    
    public static void main(String[] args) {
        String gname = "hg17";
        try {
            Genome g = Organism.findGenome(gname);
            BindingScanLoader loader = new BindingScanLoader();
            
            Vector<BindingScan> scans = new Vector<BindingScan>(loader.loadScans(g));
            Random rand = new Random();
            BindingScan scan = scans.get(rand.nextInt(scans.size()));
            System.out.println("Selected scan: " + scan.toString());
            
            Expander<Region,BindingEvent> caller = new BindingScanGenerator(loader, scan);
            
            BindingEventAnnotations annots = new BindingEventAnnotations(g, caller);
            AnnotatedBindingTableModel model = new AnnotatedBindingTableModel(annots);
            
            JTable table = new JTable(model);
            JPanel panel = new JPanel(); panel.setLayout(new BorderLayout());
            panel.add(new JScrollPane(table), BorderLayout.CENTER);
            
            JFrame f = new JFrame("Annotated Binding");
            Container c = (Container)f.getContentPane(); c.setLayout(new BorderLayout());
            c.add(panel, BorderLayout.CENTER);
            
            f.setVisible(true);
            f.pack();
            f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        } catch (NotFoundException e) {
            e.printStackTrace();
        } catch (SQLException e) {
            e.printStackTrace();
        } catch (UnknownRoleException e) {
            e.printStackTrace();
        }
    }

    // columns:
    // 0: chrom
    // 1: start
    // 2: stop
    // 3: names

    private Annotations<BindingEvent, Region> annotations;
    private LinkedList<TableModelListener> listeners;
    
    public static final int CHROMCOL = 0;
    public static final int STARTCOL = 1;
    public static final int ENDCOL = 2;
    public static final int NAMESCOL = 3;
    
    public AnnotatedBindingTableModel(Annotations<BindingEvent, Region> annots) { 
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
        case NAMESCOL: return String.class;
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
        case NAMESCOL: return "Assigned Names";
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
        case NAMESCOL: return getNames(r);
        }
        return null;
    }
    
    private String getNames(int r) { 
        BindingEvent evt = annotations.getItem(r);
        StringBuilder sb = new StringBuilder();
        boolean first = true;
        for(Region annot : annotations.getAnnotations(evt)) {
            String val = annot.toString();
            if(annot instanceof Gene) { 
                Gene g = (Gene)annot;
                val = g.getName();
            }
            sb.append((first ? "" : ", ") + val);
            first = false;
        }
        return sb.toString();
    }

    public boolean isCellEditable(int arg0, int arg1) {
        return false;
    }

    public void setValueAt(Object arg0, int arg1, int arg2) {
        throw new UnsupportedOperationException("Can't edit cell.");
    }
}

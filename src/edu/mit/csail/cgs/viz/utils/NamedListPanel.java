/*
 * Created on Aug 20, 2005
 */
package edu.mit.csail.cgs.viz.utils;

import java.util.*;

import javax.swing.*;
import java.awt.*;

import edu.mit.csail.cgs.utils.NamedValue;

/**
 * @author tdanford
 */
public class NamedListPanel<X> extends JPanel implements ListPanel<X>, ListPanelListener {
    
    private SimpleListPanel<NamedValue<X>> simpleList;
    protected ListPanelEventSource.Default source;
    
    public NamedListPanel() { 
        super();
        setLayout(new BorderLayout());
        simpleList = new SimpleListPanel<NamedValue<X>>();
        add(simpleList, BorderLayout.CENTER);
        source = new ListPanelEventSource.Default(this);
        simpleList.addListPanelListener(this);
    }
    
    public void addNamedValues(Map<String,X> valMap) { 
        for(String key : valMap.keySet()) { 
            addNamedValue(key, valMap.get(key));
        }
    }
    
    public void panelEvent(ListPanelEvent evt) { 
        int type = evt.getType();
        NamedValue<X> eData = (NamedValue<X>)evt.getData();
        source.firePanelEvent(type, eData.getData());
    }
    
    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.hyperdrive.ui.ListPanel#getNumValues()
     */
    public int getNumValues() {
        return simpleList.getNumValues();
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.hyperdrive.ui.ListPanel#getValue(int)
     */
    public X getValue(int index) {
        return simpleList.getValue(index).getData();
    }
    
    public String getName(int index) { 
        return simpleList.getValue(index).getName();
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.hyperdrive.ui.ListPanel#addValue(null)
     */
    public void addValue(X v) {
        NamedValue<X> nv = new NamedValue<X>(v);
        simpleList.addValue(nv);
    }
    
    public void addNamedValue(String n, X v) { 
        NamedValue<X> nv = new NamedValue<X>(n, v);
        simpleList.addValue(nv);        
    }
    
    public void addNamedValue(NamedValue<X> nv) { 
        simpleList.addValue(nv);
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.hyperdrive.ui.ListPanel#clear()
     */
    public void clear() {
        simpleList.clear();
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.hyperdrive.ui.ListPanel#getFirstSelectedValue()
     */
    public X getFirstSelectedValue() {
        return simpleList.getFirstSelectedValue().getData();
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.hyperdrive.ui.ListPanel#getSelectedValues()
     */
    public Collection<X> getSelectedValues() {
        LinkedList<X> dataList = new LinkedList<X>();
        Collection<NamedValue<X>> lst = simpleList.getSelectedValues();
        for(NamedValue<X> nv : lst) { 
            dataList.addLast(nv.getData());
        }
        return dataList;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.hyperdrive.ui.ListPanelEventSource#addListPanelListener(edu.mit.csail.cgs.hyperdrive.ui.ListPanelListener)
     */
    public void addListPanelListener(ListPanelListener lpl) {
        source.addListPanelListener(lpl);
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.hyperdrive.ui.ListPanelEventSource#removeListPanelListener(edu.mit.csail.cgs.hyperdrive.ui.ListPanelListener)
     */
    public void removeListPanelListener(ListPanelListener lpl) {
        source.removeListPanelListener(lpl);
    }
    
    public void removeValue(int index) { 
        simpleList.removeValue(index);
    }
    
    public int[] getSelectedIndices() { 
        return simpleList.getSelectedIndices();
    }
    
    public void removeSelectedValues() { 
        simpleList.removeSelectedValues();
    }
    
    public Collection<X> getAllValues() { 
        LinkedList<X> dataList = new LinkedList<X>();
        Collection<NamedValue<X>> lst = simpleList.getAllValues();
        for(NamedValue<X> nv : lst) { 
            dataList.addLast(nv.getData());
        }
        return dataList;
    }
}

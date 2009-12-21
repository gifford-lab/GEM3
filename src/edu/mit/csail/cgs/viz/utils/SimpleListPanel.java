/*
 * Created on Aug 19, 2005
 */
package edu.mit.csail.cgs.viz.utils;

import java.util.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseListener;

/**
 * @author tdanford
 */
public class SimpleListPanel<X> extends JPanel implements ListPanel<X> {
    
    private MutableListModel model;
    private JList list;
    protected ListPanelEventSource.Default source;
    
    public SimpleListPanel() { 
        super();
        source = new ListPanelEventSource.Default(this);
        model = new MutableListModel();
        list = new JList(model);
        setLayout(new BorderLayout());
        add(new JScrollPane(list), BorderLayout.CENTER);
    }
    
    public void addListMouseListener(MouseListener ml) { 
    	list.addMouseListener(ml);
    }
    
    public void removeListMouseListener(MouseListener ml) { 
    	list.removeMouseListener(ml);
    }
    
    public void addListPanelListener(ListPanelListener lpl) {  
        source.addListPanelListener(lpl);
    }
    public void removeListPanelListener(ListPanelListener lpl) { 
        source.removeListPanelListener(lpl);
    }
    
    public int getNumValues() { return model.getSize(); }
    public X getValue(int index) { return (X)model.elementAt(index); }
    public void addValue(X v) { model.addElement(v); }
	public void addAll(Collection<X> vals) { 
		for(X v : vals) { model.addElement(v); }
	}

    public void setValue(int index, X v) { model.setElementAt(v, index); }
    public void clear() { model.clear(); }
    public void removeValue(int index) { model.remove(index); }
    public int[] getSelectedIndices() { return list.getSelectedIndices(); }
    
    public void fireValueChanged(X v) { 
        model.fireContentsChanged(v);
    }
    
    public void fireIndexChanged(int index) { 
        model.fireContentsChanged(index);
    }

	public boolean containsValue(X v) { return model.contains(v); }
    
    public Collection<X> getAllValues() { 
        LinkedList<X> lst = new LinkedList<X>();
        for(int i = 0; i < model.getSize(); i++) { 
            lst.addLast((X)model.get(i));
        }
        return lst;
    }
    
    public void removeSelectedValues() { 
        int[] inds = getSelectedIndices();
        for(int i = inds.length-1; i >= 0; i--) { 
            removeValue(inds[i]);
        }
    }
    
    public X getFirstSelectedValue() {
        int[] inds = list.getSelectedIndices();
        if(inds.length > 0) { 
            return getValue(inds[0]);
        } else { 
            return null;
        }
    }
    
    public Collection<X> getSelectedValues() { 
        LinkedList<X> lst = new LinkedList<X>();
        int[] inds = list.getSelectedIndices();
        for(int i = 0; i < inds.length; i++) { 
            lst.addLast(getValue(inds[i]));
        }
        return lst;
    }
}

class MutableListModel extends DefaultListModel { 
    public MutableListModel() {
        super();
    }
    
    public void fireContentsChanged(Object v) { 
        int index = indexOf(v);
        fireContentsChanged(this, index, index);
    }
    
    public void fireContentsChanged(int index) { 
        fireContentsChanged(this, index, index);        
    }
}

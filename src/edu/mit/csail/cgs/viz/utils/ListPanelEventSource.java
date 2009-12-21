/*
 * Created on Aug 20, 2005
 */
package edu.mit.csail.cgs.viz.utils;

import java.util.LinkedList;

/**
 * @author tdanford
 */
public interface ListPanelEventSource {

    public void addListPanelListener(ListPanelListener lpl);
    public void removeListPanelListener(ListPanelListener lpl);
    
    public static class Default implements ListPanelEventSource {
        
        private Object source;
        private LinkedList<ListPanelListener> listeners;
        
        public Default() {
            source = this;
            listeners = new LinkedList<ListPanelListener>();
        }

        public Default(Object src) {
            source = src;
            listeners = new LinkedList<ListPanelListener>();
        }

        public void firePanelEvent(int type, Object data) { 
            ListPanelEvent evt = new ListPanelEvent(source, type, data);
            for(ListPanelListener lpl : listeners) { 
                lpl.panelEvent(evt);
            }
        }
        
        public void addListPanelListener(ListPanelListener lpl) { 
            listeners.addLast(lpl);
        }
        
        public void removeListPanelListener(ListPanelListener lpl) { 
            listeners.remove(lpl);
        }
    }
}

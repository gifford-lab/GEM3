/*
 * Created on Aug 22, 2005
 */
package edu.mit.csail.cgs.utils;

import java.util.*;

/**
 * @author tdanford
 */
public interface EventSource<E> {
    
    public void addEventListener(Listener<E> el);
    public void removeEventListener(Listener<E> el);
    public boolean hasListeners();
    
    public static class Default<Event> implements EventSource<Event> { 
        protected Object source;
        protected LinkedList<Listener<Event>> listeners;
        
        public Default() {
            source = this;
            listeners = new LinkedList<Listener<Event>>();
        }
        
        public Default(Object src) { 
            source = src;
            listeners = new LinkedList<Listener<Event>>();
        }
        
        public void fireEvent(Event e) { 
            for(Listener<Event> el : listeners) { 
                el.eventRegistered(e);
            }
        }               
        public void addEventListener(Listener<Event> el) { listeners.addLast(el); }
        public void removeEventListener(Listener<Event> el) { listeners.remove(el); }
        public boolean hasListeners() {return listeners.size() > 0;}
        
        public void clearListeners() { listeners.clear(); }
        public Collection<Listener<Event>> getListeners() { return listeners; }
    }
}

/*
 * Created on Aug 22, 2005
 */
package edu.mit.csail.cgs.utils;

/**
 * @author tdanford
 */
public interface Listener<Event> {
    public void eventRegistered(Event e);
}

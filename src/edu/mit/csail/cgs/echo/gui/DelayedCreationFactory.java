/*
 * Created on Apr 12, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.echo.gui;

import java.util.Collection;
import java.util.LinkedList;

import edu.mit.csail.cgs.utils.*;

public class DelayedCreationFactory<X> implements Runnable {

    private Collection<Listener<CreationEvent>> listeners;
    private Factory<X> factory;
    private X value;
    
    public DelayedCreationFactory(Collection<Listener<CreationEvent>> listeners, Factory<X> f) {
        factory = f;
        value = null;
        this.listeners = new LinkedList<Listener<CreationEvent>>(listeners);
    }
    
    public void run() { 
        value = factory.createObject();
        CreationEvent evt = new CreationEvent(this, value);
        for(Listener<CreationEvent> el : listeners) { 
            el.eventRegistered(evt);
        }
    }

}

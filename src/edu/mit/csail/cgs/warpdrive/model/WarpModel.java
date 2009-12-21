package edu.mit.csail.cgs.warpdrive.model;

import java.util.*;
import edu.mit.csail.cgs.utils.*;



public abstract class WarpModel implements Model {
    private HashSet<Listener<EventObject>> listeners;
    private boolean keepRunning;
    private ModelProperties props;
    
    public WarpModel() {
        listeners = new HashSet<Listener<EventObject>>();
        keepRunning = true;
        props = new ModelProperties();
    }

    public ModelProperties getProperties() {return props;}

    public boolean keepRunning() {return keepRunning;}
    public void stopRunning() {
        keepRunning = false;
        synchronized(this) {
            this.notifyAll();
        }
    }
    public synchronized void notifyListeners() {
        for (Listener<EventObject> l : listeners) {
            l.eventRegistered(new EventObject(this));
        }
    }
    public void addEventListener(Listener<EventObject> l) {
        listeners.add(l);
    }
    public void removeEventListener(Listener<EventObject> l) {
        listeners.remove(l);
        if (!hasListeners()) {
            stopRunning();
        }
    }
    public boolean hasListeners() {return listeners.size() > 0;}
}

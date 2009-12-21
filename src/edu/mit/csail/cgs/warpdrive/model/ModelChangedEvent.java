/*
 * Created on Mar 14, 2006
 */
package edu.mit.csail.cgs.warpdrive.model;

import java.util.*;

/**
 * @author tdanford
 */
public class ModelChangedEvent extends EventObject {
    
    private boolean state;

    /**
     * @param src The source object for the event.
     */
    public ModelChangedEvent(Object src, boolean s) {
        super(src);
        state = s;
    }
    
    public boolean isReady() { return state; }

}

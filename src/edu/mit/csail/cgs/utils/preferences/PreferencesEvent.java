/*
 * Created on Aug 24, 2005
 */
package edu.mit.csail.cgs.utils.preferences;

import java.util.EventObject;

/**
 * @author tdanford
 */
public class PreferencesEvent extends EventObject {
    
    public static final int OK = 0;
    public static final int CANCEL = 1;
    
    private int type;

    /**
     * @param arg0
     */
    public PreferencesEvent(Object arg0, int t) {
        super(arg0);
        type = t;
    }
    
    public int getType() { return type; }

}

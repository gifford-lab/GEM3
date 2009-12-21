/*
 * Created on Sep 29, 2005
 */
package edu.mit.csail.cgs.viz.utils;

import java.util.*;

/**
 * @author tdanford
 */
public class SelectionEvent extends EventObject {
    
    private int type;
    private Object value;
    
    public static final int OK = 0;
    public static final int CANCEL = 1;

    public SelectionEvent(Object src, int type, Object val) {
        super(src);
        this.type = type;
        value = val;
    }
    
    public int getType() { return type; }
    public boolean isOkEvent() { return type==OK; }
    public boolean isCancelEvent() { return type==CANCEL; }
    public Object getValue() { return value; }
}

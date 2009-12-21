/*
 * Created on Aug 17, 2005
 */
package edu.mit.csail.cgs.viz.utils;

import java.util.EventObject;

/**
 * @author tdanford
 */
public class ListPanelEvent extends EventObject {
    
    public static final int DISPLAY = 0;
    public static final int ENTREZ_SAVE = 1;
    public static final int EXPLORE_MOTIFS = 2;
    public static final int SAVE = 4;
    public static final int LOAD = 5;
    public static final int ADD = 6;
    public static final int REMOVE = 7;
    
    private int type;
    private Object data;

    /**
     * @param arg0
     */
    public ListPanelEvent(Object src, int type, Object data) {
        super(src);
        this.type = type;
        this.data = data;
    }
    
    public int getType() { return type; }
    public Object getData() { return data; }
}

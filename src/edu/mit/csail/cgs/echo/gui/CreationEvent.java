/*
 * Created on Apr 12, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.echo.gui;

import java.util.*;

public class CreationEvent<X> extends EventObject {
    
    private X value;
    
    public CreationEvent(Object src, X val) {
        super(src);
        value = val;
    }
    
    public X getValue() { return value; }
    
}

/*
 * Created on Aug 24, 2005
 */
package edu.mit.csail.cgs.utils.preferences;

import java.util.*;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

/**
 * @author tdanford
 */
public abstract class PreferencesPanel extends JPanel {
    
    protected Map<String,Object> values;

    public PreferencesPanel() { 
        super();
        values = new HashMap<String,Object>();
    }
    
    public Map<String,Object> getValues() { return values; }
    public void saveValues() {}
}

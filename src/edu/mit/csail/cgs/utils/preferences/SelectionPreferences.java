/*
 * Created on Sep 5, 2005
 */
package edu.mit.csail.cgs.utils.preferences;

import java.io.*;
import java.util.*;

/**
 * @author tdanford
 */
public class SelectionPreferences<X> implements Preferences<Preferences<X>> {
    
    private Vector<Preferences<X>> prefs;
    private int selected;

    public SelectionPreferences(Collection<Preferences<X>> prefs) {
        this.prefs = new Vector<Preferences<X>>(prefs);
        selected = 0;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.preferences.Preferences#getName()
     */
    public String getName() {
        return "Selection Preferences";
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.preferences.Preferences#createPanel()
     */
    public PreferencesPanel createPanel() {
        return new SelectionPreferencesPanel<X>(selected, prefs);
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.preferences.Preferences#saveFromPanel(edu.mit.csail.cgs.utils.preferences.PreferencesPanel)
     */
    public void saveFromPanel(PreferencesPanel pp) {
        Map<String,Object> values = pp.getValues();
        selected = ((Integer)values.get("selected")).intValue();
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.Factory#createObject()
     */
    public Preferences<X> createObject() {
        return prefs.get(selected);
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.Saveable#save(java.io.DataOutputStream)
     */
    public void save(DataOutputStream dos) throws IOException {
        dos.writeInt(selected);
        dos.writeInt(prefs.size());
        for(Preferences<X> p : prefs) { p.save(dos); }
    }
}

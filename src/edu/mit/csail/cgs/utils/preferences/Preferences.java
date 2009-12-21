/*
 * Created on Aug 24, 2005
 */
package edu.mit.csail.cgs.utils.preferences;

import edu.mit.csail.cgs.utils.Factory;
import edu.mit.csail.cgs.utils.Saveable;

import java.io.*;

/**
 * @author tdanford
 */
public interface Preferences<X> extends Factory<X>, Saveable {
    public String getName();
    public PreferencesPanel createPanel();
    public void saveFromPanel(PreferencesPanel pp);
}

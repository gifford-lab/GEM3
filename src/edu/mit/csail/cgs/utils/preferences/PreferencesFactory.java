/*
 * Created on Aug 27, 2005
 */
package edu.mit.csail.cgs.utils.preferences;

import edu.mit.csail.cgs.utils.Factory;

/**
 * @author tdanford
 */
public interface PreferencesFactory<X> extends Factory<Preferences<X>> {
    public Preferences<X> createObject(X v);
}

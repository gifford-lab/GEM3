package edu.mit.csail.cgs.viz.preferences;

public interface PreferencesListener {

	public void preferencesUpdated(PreferencesEvent evt);
	public void preferencesUpdateCanceled(PreferencesEvent evt);
}

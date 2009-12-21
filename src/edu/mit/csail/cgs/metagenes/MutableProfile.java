/*
 * Author: tdanford
 * Date: Aug 12, 2008
 */
package edu.mit.csail.cgs.metagenes;

import java.util.LinkedList;

public abstract class MutableProfile implements Profile {

	private LinkedList<ProfileListener> listeners;

	public MutableProfile() { 
		listeners = new LinkedList<ProfileListener>();
	}
	
	public void addProfileListener(ProfileListener pl) { 
		listeners.add(pl);
	}
	
	public void removeProfileListener(ProfileListener pl) {  
		listeners.remove(pl);
	}
	
	protected void dispatchChange(ProfileEvent e) { 
		for(ProfileListener pl : listeners) { 
			pl.profileChanged(e);
		}
	}
}

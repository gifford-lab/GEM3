package edu.mit.csail.cgs.metagenes;

import java.util.*;

public class ProfileEvent { 
	
	public static enum EventType { CHANGED, ADDED };

	private EventType type;
	private Profile changedProfile, addedProfile;
	
	public ProfileEvent(Profile c) { 
		changedProfile = c;
		addedProfile = null;
		type = EventType.CHANGED;
	}
	
	public ProfileEvent(Profile c, Profile a) { 
		changedProfile = c;
		addedProfile = a;
		type = EventType.ADDED;
	}
	
	public EventType getType() { return type; }
	public Profile changedProfile() { return changedProfile; }
	public Profile addedProfile() { return addedProfile; }
}

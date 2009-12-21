package edu.mit.csail.cgs.ewok.types;

import java.util.Collection;

public interface SelfDescribingOutput<X> extends SelfDescribingVerb { 
	public Collection<X> getValues();
}

package edu.mit.csail.cgs.conservation;

import java.util.*;

public interface ConservationDataset {
	public Vector<ExptDescriptor> getExpts();
	public Set<String> getIDs();
	
	public boolean isBound(String id, ExptDescriptor ed, BindingOptions bo);
	public Set<String> getBound(ExptDescriptor ed, BindingOptions bo);
}

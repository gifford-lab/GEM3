package edu.mit.csail.cgs.tools.hypotheses.utils;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.binding.BindingExtent;
import edu.mit.csail.cgs.ewok.verbs.Mapper;

public class ExtentToEventMapper implements Mapper<BindingExtent,BindingEvent> {
	
	public BindingEvent execute(BindingExtent ext) { 
		return new BindingEvent(ext);
	}
}

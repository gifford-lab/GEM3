package edu.mit.csail.cgs.tools.hypotheses;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.locators.ChipChipLocator;
import edu.mit.csail.cgs.datasets.locators.ExptLocator;
import edu.mit.csail.cgs.ewok.verbs.Expander;

public interface BindingExpanderFactory<LOC extends ExptLocator> {
	public Expander<Region,BindingEvent> getExpander(LOC loc);
}


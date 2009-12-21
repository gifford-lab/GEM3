/**
 * 
 */
package edu.mit.csail.cgs.warpdrive.model;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.verbs.Expander;

/**
 * @author tdanford
 */
public class BindingEventModel extends RegionExpanderModel<BindingEvent> {
    public BindingEventModel(Expander<Region,BindingEvent> exp) { 
        super(exp);
    }
}

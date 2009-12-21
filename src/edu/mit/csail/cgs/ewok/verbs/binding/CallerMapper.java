/*
 * Created on Nov 27, 2006
 */
package edu.mit.csail.cgs.ewok.verbs.binding;

import java.util.Map;

import edu.mit.csail.cgs.datasets.binding.BindingExtent;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.locators.ExptLocator;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.ewok.verbs.Mapper;

public interface CallerMapper extends Mapper<ExptLocator,Expander<Region,BindingExtent>> {
    
    public Map<String,String> getLastParams();
}
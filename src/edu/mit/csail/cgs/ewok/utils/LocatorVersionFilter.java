/*
 * Created on Apr 3, 2006
 */
package edu.mit.csail.cgs.ewok.utils;

import edu.mit.csail.cgs.datasets.locators.ExptLocator;
import edu.mit.csail.cgs.ewok.verbs.Filter;

/**
 * @author tdanford
 */
public class LocatorVersionFilter implements Filter<ExptLocator,ExptLocator> {
    
    private String version;

    public LocatorVersionFilter(String v) {
        version = v;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.ewok.verbs.Filter#execute(null)
     */
    public ExptLocator execute(ExptLocator a) {
        if(a.getNameVersion().version.equals(version)) { 
            return a;
        } else { 
            return null;
        }
    }

}

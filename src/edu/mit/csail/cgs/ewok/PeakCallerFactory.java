/*
 * Created on Sep 29, 2006
 */
package edu.mit.csail.cgs.ewok;

import edu.mit.csail.cgs.datasets.locators.ExptLocator;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.verbs.binding.PeakCaller;

/**
 * @author tdanford
 */
public interface PeakCallerFactory {
    public PeakCaller createCaller(Genome g, ExptLocator loc, Object args);
}

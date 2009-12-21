/*
 * Created on Apr 2, 2006
 */
package edu.mit.csail.cgs.ewok.utils;

import java.sql.SQLException;
import java.util.*;

import edu.mit.csail.cgs.datasets.binding.BindingExtent;
import edu.mit.csail.cgs.datasets.chipchip.Probe;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.locators.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.verbs.binding.PeakCaller;
import edu.mit.csail.cgs.ewok.verbs.binding.RegionProber;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

/**
 * @author tdanford
 */
public interface EwokBinding<Locator extends ExptLocator, P extends Probe> {
    public RegionProber<P> getProber(Locator loc);
    public Expander<Region,BindingExtent> getPeakCaller(Locator loc, double[] params);
    public Collection<Locator> getAllLocators() throws SQLException, UnknownRoleException;
    public EwokBase getBase();
}
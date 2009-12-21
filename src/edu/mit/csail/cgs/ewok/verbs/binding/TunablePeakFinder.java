/*
 * Created on Mar 10, 2006
 */
package edu.mit.csail.cgs.ewok.verbs.binding;

import edu.mit.csail.cgs.datasets.chipchip.Probe;

/**
 * @author tdanford
 */
public interface TunablePeakFinder<X extends Probe> extends PeakFinder<X> {
    public double getMaxParameter();
    public double getMinParameter();
    public double getCurrentParameter();
    public void setParameter(double p);
}

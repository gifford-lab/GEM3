/*
 * Created on Mar 3, 2006
 */
package edu.mit.csail.cgs.ewok.verbs.binding;

import java.util.*;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.chipchip.Probe;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.*;

/**
 * @author tdanford
 */
public class ThresholdPeakFinder<X extends Probe> implements TunablePeakFinder<X> {
    
    private double thresh;
    private int maxDist;
    private String key;
    private int valueIndex;
    
    private Probe current;
    private int leftDist, rightDist;

    public ThresholdPeakFinder(String exptName, double t) {
        key = exptName;
        thresh = t;
        maxDist = 1000;
        valueIndex = 0;

        current = null;
        leftDist = rightDist = -1;
    }
    
    public ThresholdPeakFinder(String exptName, double t, int vi) { 
        key = exptName;
        thresh = t;
        maxDist = 1000;
        valueIndex = vi;

        current = null;
        leftDist = rightDist = -1;
    }
    
    public void reset() { 
        current = null;
        leftDist = rightDist = -1;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.ewok.verbs.Filter#execute(null)
     */
    public BindingEvent execute(X p) {
        if(current == null) { 
            current = p;
            leftDist = p.getLocation() - maxDist; 
            rightDist = p.getLocation() + maxDist; 
            return null;
        }
        
        rightDist = p.getLocation() - current.getLocation();
        BindingEvent ret = getCurrent();
        
        leftDist = rightDist;
        rightDist = p.getLocation() + maxDist;
        current = p;
        
        return ret;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.ewok.verbs.Distiller#getCurrent()
     */
    public BindingEvent getCurrent() {
        if(current != null) { 
            double val = current.getValue(key)[valueIndex];
            if(val >= thresh) {
                int s = current.getLocation() - leftDist;
                int e = current.getLocation() + rightDist;
                if(e < s) { e = s + 1; }
                return new BindingEvent(current.getGenome(), current.getChrom(), s, e, val, 1.0, key); 
            }
        }        
        return null;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.ewok.verbs.binding.TunablePeakFinder#getMaxParameter()
     */
    public double getMaxParameter() {
        return Double.MAX_VALUE;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.ewok.verbs.binding.TunablePeakFinder#getMinParameter()
     */
    public double getMinParameter() {
        return 0.0;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.ewok.verbs.binding.TunablePeakFinder#getCurrentParameter()
     */
    public double getCurrentParameter() {
        return thresh;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.ewok.verbs.binding.TunablePeakFinder#setParameter(double)
     */
    public void setParameter(double p) {
        thresh = p;
    }

}

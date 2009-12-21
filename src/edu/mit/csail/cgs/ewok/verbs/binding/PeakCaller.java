/*
 * Created on Mar 10, 2006
 */
package edu.mit.csail.cgs.ewok.verbs.binding;

import java.util.*;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.chipchip.Probe;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.utils.Closeable;

/**
 * @author tdanford
 * 
 * PeakCaller is meant to be the basic unit of any Hyperdrive-reimplementation in 
 * ewok.  It is an Expander that turns Regions into BindingEvents, but might (as 
 * in RegionProber) hold PeakCaller-specific methods in the future (for determining
 * thresholds, for instance?)
 * 
 */
public interface PeakCaller extends Expander<Region,BindingEvent> {
    
    public static class FromTupler<X extends Probe> implements PeakCaller { 

        private RegionProber<X> prober;
        private Filter<Vector<X>,BindingEvent> finder;
        private Mapper<Iterator<X>,Iterator<Vector<X>>> tupler;
        
        public FromTupler(RegionProber<X> p, Filter<Vector<X>,BindingEvent> f, int tupleLength) { 
            prober = p;
            finder = f;
            tupler = new Tupler<X>(tupleLength);
        }

        /* (non-Javadoc)
         * @see edu.mit.csail.cgs.ewok.verbs.Expander#execute(java.lang.Object)
         */
        public Iterator<BindingEvent> execute(Region a) {
            Iterator<X> probes = prober.execute(a);
            Iterator<Vector<X>> tuples = tupler.execute(probes);
            return new FilterIterator<Vector<X>,BindingEvent>(finder, tuples);
        }
    }
    
    public static class FromFinder<X extends Probe> implements PeakCaller { 

        private RegionProber<X> prober;
        private PeakFinder<X> finder;
        
        public FromFinder(RegionProber<X> p, PeakFinder<X> f) { 
            prober = p;
            finder = f;
        }

        /* (non-Javadoc)
         * @see edu.mit.csail.cgs.ewok.verbs.Expander#execute(java.lang.Object)
         */
        public Iterator<BindingEvent> execute(Region a) {
            return new DistillerIterator<X,BindingEvent>(finder, prober.execute(a));
        }
    }
    
    public static class Wrapper implements PeakCaller {
        private Expander<Region,BindingEvent> internal;
        
        public Wrapper(Expander<Region,BindingEvent> e) { 
            internal = e;
        }
        /* (non-Javadoc)
         * @see edu.mit.csail.cgs.ewok.verbs.Expander#execute(java.lang.Object)
         */
        public Iterator<BindingEvent> execute(Region a) {
            return internal.execute(a);
        }
    }
}

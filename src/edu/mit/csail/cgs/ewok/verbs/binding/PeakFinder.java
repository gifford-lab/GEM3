/*
 * Created on Mar 10, 2006
 */
package edu.mit.csail.cgs.ewok.verbs.binding;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.chipchip.Probe;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.nouns.*;

/**
 * @author tdanford
 */
public interface PeakFinder<X extends Probe> extends Distiller<X,BindingEvent> {
    
    public static class Wrapper<A extends Probe> implements PeakFinder<A> {
        
        private Distiller<A,BindingEvent> dist;
        
        public Wrapper(Distiller<A,BindingEvent> d) { 
            dist = d;
        }

        /* (non-Javadoc)
         * @see edu.mit.csail.cgs.ewok.verbs.Distiller#execute(java.lang.Object)
         */
        public BindingEvent execute(A a) {
            return dist.execute(a);
        }

        /* (non-Javadoc)
         * @see edu.mit.csail.cgs.ewok.verbs.Distiller#getCurrent()
         */
        public BindingEvent getCurrent() {
            return dist.getCurrent();
        }

        public void reset() {
            dist.reset();
        } 
        
    }
}

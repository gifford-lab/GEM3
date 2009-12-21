/*
 * Created on Mar 10, 2006
 */
package edu.mit.csail.cgs.ewok.verbs.binding;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.nouns.*;

/**
 * @author tdanford
 */
public interface BindingEventAnalyzer extends Mapper<BindingEvent,Double> {
    
    public static class Size implements BindingEventAnalyzer {
        
        public Size() {}

        /* (non-Javadoc)
         * @see edu.mit.csail.cgs.ewok.verbs.Mapper#execute(java.lang.Object)
         */
        public Double execute(BindingEvent a) {
            return a.getSize();
        } 
        
    }
    
    public static class Confidence implements BindingEventAnalyzer {
        
        public Confidence() {}

        /* (non-Javadoc)
         * @see edu.mit.csail.cgs.ewok.verbs.Mapper#execute(java.lang.Object)
         */
        public Double execute(BindingEvent a) {
            return a.getConf();
        } 
        
    }
}

/*
 * Created on Mar 10, 2006
 */
package edu.mit.csail.cgs.ewok.verbs.binding;

import edu.mit.csail.cgs.datasets.chipchip.Probe;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.*;

/**
 * @author tdanford
 */
public interface ProbeAnalyzer extends Filter<Probe,Double> {
    
    public static class Index implements ProbeAnalyzer {
        
        private String key;
        private int index;
        
        public Index(String k, int i) { 
            key = k; 
            index = i;
        }

        /* (non-Javadoc)
         * @see edu.mit.csail.cgs.ewok.verbs.Mapper#execute(java.lang.Object)
         */
        public Double execute(Probe a) {
            if(!a.containsKey(key) || index > a.getValue(key).length) { 
                return null;
            }
            return a.getValue(key)[index]; 
        } 
        
    }
    
    public static class BayesPosterior extends Index {
        public BayesPosterior(String k) { 
            super(k, 1);
        }
    }
    
    public static class BayesStrength extends Index {
        public BayesStrength(String k) { 
            super(k, 0);
        }
    }

}

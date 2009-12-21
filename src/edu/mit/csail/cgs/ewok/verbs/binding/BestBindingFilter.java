/*
 * Created on Mar 10, 2006
 */
package edu.mit.csail.cgs.ewok.verbs.binding;

import java.util.*;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.nouns.*;

/**
 * @author tdanford
 */
public class BestBindingFilter implements Filter<Region,BindingEvent> {
    
    private PeakCaller caller;
    private BindingEventAnalyzer analyzer;
    
    public BestBindingFilter(PeakCaller pc, BindingEventAnalyzer a) { 
        caller = pc;
        analyzer = a;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.ewok.verbs.Filter#execute(java.lang.Object)
     */
    public BindingEvent execute(Region a) {
        Iterator<BindingEvent> evtItr = caller.execute(a);
		BindingEvent ret = null;
        double maxScore = -Double.MAX_VALUE;
        
        while(evtItr.hasNext()) { 
            BindingEvent e = evtItr.next();
            double score = analyzer.execute(e);
            if(ret == null || score > maxScore) { 
                maxScore = score;
                ret = e;
            }
        }

		return ret;
    }
}

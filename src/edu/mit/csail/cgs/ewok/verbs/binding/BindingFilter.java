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
public class BindingFilter implements Filter<Region,Region> {
    
    private PeakCaller caller;
    
    public BindingFilter(PeakCaller pc) { 
        caller = pc;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.ewok.verbs.Filter#execute(java.lang.Object)
     */
    public Region execute(Region a) {
        Iterator<BindingEvent> evtItr = caller.execute(a);
		Region ret = null;
        if(evtItr.hasNext()) { ret = a; }

		/*
		System.out.println("\n" + a.toString() + " --> {");
		while(evtItr.hasNext()) { 
			System.out.println(evtItr.next());
		}
		System.out.println("}");
		*/

		return ret;
    }
}

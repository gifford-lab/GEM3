/**
 * 
 */
package edu.mit.csail.cgs.warpdrive.model;

import java.util.*;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.locators.ExptLocator;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.verbs.binding.PeakCaller;
import edu.mit.csail.cgs.utils.Closeable;

/**
 * @author tdanford
 *
 */
public class BindingEventRepository implements Expander<Region,BindingEvent>, Runnable {
	
	private LinkedList<Region> toSearch;
	private PeakCaller caller;
	private ExptLocator loc;
	private boolean isFinished;
	
	private LinkedList<BindingEvent> pending, discovered;
	
	public BindingEventRepository(ExptLocator el, PeakCaller call, Collection<Region> search) { 
		pending = new LinkedList<BindingEvent>();
		discovered = new LinkedList<BindingEvent>();
		
		toSearch = new LinkedList<Region>(search);
		caller = call;
		loc = el;
		
		isFinished = false;
	}
	
	public ExptLocator getLocator() { return loc; }
	public boolean isFinished() { return isFinished; }
	
	public synchronized void addEvent(BindingEvent evt) { 
		pending.addLast(evt);
	}

	public Iterator<BindingEvent> execute(Region a) {
		synchronized(this) { 
			if(!pending.isEmpty()) { 
				discovered.addAll(pending);
				pending.clear();
			}
		}
		LinkedList<BindingEvent> overlaps = new LinkedList<BindingEvent>();
		for(BindingEvent e : discovered) { 
			if(e.overlaps(a)) { 
				overlaps.addLast(e);
			}
		}
		return overlaps.iterator();
	}
	
	public void run() { 
		while(!toSearch.isEmpty()) { 
			Region r = toSearch.removeFirst();
			Iterator<BindingEvent> itr = caller.execute(r);
			while(itr.hasNext()) { 
				addEvent(itr.next());
			}
		}
		
		if(caller instanceof Closeable) { 
			Closeable c = (Closeable)caller;
			if(!c.isClosed()) { 
				c.close();
			}
		}
		
		isFinished = true;
	}	
}


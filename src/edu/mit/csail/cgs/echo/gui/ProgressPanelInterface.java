package edu.mit.csail.cgs.echo.gui;

import java.util.*;

import edu.mit.csail.cgs.utils.EventSource;
import edu.mit.csail.cgs.utils.Listener;

public interface ProgressPanelInterface {
	
	public void registerProcessor(String id, String description);
	public void unregisterProcessor(String id);
	public void markProcessorStart(String id);
	public void markProcessorEnd(String id);
	
	public static class ProcessorPanelImpl 
		implements ProgressPanelInterface, EventSource<ChangedEvent> {
		
		private TreeMap<String,String> totalIDs;
		private Map<String,Integer> markedCounts;
		private TreeSet<String> currentIDs;
		private EventSource.Default<ChangedEvent> src;
		
		public ProcessorPanelImpl() { 
			totalIDs = new TreeMap<String,String>();
			currentIDs = new TreeSet<String>();
			markedCounts = new HashMap<String,Integer>();
			src = new EventSource.Default<ChangedEvent>();
		}
		
		public synchronized Vector<String> getTotalIDs() { 
			return new Vector<String>(totalIDs.keySet()); 
		}
		
		public synchronized Vector<String> getCurrentIDs() { 
			return new Vector<String>(currentIDs); 
		}
		
		public synchronized String getDescription(String id) { 
			return totalIDs.get(id);
		}
		
		public synchronized int getCount(String id) { 
			return markedCounts.get(id);
		}

		public synchronized void markProcessorEnd(String id) {
			currentIDs.remove(id);
			src.fireEvent(new ChangedEvent(this));
		}

		public synchronized void markProcessorStart(String id) {
			currentIDs.add(id);
			markedCounts.put(id, markedCounts.get(id)+1);
			src.fireEvent(new ChangedEvent(this));
		}

		public synchronized void registerProcessor(String id, String description) {
			totalIDs.put(id, description);
			markedCounts.put(id, 0);
		}

		public synchronized void unregisterProcessor(String id) {
			totalIDs.remove(id);
			markedCounts.remove(id);
		}

		public void addEventListener(Listener<ChangedEvent> el) {
			src.addEventListener(el);
		}

		public boolean hasListeners() {
			return src.hasListeners();
		}

		public void removeEventListener(Listener<ChangedEvent> el) {
			src.removeEventListener(el);
		} 
		
	}
}

package edu.mit.csail.cgs.metagenes;

import java.util.*;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.chipseq.*;

public class MetaProfileHandler<T extends Point, ProfileClass extends Profile> {

	private MetaProfile profile;
	private PointProfiler<T, ProfileClass> profiler, threadSafe;
	private Vector<PointAddingThread> currentlyAdding;
	
	public MetaProfileHandler(String name, BinningParameters bps, PointProfiler<T,ProfileClass> pp, boolean normalizedMeta) { 
		if(normalizedMeta)
			profile = new NormalizedMetaProfile(name, bps);
		else
			profile = new MetaProfile(name, bps);		
		profiler = pp;
		threadSafe = new ThreadSafeProfiler();
		currentlyAdding = new Vector<PointAddingThread>();
	}
	
	public MetaProfile getProfile() { return profile; }
	
	public boolean addingPoints(){
		synchronized(currentlyAdding) { 
			for(PointAddingThread pat : currentlyAdding) { 
				if(pat.running){return true;}
			}
		}return(false);
	}
	
	public void addPoints(Collection<T> points) { 
		addPoints(points.iterator());
	}
	
	public void addPoints(Iterator<T> points) { 
		PointAddingThread pat = new PointAddingThread(points);
		startAddingThread(pat);
	}
	
	private void startAddingThread(PointAddingThread pat) { 
		synchronized(currentlyAdding) { 
			currentlyAdding.add(pat);
			Thread t = new Thread(pat);
			t.start();		
		}
	}
	
	private void addingThreadFinished(PointAddingThread pat) { 
		synchronized(currentlyAdding) { 
			currentlyAdding.remove(pat);
		}
	}
	
	public void stopAllAddingThreads() { 
		synchronized(currentlyAdding) { 
			for(PointAddingThread pat : currentlyAdding) { 
				pat.stopAdding();
			}
		}
	}
	
	private class ThreadSafeProfiler implements PointProfiler<T,ProfileClass> {

		public BinningParameters getBinningParameters() {
			return profiler.getBinningParameters();
		}

		public ProfileClass execute(T a) {
			synchronized(this) { 
				return profiler.execute(a);
			}
		} 	
		public void cleanup() {}
	}
	
	private class PointAddingThread implements Runnable { 
		
		public boolean running;
		private Iterator<T> points;

		public PointAddingThread(Iterator<T> pts) { 
			running = true;
			points = pts;
		}
		
		public void stopAdding() { 
			running = false;
		}
		
		public void run() { 
			while(running && points.hasNext()) { 
				T pt = points.next();
				profile.addProfile(threadSafe.execute(pt));
			}
			running=false;
			addingThreadFinished(this);
		}
	}
}

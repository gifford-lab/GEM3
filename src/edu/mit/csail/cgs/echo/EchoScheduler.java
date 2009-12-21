/*
 * 
 */
package edu.mit.csail.cgs.echo;

import java.util.*;

import edu.mit.csail.cgs.echo.gui.*;

/**
 * @author tdanford
 * 
 * EchoScheduler is another core class of Echo -- this is an object, which holds
 * several EchoProcessors and runs them sequentially (until all are finished) 
 * in a separate thread.  
 * 
 * This is the point of entry, for adding support on multiple machines or across
 * multiple threads.
 */
public class EchoScheduler implements Runnable {

	private Vector<EchoProcessor> procs;
	private Map<EchoProcessor,String> procIDs;
	
	private GUIProgressPanel.Frame progressFrame;
	private ProgressPanelInterface progress;
	private int lastID;
	
	public EchoScheduler(Collection<EchoProcessor> p) { 
		procs = new Vector<EchoProcessor>(p);
		progressFrame = new GUIProgressPanel.Frame();
		progress = progressFrame.getProgressPanel().getInterface();
		
		lastID = 0;
		procIDs = new HashMap<EchoProcessor,String>();
		for(EchoProcessor proc : procs) { 
			procIDs.put(proc, getNextID());
		}
	}
	
	private String getNextID() { 
		return "proc" + (++lastID);
	}
	
	public void run() {
		int index = 0;
        long iter = 0;

        for(EchoProcessor proc : procIDs.keySet()) {
        	progress.registerProcessor(procIDs.get(proc), proc.toString());
        }
        
		while(!procs.isEmpty()) { 
			EchoProcessor proc = procs.get(index);
			
			if(proc.isFinished()) { 
                proc.finish();
				procs.remove(index);
				progress.unregisterProcessor(procIDs.get(proc));
			} else {
				if(proc.isReady()) {
					progress.markProcessorStart(procIDs.get(proc));
					proc.process();
					progress.markProcessorEnd(procIDs.get(proc));
				}
                
				index += 1;
			}
			
			if(index >= procs.size()) {
                index = 0; 
			}
			
            iter++;
            
            /*
            try {
                Thread.sleep((long)20);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            */
		}		
	}
}

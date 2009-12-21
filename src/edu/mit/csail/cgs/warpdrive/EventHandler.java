/*
 * Created on Mar 24, 2006
 */
package edu.mit.csail.cgs.warpdrive;

import java.util.*;
import edu.mit.csail.cgs.utils.*;

/**
 * @author tdanford
 */
public class EventHandler<X extends EventObject> implements Listener<X> {
    
    private LinkedList<X> eventQueue;
    private LinkedList<Listener<X>> listeners;
    private HandlerRunner runner;

    public EventHandler() {
        eventQueue = new LinkedList<X>();
        listeners = new LinkedList<Listener<X>>();
        
        Thread t = new Thread((runner = new HandlerRunner()));
        t.start();
    }
    
    public synchronized void addEventListener(Listener<X> m) { 
        listeners.addLast(m);
    }
    
    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.Listener#eventRegistered(java.lang.Object)
     */
    public synchronized void eventRegistered(X e) {
        System.out.println("+++ Adding HAE to Event Queue");
        eventQueue.addLast(e);
        //eventQueue.notifyAll();
        handlerNotify();
    }
    
    private synchronized void handlerWait() { 
        try {
            wait();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }
    
    private synchronized void handlerNotify() { 
        notifyAll();
    }
    
    private X getNextEvent() {
        X evt = null;
        synchronized(this) {
            while(evt == null) { 
                handlerWait();
                if(!eventQueue.isEmpty()) { 
                    System.out.println("--- Removing HAE From Event Queue");
                    evt = eventQueue.removeFirst(); 
                }
            }
        }
        return evt;
    }

    private class HandlerRunner implements Runnable {
        
        private boolean continueRunning; 
        
        public HandlerRunner() { 
            continueRunning = true;
        }
        
        public void stop() { continueRunning = false; }

        /* (non-Javadoc)
         * @see java.lang.Runnable#run()
         */
        public void run() {
            while(continueRunning) { 
                X evt = getNextEvent();
                
                synchronized(EventHandler.this) {
                    for(Listener<X> m : listeners) {
                        m.eventRegistered(evt); 
                    }
                }
            }
        }
    }
}

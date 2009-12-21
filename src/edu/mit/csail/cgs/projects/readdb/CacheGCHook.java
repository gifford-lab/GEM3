package edu.mit.csail.cgs.projects.readdb;

import java.util.logging.*;

/**
 * wait around and run garbage collection
 */
public class CacheGCHook implements Runnable {
	private Logger logger;
    
    public CacheGCHook(Logger l) {
        logger = l;
    }

    public void run() {
        while (true) {
            if (LRUCache.removed() > 200) {
                logger.log(Level.INFO,"running GC");
                LRUCache.resetRemoved();
                System.gc();
                System.runFinalization();
            } 
            try {
                Thread.sleep(500);
            } catch (InterruptedException e) {

            }
        }
    }

}
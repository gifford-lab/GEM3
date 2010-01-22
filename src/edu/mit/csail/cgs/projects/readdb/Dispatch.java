package edu.mit.csail.cgs.projects.readdb;

import java.util.*;
import java.util.logging.*;

/**
 * Dispatch is a thread that manages the WorkerThreads and assigns
 * WorkerThreads to ServerTasks.
 */

public class Dispatch implements Runnable {

    private Vector<WorkerThread> freePool, allThreads;
    private Vector<ServerTask> workQueue;
    private Server server;
    private int maxConnections;
    private int warnedMaxConn = 0;
    public Dispatch (Server s, int numThreads, int maxC) {
        server = s;
        workQueue = new Vector<ServerTask>();
        freePool = new Vector<WorkerThread>();
        allThreads = new Vector<WorkerThread>();
        for (int i = 0; i < numThreads; i++) {
            WorkerThread servthread = new WorkerThread(this);
            Thread t = new Thread(servthread);
            t.start();
            freePool.add(servthread);
            allThreads.add(servthread);
        }
        maxConnections = maxC;
    }

    /**
     * Add a new ServerTask to the set of tasks that will
     * be run.  Called by Server when it accepts a new connection.
     */
    public void addWork(ServerTask s) {
        while (workQueue.size() > maxConnections) {
            try {
                if (warnedMaxConn++ % 100 == 0) {
                    server.getLogger().log(Level.WARNING,(String.format("Hit maxconnections (%d)",maxConnections)));
                }
                Thread.sleep(10);
            } catch (InterruptedException e) {

            }
        }
        warnedMaxConn = 0;
        workQueue.add(s);
    }
    /**
     * called by WorkerThread when it's finished with a ServerTask.
     * WorkerThread has called the run method(), so getting to this
     * point means that the run method returned.  The run() method
     * should not throw any exceptions.
     *
     */
    public void freeThread(WorkerThread t, ServerTask s) {
        if (s.shouldClose()) {
            s.close();
            System.err.println("freeThread closing task " + s);
        } else {
            workQueue.add(s);
        }
        freePool.add(t);
    }
    /**
     * our main loop.  work through the tasks, seeing who appears
     * to have input
     */
    public void run() {
        int i = 0;
        while (server.keepRunning()) {            
            if (workQueue.size() > 0) {
                ServerTask s = workQueue.remove(0);
                if (s.shouldClose()) {
//                     if (server.debug()) {
//                         System.err.println("dispatch should close " + s);
//                     }
                    System.err.println("run loop closing task " + s);
                    s.close();
                } else if (s.inputAvailable()) {
                    if (freePool.size() == 0) {
                        workQueue.add(s);
//                         if (server.debug()) {
//                             System.err.println("dispatch would like to handle " + s + " but no threads");
//                         }
                        continue;
                    } else {
                        WorkerThread w = freePool.remove(0);
//                         if (server.debug()) {
//                             System.err.println("dispatch is handling " + s);
//                         }
                        w.handle(s);
                        i = 0;
                    }
                } else {
                    i++;
//                     if (server.debug()) {
//                         System.err.println("No input for " + s + " so putting it back");
//                     }
                    workQueue.add(s);
                }
            }
            try {
//                 if (server.debug()) {
//                     Thread.sleep(1000); // debugging only
//                 }

                if (workQueue.size() == 0) {
                    Thread.sleep(50);
                } else if (i > workQueue.size() * 10000) {
                    Thread.sleep(500);
                    i = 0;
                } else if (i > workQueue.size() * 1000) {
                    Thread.sleep(100);
                    HashSet<ServerTask> set = new HashSet<ServerTask>();
                    set.addAll(workQueue);
                    if (set.size() < workQueue.size()) {
                        System.err.println("Work Queue size is " + workQueue.size() + " but only " + set.size() + " unique elements");
                    }
                } else if (i > workQueue.size() * 100) {
                    Thread.sleep(10);
                } else {
                    //                    Thread.yield();
                }
            } catch (InterruptedException e) {  }
        }
        while (freePool.size() < allThreads.size()) {
            try {
                Thread.sleep(100);
            } catch (InterruptedException e) {
                
            }
        }
        for (WorkerThread t : freePool) {
            t.stopRunning();
        }
    }



}
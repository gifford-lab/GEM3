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
        int noInputAvailable = 0;
        int ranNothing = 0;
        while (server.keepRunning()) {            
            if (workQueue.size() > 0) {
                ranNothing = 0;
                ServerTask s = workQueue.remove(0);
                if (s.shouldClose()) {
                    System.err.println("run loop closing task " + s);
                    s.close();
                    noInputAvailable = 0;
                } else if (s.inputAvailable()) {
                    if (freePool.size() == 0) {
                        workQueue.add(s);
                        noInputAvailable = 0;
                        continue;
                    } else {
                        WorkerThread w = freePool.remove(0);
                        w.handle(s);
                        noInputAvailable = 0;
                    }
                } else {
                    noInputAvailable++;
                    workQueue.add(s);
                }
            } else {
                ranNothing++;
            }
            try {
                if (ranNothing > 1) {
                    Thread.sleep(2);
                } else if (noInputAvailable > workQueue.size() * 10000) {
                    Thread.sleep(500);
                    noInputAvailable = 0;
                } else {
                    Thread.yield();
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
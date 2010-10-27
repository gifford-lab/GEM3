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
    private Vector<Thread> threads;
    private Server server;
    private int maxConnections;
    private int warnedMaxConn = 0;
    public Dispatch (Server s, int numThreads, int maxC) {
        server = s;
        workQueue = new Vector<ServerTask>();
        freePool = new Vector<WorkerThread>();
        allThreads = new Vector<WorkerThread>();
        threads = new Vector<Thread>();
        for (int i = 0; i < numThreads; i++) {
            WorkerThread servthread = new WorkerThread(this);
            Thread t = new Thread(servthread);
            t.start();
            threads.add(t);
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
        System.err.println("workqueue size is " + workQueue.size() + "\n");
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
        synchronized(this) {
            notifyAll();
        }
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
        } else {
            workQueue.add(s);
        }
        freePool.add(t);
        synchronized(this) {
            notifyAll();
        }
    }
    /**
     * our main loop.  work through the tasks, seeing who appears
     * to have input
     */
    public void run() {
        int noTasksWaiting = 0;
        int noInputAvailable = 0;
        int threadCheck = 0;
        int sleep = server.getSleepiness();
        while (server.keepRunning()) {            
            if (workQueue.size() > 0) {
                noTasksWaiting = 0;
                ServerTask s = workQueue.remove(0);
                if (s.shouldClose()) {
                    s.close();
                } else if (s.inputAvailable()) {
                    while (freePool.size() == 0) {                        
                        try {
                            synchronized(this) {
                                wait(2);
                            }
                        } catch (InterruptedException e) {}
                    }                        
                    WorkerThread w = freePool.remove(0);
                    w.handle(s);
                    noInputAvailable = 0;
                } else {
                    noInputAvailable++;
                    workQueue.add(s);
                    if (noInputAvailable > workQueue.size() * 2) {
                        try {
                            synchronized(this) {
                                wait(sleep);
                            }
                        } catch (InterruptedException e) {}                
                        noInputAvailable = 0;
                    }
                }
            } else {
                noTasksWaiting++;
                try {
                    synchronized(this) {
                        wait(sleep*2);
                    }
                } catch (InterruptedException e) {}                
            } 
            if (threadCheck++ > 1000) {
                threadCheck = 0;
                for (int i = 0; i < threads.size(); i++) {
                    if (!threads.get(i).isAlive()) {
                        server.getLogger().log(Level.INFO,"Dispatch","run: DEAD THREAD.  Adding a new one");
                        WorkerThread servthread = new WorkerThread(this);
                        Thread t = new Thread(servthread);
                        t.start();
                        threads.set(i,t);
                        freePool.add(servthread);
                        try {
                            allThreads.get(i).stopRunning();
                        } catch (Exception e) {
                            server.getLogger().logp(Level.INFO,"Dispatch","run: trying to stop old thread",e.toString(),e);
                        }
                        allThreads.set(i,servthread);
                    }
                }
            }
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
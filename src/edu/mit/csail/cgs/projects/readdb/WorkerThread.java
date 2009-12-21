package edu.mit.csail.cgs.projects.readdb;

/*
 * WorkerThread is an actual thread.  It runs tasks (ServerTask) by
 * calling their run() method.  When run() returns, WorkerThread
 * returns the task to Dispatch.
 */

public class WorkerThread implements Runnable {

    private ServerTask task;
    private boolean keepRunning;
    private Dispatch dispatch;

    public WorkerThread(Dispatch d) {
        task = null;
        keepRunning = true;
        dispatch = d;
    }

    public synchronized void stopRunning() {keepRunning = false;}

    public synchronized void handle(ServerTask t) {
        task = t;
        notifyAll();
    }

    public synchronized void run() {
        while (keepRunning) {
            if (task == null) {
                try {
                    wait();
                } catch (InterruptedException e) {
                    // swallow it and go back to waiting.
                }
            } else {
                try {
                    task.run();
                } catch (Exception e) {
                    e.printStackTrace();
                }
                ServerTask t = task;
                task = null;
                dispatch.freeThread(this, t);
            }
        }

    }



}
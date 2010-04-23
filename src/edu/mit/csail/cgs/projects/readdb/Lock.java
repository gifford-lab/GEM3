package edu.mit.csail.cgs.projects.readdb;

import java.util.*;
import java.util.concurrent.locks .*;

/**
 * Static class to maintain locks on objects within the JVM.  We need this because...
 * - multiple threads can have a read lock on an object at once; synchronized doesn't
 *   do this
 * - the java.io locking stuff only locks between processes, not within the JVM
 */
public class Lock {

    private static Map<String,ReentrantReadWriteLock> locks = Collections.synchronizedMap(new HashMap<String,ReentrantReadWriteLock>());
    private static Map<Thread,Set<java.util.concurrent.locks.Lock>> threadlocks = Collections.synchronizedMap(new HashMap<Thread,Set<java.util.concurrent.locks.Lock>>());
    private static int rlcount = 0;

    /**
     * blocks to acquire a shared lock to the specified file.
     * The locking is implemented in java, so it's only
     * good for keeping out other threads.  It's keyed off the file name
     * you provide, so make sure you always generate the filename in the same way.
     */
    protected static  java.util.concurrent.locks.Lock readLock(String fname) {
        synchronized(locks) {
            if (!locks.containsKey(fname)) {
                locks.put(fname, new ReentrantReadWriteLock());
            }
        }
        Thread t = Thread.currentThread();
        synchronized(threadlocks) {
            if (!threadlocks.containsKey(t)) {
                threadlocks.put(t, new HashSet<java.util.concurrent.locks.Lock>());
            }
        }
        java.util.concurrent.locks.Lock lock = locks.get(fname).readLock();
        lock.lock();
        threadlocks.get(t).add(lock);
        //        System.err.println("READLOCK by " + t + " of " + fname + " as " + lock);
        return lock;
    }
    protected static java.util.concurrent.locks.Lock writeLock(String fname) {
        synchronized(locks) {
            if (!locks.containsKey(fname)) {
                locks.put(fname, new ReentrantReadWriteLock());
            }
        }
        Thread t = Thread.currentThread();
        synchronized(threadlocks) {
            if (!threadlocks.containsKey(t)) {
                threadlocks.put(t, new HashSet<java.util.concurrent.locks.Lock>());
            }
        }
        java.util.concurrent.locks.Lock rl = locks.get(fname).readLock();
        rl.unlock();
        threadlocks.get(t).remove(rl);
        java.util.concurrent.locks.Lock lock = locks.get(fname).writeLock();
        lock.lock();
        threadlocks.get(t).add(lock);
        //        System.err.println("WRITELOCK by " + t + " of " + fname + " as " + lock);
        return lock;
    }
    /* call to ensure that all a thread's locks have been released */
    protected static void releaseLocks() {
        Thread t = Thread.currentThread();
        if (threadlocks.containsKey(t)) {
            for (java.util.concurrent.locks.Lock l : threadlocks.get(t)) {
                //                System.err.println("UNLOCK of " + l + " by " + t);
                l.unlock();
            }
            threadlocks.get(t).clear();
        }
        if (rlcount++ > 100) {
            rlcount = 0;
            if (locks.size() > 1000) {
                /* cleanup loop so that we don't have an ever-expanding data structure */
                Collection<String> keys = locks.keySet();
                for (String k : keys) {
                    ReentrantReadWriteLock l = locks.get(k);
                    /* if there are no outstanding locks, acquire a write lock to make sure that
                       nobody else can get it, remove the ReentrantRWLock, then release
                    */
                    if (l.getReadLockCount() == 0 && l.isWriteLocked()) {
                        java.util.concurrent.locks.Lock lock = l.writeLock();
                        locks.remove(k);
                        lock.unlock();
                        lock = null;
                    }
                }
            }
        }
    }



}
package edu.mit.csail.cgs.projects.readdb;

import java.util.*;

/**
 * Static class to maintain locks on objects within the JVM.  We need this because...
 * - multiple threads can have a read lock on an object at once; synchronized doesn't
 *   do this
 * - the java.io locking stuff only locks between processes, not within the JVM
 */
public class Lock {

    private static Map<String,Set<ServerTask>> readLocks = Collections.synchronizedMap(new HashMap<String,Set<ServerTask>>());
    private static Map<String,Set<ServerTask>> acquiringWrite = Collections.synchronizedMap(new HashMap<String,Set<ServerTask>>());
    private static Map<String,ServerTask> writeLocks = Collections.synchronizedMap(new HashMap<String,ServerTask>());

    /**
     * blocks to acquire a shared lock to the specified file.
     * The locking is implemented in java, so it's only
     * good for keeping out other threads.  It's keyed off the file name
     * you provide, so make sure you always generate the filename in the same way.
     */
    protected static void readLock(ServerTask locker, String fname) {
        boolean locked = false;
        /* return now if we already have a read lock.
           this is ok outside the synchronized since
           - a thread only acquires locks for itself
           - a thread therefore won't release its lock at the same
             time that it's acquiring it
        */
        try {
            if (readLocks.containsKey(fname) &&
                readLocks.get(fname).contains(locker)) {
                return;
            }
        } catch (NullPointerException e) {
            /* ignore it.  This happens if readLocks was modified to remove
               the fname key between the two and clauses above.
            */
        }
        while (!locked) {
            if (writeLocks.containsKey(fname) || acquiringWrite.containsKey(fname)) {
                Thread.yield();
                continue;
            }
            synchronized(writeLocks) {
                synchronized(acquiringWrite) {
                    if (writeLocks.containsKey(fname) || acquiringWrite.containsKey(fname)) {
                        Thread.yield();
                        continue;
                    }
                    synchronized(readLocks) {
                        if (!readLocks.containsKey(fname)) {
                            readLocks.put(fname, new HashSet<ServerTask>());
                        }
                        readLocks.get(fname).add(locker);
                        locked = true;
                    }
                }
            }
        }
    }
    protected static void readUnLock(ServerTask locker, String fname) {
        synchronized(readLocks) {
            if (readLocks.containsKey(fname)) {
                Set<ServerTask> set = readLocks.get(fname);
                set.remove(locker);
                if (set.size() == 0) {
                    readLocks.remove(fname);
                }
            }
        }
    }
    protected static void writeLock(ServerTask locker, String fname) {
        boolean locked = false;
        try {
            if (writeLocks.containsKey(fname) &&
                writeLocks.get(fname) == locker) {
                return;
            }
        } catch (NullPointerException e) {

        }
        synchronized(acquiringWrite) {                
            if (!acquiringWrite.containsKey(fname)) {
                acquiringWrite.put(fname,new HashSet<ServerTask>());
            }
            acquiringWrite.get(fname).add(locker);
        }
        while (!locked) {
            if (writeLocks.containsKey(fname)) {
                Thread.yield();
                continue;
            }
            synchronized(writeLocks) {                
                    synchronized(readLocks) {
                    boolean justus = !readLocks.containsKey(fname) ||
                        (readLocks.containsKey(fname) &&
                         readLocks.get(fname).size() == 1 &&
                         readLocks.get(fname).contains(locker));
                    
                    // can't acquire write lock if other people haveit read locked
                    if (writeLocks.containsKey(fname) ||
                        (readLocks.containsKey(fname) && !justus)) {
                        Thread.yield();
                        continue;
                    }
                    
                    writeLocks.put(fname, locker);
                    locked = true;
                }
            }
        }
        synchronized(acquiringWrite) {                
            acquiringWrite.get(fname).remove(locker);
            if (acquiringWrite.get(fname).size() == 0) {
                acquiringWrite.remove(fname);
            }
        }
    }
    protected static void writeUnLock(ServerTask locker, String fname) {
        synchronized(writeLocks) {
            if (writeLocks.containsKey(fname) && writeLocks.get(fname) == locker) {
                writeLocks.remove(fname);
            }
        }
    }
    /* call to ensure that all a thread's locks have been released */
    protected static void releaseLocks(ServerTask t) {
        synchronized(writeLocks) {
            Set<String> s = writeLocks.keySet();
            for (String k : s) {
                if (writeLocks.get(k) == t) {
                    writeLocks.remove(k);
                }
            }            
        }
        synchronized(acquiringWrite) {
            Set<String> s = acquiringWrite.keySet();
            for (String k : s) {
                if (acquiringWrite.get(k) == t) {
                    acquiringWrite.remove(k);
                }
            }            
        }
        synchronized(readLocks) {
            Set<String> s = readLocks.keySet();
            for (String k : s) {
                Set<ServerTask> l = readLocks.get(k);
                l.remove(t);
                if (l.size() == 0) {
                    readLocks.remove(l);
                }
            }
        }
    }



}
package edu.mit.csail.cgs.projects.readdb;

import java.util.*;

/**
 * Cache of Closeable objects.   Each object
 * is associated with a key.  When an object falls
 * out of the cache, its close() method is called.
 */

public class LRUCache<X extends Closeable> {

    private static int removed = 0;

    private List<String> ordered;
    private Map<String,X> map;
    private int size;

    public LRUCache(int size) {
        ordered = Collections.synchronizedList(new ArrayList<String>());
        map = Collections.synchronizedMap(new HashMap<String,X>());
        this.size = size;
    }

    public boolean contains(String k) {
        return map.containsKey(k);
    }
    public void printKeys() {
        System.err.println(ordered.toString());
    }
    public X get(String k) {
        synchronized(map) {
            if (map.containsKey(k)) {
                //                System.err.println("GETTING " +k);
                ordered.remove(k);
                ordered.add(k);
                return map.get(k);                
            } else {
                return null;
            }
        }
    }
    public void add(String k, X o) {
        synchronized(map) {
            //            System.err.println("ADDING " +k);
            if (map.containsKey(k)) {
                remove(k);
            }
            if (ordered.size() >= size) {
                String toRemove = ordered.get(0);
                remove(toRemove);
            }
            ordered.add(k);
            map.put(k,o);
        }
    }
    public void remove(String k) {
        synchronized(map) {
            if (map.containsKey(k)) {
                //                System.err.println("REMOVING " +k);
                ordered.remove(k);
                map.remove(k);
                removed++;
            }
        }
    }
    public static int removed() {return removed;}
    public static void resetRemoved() {removed = 0;}

}
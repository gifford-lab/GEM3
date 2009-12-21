/*
 * Created on Dec 8, 2005
 */
package edu.mit.csail.cgs.utils;

import java.util.LinkedList;

/**
 * @author tdanford
 */
public interface Predicate<X> {
    public boolean accepts(X v);

    public static class All<K> implements Predicate<K> {
		public boolean accepts(K v) {
			return true;
		} 
    }
    
    public static class Not<K> implements Predicate<K> {
    	
    	protected Predicate<K> pred;
    	
    	public Not(Predicate<K> p1) { 
    		pred = p1;
    	}

		public boolean accepts(K v) {
			return !pred.accepts(v);
		} 
    }

    public static class And<K> implements Predicate<K> {
    	
    	protected LinkedList<Predicate<K>> preds;
    	
    	public And(Predicate<K> p1, Predicate<K> p2) { 
    		preds = new LinkedList<Predicate<K>>();
    		preds.addLast(p1);
    		preds.addLast(p2);
    	}
    	
    	public And() { 
    		preds = new LinkedList<Predicate<K>>();
    	}

		public boolean accepts(K v) {
			for(Predicate<K> p : preds) { 
				if(!p.accepts(v)) { 
					return false;
				}
			}
			return true;
		} 
    }

    public static class Or<K> implements Predicate<K> {
    	
    	protected LinkedList<Predicate<K>> preds;
    	
    	public Or(Predicate<K> p1, Predicate<K> p2) { 
    		preds = new LinkedList<Predicate<K>>();
    		preds.addLast(p1);
    		preds.addLast(p2);
    	}

    	public Or() {  
    		preds = new LinkedList<Predicate<K>>();
    	}
    	
		public boolean accepts(K v) {
			for(Predicate<K> p : preds) { 
				if(p.accepts(v)) { 
					return true;
				}
			}
			return false;
		} 
    }
}

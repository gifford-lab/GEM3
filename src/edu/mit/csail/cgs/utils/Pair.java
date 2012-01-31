/*
 * Created on Jun 9, 2005
 */
package edu.mit.csail.cgs.utils;

/**
 * @author tdanford
 * 
 * Typed container for a pair of objects.  
 */
public class Pair<FirstClass, LastClass> {
    
    private FirstClass first;
    private LastClass last;
    
    public Pair(FirstClass first, LastClass last) {
        this.first = first;
        this.last = last;
    }

    public FirstClass getFirst() { return first; }
    public LastClass getLast() { return last; }
    public FirstClass car() {return first;}
    public LastClass cdr() {return last;}
    
    public String toString() { return "<" + first.toString() + "," + last.toString() + ">"; }
 
    public int hashCode() { 
        int code = 17; 
        code += first.hashCode(); code *= 37;
        code += last.hashCode(); code *= 37;
        return code;
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof Pair)) { return false; }
        Pair p = (Pair)o;
        return p.first.equals(first) && p.last.equals(last);
    }
    
    public Pair<LastClass, FirstClass> swap() {
    	return new Pair<LastClass,FirstClass>(this.cdr(),this.car());
    }
    
}

package edu.mit.csail.cgs.clustering;

import java.util.Set;
import java.util.HashSet;

public class SingletonCluster<X> implements Cluster<X> {
		
    private X value;
		
    public SingletonCluster(X v) { value = v; }
		
    public X getValue() { return value; }

    public Set<X> getElements() {
        HashSet<X> set = new HashSet<X>();
        set.add(value);
        return set;
    }

    public int size() {
        return 1;
    } 
    
    public String toString() {
    	return value.toString();
    }
		
}

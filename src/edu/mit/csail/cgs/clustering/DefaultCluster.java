package edu.mit.csail.cgs.clustering;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

/**
 * @author tdanford
 *
 */
public class DefaultCluster<X> implements Cluster<X> {
	
	private HashSet<X> elmts;
	
	public DefaultCluster() { 
		elmts = new HashSet<X>(); 
	}
	
	public DefaultCluster(X v) { 
		elmts = new HashSet<X>();
		elmts.add(v);
	}
	
	public DefaultCluster(Collection<X> v) { 
		elmts = new HashSet<X>(v);
	}
	
	public void addElement(X v) { elmts.add(v); }
	public void clear() { elmts.clear(); }

	/* (non-Javadoc)
	 * @see edu.mit.csail.cgs.clustering.Cluster#getElements()
	 */
	public Set<X> getElements() {
		return new HashSet<X>(elmts);
	}

	/* (non-Javadoc)
	 * @see edu.mit.csail.cgs.clustering.Cluster#size()
	 */
	public int size() {
		return elmts.size();
	}
	
	public String toString() { 
		StringBuilder sb = new StringBuilder();
		sb.append("cluster[");
		boolean first = true;
		for(X elmt : elmts) { 
			if(!first) { sb.append(","); }
			first = false;
			sb.append(elmt.toString());
		}
		sb.append("]");
		return sb.toString();
	}
}

/**
 * 
 */
package edu.mit.csail.cgs.utils.datastructures;

import java.util.*;

/**
 * @author Timothy Danford
 */
public class TaxonomyEvent extends EventObject {
	
	public static final int ELEMENT_ADDED = 0;
	public static final int ELEMENT_REMOVED = 1;
	public static final int ELEMENT_CHANGED = 2;
	public static final int TAXON_ADDED = 3;
	public static final int TAXON_REMOVED = 4;

	private LinkedList<String> addr;
	private Object leafValue;
	private int type;

	/**
	 * @param arg0
	 */
	public TaxonomyEvent(Object src, int type, Object leaf) {
		super(src);
		addr = new LinkedList<String>();
		leafValue = leaf;
		this.type = type;
	}

	public TaxonomyEvent(Object src, TaxonomyEvent te, String a) { 
		super(src);
		addr = new LinkedList<String>(te.addr); 
		addr.addFirst(a);
		leafValue = te.leafValue;
		type = te.type;
	}
	
	public Collection<String> getAddr() { return addr; }
	public Object getLeaf() { return leafValue; }
	public int getType() { return type; }
	public void addAddr(String a) { addr.addFirst(a); }
	public String removeAddr() { return addr.removeFirst(); }
	public int getAddrLength() { return addr.size(); }
	
	public String toString() { 
		return "(" + type + ":" + getAddrString() + ") --> " + leafValue.toString();
	}
	
	public String getAddrString() { 
		StringBuilder sb = new StringBuilder();
		boolean f = true;
		for(String a : addr) { 
			if(!f) { sb.append(","); }
			f = false;
			sb.append(a);
		}
		return sb.toString();
	}
}

/**
 * 
 */
package edu.mit.csail.cgs.utils.datastructures;

import java.util.*;
import java.io.PrintStream;

import edu.mit.csail.cgs.utils.EventSource;


/**
 * @author Timothy Danford
 *
 */
public interface Taxonomy<X> extends Cloneable, EventSource<TaxonomyEvent> {
	public Collection<X> getAllElements();
	public Collection<X> getImmediateElements();
	public Taxonomy<X> getSubTaxonomy(String subAddr);
	public Taxonomy<X> getSubTaxonomy(Collection<String> addr);
	public Set<String> getAddrs();
	public boolean hasAddr(String a);
	public void addElement(X v);
	public void addElement(Collection<String> addr, X v);
	public void recursiveRemoveElement(X v);
	public void addTaxonomy(String addr, Taxonomy<X> tax);
	public int size();
	public int getNumSubTaxonomies();
	public int getNumElements();
	public Taxonomy<X> combine(Taxonomy<X> t);
	public void print(PrintStream ps);
}

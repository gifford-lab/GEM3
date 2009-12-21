/**
 * 
 */
package edu.mit.csail.cgs.utils.datastructures;

import javax.swing.event.*;
import javax.swing.tree.TreeModel;
import javax.swing.tree.TreePath;

import java.util.*;
import edu.mit.csail.cgs.utils.Listener;
import edu.mit.csail.cgs.utils.datastructures.Taxonomy;
import edu.mit.csail.cgs.utils.datastructures.TaxonomyEvent;

/**
 * @author Timothy Danford
 */
public class TaxonomyTreeModel<X> implements TreeModel, Listener<TaxonomyEvent> {
	
	private Taxonomy<X> taxonomy;
	private TaxonomyAddr<X> root;
	private LinkedList<TreeModelListener> listeners;

	public TaxonomyTreeModel(Taxonomy<X> tax) {
		root = new TaxonomyAddr<X>(tax);
		listeners = new LinkedList<TreeModelListener>();
		taxonomy = tax;
		taxonomy.addEventListener(this);
	}
	
	public void eventRegistered(TaxonomyEvent te) {
		System.out.println("Taxonomy Event Registered: " + te.toString());
		root.handleTaxonomyEvent(this, te);		
	}
	
    public void addElement(Collection<String> addr, X e) { 
        taxonomy.addElement(addr, e);
    }
    
    public void removeElement(X e) { 
    	taxonomy.recursiveRemoveElement(e);
    }
    
    protected void fireTreeNodesInserted(Object[] path, int[] childInds, Object[] children) {
        TreeModelEvent e = new TreeModelEvent(this, path, childInds, children);
        for(TreeModelListener l : listeners) { 
            l.treeNodesInserted(e);
        }
    }
    
    protected void fireTreeNodesChanged(Object[] path, int[] childInds, Object[] children) { 
        TreeModelEvent e = new TreeModelEvent(this, path, childInds, children);
        for(TreeModelListener l : listeners) { 
            l.treeNodesChanged(e);
        }
    }
    
    protected void fireTreeNodesRemoved(Object[] path, int[] childInds, Object[] children) { 
        TreeModelEvent e = new TreeModelEvent(this, path, childInds, children);
        for(TreeModelListener l : listeners) { 
            l.treeNodesRemoved(e);
        }    	
    }

	/* (non-Javadoc)
	 * @see javax.swing.tree.TreeModel#getRoot()
	 */
	public Object getRoot() {
		return root;
	}

	/* (non-Javadoc)
	 * @see javax.swing.tree.TreeModel#getChild(java.lang.Object, int)
	 */
	public Object getChild(Object parent, int index) {
		TaxonomyAddr parentAddr = (TaxonomyAddr)parent;
		if(index < parentAddr.getNumElements()) { 
			return parentAddr.getElement(index);
		} else { 
			return parentAddr.getNextAddr(index-parentAddr.getNumElements());
		}
	}

	/* (non-Javadoc)
	 * @see javax.swing.tree.TreeModel#getChildCount(java.lang.Object)
	 */
	public int getChildCount(Object parent) {
		if(parent instanceof TaxonomyAddr) { 
			TaxonomyAddr ta = (TaxonomyAddr)parent;
			return ta.getNumElements() + ta.getNumSubTaxa();
		} else { 
			return 0;
		}
	}

	/* (non-Javadoc)
	 * @see javax.swing.tree.TreeModel#isLeaf(java.lang.Object)
	 */
	public boolean isLeaf(Object value) {
		if(value instanceof TaxonomyAddr) { 
			TaxonomyAddr ta = (TaxonomyAddr)value;
			return ta.getNumElements() == 0 && ta.getNumSubTaxa() == 0;
		} else { 
			return true;
		}
	}

	/* (non-Javadoc)
	 * @see javax.swing.tree.TreeModel#valueForPathChanged(javax.swing.tree.TreePath, java.lang.Object)
	 */
	public void valueForPathChanged(TreePath path, Object value) {
	}

	/* (non-Javadoc)
	 * @see javax.swing.tree.TreeModel#getIndexOfChild(java.lang.Object, java.lang.Object)
	 */
	public int getIndexOfChild(Object parent, Object child) {
		TaxonomyAddr parentAddr = (TaxonomyAddr)parent;
		int index = parentAddr.indexOfElement(child);
		if(index != -1) { return index; }
		return parentAddr.indexOfAddr((TaxonomyAddr)child) + parentAddr.getNumElements();
	}

	/* (non-Javadoc)
	 * @see javax.swing.tree.TreeModel#addTreeModelListener(javax.swing.event.TreeModelListener)
	 */
	public void addTreeModelListener(TreeModelListener listener) {
		listeners.addLast(listener);
	}

	/* (non-Javadoc)
	 * @see javax.swing.tree.TreeModel#removeTreeModelListener(javax.swing.event.TreeModelListener)
	 */
	public void removeTreeModelListener(TreeModelListener listener) {
		listeners.remove(listener);
	}

}

class TaxonomyAddr<X> { 
	
	private Object[] path;
	private TaxonomyAddr parent;
	private LinkedList<String> addr;
	private Vector<TaxonomyAddr> nextAddrs;
	private Vector<X> elmts;
	
	public TaxonomyAddr(Taxonomy<X> base) {
		parent = null;
		path = new Object[1]; path[0] = this;
		addr = new LinkedList<String>();
		SortedSet<String> ss = new TreeSet<String>(base.getAddrs());
		nextAddrs = new Vector<TaxonomyAddr>();
		for(String s : ss) { 
			nextAddrs.add(new TaxonomyAddr<X>(this, addr, s, base));
		}
		elmts = new Vector<X>(new TreeSet<X>(base.getImmediateElements()));
	}
	
	public TaxonomyAddr(TaxonomyAddr p, Collection<String> a, String next, Taxonomy<X> base) {
		parent = p;
		path = new Object[p.path.length+1];
		for(int i = 0; i < p.path.length; i++) { path[i] = p.path[i]; }
		path[path.length-1] = this;
		
		addr = new LinkedList<String>(a);
		addr.addLast(next);
		Taxonomy<X> sub = getTaxonomy(base);
		SortedSet<String> ss = new TreeSet<String>(sub.getAddrs());
		nextAddrs = new Vector<TaxonomyAddr>();
		for(String s : ss) { 
			nextAddrs.add(new TaxonomyAddr<X>(this, addr, s, base));
		}
		elmts = new Vector<X>(new TreeSet<X>(sub.getImmediateElements()));
	}
	
	public TaxonomyAddr(TaxonomyAddr p, Collection<String> a, String next) {
		parent = p;
		path = new Object[p.path.length+1];
		for(int i = 0; i < p.path.length; i++) { path[i] = p.path[i]; }
		path[path.length-1] = this;
		
		addr = new LinkedList<String>(a);
		addr.addLast(next);
		elmts = new Vector<X>();
		nextAddrs = new Vector<TaxonomyAddr>();
	}
	
	public void handleTaxonomyEvent(TaxonomyTreeModel<X> ttm, TaxonomyEvent te) {
		System.out.println(getAddrString() + " [" + addr.size() + "] :: " + te.toString());
		int[] childinds = new int[1];
		Object[] children = new Object[1];

		if(te.getAddrLength() == 0) {
			children[0] = te.getLeaf();
			X v = (X)te.getLeaf();
			
			switch(te.getType()) { 
			case TaxonomyEvent.ELEMENT_ADDED:
				elmts.insertElementAt(v, 0);
				childinds[0] = 0;
				ttm.fireTreeNodesInserted(path, childinds, children);
				break;
			case TaxonomyEvent.ELEMENT_REMOVED:
				childinds[0] = indexOfElement(v);
				elmts.removeElementAt(childinds[0]);
				ttm.fireTreeNodesRemoved(path, childinds, children);
				break;
			case TaxonomyEvent.ELEMENT_CHANGED:
				childinds[0] = indexOfElement(v);
				ttm.fireTreeNodesChanged(path, childinds, children);
				break;

			case TaxonomyEvent.TAXON_ADDED:
			case TaxonomyEvent.TAXON_REMOVED:
			default:
				throw new IllegalArgumentException("type: " + te.getType());
			}
		} else { 
			String addr_key = te.removeAddr();
			System.out.println(">>> Addr key: " + addr_key);
			Iterator<TaxonomyAddr> taItr = nextAddrs.iterator();
			boolean found = false;
			while(taItr.hasNext() && !found) { 
				TaxonomyAddr nextA = taItr.next();
				if(found = nextA.addr.getLast().equals(addr_key)) { 
					nextA.handleTaxonomyEvent(ttm, te);
				}
			}

			if(!found) { 
				System.out.println(">> No match found to addr_key: {" + addr_key + "}");
				TaxonomyAddr newAddr = new TaxonomyAddr<X>(this, addr, addr_key);
				childinds[0] = elmts.size(); 
				System.out.println("new index: " + childinds[0]);
				children[0] = newAddr;
				nextAddrs.insertElementAt(newAddr, 0);
				ttm.fireTreeNodesInserted(path, childinds, children);
				newAddr.handleTaxonomyEvent(ttm, te);
			}
		}
	}
	
	public int indexOfAddr(TaxonomyAddr next) { 
		for(int i = 0; i < nextAddrs.size(); i++) { 
			if(nextAddrs.get(i).equals(next)) { return i; }
		}
		return -1;
	}
	
	public int indexOfElement(X v) { 
		for(int i = 0; i < elmts.size(); i++) { 
			if(elmts.get(i).equals(v)) { return i; }
		}
		return -1;		
	}
	
	public TaxonomyAddr getNextAddr(int i) { return nextAddrs.get(i); }
	public X getElement(int i) { return elmts.get(i); }
	public int getNumSubTaxa() { return nextAddrs.size(); }
	public int getNumElements() { return elmts.size(); }
	public int length() { return addr.size(); }
    public LinkedList<String> getAddr() { return addr; }
	
	public Taxonomy<X> getTaxonomy(Taxonomy<X> base) { 
		Taxonomy<X> current = base;
		for(String a : addr) { 
			if(!current.hasAddr(a)) { 
				throw new IllegalArgumentException(a);
			}
			current = current.getSubTaxonomy(a);
		}
		return current;
	}
	
	public String getLast() { return addr.getLast(); }
	
	public boolean equals(Object o) { 
		if(!(o instanceof TaxonomyAddr)) { return false; }
		TaxonomyAddr ta = (TaxonomyAddr)o;
		if(addr.size() != ta.addr.size()) { return false; }
		for(int i = 0; i < addr.size(); i++) { 
			if(!(addr.get(i).equals(ta.addr.get(i)))) { return false; }
		}
		return true;
	}
	
	public int hashCode() { 
		int code = 17;
		for(String a : addr) { 
			code += a.hashCode(); code *= 37;
		}
		return code;
	}
	
	public String getAddrString() { 
		StringBuilder sb = new StringBuilder();
		boolean first = true;
		for(String a : addr) { 
			if(!first) { sb.append(","); }
			sb.append(a);
			first = false;
		}
		return sb.toString();
	}
	
	public String toString() { 
		if(addr.size() == 0) { return "Root"; }
		return addr.getLast();
	}
}
